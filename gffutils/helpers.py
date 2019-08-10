import copy
import sys
import os
import simplejson as json
import time
import tempfile
import six
from gffutils import constants
from gffutils import bins
import gffutils
from gffutils import gffwriter
from gffutils import parser
from gffutils.attributes import dict_class

HERE = os.path.dirname(os.path.abspath(__file__))


def example_filename(fn):
    """
    Return the full path of a data file that ships with gffutils.
    """
    return os.path.join(HERE, 'test', 'data', fn)


def infer_dialect(attributes):
    """
    Infer the dialect based on the attributes.

    Parameters
    ----------
    attributes : str or iterable
        A single attributes string from a GTF or GFF line, or an iterable of
        such strings.

    Returns
    -------
    Dictionary representing the inferred dialect
    """
    if isinstance(attributes, six.string_types):
        attributes = [attributes]
    dialects = [parser._split_keyvals(i)[1] for i in attributes]
    return _choose_dialect(dialects)


def _choose_dialect(dialects):
    """
    Given a list of dialects, choose the one to use as the "canonical" version.

    If `dialects` is an empty list, then use the default GFF3 dialect

    Parameters
    ----------
    dialects : iterable
        iterable of dialect dictionaries

    Returns
    -------
    dict
    """
    # NOTE: can use helpers.dialect_compare if you need to make this more
    # complex....

    # For now, this function favors the first dialect, and then appends the
    # order of additional fields seen in the attributes of other lines giving
    # priority to dialects that come first in the iterable.
    if len(dialects) == 0:
        return constants.dialect
    final_order = []
    for dialect in dialects:
        for o in dialect['order']:
            if o not in final_order:
                final_order.append(o)
    dialect = dialects[0]
    dialect['order'] = final_order
    return dialect


def make_query(args, other=None, limit=None, strand=None, featuretype=None,
               extra=None, order_by=None, reverse=False,
               completely_within=False):
    """
    Multi-purpose, bare-bones ORM function.

    This function composes queries given some commonly-used kwargs that can be
    passed to FeatureDB methods (like .parents(), .children(), .all_features(),
    .features_of_type()).  It handles, in one place, things like restricting to
    featuretype, limiting to a genomic range, limiting to one strand, or
    returning results ordered by different criteria.

    Additional filtering/subsetting/sorting behavior should be added here.

    (Note: this ended up having better performance (and flexibility) than
    sqlalchemy)

    This function also provides support for additional JOINs etc (supplied via
    the `other` kwarg) and extra conditional clauses (`extra` kwarg).  See the
    `_QUERY` var below for the order in which they are used.

    For example, FeatureDB._relation uses `other` to supply the JOIN
    substatment, and that same method also uses `extra` to supply the
    "relations.level = ?" substatment (see the source for FeatureDB._relation
    for more details).

    `args` contains the arguments that will ultimately be supplied to the
    sqlite3.connection.execute function.  It may be further populated below --
    for example, if strand="+", then the query will include a strand clause,
    and the strand will be appended to the args.

    `args` can be pre-filled with args that are passed to `other` and `extra`.
    """

    _QUERY = ("{_SELECT} {OTHER} {EXTRA} {FEATURETYPE} "
              "{LIMIT} {STRAND} {ORDER_BY}")

    # Construct a dictionary `d` that will be used later as _QUERY.format(**d).
    # Default is just _SELECT, which returns all records in the features table.
    # (Recall that constants._SELECT gets the fields in the order needed to
    # reconstruct a Feature)
    d = dict(_SELECT=constants._SELECT, OTHER="", FEATURETYPE="", LIMIT="",
             STRAND="", ORDER_BY="", EXTRA="")

    if other:
        d['OTHER'] = other
    if extra:
        d['EXTRA'] = extra

    # If `other` and `extra` take args (that is, they have "?" in them), then
    # they should have been provided in `args`.
    required_args = (d['EXTRA'] + d['OTHER']).count('?')
    if len(args) != required_args:
        raise ValueError('Not enough args (%s) for subquery' % args)

    # Below, if a kwarg is specified, then we create sections of the query --
    # appending to args as necessary.
    #
    # IMPORTANT: the order in which things are processed here is the same as
    # the order of the placeholders in _QUERY.  That is, we need to build the
    # args in parallel with the query to avoid putting the wrong args in the
    # wrong place.

    if featuretype:
        # Handle single or iterables of featuretypes.
        #
        # e.g., "featuretype = 'exon'"
        #
        # or, "featuretype IN ('exon', 'CDS')"
        if isinstance(featuretype, six.string_types):
            d['FEATURETYPE'] = "features.featuretype = ?"
            args.append(featuretype)
        else:
            d['FEATURETYPE'] = (
                "features.featuretype IN  (%s)"
                % (','.join(["?" for _ in featuretype]))
            )
            args.extend(featuretype)

    if limit:
        # Restrict to a genomic region.  Makes use of the UCSC binning strategy
        # for performance.
        #
        # `limit` is a string or a tuple of (chrom, start, stop)
        #
        # e.g., "seqid = 'chr2L' AND start > 1000 AND end < 5000"
        if isinstance(limit, six.string_types):
            seqid, startstop = limit.split(':')
            start, end = startstop.split('-')
        else:
            seqid, start, end = limit

        # Identify possible bins
        _bins = bins.bins(int(start), int(end), one=False)

        # Use different overlap conditions
        if completely_within:
            d['LIMIT'] = (
                "features.seqid = ? AND features.start >= ? "
                "AND features.end <= ?"
            )
            args.extend([seqid, start, end])

        else:
            d['LIMIT'] = (
                "features.seqid = ? AND features.start <= ? "
                "AND features.end >= ?"
            )
            # Note order (end, start)
            args.extend([seqid, end, start])

        # Add bin clause. See issue #45.
        if len(_bins) < 900:
            d['LIMIT'] += " AND features.bin IN (%s)" % (','.join(map(str, _bins)))

    if strand:
        # e.g., "strand = '+'"
        d['STRAND'] = "features.strand = ?"
        args.append(strand)

    # TODO: implement file_order!
    valid_order_by = constants._gffkeys_extra + ['file_order', 'length']
    _order_by = []
    if order_by:
        # Default is essentially random order.
        #
        # e.g. "ORDER BY seqid, start DESC"
        if isinstance(order_by, six.string_types):
            _order_by.append(order_by)

        else:
            for k in order_by:
                if k not in valid_order_by:
                    raise ValueError("%s not a valid order-by value in %s"
                                     % (k, valid_order_by))

                # There's no length field, so order by end - start
                if k == 'length':
                    k = '(end - start)'

                _order_by.append(k)

        _order_by = ','.join(_order_by)
        if reverse:
            direction = 'DESC'
        else:
            direction = 'ASC'
        d['ORDER_BY'] = 'ORDER BY %s %s' % (_order_by, direction)

    # Ensure only one "WHERE" is included; the rest get "AND ".  This is ugly.
    where = False
    if "where" in d['OTHER'].lower():
        where = True
    for i in ['EXTRA', 'FEATURETYPE', 'LIMIT', 'STRAND']:
        if d[i]:
            if not where:
                d[i] = "WHERE " + d[i]
                where = True
            else:
                d[i] = "AND " + d[i]

    return _QUERY.format(**d), args


def _bin_from_dict(d):
    """
    Given a dictionary yielded by the parser, return the genomic "UCSC" bin
    """
    try:
        start = int(d['start'])
        end = int(d['end'])
        return bins.bins(start, end, one=True)

    # e.g., if "."
    except ValueError:
        return None


def _jsonify(x):
    """Use most compact form of JSON"""
    if isinstance(x, dict_class):
        return json.dumps(x._d, separators=(',', ':'))
    return json.dumps(x, separators=(',', ':'))


def _unjsonify(x, isattributes=False):
    """Convert JSON string to an ordered defaultdict."""
    if isattributes:
        obj = json.loads(x)
        return dict_class(obj)
    return json.loads(x)


def _feature_to_fields(f, jsonify=True):
    """
    Convert feature to tuple, for faster sqlite3 import
    """
    x = []
    for k in constants._keys:
        v = getattr(f, k)
        if jsonify and (k in ('attributes', 'extra')):
            x.append(_jsonify(v))
        else:
            x.append(v)
    return tuple(x)


def _dict_to_fields(d, jsonify=True):
    """
    Convert dict to tuple, for faster sqlite3 import
    """
    x = []
    for k in constants._keys:
        v = d[k]
        if jsonify and (k in ('attributes', 'extra')):
            x.append(_jsonify(v))
        else:
            x.append(v)
    return tuple(x)


def asinterval(feature):
    """
    Converts a gffutils.Feature to a pybedtools.Interval
    """
    import pybedtools
    return pybedtools.create_interval_from_list(str(feature).split('\t'))


def merge_attributes(attr1, attr2):
    """
    Merges two attribute dictionaries into a single dictionary.

    Parameters
    ----------
    `attr1`, `attr2` : dict

    Returns
    -------
    dict
    """

    new_d = copy.deepcopy(attr1)
    new_d.update(copy.deepcopy(attr2))

    #all of attr2 key : values just overwrote attr1, fix it
    for k, v in new_d.items():
        if not isinstance(v, list):
            new_d[k] = [v]

    for k, v in six.iteritems(attr1):
        if k in attr2:
            if not isinstance(v, list):
                v = [v]
            new_d[k].extend(v)
    return dict((k, sorted(set(v))) for k, v in new_d.items())


def dialect_compare(dialect1, dialect2):
    """
    Compares two dialects.
    """
    orig = set(dialect1.items())
    new = set(dialect2.items())
    return dict(
        added=dict(list(new.difference(orig))),
        removed=dict(list(orig.difference(new)))
    )


def sanitize_gff_db(db, gid_field="gid"):
    """
    Sanitize given GFF db. Returns a sanitized GFF db.

    Sanitizing means:

    - Ensuring that start < stop for all features
    - Standardizing gene units by adding a 'gid' attribute
      that makes the file grep-able

    TODO: Do something with negative coordinates?
    """
    def sanitized_iterator():
        # Iterate through the database by each gene's records
        for gene_recs in db.iter_by_parent_childs():
            # The gene's ID
            gene_id = gene_recs[0].id
            for rec in gene_recs:
                # Fixup coordinates if necessary
                if rec.start > rec.stop:
                    rec.start, rec.stop = rec.stop, rec.start
                # Add a gene id field to each gene's records
                rec.attributes[gid_field] = [gene_id]
                yield rec
    # Return sanitized GFF database
    sanitized_db = \
        gffutils.create_db(sanitized_iterator(), ":memory:",
                           verbose=False)
    return sanitized_db


def sanitize_gff_file(gff_fname,
                      in_memory=True,
                      in_place=False):
    """
    Sanitize a GFF file.
    """
    db = None
    if is_gff_db(gff_fname):
        # It's a database filename, so load it
        db = gffutils.FeatureDB(gff_fname)
    else:
        # Need to create a database for file
        if in_memory:
            db = gffutils.create_db(gff_fname, ":memory:",
                                    verbose=False)
        else:
            db = get_gff_db(gff_fname)
    if in_place:
        gff_out = gffwriter.GFFWriter(gff_fname,
                                      in_place=in_place)
    else:
        gff_out = gffwriter.GFFWriter(sys.stdout)
    sanitized_db = sanitize_gff_db(db)
    for gene_rec in sanitized_db.all_features(featuretype="gene"):
        gff_out.write_gene_recs(sanitized_db, gene_rec.id)
    gff_out.close()


def annotate_gff_db(db):
    """
    Annotate a GFF file by cross-referencing it with another GFF
    file, e.g. one containing gene models.
    """
    pass


def is_gff_db(db_fname):
    """
    Return True if the given filename is a GFF database.

    For now, rely on .db extension.
    """
    if not os.path.isfile(db_fname):
        return False
    if db_fname.endswith(".db"):
        return True
    return False


def to_unicode(obj, encoding='utf-8'):
    if isinstance(obj, six.string_types):
        if not isinstance(obj, six.text_type):
            obj = six.text_type(obj, encoding)
    return obj


def canonical_transcripts(db, fasta_filename):
    import pyfaidx
    fasta = pyfaidx.Fasta(fasta_filename, as_raw=True)
    for gene in db.features_of_type('gene'):

        # exons_list will contain (CDS_length, total_length, transcript, [exons]) tuples.
        exon_list = []
        for ti, transcript in enumerate(db.children(gene, level=1)):
            cds_len = 0
            total_len = 0
            exons = list(db.children(transcript, level=1))
            for exon in exons:
                exon_length = len(exon)
                if exon.featuretype == 'CDS':
                    cds_len += exon_length
                total_len += exon_length

            exon_list.append((cds_len, total_len, transcript, exons))

        # If we have CDS, then use the longest coding transcript
        if max(i[0] for i in exon_list) > 0:
            best = sorted(exon_list)[0]
        # Otherwise, just choose the longest
        else:
            best = sorted(exon_list, lambda x: x[1])[0]

        print(best)

        canonical_exons = best[-1]
        transcript = best[-2]
        seqs = [i.sequence(fasta) for i in canonical_exons]
        yield transcript, ''.join(seqs)



##
## Helpers for gffutils-cli
##
## TODO: move clean_gff here?
##
def get_gff_db(gff_fname,
               ext=".db"):
    """
    Get db for GFF file. If the database has a .db file,
    load that. Otherwise, create a named temporary file,
    serialize the db to that, and return the loaded database.
    """
    if not os.path.isfile(gff_fname):
        # Not sure how we should deal with errors normally in
        # gffutils -- Ryan?
        raise ValueError("GFF %s does not exist." % (gff_fname))
    candidate_db_fname = "%s.%s" % (gff_fname, ext)
    if os.path.isfile(candidate_db_fname):
        # Standard .db file found, so return it
        return candidate_db_fname
    # Otherwise, we need to create a temporary but non-deleted
    # file to store the db in. It'll be up to the user
    # of the function the delete the file when done.
    ## NOTE: Ryan must have a good scheme for dealing with this
    ## since pybedtools does something similar under the hood, i.e.
    ## creating temporary files as needed without over proliferation
    db_fname = tempfile.NamedTemporaryFile(delete=False)
    # Create the database for the gff file (suppress output
    # when using function internally)
    print("Creating db for %s" % (gff_fname))
    t1 = time.time()
    db = gffutils.create_db(gff_fname, db_fname.name,
                            merge_strategy="merge",
                            verbose=False)
    t2 = time.time()
    print("  - Took %.2f seconds" % (t2 - t1))
    return db
