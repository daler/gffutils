import copy
import os
import simplejson
from collections import OrderedDict
import constants
import bins
import gzip
import time
import tempfile
import gffutils

HERE = os.path.dirname(os.path.abspath(__file__))

def example_filename(fn):
    return os.path.join(HERE, 'test', 'data', fn)


def make_query(args, other=None, limit=None, strand=None, featuretype=None,
               extra=None, order_by=None, reverse=False,
               completely_within=False):
    """
    This function composes queries given some commonly-used kwargs that can be
    passed to FeatureDB methods (like .parents(), .children(), .all_features(),
    .features_of_type()).


    It also provides support for additional JOINs etc (`other`) and extra
    conditional clauses (`extra`).

    For example, `other` is used in FeatureDB._relation, where `other` is the
    JOIN substatment.

    `extra` is also used in FeatureDB._relation, where `extra` is the
    "relations.level = ?" substatment.

    `args` can be pre-filled with args that are passed to `other` and `extra`.
    """

    _QUERY = ("{_SELECT} {OTHER} {EXTRA} {FEATURETYPE} "
              "{LIMIT} {STRAND} {ORDER_BY}")

    # Will be used later as _QUERY.format(**d).  Default is just _SELECT, which
    # returns all records in the features table.
    # (Recall that constants._SELECT gets the fields in the order needed to
    # reconstruct a Feature)
    d = dict(_SELECT=constants._SELECT, OTHER="", FEATURETYPE="", LIMIT="",
             STRAND="", ORDER_BY="", EXTRA="")

    if other:
        d['OTHER'] = other
    if extra:
        d['EXTRA'] = extra

    # Quick check for args provided
    required_args = (d['EXTRA'] + d['OTHER']).count('?')
    if len(args) != required_args:
        raise ValueError('Not enough args (%s) for subquery' % args)

    # Below, if a kwarg is specified, then we create sections of the query --
    # appending to args as necessary.  Importantly, the order in which things
    # are processed here is the same as the order of the placeholders in
    # _QUERY.

    # e.g., "featuretype = 'exon'"
    #
    # or, "featuretype IN ('exon', 'CDS')"
    if featuretype:
        if isinstance(featuretype, basestring):
            d['FEATURETYPE'] = "features.featuretype = ?"
            args.append(featuretype)
        else:
            d['FEATURETYPE'] = (
                "features.featuretype IN  (%s)"
                % (','.join(["?" for _ in featuretype]))
            )
            args.extend(featuretype)

    # e.g., "seqid = 'chr2L' AND start > 1000 AND end < 5000"
    if limit:
        if isinstance(limit, basestring):
            seqid, startstop = limit.split(':')
            start, end = startstop.split('-')
        else:
            seqid, start, end = limit

        # Identify bins
        _bins = bins.bins(int(start), int(end), one=False)

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

        # add bin clause
        d['LIMIT'] += " AND features.bin IN (%s)" % (','.join(map(str, _bins)))

    # e.g., "strand = '+'"
    if strand:
        d['STRAND'] = "features.strand = ?"
        args.append(strand)

    # e.g. "ORDER BY seqid, start DESC"
    if order_by:
        if isinstance(order_by, basestring):
            order_by = [order_by]
        for k in order_by:
            if k not in constants._gffkeys_extra:
                raise ValueError("%s not a valid column" % (k))
        order_by = ','.join(order_by)
        if reverse:
            direction = 'DESC'
        else:
            direction = 'ASC'
        d['ORDER_BY'] = 'ORDER BY %s %s' % (order_by, direction)

    # Ensure only one "WHERE" is included; the rest get "AND "
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


class FeatureNotFoundError(Exception):
    """
    Error to be raised when an ID is not in the database.
    """
    def __init__(self, feature_id):
        Exception.__init__(self)
        self.feature_id = feature_id

    def __str__(self):
        return self.feature_id


class DuplicateIDError(Exception):
    pass


class AttributeStringError(Exception):
    pass


def _jsonify(x):
    """Use most compact form of JSON"""
    return simplejson.dumps(x, separators=(',', ':'))


def _unjsonify(x):
    """Convert JSON string to an ordered defaultdict."""
    return simplejson.loads(
        x, object_pairs_hook=DefaultListOrderedDict)


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
        Attribute dictionaries, assumed to be at least DefaultDict of lists,
        possibly DefaultOrderedDict.  If ordered, the first dictionary's key
        order takes precedence.
    """
    new_d = copy.deepcopy(attr1)
    for k in attr1.keys():
        if k in attr2:
            new_d[k].extend(attr2[k])
    for k in attr2.keys():
        if k not in attr1:
            new_d[k].extend(attr2[k])

    for k, v in new_d.items():
        new_d[k] = list(set(v))
    return new_d


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


class DefaultOrderedDict(OrderedDict):
    """
    OrderedDict that is also a defaultdict.

    From http://stackoverflow.com/a/6190500.
    (using callable() instead of Callable as per the comments to that answer)
    """
    def __init__(self, default_factory=None, *a, **kw):
        if (
            default_factory is not None and
            not callable(default_factory)
        ):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.items())


class DefaultListOrderedDict(DefaultOrderedDict):
    """
    defaultdict(list), but ordered.
    """
    def __init__(self, *a, **kw):
        super(DefaultListOrderedDict, self).__init__(list, *a, **kw)
        self.default_factory = list

    def __copy__(self):
        return type(self)(self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(copy.deepcopy(self.items()))


def sanitize_gff(db_fname):
    """
    Sanitize given GFF db. Return a generator of sanitized
    records.

    TODO: Do something with negative coordinates?
    """
    db = gffutils.FeatureDB(db_fname)
    for rec in db.all_features():
        if rec.start > rec.stop:
            rec.start, rec.stop = \
                rec.stop, rec.start
        yield rec

##
## Helpers for gffutils-cli
##
## TODO: move clean_gff here?
##
def get_db_fname(gff_fname,
                 ext=".db"):
    """
    Get db fname for GFF. If the database has a .db file,
    return that. Otherwise, create a named temporary file,
    serialize the db to that, and return the name of the file.

    TODO: Add parameter to control whether or not the .db
    is created or not.
    """
    if not os.path.isfile(gff_fname):
        # Not sure how we should deal with errors normally in
        # gffutils -- Ryan?
        raise Exception, "GFF %s does not exist." %(gff_fname)
        return None
    candidate_db_fname = "%s.%s" %(gff_fname, ext)
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
    print "Creating db for %s" %(gff_fname)
    t1 = time.time()
    gffutils.create_db(gff_fname, db_fname.name,
                       merge_strategy="merge",
                       verbose=False)
    t2 = time.time()
    print "  - Took %.2f seconds" %(t2 - t1)
    return db_fname.name


if __name__ == "__main__":
    d = DefaultListOrderedDict([('a', 1), ('b', 2)])
