import copy
import warnings
import collections
import tempfile
import sys
import os
import sqlite3
import six
from textwrap import dedent
from gffutils import constants
from gffutils import version
from gffutils import bins
from gffutils import helpers
from gffutils import feature
from gffutils import interface
from gffutils import iterators

import logging

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)


def deprecation_handler(kwargs):
    """
    As things change from version to version, deal with them here.
    """
    # After reconsidering, let's leave `infer_gene_extent` for another release.
    # But when it's time to deprecate it, use this code:
    if 0:
        if 'infer_gene_extent' in kwargs:
            raise ValueError(
                "'infer_gene_extent' is deprecated as of version 0.8.4 in "
                "favor of more granular control over inferring genes and/or "
                "transcripts.  The previous default was "
                "'infer_gene_extent=True`, which corresponds to the new "
                "defaults "
                "'disable_infer_genes=False' and "
                "'disable_infer_transcripts=False'. Please see the docstring "
                "for gffutils.create_db for details.")
    if len(kwargs) > 0:
        raise TypeError("unhandled kwarg in %s" % kwargs)


class _DBCreator(object):
    def __init__(self, data, dbfn, force=False, verbose=False, id_spec=None,
                 merge_strategy='error', checklines=10, transform=None,
                 force_dialect_check=False, from_string=False, dialect=None,
                 default_encoding='utf-8',
                 disable_infer_genes=False,
                 disable_infer_transcripts=False,
                 infer_gene_extent=True,
                 force_merge_fields=None,
                 text_factory=sqlite3.OptimizedUnicode,
                 pragmas=constants.default_pragmas, _keep_tempfiles=False,
                 directives=None,
                 **kwargs):
        """
        Base class for _GFFDBCreator and _GTFDBCreator; see create_db()
        function for docs
        """
        self._keep_tempfiles = _keep_tempfiles
        if force_merge_fields is None:
            force_merge_fields = []
        if merge_strategy == 'merge':
            if set(['start', 'end']).intersection(force_merge_fields):
                raise ValueError("Can't merge start/end fields since "
                                 "they must be integers")
            warn = set(force_merge_fields)\
                .intersection(['frame', 'strand'])
            for w in warn:
                warnings.warn(
                    "%s field will be merged for features with the same ID; "
                    "this may result in unusable features." % w)

        self.force_merge_fields = force_merge_fields
        self.pragmas = pragmas
        self.merge_strategy = merge_strategy
        self.default_encoding = default_encoding
        if directives is None:
            directives = []
        self.directives = directives

        if not infer_gene_extent:
            warnings.warn("'infer_gene_extent' will be deprecated. For now, "
                          "the following equivalent values were automatically "
                          "set: 'disable_infer_genes=True', "
                          "'disable_infer_transcripts=True'. Please use these "
                          "instead in the future.")
            disable_infer_genes = True
            disable_infer_transcripts = True

        self.disable_infer_genes = disable_infer_genes
        self.disable_infer_transcripts = disable_infer_transcripts

        if force:
            if os.path.exists(dbfn):
                os.unlink(dbfn)
        self.dbfn = dbfn
        self.id_spec = id_spec
        if isinstance(dbfn, six.string_types):
            conn = sqlite3.connect(dbfn)
        else:
            conn = dbfn
        self.conn = conn
        self.conn.row_factory = sqlite3.Row
        self.set_verbose(verbose)

        if text_factory is not None:
            if self.verbose == 'debug':
                logger.debug('setting text factory to %s' % text_factory)
            self.conn.text_factory = text_factory
        self._data = data

        self._orig_logger_level = logger.level

        self.iterator = iterators.DataIterator(
            data=data, checklines=checklines, transform=transform,
            force_dialect_check=force_dialect_check, from_string=from_string,
            dialect=dialect
        )
        if '_autoincrements' in kwargs:
            self._autoincrements = kwargs['_autoincrements']
        else:
            self._autoincrements = collections.defaultdict(int)

    def set_verbose(self, verbose=None):
        if verbose == 'debug':
            logger.setLevel(logging.DEBUG)
        elif verbose:
            logger.setLevel(logging.INFO)
        else:
            logger.setLevel(logging.ERROR)
        self.verbose = verbose

    def _increment_featuretype_autoid(self, key):
        self._autoincrements[key] += 1
        return '%s_%s' % (key, self._autoincrements[key])

    def _id_handler(self, f):
        """
        Given a Feature from self.iterator, figure out what the ID should be.

        This uses `self.id_spec` identify the ID.
        """

        # If id_spec is a string, convert to iterable for later
        if isinstance(self.id_spec, six.string_types):
            id_key = [self.id_spec]

        elif hasattr(self.id_spec, '__call__'):
            id_key = [self.id_spec]

        # If dict, then assume it's a feature -> attribute mapping, e.g.,
        # {'gene': 'gene_id'} for GTF
        elif isinstance(self.id_spec, dict):
            try:
                id_key = self.id_spec[f.featuretype]
                if isinstance(id_key, six.string_types):
                    id_key = [id_key]

            # Otherwise, use default auto-increment.
            except KeyError:
                return self._increment_featuretype_autoid(f.featuretype)

        # Otherwise assume it's an iterable.
        else:
            id_key = self.id_spec

        # Then try them in order, returning the first one that works:
        for k in id_key:

            if hasattr(k, '__call__'):
                _id = k(f)
                if _id:
                    if _id.startswith('autoincrement:'):
                        return self._increment_featuretype_autoid(_id[14:])
                    return _id
            else:
                # use GFF fields rather than attributes for cases like :seqid:
                # or :strand:
                if (len(k) > 3) and (k[0] == ':') and (k[-1] == ':'):
                    # No [0] here -- only attributes key/vals are forced into
                    # lists, not standard GFF fields.
                    return getattr(f, k[1:-1])
                else:
                    try:
                        return f.attributes[k][0]
                    except (KeyError, IndexError):
                        pass
        # If we get here, then default autoincrement
        return self._increment_featuretype_autoid(f.featuretype)

    def _get_feature(self, ID):
        c = self.conn.cursor()
        results = c.execute(
            constants._SELECT + ' WHERE id = ?', (ID,)).fetchone()
        return feature.Feature(dialect=self.iterator.dialect, **results)

    def _do_merge(self, f, merge_strategy, add_duplicate=False):
        """
        Different merge strategies upon name conflicts.

        "error":
            Raise error

        "warning"
            Log a warning, which indicates that all future instances of the
            same ID will be ignored

        "merge":
            Combine old and new attributes -- but only if everything else
            matches; otherwise error. This can be slow, but is thorough.

        "create_unique":
            Autoincrement based on the ID, always creating a new ID.

        "replace":
            Replaces existing database feature with `f`.
        """
        if merge_strategy == 'error':
            raise ValueError("Duplicate ID {0.id}".format(f))

        if merge_strategy == 'warning':
            logger.warning(
                "Duplicate lines in file for id '{0.id}'; "
                "ignoring all but the first".format(f))
            return None, merge_strategy

        elif merge_strategy == 'replace':
            return f, merge_strategy

        # This is by far the most complicated strategy.
        elif merge_strategy == 'merge':

            # Recall that if we made it to this method, there was at least one
            # ID collision.

            # This will eventually contain the features that match ID AND that
            # match non-attribute fields like start, stop, strand, etc.
            features_to_merge = []

            # Iterate through all features that have the same ID according to
            # the id_spec provided.
            if self.verbose == "debug":
                logger.debug('candidates with same idspec: %s'
                             % ([i.id for i in self._candidate_merges(f)]))

            # If force_merge_fields was provided, don't pay attention to these
            # fields if they're different. We are assuming the attributes field
            # will be different, hence the [:-1]
            _gffkeys_to_check = list(
                set(constants._gffkeys[:-1])
                .difference(self.force_merge_fields))

            for existing_feature in self._candidate_merges(f):
                # Check other GFF fields (if not specified in
                # self.force_merge_fields) to make sure they match.
                other_attributes_same = True
                for k in _gffkeys_to_check:
                    if getattr(existing_feature, k) != getattr(f, k):
                        other_attributes_same = False
                        break

                if other_attributes_same:
                    # All the other GFF fields match.  So this existing feature
                    # should be merged.
                    features_to_merge.append(existing_feature)
                    if self.verbose == 'debug':
                        logger.debug(
                            'same attributes between:\nexisting: %s'
                            '\nthis    : %s'
                            % (existing_feature, f))
                else:
                    # The existing feature's GFF fields don't match, so don't
                    # append anything.
                    if self.verbose == 'debug':
                        logger.debug(
                            'different attributes between:\nexisting: %s\n'
                            'this    : %s'
                            % (existing_feature, f))

            if (len(features_to_merge) == 0):
                # No merge candidates found, so we should make a new ID for
                # this feature. This can happen when idspecs match, but other
                # fields (like start/stop) are different.  Call this method
                # again, but using the "create_unique" strategy, and then
                # record the newly-created ID in the duplicates table.
                orig_id = f.id
                uniqued_feature, merge_strategy = self._do_merge(
                    f, merge_strategy='create_unique')
                self._add_duplicate(orig_id, uniqued_feature.id)
                return uniqued_feature, merge_strategy

            # Whoo! Found some candidates to merge.
            else:
                if self.verbose == 'debug':
                    logger.debug('num candidates: %s' % len(features_to_merge))

                # This is the attributes dictionary we'll be modifying.
                merged_attributes = copy.deepcopy(f.attributes)

                # Keep track of non-attribute fields (this will be an empty
                # dict if no force_merge_fields)
                final_fields = dict(
                    [(field, set([getattr(f, field)]))
                     for field in self.force_merge_fields])

                # Update the attributes
                for existing_feature in features_to_merge:
                    if self.verbose == 'debug':
                        logger.debug(
                            '\nmerging\n\n%s\n%s\n' % (f, existing_feature))
                    for k in existing_feature.attributes.keys():
                        v = merged_attributes.setdefault(k, [])
                        v.extend(existing_feature[k])
                        merged_attributes[k] = v

                    # Update the set of non-attribute fields found so far
                    for field in self.force_merge_fields:
                        final_fields[field].update(
                            [getattr(existing_feature, field)])

                # Set the merged attributes
                for k, v in merged_attributes.items():
                    merged_attributes[k] = list(set(v))
                existing_feature.attributes = merged_attributes

                # Set the final merged non-attributes
                for k, v in final_fields.items():
                    setattr(existing_feature, k, ','.join(sorted(map(str, v))))

                if self.verbose == 'debug':
                    logger.debug('\nMERGED:\n%s' % existing_feature)
                return existing_feature, merge_strategy

        elif merge_strategy == 'create_unique':
            f.id = self._increment_featuretype_autoid(f.id)
            return f, merge_strategy
        else:
            raise ValueError("Invalid merge strategy '%s'"
                             % (merge_strategy))

    def _add_duplicate(self, idspecid, newid):
        """
        Adds a duplicate ID (as identified by id_spec) and its new ID to the
        duplicates table so that they can be later searched for merging.

        Parameters
        ----------
        newid : str
            The primary key used in the features table

        idspecid : str
            The ID identified by id_spec
        """
        c = self.conn.cursor()
        try:
            c.execute(
                '''
                INSERT INTO duplicates
                (idspecid, newid)
                VALUES (?, ?)''',
                (idspecid, newid))
        except sqlite3.ProgrammingError:
            c.execute(
                '''
                INSERT INTO duplicates
                (idspecid, newid)
                VALUES (?, ?)''',
                (idspecid.decode(self.default_encoding),
                 newid.decode(self.default_encoding))
            )
        if self.verbose == 'debug':
            logger.debug('added id=%s; new=%s' % (idspecid, newid))
        self.conn.commit()

    def _candidate_merges(self, f):
        """
        Identifies those features that originally had the same ID as `f`
        (according to the id_spec), but were modified because of duplicate
        IDs.
        """
        candidates = [self._get_feature(f.id)]
        c = self.conn.cursor()
        results = c.execute(
            constants._SELECT + '''
            JOIN duplicates ON
            duplicates.newid = features.id WHERE duplicates.idspecid = ?''',
            (f.id,)
        )
        for i in results:
            candidates.append(
                feature.Feature(dialect=self.iterator.dialect, **i))
        return list(set(candidates))

    def _populate_from_lines(self, lines):
        raise NotImplementedError

    def _update_relations(self):
        raise NotImplementedError

    def _drop_indexes(self):
        c = self.conn.cursor()
        for index in constants.INDEXES:
            c.execute("DROP INDEX IF EXISTS ?", (index,))
        self.conn.commit()

    def set_pragmas(self, pragmas):
        """
        Set pragmas for the current database connection.

        Parameters
        ----------
        pragmas : dict
            Dictionary of pragmas; see constants.default_pragmas for a template
            and http://www.sqlite.org/pragma.html for a full list.
        """
        self.pragmas = pragmas
        c = self.conn.cursor()
        c.executescript(
            ';\n'.join(
                ['PRAGMA %s=%s' % i for i in self.pragmas.items()]
            )
        )
        self.conn.commit()

    def _init_tables(self):
        """
        Table creation
        """
        c = self.conn.cursor()
        v = sqlite3.sqlite_version_info
        self.set_pragmas(self.pragmas)
        c.executescript(constants.SCHEMA)
        self.conn.commit()

    def _finalize(self):
        """
        Various last-minute stuff to perform after file has been parsed and
        imported.

        In general, if you'll be adding stuff to the meta table, do it here.
        """
        c = self.conn.cursor()
        directives = self.directives + self.iterator.directives
        c.executemany('''
                      INSERT INTO directives VALUES (?)
                      ''', ((i,) for i in directives))
        c.execute(
            '''
            INSERT INTO meta (version, dialect)
            VALUES (:version, :dialect)''',
            dict(version=version.version,
                 dialect=helpers._jsonify(self.iterator.dialect))
        )

        c.executemany(
            '''
            INSERT OR REPLACE INTO autoincrements VALUES (?, ?)
            ''', list(self._autoincrements.items()))

        # These indexes are *well* worth the effort and extra storage: over
        # 500x speedup on code like this:
        #
        #   genes = []
        #   for i in db.features_of_type('snoRNA'):
        #       for k in db.parents(i, level=1, featuretype='gene'):
        #           genes.append(k.id)
        #
        logger.info("Creating relations(parent) index")
        c.execute('DROP INDEX IF EXISTS relationsparent')
        c.execute('CREATE INDEX relationsparent ON relations (parent)')
        logger.info("Creating relations(child) index")
        c.execute('DROP INDEX IF EXISTS relationschild')
        c.execute('CREATE INDEX relationschild ON relations (child)')
        logger.info("Creating features(featuretype) index")
        c.execute('DROP INDEX IF EXISTS featuretype')
        c.execute('CREATE INDEX featuretype ON features (featuretype)')
        logger.info("Creating features (seqid, start, end) index")
        c.execute('DROP INDEX IF EXISTS seqidstartend')
        c.execute('CREATE INDEX seqidstartend ON features (seqid, start, end)')
        logger.info("Creating features (seqid, start, end, strand) index")
        c.execute('DROP INDEX IF EXISTS seqidstartendstrand')
        c.execute('CREATE INDEX seqidstartendstrand ON features (seqid, start, end, strand)')

        # speeds computation 1000x in some cases
        logger.info("Running ANALYZE features")
        c.execute('ANALYZE features')

        self.conn.commit()

        self.warnings = self.iterator.warnings

    def create(self):
        """
        Calls various methods sequentially in order to fully build the
        database.
        """
        # Calls each of these methods in order.  _populate_from_lines and
        # _update_relations must be implemented in subclasses.
        self._init_tables()
        self._populate_from_lines(self.iterator)
        self._update_relations()
        self._finalize()

    # TODO: not sure this is used anywhere
    def update(self, iterator):
        self._populate_from_lines(iterator)
        self._update_relations()

    def execute(self, query):
        """
        Execute a query directly on the database.
        """
        c = self.conn.cursor()
        result = c.execute(query)
        for i in result:
            yield i

    def _insert(self, feature, cursor):
        """
        Insert a feature into the database.
        """
        try:
            cursor.execute(constants._INSERT, feature.astuple())
        except sqlite3.ProgrammingError:
            cursor.execute(
                constants._INSERT, feature.astuple(self.default_encoding))

    def _replace(self, feature, cursor):
        """
        Insert a feature into the database.
        """
        try:
            cursor.execute(
                constants._UPDATE,
                list(feature.astuple()) + [feature.id])
        except sqlite3.ProgrammingError:
            cursor.execute(
                constants._INSERT,
                list(feature.astuple(self.default_encoding)) + [feature.id])


class _GFFDBCreator(_DBCreator):
    def __init__(self, *args, **kwargs):
        """
        _DBCreator subclass specifically for working with GFF files.

        create_db() delegates to this class -- see that function for docs
        """
        super(_GFFDBCreator, self).__init__(*args, **kwargs)

    def _populate_from_lines(self, lines):
        c = self.conn.cursor()
        self._drop_indexes()
        last_perc = 0
        logger.info("Populating features")
        msg = ("Populating features table and first-order relations: "
               "%d features\r")

        # c.executemany() was not as much of an improvement as I had expected.
        #
        # Compared to a benchmark of doing each insert separately:
        # executemany using a list of dicts to iterate over is ~15% slower
        # executemany using a list of tuples to iterate over is ~8% faster
        features_seen = None
        _features, _relations = [], []
        for i, f in enumerate(lines):
            features_seen = i

            # Percent complete

            if self.verbose:
                if i % 1000 == 0:
                    sys.stderr.write(msg % i)
                    sys.stderr.flush()

            # TODO: handle ID creation here...should be combined with the
            # INSERT below (that is, don't IGNORE below but catch the error and
            # re-try with a new ID).  However, is this doable with an
            # execute-many?
            f.id = self._id_handler(f)
            try:
                self._insert(f, c)
            except sqlite3.IntegrityError:
                fixed, final_strategy = self._do_merge(f, self.merge_strategy)
                if final_strategy == 'merge':
                    c.execute(
                        '''
                        UPDATE features SET attributes = ?
                        WHERE id = ?
                        ''', (helpers._jsonify(fixed.attributes),
                              fixed.id))

                    # For any additional fields we're merging, update those as
                    # well.
                    if self.force_merge_fields:
                        _set_clause = ', '.join(
                            ['%s = ?' % field
                             for field in self.force_merge_fields])
                        values = [
                            getattr(fixed, field)
                            for field in self.force_merge_fields] + [fixed.id]
                        c.execute(
                            '''
                            UPDATE features SET %s
                            WHERE id = ?
                            ''' % _set_clause, tuple(values))

                elif final_strategy == 'replace':
                    self._replace(f, c)

                elif final_strategy == 'create_unique':
                    self._insert(f, c)

            if 'Parent' in f.attributes:
                for parent in f.attributes['Parent']:
                    c.execute(
                        '''
                        INSERT OR IGNORE INTO relations VALUES
                        (?, ?, 1)
                        ''', (parent, f.id))
        if features_seen is None:
            raise ValueError("No lines parsed -- was an empty file provided?")

        self.conn.commit()
        if self.verbose:
            logger.info(msg % i)

    def _update_relations(self):
        logger.info("Updating relations")
        c = self.conn.cursor()
        c2 = self.conn.cursor()
        c3 = self.conn.cursor()

        # TODO: pre-compute indexes?
        # c.execute('CREATE INDEX ids ON features (id)')
        # c.execute('CREATE INDEX parentindex ON relations (parent)')
        # c.execute('CREATE INDEX childindex ON relations (child)')
        # self.conn.commit()

        if isinstance(self._keep_tempfiles, six.string_types):
            suffix = self._keep_tempfiles
        else:
            suffix = '.gffutils'
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix).name
        with open(tmp, 'w') as fout:

            # Here we look for "grandchildren" -- for each ID, get the child
            # (parenthetical subquery below); then for each of those get *its*
            # child (main query below).
            #
            # Results are written to temp file so that we don't read and write at
            # the same time, which would slow things down considerably.

            c.execute('SELECT id FROM features')
            for parent in c:
                c2.execute('''
                           SELECT child FROM relations WHERE parent IN
                           (SELECT child FROM relations WHERE parent = ?)
                           ''', tuple(parent))
                for grandchild in c2:
                    fout.write('\t'.join((parent[0], grandchild[0])) + '\n')

        def relations_generator():
            with open(fout.name) as fin:
                for line in fin:
                    parent, child = line.strip().split('\t')
                    yield dict(parent=parent, child=child, level=2)

        c.executemany(
            '''
            INSERT OR IGNORE INTO relations VALUES
            (:parent, :child, :level)
            ''', relations_generator())

        # TODO: Index creation.  Which ones affect performance?
        c.execute("DROP INDEX IF EXISTS binindex")
        c.execute("CREATE INDEX binindex ON features (bin)")

        self.conn.commit()

        if not self._keep_tempfiles:
            os.unlink(fout.name)


class _GTFDBCreator(_DBCreator):
    def __init__(self, *args, **kwargs):
        """
        create_db() delegates to this class -- see that function for docs
        """
        self.transcript_key = kwargs.pop('transcript_key', 'transcript_id')
        self.gene_key = kwargs.pop('gene_key', 'gene_id')
        self.subfeature = kwargs.pop('subfeature', 'exon')
        super(_GTFDBCreator, self).__init__(*args, **kwargs)

    def _populate_from_lines(self, lines):
        msg = (
            "Populating features table and first-order relations: %d "
            "features\r"
        )

        c = self.conn.cursor()

        # Only check this many features to see if it's a gene or transcript and
        # issue the appropriate warning.
        gene_and_transcript_check_limit = 1000

        last_perc = 0
        lines_seen = 0
        for i, f in enumerate(lines):

            # See issues #48 and #20.
            if lines_seen < gene_and_transcript_check_limit:
                if (
                    f.featuretype == 'transcript' and
                    not self.disable_infer_transcripts
                ):
                    warnings.warn(
                        "It appears you have a transcript feature in your GTF "
                        "file. You may want to use the "
                        "`disable_infer_transcripts=True` "
                        "option to speed up database creation")
                elif (
                    f.featuretype == 'gene' and
                    not self.disable_infer_genes
                ):
                    warnings.warn(
                        "It appears you have a gene feature in your GTF "
                        "file. You may want to use the "
                        "`disable_infer_genes=True` "
                        "option to speed up database creation")

            lines_seen = i + 1

            # Percent complete
            if self.verbose:

                if i % 1000 == 0:
                    sys.stderr.write(msg % i)
                    sys.stderr.flush()

            f.id = self._id_handler(f)

            # Insert the feature itself...
            try:
                self._insert(f, c)
            except sqlite3.IntegrityError:
                fixed, final_strategy = self._do_merge(f, self.merge_strategy)
                if final_strategy == 'merge':
                    c.execute(
                        '''
                        UPDATE features SET attributes = ?
                        WHERE id = ?
                        ''', (helpers._jsonify(fixed.attributes),
                              fixed.id))
                    # For any additional fields we're merging, update those as
                    # well.
                    if self.force_merge_fields:
                        _set_clause = ', '.join(
                            ['%s = ?' % field
                             for field in self.force_merge_fields])
                        values = [getattr(fixed, field)
                                  for field in self.force_merge_fields]\
                            + [fixed.id]
                        c.execute(
                            '''
                            UPDATE features SET %s
                            WHERE id = ?
                            ''' % _set_clause, values)
                elif final_strategy == 'replace':
                    self._replace(f, c)
                elif final_strategy == 'create_unique':
                    self._insert(f, c)

            # For an on-spec GTF file,
            # self.transcript_key = "transcript_id"
            # self.gene_key = "gene_id"
            relations = []
            parent = None
            grandparent = None
            if self.transcript_key in f.attributes:
                parent = f.attributes[self.transcript_key][0]
                relations.append((parent, f.id, 1))

            if self.gene_key in f.attributes:
                grandparent = f.attributes[self.gene_key]
                if len(grandparent) > 0:
                    grandparent = grandparent[0]
                    relations.append((grandparent, f.id, 2))
                    if parent is not None:
                        relations.append((grandparent, parent, 1))

            # Note the IGNORE, so relationships defined many times in the file
            # (e.g., the transcript-gene relation on pretty much every line in
            # a GTF) will only be included once.
            c.executemany(
                '''
                INSERT OR IGNORE INTO relations (parent, child, level)
                VALUES (?, ?, ?)
                ''', relations
            )

        if lines_seen == 0:
            raise ValueError("No lines parsed -- was an empty file provided?")
        logger.info('Committing changes')
        self.conn.commit()
        if self.verbose:
            logger.info(msg % i)

    def _update_relations(self):

        if self.disable_infer_genes and self.disable_infer_transcripts:
            return

        # TODO: do any indexes speed this up?
        c = self.conn.cursor()
        c2 = self.conn.cursor()

        logger.info("Creating relations(parent) index")
        c.execute('DROP INDEX IF EXISTS relationsparent')
        c.execute('CREATE INDEX relationsparent ON relations (parent)')
        logger.info("Creating relations(child) index")
        c.execute('DROP INDEX IF EXISTS relationschild')
        c.execute('CREATE INDEX relationschild ON relations (child)')

        if not (self.disable_infer_genes or self.disable_infer_transcripts):
            msg = 'gene and transcript'
        elif self.disable_infer_transcripts:
            msg = 'gene'
        elif self.disable_infer_genes:
            msg = 'transcript'
        logger.info('Inferring %s extents '
                    'and writing to tempfile' % msg)

        if isinstance(self._keep_tempfiles, six.string_types):
            suffix = self._keep_tempfiles
        else:
            suffix = '.gffutils'

        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix).name
        with open(tmp, 'w') as fout:
            self._tmpfile = tmp

            # This takes some explanation...
            #
            # First, the nested subquery gets the level-1 parents of
            # self.subfeature featuretypes.  For an on-spec GTF file,
            # self.subfeature = "exon". So this subquery translates to getting the
            # distinct level-1 parents of exons -- which are transcripts.
            #
            # OK, so this first subquery is now a list of transcripts; call it
            # "firstlevel".
            #
            # Then join firstlevel on relations, but the trick is to now consider
            # each transcript a *child* -- so that relations.parent (on the first
            # line of the query) will be the first-level parent of the transcript
            # (the gene).
            #
            #
            # The result is something like:
            #
            #   transcript1     gene1
            #   transcript2     gene1
            #   transcript3     gene2
            #
            # Note that genes are repeated; below we need to ensure that only one
            # is added.  To ensure this, the results are ordered by the gene ID.
            #
            # By the way, we do this even if we're only looking for transcripts or
            # only looking for genes.

            c.execute(
                '''
                SELECT DISTINCT firstlevel.parent, relations.parent
                FROM (
                    SELECT DISTINCT parent
                    FROM relations
                    JOIN features ON features.id = relations.child
                    WHERE features.featuretype = ?
                    AND relations.level = 1
                )
                AS firstlevel
                JOIN relations ON firstlevel.parent = child
                WHERE relations.level = 1
                ORDER BY relations.parent
                ''', (self.subfeature,))

            # Now we iterate through those results (using a new cursor) to infer
            # the extent of transcripts and/or genes.

            last_gene_id = None
            n_features = 0
            for transcript_id, gene_id in c:

                if not self.disable_infer_transcripts:
                    # transcript extent
                    c2.execute(
                        '''
                        SELECT MIN(start), MAX(end), strand, seqid
                        FROM features
                        JOIN relations ON
                        features.id = relations.child
                        WHERE parent = ? AND featuretype == ?
                        ''', (transcript_id, self.subfeature))
                    transcript_start, transcript_end, strand, seqid = c2.fetchone()
                    transcript_attributes = {
                        self.transcript_key: [transcript_id],
                        self.gene_key: [gene_id]
                    }
                    transcript_bin = bins.bins(
                        transcript_start, transcript_end, one=True)

                    # Write out to file; we'll be reading it back in shortly.  Omit
                    # score, frame, source, and extra since they will always have
                    # the same default values (".", ".", "gffutils_derived", and []
                    # respectively)

                    fout.write('\t'.join(map(str, [
                        transcript_id,
                        seqid,
                        transcript_start,
                        transcript_end,
                        strand,
                        'transcript',
                        transcript_bin,
                        helpers._jsonify(transcript_attributes)
                    ])) + '\n')

                    n_features += 1

                if not self.disable_infer_genes:
                    # Infer gene extent, but only if we haven't done so already
                    if gene_id != last_gene_id:
                        c2.execute(
                            '''
                            SELECT MIN(start), MAX(end), strand, seqid
                            FROM features
                            JOIN relations ON
                            features.id = relations.child
                            WHERE parent = ? AND featuretype == ?
                            ''', (gene_id, self.subfeature))
                        gene_start, gene_end, strand, seqid = c2.fetchone()
                        gene_attributes = {self.gene_key: [gene_id]}
                        gene_bin = bins.bins(gene_start, gene_end, one=True)

                        fout.write('\t'.join(map(str, [
                            gene_id,
                            seqid,
                            gene_start,
                            gene_end,
                            strand,
                            'gene',
                            gene_bin,
                            helpers._jsonify(gene_attributes)
                        ])) + '\n')

                    last_gene_id = gene_id
                    n_features += 1

        def derived_feature_generator():
            """
            Generator of items from the file that was just created...
            """
            keys = ['parent', 'seqid', 'start', 'end', 'strand',
                    'featuretype', 'bin', 'attributes']
            with open(fout.name) as fin:
                for line in fin:
                    d = dict(list(zip(keys, line.strip().split('\t'))))
                    d.pop('parent')
                    d['score'] = '.'
                    d['source'] = 'gffutils_derived'
                    d['frame'] = '.'
                    d['extra'] = []
                    d['attributes'] = helpers._unjsonify(d['attributes'])
                    f = feature.Feature(**d)
                    f.id = self._id_handler(f)
                    yield f

        # Drop the indexes so the inserts are faster
        c.execute('DROP INDEX IF EXISTS relationsparent')
        c.execute('DROP INDEX IF EXISTS relationschild')

        # Insert the just-inferred transcripts and genes.  TODO: should we
        # *always* use "merge" here for the merge_strategy?
        logger.info("Importing inferred features into db")
        last_perc = None
        for i, f in enumerate(derived_feature_generator()):
            perc = int(i / float(n_features) * 100)
            if perc != last_perc:
                sys.stderr.write('%s of %s (%s%%)\r' % (i, n_features, perc))
                sys.stderr.flush()
            last_perc = perc
            try:
                self._insert(f, c)
            except sqlite3.IntegrityError:
                fixed, final_strategy = self._do_merge(f, 'merge')
                c.execute(
                    '''
                    UPDATE features SET attributes = ?
                    WHERE id = ?
                    ''', (helpers._jsonify(fixed.attributes),
                          fixed.id))

        logger.info("Committing changes")
        self.conn.commit()
        if not self._keep_tempfiles:
            os.unlink(fout.name)

        # TODO: recreate indexes? Typically the _finalize() method will be
        # called after this one, and indexes are created in _finalize().


def create_db(data, dbfn, id_spec=None, force=False, verbose=False,
              checklines=10, merge_strategy='error', transform=None,
              gtf_transcript_key='transcript_id', gtf_gene_key='gene_id',
              gtf_subfeature='exon', force_gff=False,
              force_dialect_check=False, from_string=False, keep_order=False,
              text_factory=sqlite3.OptimizedUnicode, force_merge_fields=None,
              pragmas=constants.default_pragmas, sort_attribute_values=False,
              dialect=None, _keep_tempfiles=False, infer_gene_extent=True,
              disable_infer_genes=False, disable_infer_transcripts=False,
              **kwargs):
    """
    Create a database from a GFF or GTF file.

    For more details on when and how to use the kwargs below, see the examples
    in the online documentation (:ref:`examples`).

    Parameters
    ----------
    data : string or iterable

        If a string (and `from_string` is False), then `data` is the path to
        the original GFF or GTF file.

        If a string and `from_string` is True, then assume `data` is the actual
        data to use.

        Otherwise, it's an iterable of Feature objects.

    dbfn : string

        Path to the database that will be created.  Can be the special string
        ":memory:" to create an in-memory database.

    id_spec : string, list, dict, callable, or None

        This parameter guides what will be used as the primary key for the
        database, which in turn determines how you will access individual
        features by name from the database.

        If `id_spec=None`, then auto-increment primary keys based on the
        feature type (e.g., "gene_1", "gene_2").  This is also the fallback
        behavior for the other values below.

        If `id_spec` is a string, then look for this key in the attributes.  If
        it exists, then use its value as the primary key, otherwise
        autoincrement based on the feature type.  For many GFF3 files, "ID"
        usually works well.

        If `id_spec` is a list or tuple of keys, then check for each one in
        order, using the first one found.  For GFF3, this might be ["ID",
        "Name"], which would use the ID if it exists, otherwise the Name,
        otherwise autoincrement based on the feature type.

        If `id_spec` is a dictionary, then it is a mapping of feature types to
        what should be used as the ID.  For example, for GTF files, `{'gene':
        'gene_id', 'transcript': 'transcript_id'}` may be useful.  The values
        of this dictionary can also be a list, e.g., `{'gene': ['gene_id',
        'geneID']}`

        If `id_spec` is a callable object, then it accepts a dictionary from
        the iterator and returns one of the following:

            * None (in which case the feature type will be auto-incremented)
            * string (which will be used as the primary key)
            * special string starting with "autoincrement:X", where "X" is
              a string that will be used for auto-incrementing.  For example,
              if "autoincrement:chr10", then the first feature will be
              "chr10_1", the second "chr10_2", and so on.

    force : bool

        If `False` (default), then raise an exception if `dbfn` already exists.
        Use `force=True` to overwrite any existing databases.

    verbose : bool

        Report percent complete and other feedback on how the db creation is
        progressing.

        In order to report percent complete, the entire file needs to be read
        once to see how many items there are; for large files you may want to
        use `verbose=False` to avoid this.

    checklines : int

        Number of lines to check the dialect.

    merge_strategy : str
        One of {merge, create_unique, error, warning, replace}.

        This parameter specifies the behavior when two items have an identical
        primary key.

        Using `merge_strategy="merge"`, then there will be a single entry in
        the database, but the attributes of all features with the same primary
        key will be merged.

        Using `merge_strategy="create_unique"`, then the first entry will use
        the original primary key, but the second entry will have a unique,
        autoincremented primary key assigned to it

        Using `merge_strategy="error"`, a :class:`gffutils.DuplicateID`
        exception will be raised.  This means you will have to edit the file
        yourself to fix the duplicated IDs.

        Using `merge_strategy="warning"`, a warning will be printed to the
        logger, and the duplicate feature will be skipped.

        Using `merge_strategy="replace"` will replace the entire existing
        feature with the new feature.

    transform : callable

        Function (or other callable object) that accepts a `Feature` object and
        returns a (possibly modified) `Feature` object.

    gtf_transcript_key, gtf_gene_key : string

        Which attribute to use as the transcript ID and gene ID respectively
        for GTF files.  Default is `transcript_id` and `gene_id` according to
        the GTF spec.

    gtf_subfeature : string

        Feature type to use as a "gene component" when inferring gene and
        transcript extents for GTF files.  Default is `exon` according to the
        GTF spec.

    force_gff : bool
        If True, do not do automatic format detection -- only use GFF.

    force_dialect_check : bool
        If True, the dialect will be checkef for every feature (instead of just
        `checklines` features).  This can be slow, but may be necessary for
        inconsistently-formatted input files.

    from_string : bool
        If True, then treat `data` as actual data (rather than the path to
        a file).

    keep_order : bool

        If True, all features returned from this instance will have the
        order of their attributes maintained.  This can be turned on or off
        database-wide by setting the `keep_order` attribute or with this
        kwarg, or on a feature-by-feature basis by setting the `keep_order`
        attribute of an individual feature.

        Note that a single order of attributes will be used for all features.
        Specifically, the order will be determined by the order of attribute
        keys in the first `checklines` of the input data. See
        helpers._choose_dialect for more information on this.

        Default is False, since this includes a sorting step that can get
        time-consuming for many features.

    infer_gene_extent : bool
        DEPRECATED in version 0.8.4. See `disable_infer_transcripts` and
        `disable_infer_genes` for more granular control.

    disable_infer_transcripts, disable_infer_genes : bool
        Only used for GTF files. By default -- and according to the GTF spec --
        we assume that there are no transcript or gene features in the file.
        gffutils then infers the extent of each transcript based on its
        constituent exons and infers the extent of each gene bases on its
        constituent transcripts.

        This default behavior is problematic if the input file already contains
        transcript or gene features (like recent GENCODE GTF files for human),
        since 1) the work to infer extents is unnecessary, and 2)
        trying to insert an inferred feature back into the database triggers
        gffutils' feature-merging routines, which can get time consuming.

        The solution is to use `disable_infer_transcripts=True` if your GTF
        already has transcripts in it, and/or `disable_infer_genes=True` if it
        already has genes in it. This can result in dramatic (100x) speedup.

        Prior to version 0.8.4, setting `infer_gene_extents=False` would
        disable both transcript and gene inference simultaneously. As of
        version 0.8.4, these argument allow more granular control.

    force_merge_fields : list
        If merge_strategy="merge", then features will only be merged if their
        non-attribute values are identical (same chrom, source, start, stop,
        score, strand, phase).  Using `force_merge_fields`, you can override
        this behavior to allow merges even when fields are different.  This
        list can contain one or more of ['seqid', 'source', 'featuretype',
        'score', 'strand', 'frame'].  The resulting merged fields will be
        strings of comma-separated values.  Note that 'start' and 'end' are not
        available, since these fields need to be integers.

    text_factory : callable
        Text factory to use for the sqlite3 database.  See
        https://docs.python.org/2/library/\
                sqlite3.html#sqlite3.Connection.text_factory
        for details. The default sqlite3.OptimizedUnicode will return Unicode
        objects only for non-ASCII data, and bytestrings otherwise.

    pragmas : dict
        Dictionary of pragmas used when creating the sqlite3 database. See
        http://www.sqlite.org/pragma.html for a list of available pragmas.  The
        defaults are stored in constants.default_pragmas, which can be used as
        a template for supplying a custom dictionary.

    sort_attribute_values : bool
        All features returned from the database will have their attribute
        values sorted.  Typically this is only useful for testing, since this
        can get time-consuming for large numbers of features.

    _keep_tempfiles : bool or string
        False by default to clean up intermediate tempfiles created during GTF
        import.  If True, then keep these tempfile for testing or debugging.
        If string, then keep the tempfile for testing, but also use the string
        as the suffix fo the tempfile. This can be useful for testing in
        parallel environments.

    Returns
    -------
    New :class:`FeatureDB` object.
    """
    _locals = locals()

    # Check if any older kwargs made it in
    deprecation_handler(kwargs)

    kwargs = dict((i, _locals[i]) for i in constants._iterator_kwargs)

    # First construct an iterator so that we can identify the file format.
    # DataIterator figures out what kind of data was provided (string of lines,
    # filename, or iterable of Features) and checks `checklines` lines to
    # identify the dialect.
    iterator = iterators.DataIterator(**kwargs)

    kwargs.update(**_locals)

    if dialect is None:
        dialect = iterator.dialect

    # However, a side-effect of this is that  if `data` was a generator, then
    # we've just consumed `checklines` items (see
    # iterators.BaseIterator.__init__, which calls iterators.peek).
    #
    # But it also chains those consumed items back onto the beginning, and the
    # result is available as as iterator._iter.
    #
    # That's what we should be using now for `data:
    kwargs['data'] = iterator._iter
    kwargs['directives'] = iterator.directives

    # Since we've already checked lines, we don't want to do it again
    kwargs['checklines'] = 0

    if force_gff or (dialect['fmt'] == 'gff3'):
        cls = _GFFDBCreator
        id_spec = id_spec or 'ID'
        add_kwargs = dict(
            id_spec=id_spec,
        )

    elif dialect['fmt'] == 'gtf':
        cls = _GTFDBCreator
        id_spec = id_spec or {'gene': 'gene_id', 'transcript': 'transcript_id'}
        add_kwargs = dict(
            transcript_key=gtf_transcript_key,
            gene_key=gtf_gene_key,
            subfeature=gtf_subfeature,
            id_spec=id_spec,
        )

    kwargs.update(**add_kwargs)
    kwargs['dialect'] = dialect
    c = cls(**kwargs)

    c.create()
    if dbfn == ':memory:':
        db = interface.FeatureDB(c.conn,
                                 keep_order=keep_order,
                                 pragmas=pragmas,
                                 sort_attribute_values=sort_attribute_values,
                                 text_factory=text_factory)
    else:
        db = interface.FeatureDB(c,
                                 keep_order=keep_order,
                                 pragmas=pragmas,
                                 sort_attribute_values=sort_attribute_values,
                                 text_factory=text_factory)

    return db
