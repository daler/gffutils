import collections
import tempfile
import sys
import os
import sqlite3
import constants
import version
import parser
import bins
import helpers
import feature

import logging

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)


class _DBCreator(object):
    def __init__(self, fn, dbfn, force=False, verbose=True, id_spec=None,
                 merge_strategy='merge', checklines=10, transform=None):
        """
        Base class for _GFFDBCreator and _GTFDBCreator; see create_db()
        function for docs
        """
        self.merge_strategy = merge_strategy

        self._autoincrements = collections.defaultdict(int)
        if force:
            if os.path.exists(dbfn):
                os.unlink(dbfn)
        self.dbfn = dbfn
        self.id_spec = id_spec
        conn = sqlite3.connect(dbfn)
        self.conn = conn
        #self.conn.text_factory = sqlite3.OptimizedUnicode
        self.conn.text_factory = str
        self.fn = fn

        self.parser = parser.Parser(
            fn, checklines=checklines, transform=transform)

        self.verbose = verbose
        self._orig_logger_level = logger.level
        if not self.verbose:
            logger.setLevel(logging.ERROR)

    def _increment_featuretype_autoid(self, key):
        self._autoincrements[key] += 1
        return '%s_%s' % (key, self._autoincrements[key])

    def _id_handler(self, d):
        """
        Given a dictionary from self.parser, figure out what the ID should be.

        `_autoincrement_key` is which field to use that will be
        auto-incremented.  Typically this will be "feature" (for exon_1,
        exon_2, etc), but another useful one is "id", which is is used for
        duplicate IDs.
        """

        # If id_spec is a string, convert to iterable for later
        if isinstance(self.id_spec, basestring):
            id_key = [self.id_spec]

        # If dict, then assume it's a feature -> attribute mapping, e.g.,
        # {'gene': 'gene_id'} for GTF
        elif isinstance(self.id_spec, dict):
            try:
                id_key = self.id_spec[d['featuretype']]
                if isinstance(id_key, basestring):
                    id_key = [id_key]

            # Otherwise, use default auto-increment.
            except KeyError:
                return self._increment_featuretype_autoid(d['featuretype'])

        # Otherwise assume it's an iterable.
        else:
            id_key = self.id_spec

        # Then try them in order, returning the first one that works:
        for k in id_key:

            if hasattr(k, '__call__'):
                _id = self.id_spec(d)
                if _id:
                    if _id.startswith('autoincrement:'):
                        return self._increment_featuretype_autoid(_id[14:])
                    return _id

            # use GFF fields rather than attributes for cases like :seqid: or
            # :strand:
            if (len(k) > 3) and (k[0] == ':') and (k[-1] == ':'):
                # No [0] here -- only attributes key/vals are forced into
                # lists, not standard GFF fields.
                return d[k[1:-1]]
            else:
                v = d['attributes'][k]
            if len(v) == 0:
                del d['attributes'][k]
            else:
                return v[0]

        # If we get here, then default autoincrement
        return self._increment_featuretype_autoid(d['featuretype'])

    def _get_feature(self, ID):
        c = self.conn.cursor()
        return dict(
            zip(
                constants._keys,
                c.execute(constants._SELECT + 'WHERE id = ?', (ID,)).fetchone()
            )
        )

    def _do_merge(self, d):
        """
        Different merge strategies upon name conflicts.

        "error":
            raise error

        "warning"
            show warning

        "merge":
            combine old and new attributes -- but only if everything else
            matches; otherwise error.  This can be slow, but is thorough.

        "create_unique":
            Autoincrement based on the ID
        """
        # Note that `d` comes in un-JSONified

        if self.merge_strategy == 'error':
            raise ValueError("Duplicate ID {id}".format(**d))
        if self.merge_strategy == 'warning':
            logger.warning(
                "Duplicate lines in file for id '{id}'; "
                "ignoring all but the first".format(**d))
            return None
        elif self.merge_strategy == 'merge':
            # retrieve the existing row
            existing_d = self._get_feature(d['id'])

            # does everything besides attributes and extra match?
            for k in constants._gffkeys[:-1]:
                # Note str() here: `existing_d` came from the db (so start/end
                # are integers) but `d` came from the file, so they are still
                # strings.
                assert str(existing_d[k]) == d[k], (
                    "Same ID, but differing info for %s field. "
                    "Line %s: \n%s" % (
                        d['id'],
                        self.parser.current_line_number,
                        self.parser.current_line))

            attributes = helpers._unjsonify(existing_d['attributes'])

            # update the attributes (using sets for de-duping)
            for k, v in d['attributes'].items():
                attributes[k] = list(set(attributes[k]).union(v))
            existing_d['attributes'] = attributes
            return existing_d
        elif self.merge_strategy == 'create_unique':
            d['id'] = self._increment_featuretype_autoid(d['id'])
            return d

    def _populate_from_lines(self, lines):
        raise NotImplementedError

    def _update_relations(self):
        raise NotImplementedError

    def _drop_indexes(self):
        c = self.conn.cursor()
        for index in constants.INDEXES:
            c.execute("DROP INDEX IF EXISTS ?", (index,))
        self.conn.commit()

    def _init_tables(self):
        """
        Table creation
        """
        c = self.conn.cursor()
        c.executescript(constants.SCHEMA)
        self.conn.commit()

    def _finalize(self):
        """
        Various last-minute stuff to perform after file has been parsed and
        imported.

        In general, if you'll be adding stuff to the meta table, do it here.
        """
        c = self.conn.cursor()

        c.executemany('''
                      INSERT INTO directives VALUES (?)
                      ''', ((i,) for i in self.parser.directives))
        c.execute(
            '''
            INSERT INTO meta (version, dialect)
            VALUES (:version, :dialect)''',
            dict(version=version.version,
                 dialect=helpers._jsonify(self.parser.dialect))
        )

        c.executemany(
            '''
            INSERT OR REPLACE INTO autoincrements VALUES (?, ?)
            ''', self._autoincrements.items())

        self.conn.commit()

        self.warnings = self.parser.warnings

    def create(self):
        """
        Calls various methods sequentially in order to fully build the
        database.
        """
        # Calls each of these methods in order.  _populate_from_lines and
        # _update_relations must be implemented in subclasses.
        self._init_tables()
        self._populate_from_lines(self.parser)
        self._update_relations()
        self._finalize()

        # reset logger to whatever it was before...
        logger.setLevel(self._orig_logger_level)


class _GFFDBCreator(_DBCreator):
    def __init__(self, *args, **kwargs):
        """
        _DBCreator subclass specifically for working with GFF files.

        create_db() delegates to this class -- see that function for docs
        """
        super(_GFFDBCreator, self).__init__(*args, **kwargs)

    def _populate_from_lines(self, lines):
        c = self.conn.cursor()
        c.execute(
            '''
            PRAGMA synchronous=NORMAL;
            ''')
        c.execute(
            '''
            PRAGMA journal_mode=WAL;
            ''')
        self._drop_indexes()
        last_perc = 0
        logger.info("Populating features")
        msg = ("Populating features table and first-order relations: "
               "%d (%d%%)\r")

        if self.verbose:
            # This is so we can provide percent-completed estimates below --
            # but it'll take time to read the file to determine how many valid
            # lines...
            self._nfeatures = float(self.parser._valid_line_count())

        # ONEBYONE is used for profiling -- how to get faster inserts?
        # ONEBYONE=False will do a single executemany
        # ONEBYONE=True will do many single execute
        #
        # c.executemany() was not as much of an improvement as I had expected.
        #
        # Compared to a benchmark of doing each insert separately:
        # executemany using a list of dicts to iterate over is ~15% slower
        # executemany using a list of tuples to iterate over is ~8% faster
        ONEBYONE = True

        _features, _relations = [], []
        for i, d in enumerate(lines):

            # Percent complete
            if self.verbose:
                perc = int(i / self._nfeatures * 100)
                if perc != last_perc:
                    sys.stderr.write(msg % (i, perc))
                    sys.stderr.flush()
                last_perc = perc

            # TODO: handle ID creation here...should be combined with the
            # INSERT below (that is, don't IGNORE below but catch the error and
            # re-try with a new ID).  However, is this doable with an
            # execute-many?
            d['id'] = self._id_handler(d)
            d['bin'] = helpers._bin_from_dict(d)

            fields = helpers._dict_to_fields(d)

            if ONEBYONE:

                # TODO: these are two, one-by-one execute statements.
                # Profiling shows that this is a slow step. Need to use
                # executemany, which probably means writing to file first.

                try:
                    c.execute(constants._INSERT, helpers._dict_to_fields(d))
                except sqlite3.IntegrityError:
                    fixed = self._do_merge(d)
                    if self.merge_strategy == 'merge':
                        c.execute(
                            '''
                            UPDATE features SET attributes = ?
                            WHERE id = ?
                            ''', (helpers._jsonify(fixed['attributes']),
                                  fixed['id']))

                    if self.merge_strategy == 'create_unique':
                        c.execute(constants._INSERT,
                                  helpers._dict_to_fields(fixed))

                # Works in all cases since attributes is a defaultdict
                for parent in d['attributes']['Parent']:
                    c.execute(
                        '''
                        INSERT OR IGNORE INTO relations VALUES
                        (?, ?, 1)
                        ''', (parent, d['id']))

            else:
                _features.append(helpers._dict_to_fields(d))

                # Works in all cases since attributes is a defaultdict
                for parent in d['attributes']['Parent']:
                    _relations.append((parent, d['id']))

        if not ONEBYONE:
            # Profiling shows that there's an extra overhead for using dict
            # syntax in sqlite3 queries.  Even though we do the lookup above
            # (when appending to _features), it's still faster to use the tuple
            # syntax.
            c.executemany(constants._INSERT, _features)

            c.executemany(
                '''
                INSERT INTO relations VALUES (?,?, 1);
                ''', _relations)

            del _relations
            del _features

        self.conn.commit()
        if self.verbose:
            sys.stderr.write('\n')

    def _update_relations(self):
        logger.info("Updating relations")
        c = self.conn.cursor()
        c2 = self.conn.cursor()
        c3 = self.conn.cursor()

        # TODO: pre-compute indexes?
        #c.execute('CREATE INDEX ids ON features (id)')
        #c.execute('CREATE INDEX parentindex ON relations (parent)')
        #c.execute('CREATE INDEX childindex ON relations (child)')
        #self.conn.commit()

        tmp = tempfile.NamedTemporaryFile(delete=False).name
        fout = open(tmp, 'w')

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
                       ''', parent)
            for grandchild in c2:
                fout.write('\t'.join((parent[0], grandchild[0])) + '\n')
        fout.close()

        def relations_generator():
            for line in open(fout.name):
                parent, child = line.strip().split('\t')
                yield dict(parent=parent, child=child, level=2)

        c.executemany(
            '''
            INSERT OR IGNORE INTO relations VALUES
            (:parent, :child, :level)
            ''', relations_generator())

        # TODO: Index creation.  Which ones affect performance?
        c.execute("CREATE INDEX binindex ON features (bin)")

        self.conn.commit()
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
        msg = ("Populating features table and first-order relations: "
               "%d (%d%%)\r")

        if self.verbose:
            # So we can provide percent-completed estimates below -- but it'll
            # take time to read the file to determine how many valid lines...
            self._nfeatures = float(self.parser._valid_line_count())

        c = self.conn.cursor()

        last_perc = 0
        for i, d in enumerate(lines):

            # Percent complete
            if self.verbose:
                perc = int(i / self._nfeatures * 100)
                if perc != last_perc:
                    sys.stderr.write(msg % (i, perc))
                    sys.stderr.flush()
                last_perc = perc

            d['id'] = self._id_handler(d)
            d['bin'] = helpers._bin_from_dict(d)

            # Insert the feature itself...
            c.execute(constants._INSERT, helpers._dict_to_fields(d))

            # TODO: This assumes an on-spec GTF file...may need a config dict
            # (or dialect?) to modify.

            relations = []
            parent = None
            grandparent = None
            if self.transcript_key in d['attributes']:
                parent = d['attributes'][self.transcript_key][0]
                relations.append((parent, d['id'], 1))

            if self.gene_key in d['attributes']:
                grandparent = d['attributes'][self.gene_key]
                if len(grandparent) > 0:
                    grandparent = grandparent[0]
                    relations.append((grandparent, d['id'], 2))
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

        self.conn.commit()

    def _update_relations(self):
        # TODO: do any indexes speed this up?

        c = self.conn.cursor()
        c2 = self.conn.cursor()

        tmp = tempfile.NamedTemporaryFile(delete=False).name
        fout = open(tmp, 'w')
        c.execute('''SELECT DISTINCT parent FROM relations''')
        for parent, in c:
            c2.execute(
                '''
                SELECT MIN(start), MAX(end), level, strand, seqid
                FROM features
                JOIN relations ON
                features.id = relations.child
                WHERE parent = ? AND featuretype == ?
                ''', (parent, self.subfeature))
            start, end, level, strand, seqid = c2.fetchone()

            # if level == 2, then it's a gene
            # if level == 1, then it's a transcript

            if level == 1:
                attributes = {self.transcript_key: [parent]}
                feature = 'transcript'

            if level == 2:
                attributes = {self.gene_key: [parent]}
                feature = 'gene'

            # We need to recalculate a bin for entire parent feature
            _bin = bins.bins(start, end, one=True)

            # Write out to file; we'll be reading it back in shortly.  Omit
            # score, frame, source, and extra since they will have the same
            # default values for all all genes.
            fout.write('\t'.join(map(str, [
                parent, seqid, start, end, strand, feature, _bin,
                helpers._jsonify(attributes)])) + '\n')
        fout.close()

        def derived_feature_generator():
            """
            Generator of items from the file that was just created...
            """
            keys = ['parent', 'seqid', 'start', 'end', 'strand',
                    'featuretype', 'bin', 'attributes']
            for line in open(fout.name):
                d = dict(zip(keys, line.strip().split('\t')))
                d['score'] = '.'
                d['source'] = 'gffutils_derived'
                d['frame'] = '.'
                d['extra'] = []
                d['attributes'] = helpers._unjsonify(d['attributes'])
                d['id'] = self._id_handler(d)
                yield helpers._dict_to_fields(d)

        # Batch-insert the derived features
        c.executemany(constants._INSERT, derived_feature_generator())

        self.conn.commit()
        os.unlink(fout.name)

        # TODO: recreate indexes?

    def execute(self, query):
        """
        Execute a query directly on the database.
        """
        c = self.conn.cursor()
        c.execute(query)
        for i in cursor:
            yield i


class FeatureDB(object):
    def __init__(self, dbfn):
        """
        Connect to a database created by :func:`gffutils.create_db`.

        Parameters
        ----------

        `dbfn` : str

            Path to a database created by :func:`gffutils.create_db`.


        .. note::

            `dbfn` can also be a subclass of :class:`_DBCreator`, useful for
            when :func:`gffutils.create_db` is provided the ``dbfn=":memory:"``
            kwarg.

        """
        # Since specifying the string ":memory:" will actually try to connect
        # to a new, separate (and empty) db in memory, we can alternatively
        # pass in a _DBCreator instance to use its db.
        if isinstance(dbfn, _DBCreator):
            if dbfn.dbfn == ':memory:':
                self.conn = dbfn.conn
            else:
                self.dbfn = dbfn.dbfn
                self.conn = sqlite3.connect(self.dbfn)
        else:
            self.dbfn = dbfn
            self.conn = sqlite3.connect(self.dbfn)

        self.conn.text_factory = str

        c = self.conn.cursor()

        # Load some meta info
        c.execute(
            '''
            SELECT version, dialect FROM meta
            ''')
        version, dialect = c.fetchone()
        self.version = version
        self.dialect = helpers._unjsonify(dialect)

        # Load directives from db
        c.execute(
            '''
            SELECT directive FROM directives
            ''')
        self.directives = [directive[0] for directive in c if directive]

        # Load autoincrements so that when we add new features, we can start
        # autoincrementing from here (instead of from 1, which would cause name
        # collisions)
        c.execute(
            '''
            SELECT base, n FROM autoincrements
            ''')
        self._autoincrements = dict(c)

    def schema(self):
        """
        Returns the database schema as a string.
        """
        c = self.conn.cursor()
        c.execute(
            '''
            SELECT sql FROM sqlite_master
            ''')
        results = []
        for i, in c:
            if i is not None:
                results.append(i)
        return '\n'.join(results)

    def __getitem__(self, key):
        c = self.conn.cursor()
        c.execute(constants._SELECT + ' WHERE id = ?', (key,))
        results = c.fetchall()
        if len(results) != 1:
            raise helpers.FeatureNotFoundError(key)
        return feature.Feature(*results, dialect=self.dialect)

    def count_features_of_type(self, featuretype):
        """
        Simple count of features.

        Can be faster than "grep", and is faster than checking the length of
        results from :meth:`gffutils.FeatureDB.features_of_type`.

        Parameters
        ----------

        `featuretype` : string

            Feature type (e.g., "gene") to count.

        Returns
        -------
        The number of features of this type, as an integer

        """
        c = self.conn.cursor()
        c.execute(
            '''
            SELECT count() FROM features
            WHERE feature = ?
            ''', (featuretype,))
        results = c.fetchone()
        if results is not None:
            results = results[0]
        return results

    def features_of_type(self, featuretype, seqid=None, start=None, end=None,
                         strand=None):
        """
        Returns an iterator of :class:`gffutils.Feature` objects.

        Parameters
        ----------
        `featuretype` : string

            Feature type to return

        `seqid` : string

            If not `None`, restrict results to this seqid

        `start` : int

            If not `None`, restrict results to only >= `start`

        `end` : int

            If not `None`, restrict results to only <= `end`

        `strand` : { "-", "+", "." }

            If not `None`, restrict to a strand.
        """
        select_clause = constants._SELECT + ' WHERE featuretype = ?'
        args = [featuretype]
        filter_clause = ''

        if seqid is not None:
            filter_clause += ' AND seqid = ?'
            args.append(seqid)

        # TODO: add bin constraints, too

        if start is not None:
            filter_clause += ' AND start >= ?'
            args.append(start)

        if end is not None:
            filter_clause += ' AND end <= ?'
            args.append(end)

        if strand is not None:
            filter_clause += ' AND strand = ?'
            args.append(strand)

        c = self.conn.cursor()
        c.execute(
            select_clause + filter_clause + ' ORDER BY start',
            tuple(args)
        )
        for i in c:
            yield feature.Feature(i, dialect=self.dialect)

    def all_features(self, seqid=None, start=None, end=None, strand=None):
        """
        Returns an iterator of all :class:`Feature` objects in the database.

        Parameters
        ----------

        `seqid` : string

            If not `None`, restrict results to this seqid

        `start` : int

            If not `None`, restrict results to only >= `start`

        `end` : int

            If not `None`, restrict results to only <= `end`

        `strand` : { "-", "+", "." }

            If not `None`, restrict to a strand.

        """
        select_clause = constants._SELECT
        args = []
        filter_clause = 'WHERE'

        if seqid is not None:
            filter_clause += ' AND seqid = ?'
            args.append(seqid)

        # TODO: add bin constraints, too

        if start is not None:
            filter_clause += ' AND start >= ?'
            args.append(start)

        if end is not None:
            filter_clause += ' AND end <= ?'
            args.append(end)

        if strand is not None:
            filter_clause += ' AND strand = ?'
            args.append(strand)

        # Nothing was added above, so reset to empty string
        if filter_clause == 'WHERE':
            filter_clause = ''

        c = self.conn.cursor()
        c.execute(
            select_clause + filter_clause + ' ORDER BY start',
            tuple(args)
        )
        for i in c:
            yield feature.Feature(i, dialect=self.dialect)

    def featuretypes(self):
        """
        Iterator of feature types found in the database.
        """
        c = self.conn.cursor()
        c.execute(
            '''
            SELECT DISTINCT feature from features
            ''')
        for i, in c:
            yield i

    def _relation(self, id, join_on, join_to, level=None, featuretype=None):
        """
        Parameters
        ----------

        `id` : string or a Feature object

        `level` : None or int

            If `level=None` (default), then return all children regardless
            of level.  If `level` is an integer, then constrain to just that
            level.

        `featuretype` : str or None

            If `featuretype` is not `None`, then constrain results to just that
            featuretype.
        """

        if isinstance(id, feature.Feature):
            id = id.id

        c = self.conn.cursor()

        args = [id]

        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = ' AND features.featuretype = ?'
            args.append(featuretype)

        level_clause = ''
        if level is not None:
            level_clause = ' AND relations.level = ?'
            args.append(level)

        # modify _SELECT so that only unique results are returned
        query = constants._SELECT.replace("SELECT", "SELECT DISTINCT")
        c.execute(
            query +
            '''
            JOIN relations
            ON relations.{join_on} = features.id
            WHERE relations.{join_to} = ?
            {featuretype_clause}
            {level_clause}
            ORDER BY start'''.format(**locals()),
            tuple(args))
        for i in c:
            yield feature.Feature(i, dialect=self.dialect)

    def children(self, id, level=None, featuretype=None):
        """
        Return children of feature `id`.

        {_relation_docstring}
        """
        return self._relation(
            id, join_on='child', join_to='parent', level=level,
            featuretype=featuretype)

    def parents(self, id, level=None, featuretype=None):
        """
        Return parents of feature `id`.

        {_relation_docstring}
        """
        return self._relation(
            id, join_on='parent', join_to='child', level=level,
            featuretype=featuretype)

    def execute(self, query):
        """
        Execute arbitrary queries on the db.

        .. seealso::

                :class:`FeatureDB.schema` may be helpful when writing your own
                queries.

        Parameters
        ----------

        `query` : str

            Query to execute -- trailing ";" optional.

        """
        c = self.conn.cursor()
        return c.execute(query)

    def region(self, region, featuretype=None, completely_within=False):
        """
        Return features with any part overlapping `region`.

        Parameters
        ----------
        `region` : string, tuple, or Feature instance
            If string, then of the form "seqid:start-end".  If tuple, then
            (seqid, start, end).  If :class:`Feature`, then use the features
            seqid, start, and end values.

        `featuretype` : None, string, or iterable
            If not None, then restrict output.  If string, then only report
            that feature type.  If iterable, then report all featuretypes in
            the iterable.

        `completely_within` : bool
            If False (default), returns features that overlap `region`, even
            partially.  If True, only return features that are completely
            within `region`.
        """
        if isinstance(region, basestring):
            seqid, coords = region.split(':')
            start, end = coords.split('-')

        elif isinstance(region, feature.Feature):
            seqid = region.seqid
            start = region.start
            end = region.end
        else:
            seqid, start, end = region

        # Get a list of all possible bins for this region
        _bins = list(bins.bins(start, end, one=False))

        if completely_within:
            position_clause = 'start >= ? AND end <= ?'
            args = [seqid, start, end]
        else:
            position_clause = 'start < ? AND end > ?'
            # note start/end swap
            args = [seqid, end, start]

        args += _bins

        _bin_clause = ' or ' .join(['bin = ?' for _ in _bins])

        query = ' '.join([
            constants._SELECT,
            'WHERE seqid = ? AND', position_clause,
            'AND', '(', _bin_clause, ')'])

        # Add the featuretype clause
        if featuretype is not None:
            if isinstance(featuretype, basestring):
                featuretype = [featuretype]
            feature_clause = ' or '.join(
                ['featuretype = ?' for _ in featuretype])
            query += ' AND (%s) ' % feature_clause
            args.extend(featuretype)

        c = self.conn.cursor()
        c.execute(query, tuple(args))
        for i in c:
            yield feature.Feature(i, dialect=self.dialect)

    # Recycle the docs for _relation so they stay consistent between parents()
    # and children()
    children.__doc__ = children.__doc__.format(
        _relation_docstring=_relation.__doc__)
    parents.__doc__ = parents.__doc__.format(
        _relation_docstring=_relation.__doc__)


def create_db(fn, dbfn, id_spec=None, force=False, verbose=True, checklines=10,
              merge_strategy='error', transform=None,
              gtf_transcript_key='transcript_id', gtf_gene_key='gene_id',
              gtf_subfeature='exon', force_gff=False):
    """
    Create a database from a GFF or GTF file.

    Parameters
    ----------
    `fn` : string

        Path to the original GFF or GTF file.

    `dbfn` : string

        Path to the database that will be created.  Can be the special string
        ":memory:" to create an in-memory database.

    `id_spec` : string, list, dict, callable, or None

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
        the parser and returns one of the following:

            * None (in which case the feature type will be auto-incremented)
            * string (which will be used as the primary key)
            * special string starting with "autoincrement:X", where "X" is
              a string that will be used for auto-incrementing.  For example,
              if "autoincrement:chr10", then the first feature will be
              "chr10_1", the second "chr10_2", and so on.

    `force` : bool

        If `False` (default), then raise an exception if `dbfn` already exists.
        Use `force=True` to overwrite any existing databases.

    `verbose` : bool

        Report percent complete and other feedback on how the db creation is
        progressing.

        In order to report percent complete, the entire file needs to be read
        once to see how many items there are; for large files you may want to
        use `verbose=False` to avoid this.

    `checklines` : int

        Number of lines to check the dialect.

    `merge_strategy` : { "merge", "create_unique", "error", "warning" }

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

    `transform` : callable

        Function (or other callable object) that accepts a dictionary and
        returns a dictionary.
    """

    # Auto-detect format
    p = parser.Parser(fn, checklines=checklines, transform=transform)
    dialect = p._sniff()
    if (dialect['fmt'] == 'gff3') or force_gff:
        cls = _GFFDBCreator
        id_spec = id_spec or 'ID'
        kwargs = {}
    elif dialect['fmt'] == 'gtf':
        cls = _GTFDBCreator
        id_spec = id_spec or {'gene': 'gene_id'}
        kwargs = dict(
            transcript_key=gtf_transcript_key,
            gene_key=gtf_gene_key,
            subfeature=gtf_subfeature)

    c = cls(fn, dbfn, id_spec=id_spec, force=force, verbose=verbose,
            merge_strategy=merge_strategy, transform=transform, **kwargs)
    c.create()
    db = FeatureDB(c)
    db.parser = c.parser
    return db
