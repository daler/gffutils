import sqlite3

import create
import bins
import helpers
import constants
import feature

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
        if isinstance(dbfn, create._DBCreator):
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
            WHERE featuretype = ?
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

