import collections
import os
import six
import sqlite3
import shutil
import warnings
from gffutils import bins
from gffutils import helpers
from gffutils import constants
from gffutils import merge_criteria as mc
from gffutils.feature import Feature
from gffutils.exceptions import FeatureNotFoundError


def assign_child(parent, child):
    """
    Helper for add_relation()
    Sets 'Parent' attribute to parent['ID']

    Parameters
    ----------

    parent : Feature
        Parent Feature

    :param parent: parent Feature

    child : Feature
        Child Feature

    Returns
    -------

    Child Feature
    """
    child.attributes['Parent'] = parent['ID']
    return child


# Reusable constant for FeatureDB.merge()
no_children = tuple()

def _finalize_merge(feature, feature_children):
    """
    Helper for FeatureDB.merge() to update source and assign children property

    Parameters
    ----------
    feature : Feature
        feature to finalise
    
    feature_children
        list of children to assign
    
    Returns
    -------
    feature, modified
    """
    if len(feature_children) > 1:
        feature.source = ','.join(set(child.source for child in feature_children))
        feature.children = feature_children
    else:
        feature.children = no_children
    return feature

class FeatureDB(object):
    # docstring to be filled in for methods that call out to
    # helpers.make_query()
    _method_doc = """
        limit : string or tuple
            Limit the results to a genomic region.  If string, then of the form
            "seqid:start-end"; if tuple, then (seqid, start, end).

        strand : "-" | "+" | "."
            Limit the results to one strand

        featuretype : string or tuple
            Limit the results to one or several featuretypes.

        order_by : string or tuple
            Order results by one or many fields; the string or tuple items must
            be in: 'seqid', 'source', 'featuretype', 'start', 'end', 'score',
            'strand', 'frame', 'attributes', 'extra'.

        reverse : bool
            Change sort order; only relevant if `order_by` is not None.  By
            default, results will be in ascending order, so use `reverse=True`
            for descending.

        completely_within : bool
            If False (default), a single bp overlap with `limit` is sufficient
            to return a feature; if True, then the feature must be completely
            within `limit`. Only relevant when `limit` is not None.
    """

    def __init__(self, dbfn, default_encoding='utf-8', keep_order=False,
                 pragmas=constants.default_pragmas,
                 sort_attribute_values=False,
                 text_factory=sqlite3.OptimizedUnicode):
        """
        Connect to a database created by :func:`gffutils.create_db`.

        Parameters
        ----------

        dbfn : str

            Path to a database created by :func:`gffutils.create_db`.

        text_factory : callable

            Optionally set the way sqlite3 handles strings.  Default is
            sqlite3.OptimizedUnicode, which returns ascii when possible,
            unicode otherwise

        encoding : str

            When non-ASCII characters are encountered, assume they are in this
            encoding.

        keep_order : bool

            If True, all features returned from this instance will have the
            order of their attributes maintained.  This can be turned on or off
            database-wide by setting the `keep_order` attribute or with this
            kwarg, or on a feature-by-feature basis by setting the `keep_order`
            attribute of an individual feature.

            Default is False, since this includes a sorting step that can get
            time-consuming for many features.

        pragmas : dict
            Dictionary of pragmas to use when connecting to the database.  See
            http://www.sqlite.org/pragma.html for the full list of
            possibilities, and constants.default_pragmas for the defaults.
            These can be changed later using the :meth:`FeatureDB.set_pragmas`
            method.

        Notes
        -----

            `dbfn` can also be a subclass of :class:`_DBCreator`, useful for
            when :func:`gffutils.create_db` is provided the ``dbfn=":memory:"``
            kwarg.

        """
        # Since specifying the string ":memory:" will actually try to connect
        # to a new, separate (and empty) db in memory, we can alternatively
        # pass in a sqlite connection instance to use its existing, in-memory
        # db.
        from gffutils import create
        if isinstance(dbfn, create._DBCreator):
            self.conn = dbfn.conn
            self.dbfn = dbfn.dbfn

        elif isinstance(dbfn, sqlite3.Connection):
            self.conn = dbfn
            self.dbfn = dbfn
        # otherwise assume dbfn is a string.
        elif dbfn == ':memory:':
            raise ValueError(
                "cannot connect to memory db; please provide the connection")
        else:
            if not os.path.exists(dbfn):
                raise ValueError("Database file %s does not exist" % dbfn)
            self.dbfn = dbfn
            self.conn = sqlite3.connect(self.dbfn)

        if text_factory is not None:
            self.conn.text_factory = text_factory
        self.conn.row_factory = sqlite3.Row

        self.default_encoding = default_encoding
        self.keep_order = keep_order
        self.sort_attribute_values = sort_attribute_values
        c = self.conn.cursor()

        # Load some meta info
        # TODO: this is a good place to check for previous versions, and offer
        # to upgrade...
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
        # autoincrementing from where we last left off (instead of from 1,
        # which would cause name collisions)
        c.execute(
            '''
            SELECT base, n FROM autoincrements
            ''')
        self._autoincrements = collections.defaultdict(int, c)

        self.set_pragmas(pragmas)

        if not self._analyzed():
            warnings.warn(
                "It appears that this database has not had the ANALYZE "
                "sqlite3 command run on it. Doing so can dramatically "
                "speed up queries, and is done by default for databases "
                "created with gffutils >0.8.7.1 (this database was "
                "created with version %s) Consider calling the analyze() "
                "method of this object." % self.version)

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

    def _feature_returner(self, **kwargs):
        """
        Returns a feature, adding additional database-specific defaults
        """
        kwargs.setdefault('dialect', self.dialect)
        kwargs.setdefault('keep_order', self.keep_order)
        kwargs.setdefault('sort_attribute_values', self.sort_attribute_values)
        return Feature(**kwargs)

    def _analyzed(self):
        res = self.execute(
            """
            SELECT name FROM sqlite_master WHERE type='table'
            AND name='sqlite_stat1';
            """)
        return len(list(res)) == 1

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
        if isinstance(key, Feature):
            key = key.id
        c = self.conn.cursor()
        try:
            c.execute(constants._SELECT + ' WHERE id = ?', (key,))
        except sqlite3.ProgrammingError:
            c.execute(
                constants._SELECT + ' WHERE id = ?',
                (key.decode(self.default_encoding),))
        results = c.fetchone()
        # TODO: raise error if more than one key is found
        if results is None:
            raise FeatureNotFoundError(key)
        return self._feature_returner(**results)

    def count_features_of_type(self, featuretype=None):
        """
        Simple count of features.

        Can be faster than "grep", and is faster than checking the length of
        results from :meth:`gffutils.FeatureDB.features_of_type`.

        Parameters
        ----------

        featuretype : string

            Feature type (e.g., "gene") to count.  If None, then count *all*
            features in the database.

        Returns
        -------
        The number of features of this type, as an integer

        """
        c = self.conn.cursor()
        if featuretype is not None:
            c.execute(
                '''
                SELECT count() FROM features
                WHERE featuretype = ?
                ''', (featuretype,))
        else:
            c.execute(
                '''
                SELECT count() FROM features
                ''')

        results = c.fetchone()
        if results is not None:
            results = results[0]
        return results

    def features_of_type(self, featuretype, limit=None, strand=None,
                         order_by=None, reverse=False,
                         completely_within=False):
        """
        Returns an iterator of :class:`gffutils.Feature` objects.

        Parameters
        ----------
        {_method_doc}
        """
        query, args = helpers.make_query(
            args=[],
            limit=limit,
            featuretype=featuretype,
            order_by=order_by,
            reverse=reverse,
            strand=strand,
            completely_within=completely_within,
        )

        for i in self._execute(query, args):
            yield self._feature_returner(**i)

    # TODO: convert this to a syntax similar to itertools.groupby
    def iter_by_parent_childs(self, featuretype="gene", level=None,
                              order_by=None, reverse=False,
                              completely_within=False):
        """
        For each parent of type `featuretype`, yield a list L of that parent
        and all of its children (`[parent] + list(children)`). The parent will
        always be L[0].

        This is useful for "sanitizing" a GFF file for downstream tools.

        Additional kwargs are passed to :meth:`FeatureDB.children`, and will
        therefore only affect items L[1:] in each yielded list.
        """
        # Get all the parent records of the requested feature type
        parent_recs = self.all_features(featuretype=featuretype)
        # Return a generator of these parent records and their
        # children
        for parent_rec in parent_recs:
            unit_records = \
                [parent_rec] + list(self.children(parent_rec.id))
            yield unit_records

    def all_features(self, limit=None, strand=None, featuretype=None,
                     order_by=None, reverse=False, completely_within=False):
        """
        Iterate through the entire database.

        Returns
        -------
        A generator object that yields :class:`Feature` objects.

        Parameters
        ----------
        {_method_doc}
        """
        query, args = helpers.make_query(
            args=[],
            limit=limit,
            strand=strand,
            featuretype=featuretype,
            order_by=order_by,
            reverse=reverse,
            completely_within=completely_within
        )
        for i in self._execute(query, args):
            yield self._feature_returner(**i)

    def featuretypes(self):
        """
        Iterate over feature types found in the database.

        Returns
        -------
        A generator object that yields featuretypes (as strings)
        """
        c = self.conn.cursor()
        c.execute(
            '''
            SELECT DISTINCT featuretype from features
            ''')
        for i, in c:
            yield i

    def _relation(self, id, join_on, join_to, level=None, featuretype=None,
                  order_by=None, reverse=False, completely_within=False,
                  limit=None):

        # The following docstring will be included in the parents() and
        # children() docstrings to maintain consistency, since they both
        # delegate to this method.
        """
        Parameters
        ----------

        id : string or a Feature object

        level : None or int

            If `level=None` (default), then return all children regardless
            of level.  If `level` is an integer, then constrain to just that
            level.
        {_method_doc}

        Returns
        -------
        A generator object that yields :class:`Feature` objects.
        """

        if isinstance(id, Feature):
            id = id.id

        other = '''
        JOIN relations
        ON relations.{join_on} = features.id
        WHERE relations.{join_to} = ?
        '''.format(**locals())
        args = [id]

        level_clause = ''
        if level is not None:
            level_clause = 'relations.level = ?'
            args.append(level)

        query, args = helpers.make_query(
            args=args,
            other=other,
            extra=level_clause,
            featuretype=featuretype,
            order_by=order_by,
            reverse=reverse,
            limit=limit,
            completely_within=completely_within,
        )

        # modify _SELECT so that only unique results are returned
        query = query.replace("SELECT", "SELECT DISTINCT")
        for i in self._execute(query, args):
            yield self._feature_returner(**i)

    def children(self, id, level=None, featuretype=None, order_by=None,
                 reverse=False, limit=None, completely_within=False):
        """
        Return children of feature `id`.
        {_relation_docstring}
        """
        return self._relation(
            id, join_on='child', join_to='parent', level=level,
            featuretype=featuretype, order_by=order_by, reverse=reverse,
            limit=limit, completely_within=completely_within)

    def parents(self, id, level=None, featuretype=None, order_by=None,
                reverse=False, completely_within=False, limit=None):
        """
        Return parents of feature `id`.
        {_relation_docstring}
        """
        return self._relation(
            id, join_on='parent', join_to='child', level=level,
            featuretype=featuretype, order_by=order_by, reverse=reverse,
            limit=limit, completely_within=completely_within)

    def _execute(self, query, args):
        self._last_query = query
        self._last_args = args
        c = self.conn.cursor()
        c.execute(query, tuple(args))
        return c

    def execute(self, query):
        """
        Execute arbitrary queries on the db.

        .. seealso::

                :class:`FeatureDB.schema` may be helpful when writing your own
                queries.

        Parameters
        ----------

        query : str

            Query to execute -- trailing ";" optional.

        Returns
        -------
        A sqlite3.Cursor object that can be iterated over.
        """
        c = self.conn.cursor()
        return c.execute(query)

    def analyze(self):
        """
        Runs the sqlite ANALYZE command to potentially speed up queries
        dramatically.
        """
        self.execute('ANALYZE features')
        self.conn.commit()

    def region(self, region=None, seqid=None, start=None, end=None,
               strand=None, featuretype=None, completely_within=False):
        """
        Return features within specified genomic coordinates.

        Specifying genomic coordinates can be done in a flexible manner

        Parameters
        ----------
        region : string, tuple, or Feature instance
            If string, then of the form "seqid:start-end".  If tuple, then
            (seqid, start, end).  If :class:`Feature`, then use the features
            seqid, start, and end values.

            This argument is mutually exclusive with start/end/seqid.

            *Note*: By design, even if a feature is provided, its strand will
            be ignored.  If you want to restrict the output by strand, use the
            separate `strand` kwarg.

        strand : + | - | . | None
            If `strand` is provided, then only those features exactly matching
            `strand` will be returned. So `strand='.'` will only return
            unstranded features. Default is `strand=None` which does not
            restrict by strand.

        seqid, start, end, strand
            Mutually exclusive with `region`.  These kwargs can be used to
            approximate slice notation; see "Details" section below.

        featuretype : None, string, or iterable
            If not None, then restrict output.  If string, then only report
            that feature type.  If iterable, then report all featuretypes in
            the iterable.

        completely_within : bool
            By default (`completely_within=False`), returns features that
            partially or completely overlap `region`.  If
            `completely_within=True`, features that are completely within
            `region` will be returned.

        Notes
        -------

        The meaning of `seqid`, `start`, and `end` is interpreted as follows:

        ====== ====== ===== ======================================
        seqid  start  end   meaning
        ====== ====== ===== ======================================
        str    int    int   equivalent to `region` kwarg
        None   int    int   features from all chroms within coords
        str    None   int   equivalent to [:end] slice notation
        str    int    None  equivalent to [start:] slice notation
        None   None   None  equivalent to FeatureDB.all_features()
        ====== ====== ===== ======================================

        If performance is a concern, use `completely_within=True`. This allows
        the query to be optimized by only looking for features that fall in the
        precise genomic bin (same strategy as UCSC Genome Browser and
        BEDTools). Otherwise all features' start/stop coords need to be
        searched to see if they partially overlap the region of interest.

        Examples
        --------

        - `region(seqid="chr1", start=1000)` returns all features on chr1 that
          start or extend past position 1000

        - `region(seqid="chr1", start=1000, completely_within=True)` returns
          all features on chr1 that start past position 1000.

        - `region("chr1:1-100", strand="+", completely_within=True)` returns
          only plus-strand features that completely fall within positions 1 to
          100 on chr1.

        Returns
        -------
        A generator object that yields :class:`Feature` objects.
        """
        # Argument handling.
        if region is not None:
            if (seqid is not None) or (start is not None) or (end is not None):
                raise ValueError(
                    "If region is supplied, do not supply seqid, "
                    "start, or end as separate kwargs")
            if isinstance(region, six.string_types):
                toks = region.split(':')
                if len(toks) == 1:
                    seqid = toks[0]
                    start, end = None, None
                else:
                    seqid, coords = toks[:2]
                    if len(toks) == 3:
                        strand = toks[2]
                    start, end = coords.split('-')

            elif isinstance(region, Feature):
                seqid = region.seqid
                start = region.start
                end = region.end
                strand = region.strand

            # otherwise assume it's a tuple
            else:
                seqid, start, end = region[:3]

        # e.g.,
        #   completely_within=True..... start >= {start} AND end <= {end}
        #   completely_within=False.... start <  {end}   AND end >  {start}
        if completely_within:
            start_op = '>='
            end_op = '<='
        else:
            start_op = '<'
            end_op = '>'
            end, start = start, end

        args = []
        position_clause = []
        if seqid is not None:
            position_clause.append('seqid = ?')
            args.append(seqid)
        if start is not None:
            start = int(start)
            position_clause.append('start %s ?' % start_op)
            args.append(start)
        if end is not None:
            end = int(end)
            position_clause.append('end %s ?' % end_op)
            args.append(end)

        position_clause = ' AND '.join(position_clause)

        # Only use bins if we have defined boundaries and completely_within is
        # True. Otherwise you can't know how far away a feature stretches
        # (which means bins are not computable ahead of time)
        _bin_clause = ''
        if (start is not None) and (end is not None) and completely_within:
            if start <= bins.MAX_CHROM_SIZE and end <= bins.MAX_CHROM_SIZE:
                _bins = list(bins.bins(start, end, one=False))
                # See issue #45
                if len(_bins) < 900:
                    _bin_clause = ' or ' .join(['bin = ?' for _ in _bins])
                    _bin_clause = 'AND ( %s )' % _bin_clause
                    args += _bins

        query = ' '.join([
            constants._SELECT,
            'WHERE ',
            position_clause,
            _bin_clause])

        # Add the featuretype clause
        if featuretype is not None:
            if isinstance(featuretype, six.string_types):
                featuretype = [featuretype]
            feature_clause = ' or '.join(
                ['featuretype = ?' for _ in featuretype])
            query += ' AND (%s) ' % feature_clause
            args.extend(featuretype)

        if strand is not None:
            strand_clause = ' and strand = ? '
            query += strand_clause
            args.append(strand)

        c = self.conn.cursor()
        self._last_query = query
        self._last_args = args
        self._context = {
            'start': start,
            'end': end,
            'seqid': seqid,
            'region': region,
        }
        c.execute(query, tuple(args))
        for i in c:
            yield self._feature_returner(**i)

    def interfeatures(self, features, new_featuretype=None,
                      merge_attributes=True, dialect=None,
                      attribute_func=None, update_attributes=None):
        """
        Construct new features representing the space between features.

        For example, if `features` is a list of exons, then this method will
        return the introns.  If `features` is a list of genes, then this method
        will return the intergenic regions.

        Providing N features will return N - 1 new features.

        This method purposefully does *not* do any merging or sorting of
        coordinates, so you may want to use :meth:`FeatureDB.merge` first, or
        when selecting features use the `order_by` kwarg, e.g.,
        `db.features_of_type('gene', order_by=('seqid', 'start'))`.

        Parameters
        ----------
        features : iterable of :class:`feature.Feature` instances
            Sorted, merged iterable

        new_featuretype : string or None
            The new features will all be of this type, or, if None (default)
            then the featuretypes will be constructed from the neighboring
            features, e.g., `inter_exon_exon`.

        merge_attributes : bool
            If True, new features' attributes will be a merge of the neighboring
            features' attributes.  This is useful if you have provided a list of
            exons; the introns will then retain the transcript and/or gene
            parents as a single item. Otherwise, if False, the attribute will
            be a comma-separated list of values, potentially listing the same
            gene ID twice.

        attribute_func : callable or None
            If None, then nothing special is done to the attributes.  If
            callable, then the callable accepts two attribute dictionaries and
            returns a single attribute dictionary.  If `merge_attributes` is
            True, then `attribute_func` is called before `merge_attributes`.
            This could be useful for manually managing IDs for the new
            features.

        update_attributes : dict
            After attributes have been modified and merged, this dictionary can
            be used to replace parts of the attributes dictionary.

        Returns
        -------
        A generator that yields :class:`Feature` objects
        """
        for i, f in enumerate(features):
            # no inter-feature for the first one
            if i == 0:
                interfeature_start = f.stop
                last_feature = f
                continue

            interfeature_stop = f.start
            if new_featuretype is None:
                new_featuretype = 'inter_%s_%s' % (
                    last_feature.featuretype, f.featuretype)
            if last_feature.strand != f.strand:
                new_strand = '.'
            else:
                new_strand = f.strand

            if last_feature.chrom != f.chrom:
                # We've moved to a new chromosome.  For example, if we're
                # getting intergenic regions from all genes, they will be on
                # different chromosomes. We still assume sorted features, but
                # don't complain if they're on different chromosomes -- just
                # move on.
                last_feature = f
                continue

            strand = new_strand
            chrom = last_feature.chrom

            # Shrink
            interfeature_start += 1
            interfeature_stop -= 1

            if merge_attributes:
                new_attributes = helpers.merge_attributes(
                    last_feature.attributes, f.attributes)
            else:
                new_attributes = {}

            if update_attributes:
                new_attributes.update(update_attributes)

            new_bin = bins.bins(
                interfeature_start, interfeature_stop, one=True)
            _id = None
            fields = dict(
                seqid=chrom,
                source='gffutils_derived',
                featuretype=new_featuretype,
                start=interfeature_start,
                end=interfeature_stop,
                score='.',
                strand=strand,
                frame='.',
                attributes=new_attributes,
                bin=new_bin)

            if dialect is None:
                # Support for @classmethod -- if calling from the class, then
                # self.dialect is not defined, so defer to Feature's default
                # (which will be constants.dialect, or GFF3).
                try:
                    dialect = self.dialect
                except AttributeError:
                    dialect = None
            yield self._feature_returner(**fields)
            interfeature_start = f.stop
            last_feature = f

    def delete(self, features, make_backup=True, **kwargs):
        """
        Delete features from database.

        features : str, iterable, FeatureDB instance
            If FeatureDB, all features will be used. If string, assume it's the
            ID of the feature to remove. Otherwise, assume it's an iterable of
            Feature objects. The classes in gffutils.iterators may be helpful
            in this case.

        make_backup : bool
            If True, and the database you're about to update is a file on disk,
            makes a copy of the existing database and saves it with a .bak
            extension.

        Returns
        -------
        FeatureDB object, with features deleted.
        """
        if make_backup:
            if isinstance(self.dbfn, six.string_types):
                shutil.copy2(self.dbfn, self.dbfn + '.bak')

        c = self.conn.cursor()
        query1 = """
        DELETE FROM features WHERE id = ?
        """
        query2 = """
        DELETE FROM relations WHERE parent = ? OR child = ?
        """
        if isinstance(features, FeatureDB):
            features = features.all_features()
        if isinstance(features, six.string_types):
            features = [features]
        if isinstance(features, Feature):
            features = [features]
        for feature in features:
            if isinstance(feature, six.string_types):
                _id = feature
            else:
                _id = feature.id
            c.execute(query1, (_id,))
            c.execute(query2, (_id, _id))
        self.conn.commit()
        return self

    def update(self, data, make_backup=True, **kwargs):
        """
        Update the on-disk database with features in `data`.

        WARNING: If you used any non-default kwargs for gffutils.create_db when
        creating the database in the first place (especially
        `disable_infer_transcripts` or `disable_infer_genes`) then you should
        use those same arguments here.

        The returned object is the same FeatureDB, but since it is pointing to
        the same database and that has been just updated, the new features can
        now be accessed.

        Parameters
        ----------

        data : str, iterable, FeatureDB instance
            If FeatureDB, all data will be used. If string, assume it's
            a filename of a GFF or GTF file.  Otherwise, assume it's an
            iterable of Feature objects.  The classes in gffutils.iterators may
            be helpful in this case.

        make_backup : bool
            If True, and the database you're about to update is a file on disk,
            makes a copy of the existing database and saves it with a .bak
            extension.

        kwargs :
            Additional kwargs are passed to gffutils.create_db; see the help
            for that function for details.

        Returns
        -------
        Returns self but with the underlying database updated.
        """

        from gffutils import create
        from gffutils import iterators
        if make_backup:
            if isinstance(self.dbfn, six.string_types):
                shutil.copy2(self.dbfn, self.dbfn + '.bak')

        # get iterator-specific kwargs
        _iterator_kwargs = {}
        for k, v in kwargs.items():
            if k in constants._iterator_kwargs:
                _iterator_kwargs[k] = v

        # Handle all sorts of input
        data = iterators.DataIterator(data, **_iterator_kwargs)
        kwargs['_autoincrements'] = self._autoincrements

        if self.dialect['fmt'] == 'gtf':
            if 'id_spec' not in kwargs:
                kwargs['id_spec'] = {
                    'gene': 'gene_id', 'transcript': 'transcript_id'}
            db = create._GTFDBCreator(
                data=data, dbfn=self.dbfn, dialect=self.dialect, **kwargs)
        elif self.dialect['fmt'] == 'gff3':
            if 'id_spec' not in kwargs:
                kwargs['id_spec'] = 'ID'
            db = create._GFFDBCreator(
                data=data, dbfn=self.dbfn, dialect=self.dialect, **kwargs)

        else:
            raise ValueError

        peek, data._iter = iterators.peek(data._iter, 1)
        if len(peek) == 0: return db  # If the file is empty then do nothing

        db._populate_from_lines(data)
        db._update_relations()

        # Note that the autoincrements gets updated here
        db._finalize()

        # Read it back in directly from the stored autoincrements table
        self._autoincrements.update(db._autoincrements)
        return self

    def add_relation(self, parent, child, level, parent_func=None,
                     child_func=None):
        """
        Manually add relations to the database.

        Parameters
        ----------
        parent : str or Feature instance
             Parent feature to add.

        child : str or Feature instance
            Child feature to add

        level : int
            Level of the relation.  For example, if parent is a gene and child
            is an mRNA, then you might want level to be 1.  But if child is an
            exon, then level would be 2.

        parent_func, child_func : callable
            These optional functions control how attributes are updated in the
            database.  They both have the signature `func(parent, child)` and
            must return a [possibly modified] Feature instance.  For example,
            we could add the child's database id as the "child" attribute in
            the parent::

                def parent_func(parent, child):
                    parent.attributes['child'] = child.id

            and add the parent's "gene_id" as the child's "Parent" attribute::

                def child_func(parent, child):
                    child.attributes['Parent'] = parent['gene_id']

        Returns
        -------
        FeatureDB object with new relations added.
        """
        if isinstance(parent, six.string_types):
            parent = self[parent]
        if isinstance(child, six.string_types):
            child = self[child]

        c = self.conn.cursor()
        c.execute('''
                  INSERT INTO relations (parent, child, level)
                  VALUES (?, ?, ?)''', (parent.id, child.id, level))

        if parent_func is not None:
            parent = parent_func(parent, child)
            self._update(parent, c)
        if child_func is not None:
            child = child_func(parent, child)
            self._update(child, c)

        self.conn.commit()
        return self

    def _update(self, feature, cursor):
        values = [list(feature.astuple()) + [feature.id]]
        cursor.execute(
            constants._UPDATE, *tuple(values))

    def _insert(self, feature, cursor):
        """
        Insert a feature into the database.
        """
        try:
            cursor.execute(constants._INSERT, feature.astuple())
        except sqlite3.ProgrammingError:
            cursor.execute(
                constants._INSERT, feature.astuple(self.default_encoding))

    def create_introns(self, exon_featuretype='exon',
                       grandparent_featuretype='gene', parent_featuretype=None,
                       new_featuretype='intron', merge_attributes=True):
        """
        Create introns from existing annotations.


        Parameters
        ----------
        exon_featuretype : string
            Feature type to use in order to infer introns.  Typically `"exon"`.

        grandparent_featuretype : string
            If `grandparent_featuretype` is not None, then group exons by
            children of this featuretype.  If `granparent_featuretype` is
            "gene" (default), then introns will be created for all first-level
            children of genes.  This may include mRNA, rRNA, ncRNA, etc.  If
            you only want to infer introns from one of these featuretypes
            (e.g., mRNA), then use the `parent_featuretype` kwarg which is
            mutually exclusive with `grandparent_featuretype`.

        parent_featuretype : string
            If `parent_featuretype` is not None, then only use this featuretype
            to infer introns.  Use this if you only want a subset of
            featuretypes to have introns (e.g., "mRNA" only, and not ncRNA or
            rRNA). Mutually exclusive with `grandparent_featuretype`.

        new_featuretype : string
            Feature type to use for the inferred introns; default is
            `"intron"`.

        merge_attributes : bool
            Whether or not to merge attributes from all exons. If False then no
            attributes will be created for the introns.

        Returns
        -------
        A generator object that yields :class:`Feature` objects representing
        new introns

        Notes
        -----
        The returned generator can be passed directly to the
        :meth:`FeatureDB.update` method to permanently add them to the
        database, e.g., ::

            db.update(db.create_introns())

        """
        if (grandparent_featuretype and parent_featuretype) or (
            grandparent_featuretype is None and parent_featuretype is None
        ):
            raise ValueError("exactly one of `grandparent_featuretype` or "
                             "`parent_featuretype` should be provided")

        if grandparent_featuretype:
            def child_gen():
                for gene in self.features_of_type(grandparent_featuretype):
                    for child in self.children(gene, level=1):
                        yield child
        elif parent_featuretype:
            def child_gen():
                for child in self.features_of_type(parent_featuretype):
                    yield child

        for child in child_gen():
            exons = self.children(child, level=1, featuretype=exon_featuretype,
                                  order_by='start')
            for intron in self.interfeatures(
                exons, new_featuretype=new_featuretype,
                merge_attributes=merge_attributes, dialect=self.dialect
            ):
                yield intron

    def _old_merge(self, features, ignore_strand=False):
        """
        DEPRECATED, only retained here for backwards compatibility. Please use
        merge().

        Merge overlapping features together.

        Parameters
        ----------

        features : iterator of Feature instances

        ignore_strand : bool
            If True, features on multiple strands will be merged, and the final
            strand will be set to '.'.  Otherwise, ValueError will be raised if
            trying to merge features on differnt strands.

        Returns
        -------
        A generator object that yields :class:`Feature` objects representing
        the newly merged features.
        """

        # Consume iterator up front...
        features = list(features)

        if len(features) == 0:
            raise StopIteration

        # Either set all strands to '+' or check for strand-consistency.
        if ignore_strand:
            strand = '.'
        else:
            strands = [i.strand for i in features]
            if len(set(strands)) > 1:
                raise ValueError('Specify ignore_strand=True to force merging '
                                 'of multiple strands')
            strand = strands[0]

        # Sanity check to make sure all features are from the same chromosome.
        chroms = [i.chrom for i in features]
        if len(set(chroms)) > 1:
            raise NotImplementedError('Merging multiple chromosomes not '
                                      'implemented')
        chrom = chroms[0]

        # To start, we create a merged feature of just the first feature.
        current_merged_start = features[0].start
        current_merged_stop = features[0].stop

        # We don't need to check the first one, so start at feature #2.
        for feature in features[1:]:
            # Does this feature start within the currently merged feature?...
            if feature.start <= current_merged_stop + 1:
                # ...It starts within, so leave current_merged_start where it
                # is.  Does it extend any farther?
                if feature.stop >= current_merged_stop:
                    # Extends further, so set a new stop position
                    current_merged_stop = feature.stop
                else:
                    # If feature.stop < current_merged_stop, it's completely
                    # within the previous feature.  Nothing more to do.
                    continue
            else:
                # The start position is outside the merged feature, so we're
                # done with the current merged feature.  Prepare for output...
                merged_feature = dict(
                    seqid=feature.chrom,
                    source='.',
                    featuretype=feature.featuretype,
                    start=current_merged_start,
                    end=current_merged_stop,
                    score='.',
                    strand=strand,
                    frame='.',
                    attributes='')
                yield self._feature_returner(**merged_feature)

                # and we start a new one, initializing with this feature's
                # start and stop.
                current_merged_start = feature.start
                current_merged_stop = feature.stop

        # need to yield the last one.
        if len(features) == 1:
            feature = features[0]
        merged_feature = dict(
            seqid=feature.chrom,
            source='.',
            featuretype=feature.featuretype,
            start=current_merged_start,
            end=current_merged_stop,
            score='.',
            strand=strand,
            frame='.',
            attributes='')
        yield self._feature_returner(**merged_feature)


    def merge(self, features,
              merge_criteria=(mc.seqid, mc.overlap_end_inclusive, mc.strand, mc.feature_type),
              multiline=False):
        """
        Merge features matching criteria together


        Returned Features have a special property called 'children' that is
        a list of the component features.  This only exists for the lifetime of
        the Feature instance.

        Parameters
        ----------
        features : iterable
            Iterable of Feature instances to merge

        merge_criteria : list
            List of merge criteria callbacks. All must evaluate to True in
            order for a feature to be merged. See notes below on callback
            signature.

        multiline : bool
            True to emit multiple features with the same ID attribute, False
            otherwise.

        Returns
        -------
        Generator yielding merged Features

        Notes
        -----

        See the `gffutils.merge_criteria` module (imported here as `mc`) for
        existing callback functions. For writing custom callbacks, functions
        must have the following signature::

            callback(
                acc: gffutils.Feature,
                cur: gffutils.Feature,
                components: [gffutils.Feature]
            ) -> bool

        Where:

            - `acc`: current accumulated feature
            - `cur`: candidate feature to merge
            - `components`: list of features that compose acc

        The function should return True to merge `cur` into `acc`, False to set
        `cur` to `acc` (that is, start a new merged feature).


        If merge criteria allows different feature types then the merged
        features' feature types should have their featuretype property
        reassigned to a more specific ontology value.
        """
        if not isinstance(merge_criteria, list):
            try:
                merge_criteria = list(merge_criteria)
            except TypeError:
                merge_criteria = [merge_criteria]

        # To start, we create a merged feature of just the first feature.
        features = iter(features)
        last_id = None
        current_merged = None
        feature_children = []

        for feature in features:
            if current_merged is None:
                if all(criteria(feature, feature, feature_children) for criteria in merge_criteria):
                    current_merged = feature
                    feature_children = [feature]
                else:
                    yield _finalize_merge(feature, no_children)
                    last_id = None
                continue

            if len(feature_children) == 0:  # current_merged is last feature and unchecked
                if all(criteria(current_merged, current_merged, feature_children) for criteria in merge_criteria):
                    feature_children.append(current_merged)
                else:
                    yield _finalize_merge(current_merged, no_children)
                    current_merged = feature
                    last_id = None
                    continue

            if all(criteria(current_merged, feature, feature_children) for criteria in merge_criteria):
                # Criteria satisfied, merge
                # TODO Test multiline records and iron out the following code
                # if multiline and (feature.start > current_merged.end + 1 or feature.end + 1 < current_merged.start):
                #    # Feature is possibly multiline (discontiguous), keep ID but start new record
                #    yield _finalize_merge(current_merged, feature_children)
                #    current_merged = feature
                #    feature_children = [feature]

                if len(feature_children) == 1:
                    # Current merged is only child and merge is going to occur, make copy
                    current_merged = vars(current_merged).copy()
                    del current_merged['attributes']
                    del current_merged['extra']
                    del current_merged['dialect']
                    del current_merged['keep_order']
                    del current_merged['sort_attribute_values']
                    current_merged = self._feature_returner(**current_merged)
                    if not last_id:
                        # Generate unique ID for new Feature
                        self._autoincrements[current_merged.featuretype] += 1
                        last_id = current_merged.featuretype + '_' + str(
                            self._autoincrements[current_merged.featuretype])
                    current_merged['ID'] = last_id
                    current_merged.id = last_id

                feature_children.append(feature)

                # Set mismatched properties to ambiguous values
                if feature.seqid not in current_merged.seqid.split(','): current_merged.seqid += ',' + feature.seqid
                if feature.strand != current_merged.strand: current_merged.strand = '.'
                if feature.frame != current_merged.frame: current_merged.frame = '.'
                if feature.featuretype != current_merged.featuretype: current_merged.featuretype = "sequence_feature"

                if feature.start < current_merged.start:
                    # Extends prior, so set a new start position
                    current_merged.start = feature.start

                if feature.end > current_merged.end:
                    # Extends further, so set a new stop position
                    current_merged.end = feature.end

            else:
                yield _finalize_merge(current_merged, feature_children)
                current_merged = feature
                feature_children = []
                last_id = None

        if current_merged:
            yield _finalize_merge(current_merged, feature_children)

    def merge_all(self,
                  merge_order=('seqid', 'featuretype', 'strand', 'start'),
                  merge_criteria=(mc.seqid, mc.overlap_end_inclusive, mc.strand, mc.feature_type),
                  featuretypes_groups=(None,),
                  exclude_components=False):
        """
        Merge all features in database according to criteria.
        Merged features will be assigned as children of the merged record.
        The resulting records are added to the database.

        Parameters
        ----------
        merge_order : list
            Ordered list of columns with which to group features before evaluating criteria
        merge_criteria : list
            List of merge criteria callbacks. See merge().
        featuretypes_groups : list
            iterable of sets of featuretypes to merge together
        exclude_components : bool
            True: child features will be discarded. False to keep them.

        Returns
        -------
            list of merge features
        """

        if not len(featuretypes_groups):
            # Can't be empty
            featuretypes_groups = (None,)

        result_features = []

        # Merge features per featuregroup
        for featuregroup in featuretypes_groups:
            for merged in self.merge(self.all_features(featuretype=featuregroup, order_by=merge_order),
                                merge_criteria=merge_criteria):
                # If feature is result of merge
                if merged.children:
                    self._insert(merged, self.conn.cursor())
                    if exclude_components:
                        # Remove child features from DB
                        self.delete(merged.children)
                    else:
                        # Add child relations to DB
                        for child in merged.children:
                            self.add_relation(merged, child, 1, child_func=assign_child)
                    result_features.append(merged)
                else:
                    pass  # Do nothing, feature is already in DB

        return result_features

    def children_bp(self, feature, child_featuretype='exon', merge=False,
                    ignore_strand=False):
        """
        Total bp of all children of a featuretype.

        Useful for getting the exonic bp of an mRNA.

        Parameters
        ----------

        feature : str or Feature instance

        child_featuretype : str
            Which featuretype to consider.  For example, to get exonic bp of an
            mRNA, use `child_featuretype='exon'`.

        merge : bool
            Whether or not to merge child features together before summing
            them.

        ignore_strand : bool
            If True, then overlapping features on different strands will be
            merged together; otherwise, merging features with different strands
            will result in a ValueError.

        Returns
        -------
        Integer representing the total number of bp.
        """

        children = self.children(feature, featuretype=child_featuretype,
                                 order_by='start')
        if merge:
            children = self.merge(children, ignore_strand=ignore_strand)

        total = 0
        for child in children:
            total += len(child)
        return total

    def bed12(self, feature, block_featuretype=['exon'],
              thick_featuretype=['CDS'], thin_featuretype=None,
              name_field='ID', color=None):
        """
        Converts `feature` into a BED12 format.

        GFF and GTF files do not necessarily define genes consistently, so this
        method provides flexiblity in specifying what to call a "transcript".

        Parameters
        ----------
        feature : str or Feature instance
            In most cases, this feature should be a transcript rather than
            a gene.

        block_featuretype : str or list
            Which featuretype to use as the exons. These are represented as
            blocks in the BED12 format.  Typically 'exon'.

            Use the `thick_featuretype` and `thin_featuretype` arguments to
            control the display of CDS as thicker blocks and UTRs as thinner
            blocks.

            Note that the features for `thick` or `thin` are *not*
            automatically included in the blocks; if you do want them included,
            then those featuretypes should be added to this `block_features`
            list.

            If no child features of type `block_featuretype` are found, then
            the full `feature` is returned in BED12 format as if it had
            a single exon.

        thick_featuretype : str or list
            Child featuretype(s) to use in order to determine the boundaries of
            the "thick" blocks. In BED12 format, these represent coding
            sequences; typically this would be set to "CDS".  This argument is
            mutually exclusive with `thin_featuretype`.

            Specifically, the BED12 thickStart will be the start coord of the
            first `thick` item and the thickEnd will be the stop coord of the
            last `thick` item.

        thin_featuretype : str or list
            Child featuretype(s) to use in order to determine the boundaries of
            the "thin" blocks.  In BED12 format, these represent untranslated
            regions.  Typically "utr" or ['three_prime_UTR', 'five_prime_UTR'].
            Mutually exclusive with `thick_featuretype`.

            Specifically, the BED12 thickStart field will be the stop coord of
            the first `thin` item and the thickEnd field will be the start
            coord of the last `thin` item.

        name_field : str
            Which attribute of `feature` to use as the feature's name.  If this
            field is not present, a "." placeholder will be used instead.

        color : None or str
            If None, then use black (0,0,0) as the RGB color; otherwise this
            should be a comma-separated string of R,G,B values each of which
            are integers in the range 0-255.
        """
        if thick_featuretype and thin_featuretype:
            raise ValueError("Can only specify one of `thick_featuertype` or "
                             "`thin_featuretype`")

        exons = list(self.children(feature, featuretype=block_featuretype,
                                   order_by='start'))
        if len(exons) == 0:
            exons = [feature]
        feature = self[feature]
        first = exons[0].start
        last = exons[-1].stop

        if first != feature.start:
            raise ValueError(
                "Start of first exon (%s) does not match start of feature (%s)"
                % (first, feature.start))
        if last != feature.stop:
            raise ValueError(
                "End of last exon (%s) does not match end of feature (%s)"
                % (last, feature.stop))

        if color is None:
            color = '0,0,0'
        color = color.replace(' ', '').strip()

        # Use field names as defined at
        # http://genome.ucsc.edu/FAQ/FAQformat.html#format1
        chrom = feature.chrom
        chromStart = feature.start - 1
        chromEnd = feature.stop
        orig = constants.always_return_list
        constants.always_return_list = True
        try:
            name = feature[name_field][0]
        except KeyError:
            name = "."

        constants.always_return_list = orig

        score = feature.score
        if score == '.':
            score = '0'
        strand = feature.strand
        itemRgb = color
        blockCount = len(exons)
        blockSizes = [len(i) for i in exons]
        blockStarts = [i.start - 1 - chromStart for i in exons]

        if thick_featuretype:
            thick = list(self.children(feature, featuretype=thick_featuretype,
                                       order_by='start'))
            if len(thick) == 0:
                thickStart = feature.start
                thickEnd = feature.stop
            else:
                thickStart = thick[0].start - 1  # BED 0-based coords
                thickEnd = thick[-1].stop

        if thin_featuretype:
            thin = list(self.children(feature, featuretype=thin_featuretype,
                                      order_by='start'))
            if len(thin) == 0:
                thickStart = feature.start
                thickEnd = feature.stop
            else:
                thickStart = thin[0].stop
                thickEnd = thin[-1].start - 1  # BED 0-based coords

        tst = chromStart + blockStarts[-1] + blockSizes[-1]
        assert tst == chromEnd, "tst=%s; chromEnd=%s" % (tst, chromEnd)

        fields = [
            chrom,
            chromStart,
            chromEnd,
            name,
            score,
            strand,
            thickStart,
            thickEnd,
            itemRgb,
            blockCount,
            ','.join(map(str, blockSizes)),
            ','.join(map(str, blockStarts))]
        return '\t'.join(map(str, fields))

    # Recycle the docs for _relation so they stay consistent between parents()
    # and children()
    children.__doc__ = children.__doc__.format(
        _relation_docstring=_relation.__doc__)
    parents.__doc__ = parents.__doc__.format(
        _relation_docstring=_relation.__doc__)

    # Add the docs for methods that call helpers.make_query()
    for method in [parents, children, features_of_type, all_features]:
        method.__doc__ = method.__doc__.format(
            _method_doc=_method_doc)
