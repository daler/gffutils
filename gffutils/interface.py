import sqlite3

import create
import bins
import helpers
import constants
from feature import Feature, feature_from_line


class FeatureDB(object):
    # docstring to be filled in for methods that call out to
    # helpers.make_query()
    _method_doc = """
        `limit`: string or tuple
            Limit the results to a genomic region.  If string, then of the form
            "seqid:start-end"; if tuple, then (seqid, start, end).

        `strand`: "-" | "+" | "."
            Limit the results to one strand

        `featuretype`: string or tuple
            Limit the results to one or several featuretypes.

        `order_by`: string or tuple
            Order results by one or many fields; the string or tuple items must
            be in: 'seqid', 'source', 'featuretype', 'start', 'end', 'score',
            'strand', 'frame', 'attributes', 'extra'.

        `reverse`: bool
            Change sort order; only relevant if `order_by` is not None.  By
            default, results will be in ascending order, so use `reverse=True`
            for descending.

        `completely_within': bool
            If False (default), a single bp overlap with `limit` is sufficient
            to return a feature; if True, then the feature must be completely
            within `limit`. Only relevant when `limit` is not None.
    """

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
            self._DBCreator_instance = dbfn
            if dbfn.dbfn == ':memory:':
                self.conn = dbfn.conn
            else:
                self.dbfn = dbfn.dbfn
                self.conn = sqlite3.connect(self.dbfn)
        else:
            self.dbfn = dbfn
            self.conn = sqlite3.connect(self.dbfn)

        self.conn.text_factory = str
        self.conn.row_factory = sqlite3.Row

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
        if isinstance(key, Feature):
            key = key.id
        c = self.conn.cursor()
        c.execute(constants._SELECT + ' WHERE id = ?', (key,))
        results = c.fetchall()
        if len(results) != 1:
            raise helpers.FeatureNotFoundError(key)
        return Feature(dialect=self.dialect, **results[0])

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
            yield Feature(dialect=self.dialect, **i)

    def iter_by_parent_childs(self, featuretype="gene"):
        """
        Iterate through GFF database by parent-child units. Return
        a generator where each item is a feature of 'featuretype'
        and all of its children. For example, when featuretype is gene,
        this returns a generator where each item is a list of records
        belonging to a gene (where the first record is the 'gene' record)
        itself.

        TODO: Maybe add option to limit this by depth?

        Alternative names:
          iter_by_parent_children
          iter_by_parent_childs
          iter_by_parent_unit
          iter_by_unit
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
        Returns an iterator of all :class:`Feature` objects in the database.

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
            yield Feature(dialect=self.dialect, **i)

    def featuretypes(self):
        """
        Iterator of feature types found in the database.
        """
        c = self.conn.cursor()
        c.execute(
            '''
            SELECT DISTINCT featuretype from features
            ''')
        for i, in c:
            yield i

    def _relation(self, id, join_on, join_to, level=None, featuretype=None,
                  order_by=None, reverse=False, completely_within=False):
        """
        Parameters
        ----------

        `id` : string or a Feature object

        `level` : None or int

            If `level=None` (default), then return all children regardless
            of level.  If `level` is an integer, then constrain to just that
            level.
        {_method_doc}
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
            completely_within=completely_within,
        )

        # modify _SELECT so that only unique results are returned
        query = query.replace("SELECT", "SELECT DISTINCT")
        for i in self._execute(query, args):
            yield Feature(dialect=self.dialect, **i)

    def children(self, id, level=None, featuretype=None, order_by=None,
                 reverse=False, completely_within=False):
        """
        Return children of feature `id`.

        {_relation_docstring}
        """
        return self._relation(
            id, join_on='child', join_to='parent', level=level,
            featuretype=featuretype, order_by=order_by, reverse=reverse,
            completely_within=completely_within)

    def parents(self, id, level=None, featuretype=None, order_by=None,
                reverse=False, completely_within=False):
        """
        Return parents of feature `id`.

        {_relation_docstring}
        """
        return self._relation(
            id, join_on='parent', join_to='child', level=level,
            featuretype=featuretype, order_by=order_by, reverse=reverse,
            completely_within=completely_within)

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

        elif isinstance(region, Feature):
            seqid = region.seqid
            start = region.start
            end = region.end
        else:
            seqid, start, end = region

        # Get a list of all possible bins for this region
        _bins = list(bins.bins(int(start), int(end), one=False))

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
            yield Feature(dialect=self.dialect, **i)

    @classmethod
    def interfeatures(self, features, new_featuretype=None,
                      merge_attributes=True, dialect=None):
        """
        Construct new features representing the space between features.

        For example, if `features` is a list of exons, then this method will
        return the introns.  If `features` is a list of genes, then this method
        will return the intergenic regions.

        Providing N features will return N - 1 new features.

        This method purposefully does *not* do any merging or sorting, so you
        may want to use :meth:`FeatureDB.merge` first.

        The new features' attributes will be a merge of the neighboring
        features' attributes.  This is useful if you have provided a list of
        exons; the introns will then retain the transcript and/or gene parents.

        Parameters
        ----------
        `features` : iterable of :class:`feature.Feature` instances
            Sorted, merged iterable

        `new_featuretype` : string or None
            The new features will all be of this type, or, if None (default)
            then the featuretypes will be constructed from the neighboring
            features, e.g., `inter_exon_exon`.

        `attribute_func` : callable or None
            If None, then nothing special is done to the attributes.  If
            callable, then the callable accepts two attribute dictionaries and
            returns a single attribute dictionary.  If `merge_attributes` is
            True, then `attribute_func` is called before `merge_attributes`.
            This could be useful for manually managing IDs for the new
            features.

        >>> features = [
        ... "chr1 . exon 1   100 . + . ID=exon1; Parent=mRNA1",
        ... "chr1 . exon 200 250 . + . ID=exon2; Parent=mRNA1",
        ... "chr1 . exon 500 600 . + . ID=exon3; Parent=mRNA1",
        ... ]
        >>> features = [feature_from_line(i) for i in features]
        >>> for i in FeatureDB.interfeatures(
        ... features, new_featuretype="intron"):
        ...     print i  # doctest: +NORMALIZE_WHITESPACE
        chr1 gffutils_derived intron 101 199 . + . ID=exon1,exon2;Parent=mRNA1
        chr1 gffutils_derived intron 251 499 . + . ID=exon1,exon3;Parent=mRNA1

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
            assert last_feature.strand == f.strand
            assert last_feature.chrom == f.chrom
            strand = last_feature.strand
            chrom = last_feature.chrom

            # Shrink
            interfeature_start += 1
            interfeature_stop -= 1

            new_attributes = helpers.merge_attributes(
                last_feature.attributes, f.attributes)

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
            yield Feature(dialect=dialect, **fields)
            interfeature_start = f.stop

    def update(self, features, merge_strategy=None, transform=None,
               id_spec=None, verbose=False):
        if self._DBCreator_instance is not None:
            db = self._DBCreator_instance

        else:
            # TODO: Currently create._DBCreator only accepts a filename, which
            # it creates a parser.Parser.  It should accept an iterable of
            # Features, but this will take some refactoring.
            #
            # So for now, create an intermediate tempfile.
            tmp = tempfile.NamedTemporaryFile()
            for f in features:
                tmp.write(str(f) + '\n')
            tmp.seek(0)
            if self.dialect['fmt'] == 'gtf':
                db = create._GTFCreator(fn=tmp.name, dbfn=self.dbfn,
                                        verbose=verbose)
            elif self.dialect['fmt'] == 'gff3':
                db = create._GFFDBCreator(fn=tmp.name, dbfn=self.dbfn,
                                          verbose=verbose)

        if merge_strategy:
            db.merge_strategy = merge_strategy
        if id_spec:
            db.id_spec = id_spec

        db.set_verbose(False)
        db._populate_from_lines(features)
        db._update_relations()
        db._finalize()

    def create_introns(self, exon_featuretype='exon',
                       grandparent_featuretype='gene', parent_featuretype=None,
                       new_featuretype='intron', merge_attributes=True):
        """
        Create introns from existing annotations.


        Parameters
        ----------
        `exon_featuretype`: string
            Feature type to use in order to infer introns.  Typically `"exon"`.

        `grandparent_featuretype` : string
            If `grandparent_featuretype` is not None, then group exons by
            children of this featuretype.  If `granparent_featuretype` is
            "gene" (default), then all genes will be retrieved, then *all*
            children of each gene (e.g., mRNA, rRNA, ncRNA, etc) will then be
            searched for exons from which to infer introns. Mutually exclusive
            with `parent_featuretype`.

        `parent_featuretype` : string
            If `parent_featuretype` is not None, then only use this featuretype
            to infer introns.  Use this if you only want a subset of
            featuretypes to have introns (e.g., "mRNA" only, and not ncRNA or
            rRNA). Mutually exclusive with `grandparent_featuretype`.

        `new_featuretype`: string
            Feature type to use for the inferred introns; default is
            `"intron"`.

        `merge_attributes`: bool
            Whether or not to merge attributes from all exons; if False then no
            attributes will be created for the introns.
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

    def merge(self, feature, ignore_strand=False):

        # Consume iterator up front...
        features = list(features)

        # Either set all strands to '+' or check for strand-consistency.
        if ignore_strand:
            strand = '+'
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
                merged_feature = Feature(chrom=feature.chrom,
                                         source='.',
                                         featuretype=featuretype,
                                         start=current_merged_start,
                                         stop=current_merged_stop,
                                         score='.',
                                         strand=strand,
                                         frame='.',
                                         attributes='')
                yield merged_feature

                # and we start a new one, initializing with this feature's
                # start and stop.
                current_merged_start = feature.start
                current_merged_stop = feature.stop

        # need to yield the last one.
        if len(features) == 1:
            feature = features[0]
        merged_feature = Feature(chrom=feature.chrom,
                                 source='.',
                                 featuretype=featuretype,
                                 start=current_merged_start,
                                 stop=current_merged_stop,
                                 score='.',
                                 strand=strand,
                                 frame='.',
                                 attributes='')
        yield merged_feature

    def children_bp(self, feature, child_featuretype='exon', merge=False):
        """
        Returns the total bp of all children of a featuretype

        Useful for getting the exonic bp of an mRNA.
        """
        children = self.children(feature, featuretype=child_featuretype,
                                 order_by='start')
        if merge:
            children = self.merge(children)

        total = 0
        for child in children:
            total += len(child)
        return total

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
