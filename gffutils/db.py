import sqlite3
import tempfile
import os
import time
import sys
from gfffeature import GFFFile, Feature
from helpers import FeatureNotFoundError, asinterval


class DBCreator(object):
    def __init__(self, fn, dbfn, force=False, verbose=True, **kwargs):
        """
        Generic feature database constructor.  Subclasses must override the
        populate_from_features method and update_relations method.
        """
        if force:
            if os.path.exists(dbfn):
                os.unlink(dbfn)
        self.dbfn = dbfn
        conn = sqlite3.connect(dbfn)
        conn.text_factory = sqlite3.OptimizedUnicode
        self.conn = conn
        self.fn = fn

        gfffile = GFFFile(fn, **kwargs)
        self.nfeatures = len(gfffile)
        self.features = gfffile
        self.verbose = verbose

    def populate_from_features(self, features):
        pass

    def update_relations(self):
        pass

    def drop_indexes(self):
        c = self.conn.cursor()
        c.execute('DROP INDEX IF EXISTS ids')
        c.execute('DROP INDEX IF EXISTS parentindex')
        c.execute('DROP INDEX IF EXISTS childindex')
        self.conn.commit()

    def init_tables(self):
        c = self.conn.cursor()
        c.executescript('''
        CREATE TABLE features (
                                id text,
                                chrom text,
                                start int,
                                stop int,
                                strand text,
                                featuretype text,
                                score float,
                                source text,
                                frame text,
                                attributes text,
                                primary key (id)
                              );
        CREATE TABLE relations (
            parent text,
            child text,
            level int,
            primary key(parent,child,level) );

        CREATE TABLE meta (
            filetype text
            );
        ''')
        self.conn.commit()

    def set_filetype(self):
        c = self.conn.cursor()
        c.execute('''
        INSERT INTO meta VALUES (?)''', (self.filetype,))
        self.conn.commit()

    def create(self):
        self.init_tables()
        self.populate_from_features(self.features)
        self.update_relations()
        self.set_filetype()


class GFFDBCreator(DBCreator):
    def __init__(self, *args, **kwargs):
        self.filetype = 'gff'
        DBCreator.__init__(self, *args, **kwargs)

    def populate_from_features(self, features):
        t0 = time.time()
        c = self.conn.cursor()
        self.drop_indexes()
        nfeatures = float(self.nfeatures)
        last_perc = 0
        for i, feature in enumerate(features):
            if self.verbose:
                perc = int(i / nfeatures * 100)
                if perc != last_perc:
                    msg = 'Populating features table and first-order '
                    msg += 'relations: %d (%d%%)' % (i, perc)

                    sys.stderr.write('\r' + msg)
                    sys.stderr.flush()
                last_perc = perc

            if not feature.id:
                new_id = '%s:%s:%s-%s:%s' % (feature.featuretype,
                        feature.chrom, feature.start, feature.stop,
                        feature.strand)
                feature.id = new_id
            c.execute('''
                      INSERT OR IGNORE INTO features VALUES
                      (?,?,?,?,?,?,?,?,?,?)
                      ''', (
                      feature.id,
                      feature.chrom,
                      feature.start,
                      feature.stop,
                      feature.strand,
                      feature.featuretype,
                      feature.score,
                      feature.source,
                      feature.frame,
                      feature._str_attributes))

            # If this feature has a parent, we have a relation to add to the
            # relations table.
            try:
                parents = feature.attributes['Parent']
                if isinstance(parents, basestring):
                    parents = [parents]
                child = feature.id
                for parent in parents:
                    c.execute('INSERT INTO relations VALUES (?,?,?)',
                            (parent, child, 1))
            except KeyError:
                pass

        self.conn.commit()

        t1 = time.time()
        if self.verbose:
            sys.stderr.write('...done in %.1fs\n' % (t1 - t0))
            sys.stderr.flush()

    def update_relations(self):
        t0 = time.time()
        if self.verbose:
            sys.stderr.write('Updating 2nd-order relations...')
        c = self.conn.cursor()

        # create 2 more cursors so we can iterate over one while querying on
        # the other's iteration
        c2 = self.conn.cursor()
        c3 = self.conn.cursor()

        # First index for fast lookup...
        c.execute('CREATE INDEX ids ON features (id)')
        c.execute('CREATE INDEX parentindex ON relations (parent)')
        c.execute('CREATE INDEX childindex ON relations (child)')
        self.conn.commit()

        c.execute('SELECT id FROM features')

        # Strategy is to get each ID in the db, then see if it has a listed
        # child (e.g., a transcript). If so, then see if that child also has a
        # child (e.g., an exon).  Then write the relation of parent gene to
        # grandchild exon to tempfile for later import.
        tmp = tempfile.mktemp()
        fout = open(tmp, 'w')
        for parent in c:
            parent = parent[0]

            # Here we get the first-level child from the initial import.   This
            # data was contained in the "Parent=" attribute of each GFF
            # feature.
            c2.execute('''
            SELECT child FROM relations
            WHERE parent = ? AND level=1''', (parent,))

            # For each of those children, we get *their* children -- so the
            # 'grandchildren' of the original parent.
            for child in c2:
                child = child[0]
                c3.execute('''
                SELECT child FROM relations
                WHERE parent = ? AND level=1''', (child,))
                for grandchild in c3:
                    grandchild = grandchild[0]
                    fout.write('%s\t%s\n' % (parent, grandchild))
        fout.close()

        # Import the "grandchild" file into the relations table
        for line in open(tmp):
            parent, child = line.strip().split('\t')
            c.execute('''
            INSERT OR IGNORE INTO relations
            VALUES (?,?,?)''', (parent, child, 2))

        t1 = time.time()
        if self.verbose:
            sys.stderr.write('done in %.1fs\n' % (t1 - t0))
            sys.stderr.flush()

        # Re-create indexes.
        # TODO: performance testing to see how much of a difference these make
        if self.verbose:
            sys.stderr.write('Creating indexes...')
            sys.stderr.flush()
        c.execute('drop index childindex')
        c.execute('drop index parentindex')
        c.execute('create index parentindex on relations (parent)')
        c.execute('create index childindex on relations (child)')
        c.execute('create index starts on features(start)')
        c.execute('create index stops on features(stop)')
        c.execute('create index startstrand on features(start, strand)')
        c.execute('create index stopstrand on features(stop,strand)')
        c.execute('create index featuretypes on features(featuretype)')

        self.conn.commit()
        os.unlink(tmp)
        t2 = time.time()
        if self.verbose:
            sys.stderr.write('done in %.1fs\n' % (t2 - t1))
            sys.stderr.flush()


class GTFDBCreator(DBCreator):
    def __init__(self, fn, dbfn, *args, **kwargs):
        DBCreator.__init__(self, fn, dbfn, *args, **kwargs)
        self.filetype = 'gtf'

    def populate_from_features(self, features):
        t0 = time.time()
        self.drop_indexes()
        c = self.conn.cursor()
        for feature in features:
            parent = feature.attributes['transcript_id']
            grandparent = feature.attributes['gene_id']

            # A database-specific ID to use
            ID = '%s:%s:%s-%s:%s' % (
                    feature.featuretype, feature.chrom, feature.start,
                    feature.stop, feature.strand)

            # If it's an exon, its attributes include its parent transcript
            # and its 'grandparent' gene.  So we can insert these
            # relationships into the relations table now.

            # Note that the table schema has (parent,child) as a primary
            # key, so the INSERT OR IGNORE won't add multiple entries for a
            # single (parent,child) relationship

            # The gene has a grandchild exon
            c.execute('''
            INSERT OR IGNORE INTO relations VALUES (?,?,?)
            ''', (grandparent, ID, 2))

            # The transcript has a child exon
            c.execute('''
            INSERT OR IGNORE INTO relations VALUES (?,?,?)
            ''', (parent, ID, 1))

            # The gene has a child transcript
            c.execute('''
            INSERT OR IGNORE INTO relations VALUES (?,?,?)
            ''', (grandparent, parent, 1))

            # Insert the feature into the features table.
            c.execute('''
                      INSERT OR IGNORE INTO features
                      VALUES (?,?,?,?,?,?,?,?,?,?)
                      ''', (ID,
                       feature.chrom,
                       feature.start,
                       feature.stop,
                       feature.strand,
                       feature.featuretype,
                       feature.score,
                       feature.source,
                       feature.frame,
                       feature._str_attributes))
        self.conn.commit()

    def update_relations(self):
        self.drop_indexes()
        c = self.conn.cursor()
        c2 = self.conn.cursor()
        c.execute('CREATE INDEX ids ON features (id)')
        c.execute('CREATE INDEX parentindex ON relations (parent)')
        c.execute('CREATE INDEX childindex ON relations (child)')

        tmp = tempfile.mktemp()
        fout = open(tmp, 'w')
        c.execute("SELECT DISTINCT parent FROM relations")
        for parent in c:
            parent = parent[0]
            c2.execute("""
                       SELECT min(start), max(stop), level, strand, chrom
                       FROM features
                       JOIN relations ON
                       features.ID = relations.child
                       WHERE
                       parent = ?
                       AND featuretype == "exon"
                       """, (parent,))
            child_limits = c2.fetchone()
            start, end, level, strand, chrom = child_limits

            # The strategy here is to write parents to file, and later read the
            # file back into the database so that we don't keep writing to the
            # database while we're in the middle of consuming a query results
            # iterator...

            # Using some assumptions we can make some shortcuts to determine
            # what featuretype this parent really is.  Since the features table
            # only contains exons, and we're joining on the features table,
            # only exon children or only grandchildren will be returned (this
            # can be verified by the commented-out test code above).  If only
            # grandchildren were returned (i.e., level=2) then the parent is a
            # gene.

            # In addition, we need to create a fake attributes string so that
            # later, upon accessing the database, the GTFDB wrapper will know
            # how to assign an ID.  This is sort of a hack in order to maintain
            # cleaner class inheritance between GFFFeatures and GTFFeatures.

            # Since the relations table only has transcript children in it, any
            # parents at level 1 are mRNAs.
            #
            # WARNING: this does NOT account for non-coding RNAs -- they will
            # still be called "mRNA"

            if level == 1:
                featuretype = 'mRNA'
                attributes = 'transcript_id "%s"; ' % parent

            # On the other hand, level 2 parents means that the parent is a
            # gene.
            if level == 2:
                featuretype = 'gene'
                attributes = 'gene_id "%s"; ' % parent

            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                parent, chrom, start, end, strand, featuretype, attributes))
        fout.close()

        # Now that the file has been created and we're done querying for the
        # parents, slurp in the file and insert everything into the db.
        fin = open(fout.name)
        for line in fin:
            score = '.'
            source = 'derived'
            frame = '.'
            c.execute("""
                      INSERT OR IGNORE INTO features (id, chrom, start, stop,
                      strand, featuretype, attributes, score, source, frame)
                      VALUES (?,?,?,?,?,?,?,?,?,?)
                      """, line.strip().split('\t') + [score, source, frame])

        self.conn.commit()
        os.remove(tmp)

        c.execute('DROP INDEX ids')
        c.execute('CREATE INDEX ids ON features (id)')
        self.conn.commit()


class FeatureDB:
    def __init__(self, db_fn):
        """
        Class to interface with a database created by the create_gffdb()
        function.

        *db_fn* is the filename of the database to use.  See the help for
        create_gffdb() for much more information on the creation and design
        ideas behind the database.
        """
        self.db_fn = db_fn
        self.conn = sqlite3.connect(db_fn)
        self.conn.text_factory = str

        # SELECT statement to use when you want all the fields of a feature,
        # plus the db ID.  Results from this are usually sent right to
        # self._newfeature()
        self.SELECT = """
        SELECT id, chrom, source, featuretype, start, stop, score, strand,
        frame, attributes
        """

        c = self.conn.cursor()
        try:
            c.execute('''
            SELECT filetype FROM meta
            ''')
            self.filetype = c.fetchone()[0]
        except sqlite3.OperationalError:
            c.execute(self.SELECT + ' FROM features LIMIT 1')
            results = c.fetchone()[0]
            feature = self._newfeature(results)
            self.filetype = feature.filetype

    def _newfeature(self, *args):
        feature = Feature(*args[1:])
        feature.dbid = args[0]
        return feature

    def __getitem__(self, id):
        c = self.conn.cursor()
        if isinstance(id, Feature):
            id = id.id
        c.execute('''
                  %s from features where id = ?
                  ''' % self.SELECT, (id,))
        results = c.fetchall()
        if not len(results) == 1:
            raise FeatureNotFoundError(id)
        return self._newfeature(*results[0])

    def attribute_search(self, text, featuretype='gene'):
        """
        Looks for *text* within the "attributes" field of features of
        *featuretype* ('gene' by default).  Useful for looking up genes based
        on symbol or name rather than accession number.  Returns a list of
        matching Features.

        Uses SQL's LIKE operator, which is case-insensitive.
        """
        text = '%' + text + '%'
        c = self.conn.cursor()
        c.execute('''
                  %s from features where attributes
                  like "%s" and featuretype = ?
                  ''' % (self.SELECT, text,), (featuretype,))
        results = []
        for result in c.fetchall():
            results.append(self._newfeature(*result))
        return results

    def features_of_type(self, featuretype, chrom=None, start=None, stop=None,
            strand=None):
        """
        Returns an iterator of GFFFeature objects that have the feature type
        *featuretype*.  You can optionally specify chrom, strand, start, and
        stop filters.

        For example, to get all exons within a range::

            self.features_of_type('exon', chrom='chr3L', start=11100,
                                  stop=12200)

        Any of the filters can be empty::

            # exons on any chromosome after the first 1kb
            self.features_of_type('exon',start=1000)

            # genes on + strand of chrX
            self.features_of_type('gene',chrom='chrX', strand='+')

        For more complicated stuff, (e.g., filtering by multiple chromosomes)
        you'll probably want to use self.execute() directly and write your own
        SQL queries . . . or just call this method iteratively.
        """
        filter_clause = ''
        if chrom is not None:
            filter_clause += ' AND chrom = "%s"' % chrom

        if start is not None:
            filter_clause += ' AND start >= %s' % start
        if stop is not None:
            filter_clause += ' AND stop <= %s' % stop
        if strand is not None:
            filter_clause += ' AND strand = "%s"' % strand

        c = self.conn.cursor()
        c.execute('''
                  %s
                  FROM
                  features
                  WHERE featuretype = ?
                  %s
                  ORDER BY start
                  ''' % (self.SELECT, filter_clause), (featuretype, ))
        for i in c:
            yield self[i[0]]

    def featuretypes(self):
        """
        Returns an iterator of the different feature types found in the
        database.
        """
        c = self.conn.cursor()
        c.execute('''
                  SELECT DISTINCT featuretype from features
                  ''')
        for i in c:
            yield i[0]

    def all(self):
        """
        Returns an iterator of ALL features in the database.
        """
        c = self.conn.cursor()
        c.execute('''
                  %s
                  FROM
                  features
                  ''' % self.SELECT)
        for i in c:
            yield self._newfeature(*i)

    def exonic_bp(self, id, ignore_strand=False):
        """
        Merges all exons of the gene *id* and sums the total bp.
        *ignore_strand* is passed to merge_features(). Useful for calculating
        RPKM for full gene models.
        """
        if isinstance(id, Feature):
            id = id.id
        exons = self.children(id, featuretype='exon', level=2)
        exons = list(exons)
        if len(exons) == 0:
            return 0
        merged_exons = self.merge_features(exons, ignore_strand)
        total_exon_bp = 0
        for exon in merged_exons:
            total_exon_bp += len(exon)
        return total_exon_bp

    def closest_feature(self, chrom, pos, strand=None, featuretype=None,
            ignore=None, end=None):
        """
        Returns the closest feature and the distance to the start position of
        that feature from *pos*.  It's up the caller to determine things like
        upstream/downstream, TSS, TTS, midpoints, etc.

        *ignore* is either None (default) or a list of IDs that you would like
        to ignore.

        If you are providing the *chrom*, *pos* of a gene's TSS, then you will
        have to add that gene's ID to the *ignore* list so you don't get that
        same gene returned as the closest (which it is).

        In this case, you'll probably want to provide *featuretype* too (e.g.,
        featuretype='gene') -- otherwise you'll also have to add the gene's
        mRNAs, exons, introns, CDSs, and any other related features to the
        *ignore* list.

        The *ignore* list also makes it possible to get the second-closest
        feature to a coord by calling this method to get the closest feature,
        then adding that closest feature's ID to the *ignore* list and calling
        this method a second time.

        Note that providing a *featuretype* and *strand* will speed things up
        considerably.

        *end* refers to which end of features you want to look for; this can be
        "start" or "stop" and its meaning depends on strand.

        For example, the kwargs  {strand:"-", end:"stop", featuretype:"gene"}
        will look for the closest minus-strand TSS; use {strand:"+",
        end:"start", featuretype:"gene"} for plus-strand closest TSS.

        For upstream/downstream work, the idiom is to use a while loop,
        ignoring consecutively close features that are close but in the wrong
        direction.

        For example, to get the closest upstream TSS that's on the plus strand
        to chr3R:50000::

            >>> ignore = []
            >>> while (closest_pos > pos):
            ...     dist,closest_plus_feature =
            G.closest_feature(chrom='chr3R', pos=50000, strand="+",
            end="start", ignore=ignore, featuretype='gene')
            ...     closest_pos = closest_plus_feature.TSS
            ...     ignore.append(closest_plus_feature.id)

        And here's the closest upstream TSS on the minus strand::

            >>> ignore = []
            >>> while (closest_pos > pos):
            ...     dist,closest_plus_feature =
            G.closest_feature(chrom='chr3R', pos=50000, strand="-", end="end",
            ignore=ignore, featuretype='gene')
            ...     closest_pos = closest_plus_feature.TSS
            ...     ignore.append(closest_plus_feature.id)

        """

        # e.g., AND id != FBgn0001 AND id != FBgn0002
        ignore_clause = ''
        if ignore is not None:
            if type(ignore) is str:
                ignore = [ignore]
            for i in ignore:
                ignore_clause += ' AND id != "%s" ' % i

        strand_clause = ''
        if strand is not None:
            strand_clause = ' AND strand="%s" ' % strand

        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = ' AND featuretype="%s" ' % featuretype

        c = self.conn.cursor()

        candidates = []
        if end is None:
            end = ['start', 'stop']

        if 'start' in end:
            c.execute('''
            SELECT abs(%(pos)s-start) as AAA,id FROM features
            WHERE
            chrom = ?
            %(ignore_clause)s
            %(strand_clause)s
            %(featuretype_clause)s
            ORDER BY AAA LIMIT 1
            ''' % locals(), (chrom,))
            candidates.append(c.fetchone())

        if 'stop' in end:
            c.execute('''
            SELECT abs(%(pos)s-stop) as AAA,id FROM features
            WHERE
            chrom = ?
            %(ignore_clause)s
            %(strand_clause)s
            %(featuretype_clause)s
            ORDER BY AAA LIMIT 1
            ''' % locals(), (chrom,))
            candidates.append(c.fetchone())

        candidates = [i for i in candidates if i is not None]
        if len(candidates) == 0:
            return None, None

        candidates.sort(key=lambda x: x[0])
        dist, feature_id = candidates[0]
        return dist, self[feature_id]

    def overlapping_features(self, chrom, start, stop, featuretype=None,
            strand=None, completely_within=False):
        """
        Returns an iterator of features of type *featuretype* that overlap the
        coords provided. If *featuretype* is None (default), then all features
        will be returned.

        If *strand* is not None, then only features of the strand specifed will
        be returned.

        If *completely_within* is True, a feature needs to be completely within
        the boundaries of *id*.  If False (default), features that extend out
        on either side will be included as well.
        """

        if strand is None:
            strand_clause = ''
        else:
            strand_clause = ' AND strand = "%s"' % strand

        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = ' AND featuretype = "%s"' % featuretype

        if completely_within:
            within_clause = '''
            AND ((start BETWEEN %(start)s AND %(stop)s)
            AND (stop BETWEEN %(start)s AND %(stop)s))''' % locals()
        else:
            within_clause = ''' AND (
                                        (start BETWEEN %(start)s AND %(stop)s)
                                        OR
                                        (stop BETWEEN %(start)s AND %(stop)s)
                                        OR
                                        (start < %(stop)s AND stop > %(start)s)
                                     ) ''' % locals()
        c = self.conn.cursor()
        c.execute('''
        %s
        FROM features WHERE
        chrom = ?
        %s
        %s
        %s
        ''' % (self.SELECT, strand_clause, featuretype_clause, within_clause),
        (chrom,))
        for i in c:
            yield self._newfeature(*i)

    def merge_features(self, features, ignore_strand=False):
        """
        Returns an iterator of merged features (overlapping features are
        combined into one, spanning from the lowest start to the highest stop
        position)

        *features* is an iterable of features that you want to merge together.

        *features* will be converted into a list because this method needs to
        sort the features.

        If *ignore_strand* is True, strand will be forced to all '+' and no
        checking will be done on the input

        If *ignore_strand* is False (default), then *features* must be all on
        the same strand.

        The new *featuretype* will be a joining of the various featuretypes
        that were provided.  If all the features were 'exon', then the new,
        merged feature will have a featuretype of 'merged_exon'.  If there were
        some introns in *features*, then the new merged feature will have a
        featuretype of 'merged_exon_intron'.

        Note that merged features are not saved to the database -- they are
        created on-the-fly and only exist in memory.
        """
        # If it quacks like an iterator...then turn it into a list.
        if hasattr(features, 'next'):
            features = list(features)

        # Quick function to sort by start position
        def start_pos(x):
            return x.start

        # Sort 'em by start position
        features.sort(key=start_pos)

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

        # If they were all exons, merged objects will have featuretype
        # 'merged_exon'.  If features included both exon and CDS, will have
        # featuretype 'merged_exon_CDS'.
        featuretypes = list(set([i.featuretype for i in features]))
        featuretypes.sort()
        featuretypes = map(str, featuretypes)
        featuretypes = '_'.join(featuretypes)
        featuretype = 'merged_%s' % featuretypes

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

    def interfeatures(self, features):
        """
        Given an iterable of *features*, returns an iterator of new features
        defined by the intervening space between each consecutive feature and
        the one before it.

        If you provide N features, you'll get back N-1 intervening features.

        So if *features* contains the exons of a transcript, this method will
        return the introns.

        This is a purposefully naive method that does NOT do sorting or merging
        -- it simply returns the intervening space between the features
        provided.

        If the features returned are not sorted, you may get overlapping
        results.

        Example for getting all the 'non-exonic' space on the plus strand of
        chrX::

            >>> # merge_features needs a single chrom and single strand
            >>> exons = G.features_of_type('exon',chrom='chrX',strand='+')
            >>> merged_exons = G.merge_features(exons)
            >>> non_exonic = G.interfeatures(merged_exons)

        Example for getting the introns of a transcript::

            >>> transcript = G.random_feature('mRNA')
            >>> exons = G.children(transcript,'exon')
            >>> introns = G.interfeatures(exons)
        """
        for i, feature in enumerate(features):
            if i == 0:
                interfeature_start = feature.stop
                last_feature = feature
                continue
            interfeature_stop = feature.start
            featuretype = 'inter_%s_%s' % (
                    last_feature.featuretype, feature.featuretype)
            assert last_feature.strand == feature.strand
            assert last_feature.chrom == feature.chrom
            strand = last_feature.strand
            chrom = last_feature.chrom

            # Shrink
            interfeature_start += 1
            interfeature_stop -= 1

            yield Feature(chrom=chrom,
                          source='.',
                          featuretype=featuretype,
                          start=interfeature_start,
                          stop=interfeature_stop,
                          score='.',
                          strand=strand,
                          frame='.',
                          attributes='')
            interfeature_start = feature.stop

    def chromosomes(self):
        """
        Returns a list of chromosomes in the database
        """
        c = self.conn.cursor()
        c.execute('''
        SELECT DISTINCT chrom FROM features
        ''')
        return [i[0] for i in c.fetchall()]

    def strands(self):
        """
        Returns a list of strands in the database
        """
        c = self.conn.cursor()
        c.execute('''
        SELECT DISTINCT strand FROM features
        ''')
        return [i[0] for i in c.fetchall()]

    def children(self, id, level=1, featuretype=None):
        """
        Returns an iterator of the children *level* levels away.  Immediate
        children are level=1; "grandchildren" are level=2.  The child
        transcripts of a gene will be at level=1; the child exons of the same
        gene are at level=2.  Returns an error if there are not enough levels
        of children for the level specified.
        """

        if level > 2:
            raise NotImplementedError('Levels > 2 not supported yet.')

        if isinstance(id, Feature):
            id = id.id

        cursor = self.conn.cursor()

        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = '''
            AND features.featuretype = "%s"''' % featuretype

        cursor.execute('''
        SELECT DISTINCT
        id, chrom, source, featuretype, start, stop, score, strand, frame,
        attributes
        FROM features JOIN relations
        ON relations.child = features.id
        WHERE relations.parent = ?
        AND relations.level = ?
        %s
        ORDER BY start''' % (
            featuretype_clause), (id, level))
        for i in cursor:
            yield self._newfeature(*i)

    def parents(self, id, level=1, featuretype=None):
        """
        Returns an iterator of the parents *level* levels away.  Immediate
        children are level=1; "grandchildren" are level=2.  The child
        transcripts of a gene will be at level=1; the child exons of the same
        gene are at level=2.  Returns an error if there are not enough levels
        of children for the level specified."""

        if level > 2:
            raise NotImplementedError('Levels > 2 not supported yet.')

        if isinstance(id, Feature):
            id = id.id

        cursor = self.conn.cursor()
        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = '''
            AND features.featuretype = "%s"''' % featuretype
        cursor.execute('''
            SELECT DISTINCT
            id, chrom, source, featuretype, start, stop, score, strand, frame,
            attributes
            FROM features JOIN relations
            ON relations.parent = features.id
            WHERE relations.child = ?
            AND relations.level = ?
            %s
            ORDER BY start''' % (featuretype_clause), (id, level))
        for i in cursor:
            yield self._newfeature(*i)

    def execute(self, query):
        """
        Returns an iterator of results from any arbitrary query passed in.

        No type conversion is done, so you won't get GFFFeature objects back,
        just tuples of strings.
        """
        cursor = self.conn.cursor()
        cursor.execute(query)
        for i in cursor:
            yield i

    def to_GTF(self, gene, include_utrs=True):
        """
        Converts the gene into a series of GTF features that can be output to a
        file.

        Uses the GTF2.2 spec from http://mblab.wustl.edu/GTF22.html

        * Featuretypes CDS, start_codon, stop_codon are required.
        * start codon included in first CDS
        * stop codon not included in last CDS; the last CDS must shrink by 3bp
        * non-contiguous start/stop codons indicated by nonzero frame field; if
          nonzero then there must be multiple features for the same codon
        * frame:
            * 0 = starts with whole codon at 5'-most base
            * 1 = one extra base (the last of a codon) before the first whole codon
            * 2 = two extra bases

        """
        # how we're sorting exons...
        def by_start(x):
            return x.start

        # first, get all the transcripts.
        for transcript in self.children(gene, level=1):
            gene_id = gene.id
            transcript_id = transcript.id

            # all exons of this transcript will have the same gene/transcript
            # parents, so just create the attributes once
            attributes = 'gene_id "%s"; transcript_id "%s";' % (
                    gene_id, transcript_id)

            # Sort exons separately from CDSs.
            exons = sorted(
                    self.children(
                        transcript, level=1, featuretype='exon'), key=by_start)
            cdss = sorted(
                    self.children(
                        transcript, level=1, featuretype='CDS'), key=by_start)

            # needs more testing before implementation
            utrs = []
            """
            if include_utrs:
                five_prime_utrs = sorted(
                        self.children(
                                      transcript,
                                      level=1, featuretype='five_prime_UTR'),
                            key=by_start)
                for five_prime_utr in five_prime_utrs:
                    five_prime_utr.featuretype = '5UTR'
                    utrs.append(five_prime_utr)

                # 3' UTRs have stop codon positions removed from them
                three_prime_utrs = sorted(
                        self.children(
                            transcript, level=1,
                            featuretype='three_prime_UTR'), key=by_start)
                for three_prime_utr in three_prime_utrs:
                    three_prime_utr.featuretype = '3UTR'
                    if three_prime_utr.strand == '-':
                        three_prime_utr.stop -= 3
                    if three_prime_utr.strand == '+':
                        three_prime_utr.start += 3
                    utrs.append(three_prime_utr)
            """
            # First the exons, which are easy:
            for feature in exons:
                yield Feature(chrom=feature.chrom,
                                 start=feature.start,
                                 stop=feature.stop,
                                 featuretype=feature.featuretype,
                                 score=feature.score,
                                 frame=feature.frame,
                                 strand=feature.strand,
                                 source=feature.source,
                                 attributes=attributes)

            # Construct new features for start_codon and stop_codon.
            # Contiguous start and stop codons always have frame=0, says the
            # GTF 2.2 spec.  TODO: how to identify from GFF when a codon is
            # split?

            # Only consider CDSs and start/stop if this is a coding gene.
            # Rather than rely on featuretype=='mRNA", here we look for the
            # presence of CDSs.
            if len(cdss) > 0:
                if transcript.strand == '+':
                    # start codon
                    yield Feature(chrom=transcript.chrom,
                                     start=cdss[0].start,
                                     stop=cdss[0].start + 2,
                                     strand=transcript.strand,
                                     featuretype='start_codon',
                                     source=transcript.source,
                                     score='.',
                                     frame='0',
                                     attributes=attributes)

                    # stop codon
                    yield Feature(chrom=transcript.chrom,
                                     start=cdss[-1].stop - 2,
                                     stop=cdss[-1].stop,
                                     strand=transcript.strand,
                                     featuretype='stop_codon',
                                     source=transcript.source,
                                     score='.',
                                     frame='0',
                                     attributes=attributes)

                    # write out the first CDS, but only if a multi-CDS gene -- since
                    # single-CDS genes need to have stop_codon coords trimmed
                    nextframe = 0
                    if len(cdss) > 1:
                        yield Feature(chrom=transcript.chrom,
                                      source=transcript.source,
                                      featuretype='CDS',
                                      strand=transcript.strand,
                                      start=cdss[0].start,
                                      stop=cdss[0].stop,
                                      score='.',
                                      frame=str(nextframe),
                                      attributes=attributes)
                        nextframe = (3 - ((len(cdss[0]) - nextframe) % 3)) % 3

                    for cds in cdss[1:-1]:
                        yield Feature(chrom=transcript.chrom,
                                      source=transcript.source,
                                      strand=transcript.strand,
                                      featuretype='CDS',
                                      start=cds.start,
                                      stop=cds.stop,
                                      score='.',
                                      frame=str(nextframe),
                                      attributes=attributes)
                        nextframe = (3 - ((len(cds) - nextframe) % 3)) % 3

                    # last CDS has stop_codon's coords removed
                    yield Feature(chrom=transcript.chrom,
                                  source=transcript.source,
                                  strand=transcript.strand,
                                  featuretype='CDS',
                                  start=cdss[-1].start,
                                  stop=cdss[-1].stop - 3,
                                  score='.',
                                  frame=str(nextframe),
                                  attributes=attributes)

                # Or, if we're on the other strand,
                if transcript.strand == '-':

                    # minus-strand has start codon included in CDS coords
                    yield Feature(chrom=transcript.chrom,
                                     start=cdss[-1].stop - 2,
                                     stop=cdss[-1].stop,
                                     strand=transcript.strand,
                                     featuretype='start_codon',
                                     source=transcript.source,
                                     score='.',
                                     frame='0',
                                     attributes=attributes)

                    # minus-strand stop codon needs to be removed from cdss[0]
                    yield Feature(chrom=transcript.chrom,
                                     start=cdss[0].start,
                                     stop=cdss[0].start + 2,
                                     strand=transcript.strand,
                                     featuretype='stop_codon',
                                     source=transcript.source,
                                     score='.',
                                     frame='0',
                                     attributes=attributes)

                    # first CDS is last sorted; only deal with it here a multi-CDS gene
                    nextframe = 0
                    if len(cdss) > 1:
                        yield Feature(chrom=transcript.chrom,
                                      source=transcript.source,
                                      featuretype='CDS',
                                      strand=transcript.strand,
                                      start=cdss[-1].start,
                                      stop=cdss[-1].stop,
                                      score='.',
                                      frame=str(nextframe),
                                      attributes=attributes)
                        nextframe = (3 - ((len(cdss[-1]) - nextframe) % 3)) % 3

                    for cds in cdss[1:-1][::-1]:
                        yield Feature(chrom=transcript.chrom,
                                      source=transcript.source,
                                      strand=transcript.strand,
                                      featuretype='CDS',
                                      start=cds.start,
                                      stop=cds.stop,
                                      score='.',
                                      frame=str(nextframe),
                                      attributes=attributes)
                        nextframe = (3 - ((len(cds) - nextframe) % 3)) % 3

                    # last CDS in gene (first in sorted list) has stop_codon's
                    # coords removed
                    yield Feature(chrom=transcript.chrom,
                                  source=transcript.source,
                                  strand=transcript.strand,
                                  featuretype=cdss[0].featuretype,
                                  start=cdss[0].start + 3,
                                  stop=cdss[0].stop,
                                  score='.',
                                  frame=str(nextframe),
                                  attributes=attributes)

    def refFlat(self, gene):
        """Writes the gene out as a RefFlat format:

            geneName altname chrom strand txStart txEnd cdsStart cdsEnd
            exonCount exonStarts exonEnds

        Caveats:

            cdsStart and cdsEnd are the min and max positions, respectively, of
            all CDSs in the _gene_.  So a particular isoform's CDS may not
            actually start on the gene-wide minimum CDS position.

            Assumes that there was an attribute in the GFF file called 'Name'.
            This will then be found in the gene.attributes.Name attribute.
        """

        geneID = gene.id
        for transcript in self.children(geneID, level=1, featuretype='mRNA'):

            txStart = transcript.start
            txEnd = transcript.stop
            exons = []
            cdss = []
            exon_count = 0
            for i in self.children(transcript, 1):
                if i.featuretype == 'CDS':
                    cdss.append(i)
                if i.featuretype == 'exon':
                    exons.append(i)

            if len(cdss) == 0:
                cdsStart = 'NA'
                cdsEnd = 'NA'
            else:
                cdsStart = min([i.start for i in cdss])
                cdsEnd = max([i.stop for i in cdss])

            def exon_sort(e):
                return e.start

            exons.sort(key=exon_sort)

            exonStarts = ','.join([str(i.start) for i in exons]) + ','
            exonEnds = ','.join([str(i.stop) for i in exons]) + ','
            try:
                name = transcript.attributes.Name[0]
            except AttributeError:
                name = transcript.id
            line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                    name,
                    gene.id,
                    gene.chrom,
                    gene.strand,
                    gene.start,
                    gene.stop,
                    cdsStart,
                    cdsEnd,
                    len(exons),
                    exonStarts,
                    exonEnds)
            yield line

    def count_features_of_type(self, featuretype):
        """
        More memory efficient than the functionally-equivalent::

            >>> len(list(self.features_of_type(featuretype)))
        """
        c = self.conn.cursor()
        c.execute('''
        SELECT count() FROM features
        WHERE featuretype = ?''', (featuretype,))
        results = c.fetchone()
        if results is not None:
            results = results[0]
        return results

    def promoter(self, id, dist=1000, truncate_at_next_feature=None,
            direction='upstream'):
        """
        Returns a new Feature of featuretype "promoter", with the definition
        of the promoter guided by the kwargs described below.

        *dist* (default 1000) is the distance in bp from TSS that you want to
        include.  The returned feature will have a legnth of *dist* + 1, since
        the TSS will be included as well.

        TSS is considered the feature start position if the feature is on the
        plus strand, or the feature stop position if on the minus strand. The
        application of *dist* to getting a promoter can be changed by the other
        kwargs below:

        *direction* can be one of 'upstream', 'downstream' or 'both'.

        If *truncate_at_next_feature* (None by default) is a featuretype (e.g.,
        'gene' or 'mRNA') then first try and go up to *dist* bp upstream from
        feature start.  If there is another feature of the specified
        featuretype within *dist* upstream on either strand, then the promoter
        region returned will be truncated to stop at this feature.  This is
        useful if you don't want your promoter regions to overlap with another
        gene.

        If *truncate_at_next_feature* is not None AND *bidirectional* is True,
        then the downstream distance will not be truncated.

        Note that without supplying the chromosome limits, it's impossible to
        know whether the returned promoter runs off the end of the chromosome.
        The lower boundary is accounted for (i.e., no negative coordinates) but
        checking the upper boundary is the responsibility of the calling
        function.
        """
        if not isinstance(id, Feature):
            feature = self[id]
        else:
            feature = id

        c = self.conn.cursor()

        if truncate_at_next_feature is not None:
            # Locate closest feature.  Different strands have different
            # queries:
            #
            # for (+) strand, look for other features' stops that are less than
            # to this features's start and find out how far away it is:
            if feature.strand == '+':
                c.execute('''
                SELECT min(%s-stop)
                FROM features WHERE
                featuretype = ?
                AND chrom = ?
                AND id != ?
                AND stop < ?
                ''' % (feature.start), (truncate_at_next_feature,
                                       feature.chrom,
                                       feature.id,
                                       feature.start))
                closest_feature_dist = c.fetchone()[0]
                if (closest_feature_dist < dist) \
                        and (closest_feature_dist is not None):
                    dist = closest_feature_dist

            # for (-) strand, look for other features' starts that are greater
            # than this feature's stop and find out how far away it is:
            if feature.strand == '-':
                c.execute('''
                SELECT min(start-%s)
                FROM features WHERE
                featuretype = ?
                AND chrom = ?
                AND id != ?
                AND start > ?
                ''' % (feature.stop), (truncate_at_next_feature, feature.chrom,
                                       feature.id, feature.stop))
                closest_feature_dist = c.fetchone()[0]
                if (closest_feature_dist < dist) \
                        and (closest_feature_dist is not None):
                    dist = closest_feature_dist

        if feature.strand == '+':
            TSS = feature.start
            upstream = TSS - dist
            downstream = TSS + dist

        elif feature.strand == '-':
            TSS = feature.stop
            upstream = TSS + dist
            downstream = TSS - dist

        else:
            raise ValueError('Feature strand is "%s" (%s) so promoter '
                             'is ambiguous' % (feature.strand, feature))

        # Negative coords not allowed; truncate to beginning of chrom
        if upstream < 1:
            upstream = 1

        if downstream < 1:
            downstream = 1

        if direction == 'both':
            coords = [upstream, downstream]
        elif direction == 'upstream':
            coords = [upstream, TSS]
        elif direction == 'downstream':
            coords = [downstream, TSS]
        else:
            raise ValueError("need to have a direction of either 'both', "
                             "'upstream', or 'downstream'")
        coords.sort()
        start, stop = coords

        promoter = Feature(
                chrom=feature.chrom,
                source='imputed',
                featuretype='promoter',
                start=start,
                stop=stop,
                strand=feature.strand,
                score='.',
                frame='.',
                attributes='')
        return promoter

    def random_feature(self, featuretype=None):
        """
        Chooses a random feature from the database.  Useful for testing or
        experimenting with the module.  Specify a *featuretype* to restrict
        results to that feature type.

        Idea from here:

            http://www.mail-archive.com/sqlite-users@sqlite.org/msg14657.html
        """

        c = self.conn.cursor()
        featuretype_clause = ''
        featuretype_subclause = ''
        if featuretype is not None:
            featuretype_clause = 'featuretype = "%s" AND ' % featuretype
            featuretype_subclause = 'WHERE featuretype = "%s"' % featuretype
            featuretype_subclause = ''
        c.execute('''
        %s
        FROM features
        WHERE
        %s
        rowid >= abs(random()) %% (SELECT count() FROM features %s) LIMIT 1
        ''' % (self.SELECT, featuretype_clause, featuretype_subclause))
        results = c.fetchone()
        return self._newfeature(*results)

    def coding_genes(self):
        """
        Returns an iterator of genes that contain CDSs. Useful for if you want
        to exclude tRNA, various ncRNAs, etc, since they are also annotated
        with featuretype "gene" and contain 'exon' features (at least in
        FlyBase GFFs)
        """
        for g in self.features_of_type('gene'):
            coding = False
            for grandchild in self.children(g.id, level=2):
                if grandchild.featuretype == 'CDS':
                    coding = True
            if coding:
                yield g

    def n_gene_isoforms(self, geneID):
        """
        Returns the number of isoforms that this gene has.
        """
        if type(geneID) is not str:
            geneID = geneID.id

        n = 0
        for i in self.children(geneID, level=1, featuretype='mRNA'):
            n += 1
        return n

    def n_exon_isoforms(self, exonID):
        """
        Returns the number of isoforms that this exon is found in.
        """

        if type(exonID) is not str:
            exonID = exonID.id

        c = self.conn.cursor()
        c.execute('''
        SELECT count() FROM relations
        JOIN features
        ON relations.parent=features.id
        WHERE relations.child=?
        AND
        relations.level=1
        AND
        features.featuretype="mRNA"
        ''', (exonID,))
        return c.fetchone()[0]

    def exons_gene(self, exonID):
        """
        Returns the ID of the exon's parent gene.  Fast, single-purpose
        method that doesn't do the type conversion or sorting of
        self.parents().
        """

        if type(exonID) is not str:
            exonID = exonID.id

        c = self.conn.cursor()
        c.execute('''
        SELECT parent FROM relations
        JOIN features
        ON relations.parent=features.id
        WHERE child=?
        AND
        relations.level=2
        AND
        features.featuretype="gene"
        ''', (exonID,))
        return c.fetchone()[0]
