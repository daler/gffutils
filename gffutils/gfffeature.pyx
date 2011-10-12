import gzip

cdef class GFFFile:
    cdef public str filename, filetype
    cdef public object fn
    cdef object ignore, only

    def __init__(self, str filename, ignore=None, only=None):
        self.filename = filename
        self.ignore = ignore
        self.only = only
        if only and ignore:
            raise ValueError('Please specify only one of "ignore" or "only"')
        if filename.endswith('.gz'):
            self.fn = gzip.open(filename)
        else:
            self.fn = open(filename)
        try:
            first_feature = self.__next__()
        except StopIteration:
            raise ValueError('No features found in file "%s"' % filename)
        self.filetype = first_feature.filetype
        self.fn.seek(0)

    cdef int _valid_line(self, list fields):
        if fields[0].startswith('#'):
            return 0
        if len(fields) != 9:
            return 0
        if self.ignore:
            if fields[2] in self.ignore:
                return 0
        if self.only:
            if fields[2] in self.only:
                return 1
        return 1

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            line = self.fn.next()

            # special case for flybase files
            if line.startswith('##FASTA'):
                raise StopIteration

            fields = line.split('\t')
            if self._valid_line(fields):
                break
        fields[3] = int(fields[3])
        fields[4] = int(fields[4])
        return Feature(*fields)

    def __len__(self):
        n = sum(1 for _ in self)
        self.reset()
        return n

    def reset(self):
        self.fn.seek(0)

cdef class Feature:
    cdef public int start, stop
    cdef public str chrom, featuretype, source, score, strand, frame
    cdef public str _id, dbid
    cdef public str _str_attributes, filetype
    cdef object _attributes

    def __init__(self, str chrom=".", str source=".", str featuretype=".",
                 int start=1, int stop=1, str score=".", str strand=".",
                 str frame=".", str attributes="", str filetype=""):
        self.chrom = chrom
        self.source = source
        self.featuretype = featuretype
        self.start = start
        self.stop = stop
        self.score = score
        self.strand = strand
        self.frame = frame
        self._str_attributes = attributes.strip()
        self._attributes = None
        self.filetype = filetype
        self.dbid = ""

        # filetype wasn't provided, so we have to guess...guess it's a GFF
        # if there are as many "=" as ";",possibly missing a ";" at the
        # end of the attributes string.
        if self.filetype == "":
            semicolons = attributes.count(';')
            if attributes.count('=') > semicolons - 1:
                self.filetype = 'gff'
            else:
                self.filetype = 'gtf'

    property attributes:
        def __get__(self):
            # Lazy evaluation
            if self._attributes:
                return self._attributes
            else:
                self._attributes = Attributes(
                        self._str_attributes, filetype=self.filetype)
                return self._attributes

        def __set__(self, value):
            if not isinstance(value, Attributes):
                raise ValueError('Feature.attributes must be an '
                                 'Attributes object')
            self._attributes = value

    def __len__(self):
        return self.stop - self.start + 1

    def __str__(self):
        fields = [self.chrom, self.source, self.featuretype,
                str(self.start), str(self.stop), self.score, self.strand,
                self.frame, str(self.attributes)]
        return '\t'.join(fields)

    def __richcmp__(Feature self, Feature other, int cmp):
        """
        <    0
        ==   2
        >    4
        <=   1
        !=   3
        >=   5
        """
        if cmp == 2:
            return str(self) == str(other)
        if cmp == 3:
            return str(self) != str(other)
        else:
            raise ValueError('Comparisons other than "==" and "!=" are not '
                             'supported for Feature objects')

    def __repr__(self):
        return '<Feature: %s, %s:%s-%s (%s)>' % (
                self.featuretype, self.chrom, self.start, self.stop, self.strand)

    property id:
        def __get__(self):
            if self._id:
                return self._id
            if self.filetype == 'gff':
                for key in ("ID", "Name", "gene_name"):
                    try:
                        self._id = self.attributes[key]
                        return self._id
                    except KeyError:
                        pass
                return self.dbid

            if self.filetype == 'gtf':
                if self.featuretype not in ('gene', 'mRNA'):
                    self._id = '%s:%s:%s-%s:%s' % (
                            self.featuretype,
                            self.chrom,
                            self.start,
                            self.stop,
                            self.strand)
                    return self._id
                return self.dbid

        def __set__(self, value):
            self._id = value

    def tostring(self):
        """
        Backwards compatibility function -- simply returns str(self).
        """
        return str(self)

cdef class Attributes:
    """
    Class that acts like a dictionary but prints attributes nicely according to
    filetype.

    Constructor:

        Attributes(attr_str="", filetype="gff")

    Example usage:

        gff_attrs = Attributes('ID=FBgn000001;')
        gtf_attrs = Attributes('gene_id "FBgn000001";')
    """
    cdef dict _attr_dict
    cdef str _attr_str, sep, field_sep, trailing_sep, filetype
    cdef list _field_order
    def __init__(self, attr_str="", filetype="gff"):
        self._attr_str = attr_str.strip()
        self._attr_dict = {}
        self._field_order = []
        # quick exit
        if attr_str == "":
            return

        self.sep = ';'

        if filetype == 'gff':
            self.filetype = 'gff'
            self.field_sep = '='

        if filetype == 'gtf':
            self.filetype = 'gtf'
            self.field_sep = ' '

        # If input had a separator on the end, so should output
        if self._attr_str[-1] == self.sep:
            self.trailing_sep = self.sep
        else:
            self.trailing_sep = ""

        for attribute in attr_str.strip().split(self.sep):
            if attribute.startswith(' '):
                prefix_space = ' '
            else:
                prefix_space = ''
            attribute = attribute.strip()
            if attribute:
                field, value = attribute.split(self.field_sep)

                # comma-separated lists turn into lists; otherwise string
                value = value.replace('"', '').split(',')
                if len(value) == 1:
                    value = value[0]

                self._attr_dict[field] = value

                # Keep track of order, and whether or not the attribute has
                # a preceding space, so that the output is identical
                self._field_order.append((field, prefix_space))

    def __setitem__(self, key, value):
        """
        Sets both the key/item in self.dict *as well as* the interval object's
        attrs field if it's a GFF Interval
        """
        self._attr_dict[key] = value
        if key not in self._field_order:
            # default to no space for GFF, but with space for GTF
            if self.filetype == 'gff':
                self._field_order.append((key, ""))
            if self.filetype == 'gtf':
                self._field_order.append((key, " "))

    def __getitem__(self, key):
        return self._attr_dict[key]

    def __str__(self):
        # stringify all items first
        if len(self._field_order) == 0:
            return ""
        attributes = []
        if self.filetype == 'gtf':
            quotes = '"'
        else:
            quotes = ""
        for field, space in self._field_order:
            value = self._attr_dict[field]
            if isinstance(value, basestring):
                attributes.append((space + field, quotes + value + quotes))
            else:
                attributes.append((space + field, quotes + ','.join(value) + quotes))
        return self.sep.join([self.field_sep.join(kvs) \
                for kvs in attributes]) + self.trailing_sep

    def __repr__(self):
        return repr(self._attr_dict)
