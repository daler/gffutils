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
        first_feature = self.__next__()
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
    cdef public str _id
    cdef public str _str_attributes, filetype
    cdef object _attributes

    def __init__(self, str chrom=".", str source=".", str featuretype=".",
                 int start=1, int stop=1, str score=".", str strand=".",
                 str frame=".", str attributes=""):
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
        self.filetype = ""

        if '=' in attributes:
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
            if self.filetype == 'gtf':
                if self.featuretype not in ('gene', 'mRNA'):
                    self._id = '%s:%s:%s:%s:%s' % (
                            self.featuretype,
                            self.chrom,
                            self.start,
                            self.stop,
                            self.strand)
                    return self._id

        def __set__(self, value):
            self._id = value


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

    def __init__(self, attr_str="", filetype="gff"):
        self._attr_str = attr_str.strip()
        self._attr_dict = {}
        # quick exit
        if attr_str == "":
            return

        if filetype == 'gff':
            self.filetype = 'gff'
            self.sep, self.field_sep = (';', '=')
        if filetype == 'gtf':
            self.filetype = 'gtf'
            self.sep, self.field_sep = (';', ' ')

        if attr_str[-1] == self.sep:
            self.trailing_sep = self.sep
        else:
            self.trailing_sep = ""

        kvs = map(str.strip, attr_str.strip().split(self.sep))
        for kv in kvs:
            if kv:
                field, value = kv.split(self.field_sep)
                value = value.replace('"', '').split(',')
                if len(value) == 1:
                    value = value[0]
                self._attr_dict[field] = value

    def __setitem__(self, key, value):
        """
        Sets both the key/item in self.dict *as well as* the interval object's
        attrs field if it's a GFF Interval
        """
        self._attr_dict[key] = value

    def __getitem__(self, key):
        return self._attr_dict[key]

    def __str__(self):
        # stringify all items first
        items = []
        if self.filetype == 'gtf':
            quotes = '"'
        else:
            quotes = ""
        for i, j in self._attr_dict.items():
            if isinstance(j, basestring):
                items.append((i, quotes + j + quotes))
        return self.sep.join([self.field_sep.join(kvs) \
                for kvs in items]) + self.trailing_sep

    def __repr__(self):
        return repr(self._attr_dict)
