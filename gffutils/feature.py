from collections import defaultdict, MutableMapping
import simplejson
import helpers
import constants
import parser
import bins

# http://stackoverflow.com/a/3387975
class Attributes(MutableMapping):
    def __init__(self, *args, **kwargs):
        self._order = []
        self._d = dict()
        self.update(dict(*args, **kwargs))

    def __setitem__(self, key, value):
        if not isinstance(value, list):
            raise ValueError(
                'Attribute value "%s" for key "%s" must be a list'
                % (repr(value), repr(key)))
        self._d[key] = value
        self._order.append(key)

    def __getitem__(self, key):
        return self._d[key]

    def __delitem__(self, key):
        del self._d[key]

    def __iter__(self):
        return iter(self.keys)

    def __len__(self):
        return len(self._d)

    def keys(self):
        return self._order

    def values(self):
        return [self._d[k] for k in self._order]

    def items(self):
        return [(k, self._d[k]) for k in self._order]

    def __str__(self):
        s = []
        for i in self.items():
            s.append("%s: %s" % i)
        return '\n'.join(s)

    def update(self, *args, **kwargs):
        for k, v in dict(*args, **kwargs).iteritems():
            self[k] = v


class Feature(object):
    def __init__(self, seqid=".", source=".", featuretype=".",
                 start=".", end=".", score=".", strand=".", frame=".",
                 attributes=None, extra=None, bin=None, id=None, dialect=None,
                 file_order=None):
        """
        Represents a feature from the database.

        When printed, reproduces the original line from the file as faithfully
        as possible using `dialect`.

        Parameters
        ----------
        `seqid` : string

        `source`: string

        `featuretype` : string

        `start`, `end`: int or "."
            If "." (the default placeholder for GFF files), then the
            corresponding attribute will be None.

        `score`: string

        `strand`: string

        `frame`: string

        `attributes`: string or dict
            If a string, first assume it is serialized JSON; if this fails then
            assume it's the original key/vals string.  If it's a dictionary
            already, then use as-is.

            The end result is that this instance's `attributes` attribute will
            always be a dictionary.

            Upon printing, the attributes will be reconstructed based on this
            dictionary and the dialect -- except if the original attributes
            string was provided, in which case that will be used directly.

        `extra`: string or list
            Additional fields after the canonical 9 fields for GFF/GTF.

            If a string, then first assume it's serialized JSON; if this fails
            then assume it's a tab-delimited string of additional fields.  If
            it's a list already, then use as-is.

        `bin`: int
            UCSC genomic bin. If None, will be created based on provided
            start/end; if start or end is "." then bin will be None.

        `id` : None or string
            Database-specific primary key for this feature.  The only time this
            should not be None is if this feature is coming from a database, in
            which case it will be filled in automatically.

        `dialect` : dict or None

            The dialect to use when reconstructing attribute strings; defaults
            to the GFF3 spec.  :class:`FeatureDB` objects will automatically
            attach the dialect from the original file.

        `file_order` : int
            This is the `rowid` special field used in a sqlite3 database; this
            is provided by FeatureDB.

        """
        # start/end can be provided as int-like, ".", or None, but will be
        # converted to int or None
        if start == ".":
            start = None
        elif start is not None:
            start = int(start)
        if end == ".":
            end = None
        elif end is not None:
            end = int(end)

        # Flexible handling of attributes:
        # If dict, then use that; otherwise assume JSON; otherwise assume
        # original string.
        self._orig_attribute_str = None
        attributes = attributes or helpers.DefaultListOrderedDict()
        if isinstance(attributes, basestring):
            try:
                attributes = helpers._unjsonify(attributes)

            # it's a string but not JSON: assume original attributes string.
            except simplejson.JSONDecodeError:
                # Saved for later printing
                self._orig_attribute_str = attributes

                # But Feature.attributes is still a dict
                attributes, _dialect = parser._split_keyvals(attributes)

                # Use this dialect if none provided.
                dialect = dialect or _dialect

        # If string, then assume tab-delimited; otherwise list
        extra = extra or []
        if isinstance(extra, basestring):
            try:
                extra = helpers._unjsonify(extra)
            except simplejson.JSONDecodeError:
                extra = extra.split('\t')

        # Calculate bin if not provided
        if bin is None:
            try:
                bin = bins.bins(start, end, one=True)
            except TypeError:
                bin = None

        self.seqid = seqid
        self.source = source
        self.featuretype = featuretype
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = attributes
        self.extra = extra
        self.bin = bin
        self.id = id
        self.dialect = dialect or constants.dialect
        self.file_order = file_order

    def __repr__(self):
        memory_loc = hex(id(self))

        # Reconstruct start/end as "."
        if self.start is None:
            start = '.'
        else:
            start = self.start
        if self.end is None:
            end = '.'
        else:
            end = self.end

        return (
            "<Feature {x.featuretype} ({x.seqid}:{start}-{end}"
            "[{x.strand}]) at {loc}>".format(x=self, start=start, end=end,
                                             loc=memory_loc))

    def __getitem__(self, key):
        if key in constants._keys:
            return getattr(self, key)
        else:
            raise KeyError

    def __setitem__(self, key, value):
        if key in constants._keys:
            self.key = value
        else:
            raise KeyError

    def __str__(self):
        # All fields but attributes (and extra).
        items = [getattr(self, k) for k in constants._gffkeys[:-1]]
        if items[3] is None:
            items[3] = "."
        if items[4] is None:
            items[4] = "."

        # Reconstruct from dict and dialect (only if original attributes
        # weren't provided)
        if self._orig_attribute_str:
            reconstructed_attributes = self._orig_attribute_str
        else:
            reconstructed_attributes = parser._reconstruct(
                self.attributes, self.dialect)

        # Final line includes reconstructed as well as any previously-added
        # "extra" fields
        items.append(reconstructed_attributes)
        if self.extra:
            items.append('\t'.join(self.extra))
        return '\t'.join(map(str, items))

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    def __len__(self):
        return self.stop - self.start + 1

    @property
    def chrom(self):
        return self.seqid

    @chrom.setter
    def chrom(self, value):
        self.seqid = value

    @property
    def stop(self):
        return self.end

    @stop.setter
    def stop(self, value):
        self.end = value

    def astuple(self):
        t = []
        for i in constants._keys:
            v = getattr(self, i)
            if i in ('attributes', 'extra'):
                v = helpers._jsonify(v)
            t.append(v)
        return tuple(t)


def feature_from_line(line, dialect=None):
    """
    Given a line from a GFF file, return a Feature object

    Parameters
    ----------
    `line`: string
        As long as there are only 9 fields (standard GFF/GTF), then it's OK to
        use spaces instead of tabs to separate fields in `line`.  But if >9
        fields are to be used, then tabs must be used.
    """
    if '\t' in line:
        fields = line.rstrip('\n\r').split('\t')
    else:
        fields = line.rstrip('\n\r').split(None, 8)
    try:
        attr_string = fields[8]
    except IndexError:
        attr_string = ""
    attrs, _dialect = parser._split_keyvals(attr_string, dialect=dialect)
    d = dict(zip(constants._gffkeys, fields))
    d['attributes'] = attrs
    d['extra'] = fields[9:]
    if dialect is None:
        dialect = _dialect
    return Feature(dialect=dialect, **d)
