import helpers
import constants
import parser
import bins


class Feature(object):
    def __init__(self, fields, dialect=None):
        """
        Represents a feature from the database.

        When printed, reproduces the original line from the file as faithfully
        as possible using `dialect`.

        Parameters
        `fields` : list

            List of fields from the database.  See :attr:`constants._keys` for
            the order.

        `dialect` : dict or None

            The dialect to use when reconstructing attribute strings; defaults
            to the GFF3 spec.  :class:`FeatureDB` objects will automatically
            attach the dialect from the original file.

        """
        self._fields = fields
        self.dialect = dialect or constants.dialect
        for k, v in zip(constants._keys, fields):
            if k in ('attributes', 'extra'):
                if isinstance(v, basestring):
                    v = helpers._unjsonify(v)
            setattr(self, k, v)

    def __repr__(self):
        memory_loc = hex(id(self))
        return (
            "<Feature {x.featuretype} ({x.seqid}:{x.start}-{x.end}"
            "[{x.strand}]) at {loc}>".format(x=self, loc=memory_loc))

    def __str__(self):
        # all but attributes
        items = [getattr(self, k) for k in constants._gffkeys[:-1]]

        # reconstruct from dict and dialect
        reconstructed_attributes = parser._reconstruct(
            self.attributes, self.dialect)

        # final line includes reconstructed as well as any previously-added
        # "extra" fields
        items.append(reconstructed_attributes)
        print self.extra
        if self.extra:
            items.append('\t'.join(self.extra))
        return '\t'.join(map(str, items))

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

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


def feature_from_line(line):
    """
    Given a line from a GFF file, return a Feature object
    """
    fields = line.rstrip('\n\r').split(None, 9)
    attrs, dialect = parser._split_keyvals(fields[8])
    fields[3] = int(fields[3])
    fields[4] = int(fields[4])
    _bin = bins.bins(int(fields[3]), int(fields[4]), one=True)
    _id = None
    items = [_id] + fields[:8] + [attrs] + [fields[9:]]

    return Feature(items, dialect=dialect)
