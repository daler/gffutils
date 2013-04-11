import helpers
import constants
import parser

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
        if self.extra:
            items.append('\t'.join(extra))
        return '\t'.join(map(str, items))

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)
