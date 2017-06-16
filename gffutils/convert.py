"""
Conversion functions that operate on :class:`FeatureDB` classes.
"""

import six

def to_bed12(f, db, child_type='exon', name_field='ID'):
    """
    Given a top-level feature (e.g., transcript), construct a BED12 entry
    Parameters
    ----------
    f : Feature object or string
        This is the top-level feature represented by one BED12 line.  For
        a canonical GFF or GTF, this will generally be a transcript.
    db : a FeatureDB object
        This is need to get the children for the feature
    child_type : str
        Featuretypes that will be represented by the BED12 "blocks".  Typically
        "exon".
    name_field : str
        Attribute to be used in the "name" field of the BED12 entry.  Usually
        "ID" for GFF; "transcript_id" for GTF.
    """
    if isinstance(f, six.string_types):
        f = db[f]
    children = list(db.children(f, featuretype=child_type, order_by='start'))
    sizes = [len(i) for i in children]
    starts = [i.start - f.start for i in children]
    fields = [
        f.chrom,
        f.start - 1,  # GTF -> BED coord system
        f.stop,
        f.attributes.get(name_field, ['.'])[0],
        f.score,
        f.strand,
        f.start,
        f.stop,
        '0,0,0',
        len(children),
        ','.join(map(str, sizes)),
        ','.join(map(str, starts))
    ]
    return '\t'.join(map(str, fields)) + '\n'
