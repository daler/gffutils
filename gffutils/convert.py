"""
Conversion functions that operate on :class:`FeatureDB` classes.
"""

def to_bed12(f, db, child_type='exon', name_field='ID'):
    """
    Given a top-level feature (e.g., transcript), construct a BED12 entry
    """
    if isinstance(f, basestring):
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

