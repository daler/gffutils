"""
Module for integration with BioPython, specifically SeqRecords and SeqFeature
objects.
"""
try:
    from Bio.SeqFeature import SeqFeature, FeatureLocation
except ImportError:
    raise ImportError(
        "BioPython must be installed to use this module")
from . import Feature

_biopython_strand = {
    '+':  1,
    '-': -1,
    '.':  0,
}
_feature_strand = dict((v, k) for k, v in _biopython_strand.items())


def to_seqfeature(f):
    """
    Converts a gffutils.Feature object to a Bio.SeqFeature object.

    The GFF fields `source`, `score`, `seqid`, and `frame` are stored as
    qualifiers.  GFF `attributes` are also stored as qualifiers.
    """
    qualifiers = {
        'source': [f.source],
        'score': [f.score],
        'seqid': [f.seqid],
        'frame': [f.frame],
    }
    qualifiers.update(f.attributes)
    return SeqFeature(
        FeatureLocation(f.start, f.stop),
        id=f.id,
        type=f.featuretype,
        strand=_biopython_strand[f.strand],
        qualifiers=qualifiers
    )


def from_seqfeature(s, **kwargs):
    """
    Converts a Bio.SeqFeature object to a gffutils.Feature object.

    The GFF fields `source`, `score`, `seqid`, and `frame` are assumed to be
    stored as qualifiers.  Any other qualifiers will be assumed to be GFF
    attributes.
    """
    source = s.qualifiers.get('source', '.')[0]
    score = s.qualifiers.get('score', '.')[0]
    seqid = s.qualifiers.get('seqid', '.')[0]
    frame = s.qualifiers.get('frame', '.')[0]
    strand = _feature_strand[s.strand]
    start = s.location.start.position
    stop = s.location.end.position
    featuretype = s.type
    id = s.id
    attributes = dict(s.qualifiers)
    attributes.pop('source', '.')
    attributes.pop('score', '.')
    attributes.pop('seqid', '.')
    attributes.pop('frame', '.')
    return Feature(seqid, source, featuretype, start, stop, score, strand,
                   frame, attributes, **kwargs)
