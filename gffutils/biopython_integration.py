"""
Module for integration with BioPython, specifically SeqRecords and SeqFeature
objects.
"""
import six
try:
    from Bio.SeqFeature import SeqFeature, FeatureLocation
except ImportError:
    import warnings
    warnings.warn("BioPython must be installed to use this module")
from .feature import Feature, feature_from_line

_biopython_strand = {
    '+':  1,
    '-': -1,
    '.':  0,
}
_feature_strand = dict((v, k) for k, v in _biopython_strand.items())


def to_seqfeature(feature):
    """
    Converts a gffutils.Feature object to a Bio.SeqFeature object.

    The GFF fields `source`, `score`, `seqid`, and `frame` are stored as
    qualifiers.  GFF `attributes` are also stored as qualifiers.

    Parameters
    ----------
    feature : Feature object, or string
        If string, assume it is a GFF or GTF-format line; otherwise just use
        the provided feature directly.
    """
    if isinstance(feature, six.string_types):
        feature = feature_from_line(feature)

    qualifiers = {
        'source': [feature.source],
        'score': [feature.score],
        'seqid': [feature.seqid],
        'frame': [feature.frame],
    }
    qualifiers.update(feature.attributes)
    return SeqFeature(
        # Convert from GFF 1-based to standard Python 0-based indexing used by
        # BioPython
        FeatureLocation(feature.start - 1, feature.stop),
        id=feature.id,
        type=feature.featuretype,
        strand=_biopython_strand[feature.strand],
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

    # BioPython parses 1-based GenBank positions into 0-based for use within
    # Python.  We need to convert back to 1-based GFF format here.
    start = s.location.start.position + 1
    stop = s.location.end.position
    featuretype = s.type
    id = s.id
    attributes = dict(s.qualifiers)
    attributes.pop('source', '.')
    attributes.pop('score', '.')
    attributes.pop('seqid', '.')
    attributes.pop('frame', '.')
    return Feature(seqid, source, featuretype, start, stop, score, strand,
                   frame, attributes, id=id, **kwargs)
