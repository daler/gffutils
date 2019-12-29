"""
Merge criteria used by FeatureDB merge() and merge_all()
"""
from gffutils.feature import Feature


def seqid(acc, cur, components):
    return cur.seqid == acc.seqid


def strand(acc, cur, components):
    return acc.strand == cur.strand


def feature_type(acc, cur, components):
    return acc.featuretype == cur.featuretype


def exact_coordinates_only(acc, cur, components):
    return cur.start == acc.start and cur.stop == acc.end


def overlap_end_inclusive(acc, cur, components):
    return acc.start <= cur.start <= acc.end + 1


def overlap_start_inclusive(acc, cur, components):
    return acc.start <= cur.end + 1 <= acc.end + 1


def overlap_any_inclusive(acc, cur, components):
    return acc.start <= cur.start <= acc.end + 1 or acc.start <= cur.end + 1 <= acc.end + 1
