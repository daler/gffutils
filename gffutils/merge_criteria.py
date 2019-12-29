"""
Merge criteria used by FeatureDB merge() and merge_all()
"""
from gffutils.feature import Feature


def seqid(acc: Feature, cur: Feature, components: [Feature]):
    return cur.seqid == acc.seqid


def strand(acc: Feature, cur: Feature, components: [Feature]):
    return acc.strand == cur.strand


def feature_type(acc: Feature, cur: Feature, components: [Feature]):
    return acc.featuretype == cur.featuretype


def exact_coordinates_only(acc: Feature, cur: Feature, components: [Feature]):
    return cur.start == acc.start and cur.stop == acc.end


def overlap_end_inclusive(acc: Feature, cur: Feature, components: [Feature]):
    return acc.start <= cur.start <= acc.end + 1


def overlap_start_inclusive(acc: Feature, cur: Feature, components: [Feature]):
    return acc.start <= cur.end + 1 <= acc.end + 1


def overlap_any_inclusive(acc: Feature, cur: Feature, components: [Feature]):
    return acc.start <= cur.start <= acc.end + 1 or acc.start <= cur.end + 1 <= acc.end + 1


def overlap_end_threshold(threshold: int):
    def partial(acc: Feature, cur: Feature, components: [Feature]):
        return acc.start <= cur.start <= acc.end + threshold
    return partial


def overlap_start_threshold(threshold: int):
    def partial(acc: Feature, cur: Feature, components: [Feature]):
        return acc.start - threshold <= cur.end + 1 <= acc.end + 1

    return partial


def overlap_any_threshold(threshold: int):
    def partial(acc: Feature, cur: Feature, components: [Feature]):
        return acc.start - threshold <= cur.end + 1 <= acc.end + 1 or acc.start <= cur.start <= acc.end + threshold
    return partial
