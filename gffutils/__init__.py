from gfffeature import Feature, GFFFile
from db import GFFDBCreator, GTFDBCreator, FeatureDB
from helpers import FeatureNotFoundError, clean_gff, inspect_featuretypes
import copy_reg


def feature_constructor(fields):
    return Feature(*fields)


def feature_reducer(feature):
    return feature_constructor, (
            (
                feature.chrom, feature.source, feature.featuretype,
                feature.start, feature.stop, feature.score, feature.strand,
                feature.frame, feature._str_attributes, feature.filetype
            ), )

copy_reg.pickle(Feature, feature_reducer, feature_constructor)
