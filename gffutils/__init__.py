import os
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


def create_db(fn, dbfn, verbose=True, force=False):
    """
    Detects input format (GFF or GTF and calls the appropriate creation code,
    resulting in a new database ready for use with FeatureDB.

    `fn` is the filename of a GFF or GTF file

    `dbfn` is the filename of the database to create.

    Use `verbose=False` to disable progress info.

    Use `force=True` to force overwrite of `dbfn` if it already exists.
    """
    first_feature = GFFFile(fn).next()
    fmt = first_feature.filetype
    if fmt == 'gff':
        creator = GFFDBCreator(fn, dbfn, verbose=verbose)
    else:
        creator = GTFDBCreator(fn, dbfn, verbose=verbose)
    creator.create()
