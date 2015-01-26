from gffutils import iterators
from gffutils import interface
from collections import Counter
import sys


def inspect(data, look_for=['featuretype', 'chrom', 'attribute_keys',
                            'feature_count'], limit=None, verbose=True):
    """
    Inspect a GFF or GTF data source.

    This function is useful for figuring out the different featuretypes found
    in a file (for potential removal before creating a FeatureDB).

    Returns a dictionary with a key for each item in `look_for` and
    a corresponding value that is a dictionary of how many of each unique item
    were found.

    There will always be a `feature_count` key, indicating how many features
    were looked at (if `limit` is provided, then `feature_count` will be the
    same as `limit`).

    For example, if `look_for` is ['chrom', 'featuretype'], then the result
    will be a dictionary like::

        {
            'chrom': {
                'chr1': 500,
                'chr2': 435,
                'chr3': 200,
                ...
                ...
            }.

            'featuretype': {
                'gene': 150,
                'exon': 324,
                ...
            },

            'feature_count': 5000

        }


    Parameters
    ----------
    data : str, FeatureDB instance, or iterator of Features
        If `data` is a string, assume it's a GFF or GTF filename.  If it's
        a FeatureDB instance, then its `all_features()` method will be
        automatically called. Otherwise, assume it's an iterable of Feature
        objects.

    look_for : list
        List of things to keep track of. Options are:

            - any attribute of a Feature object, such as chrom, source, start,
              stop, strand.

            - "attribute_keys", which will look at all the individual
              attribute keys of each feature

    limit : int
        Number of features to look at.  Default is no limit.

    verbose : bool
        Report how many features have been processed.

    Returns
    -------
    dict
    """

    results = {}
    obj_attrs = []
    for i in look_for:
        if i not in ['attribute_keys', 'feature_count']:
            obj_attrs.append(i)
        results[i] = Counter()

    attr_keys = 'attribute_keys' in look_for

    d = iterators.DataIterator(data)
    feature_count = 0
    for f in d:
        if verbose:
            sys.stderr.write('\r%s features inspected' % feature_count)
            sys.stderr.flush()

        for obj_attr in obj_attrs:
            results[obj_attr].update([getattr(f, obj_attr)])

        if attr_keys:
            results['attribute_keys'].update(f.attributes.keys())

        feature_count += 1
        if limit and feature_count == limit:
            break

    new_results = {}
    for k, v in results.items():
        new_results[k] = dict(v)

    new_results['feature_count'] = feature_count
    return new_results
