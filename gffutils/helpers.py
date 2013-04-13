import os
import simplejson
from collections import OrderedDict
import constants
import bins

HERE = os.path.dirname(os.path.abspath(__file__))


def example_filename(fn):
    return os.path.join(HERE, 'test', 'data', fn)


def _bin_from_dict(d):
    """
    Given a dictionary yielded by the parser, return the genomic "UCSC" bin
    """
    try:
        start = int(d['start'])
        end = int(d['end'])
        return bins.bins(start, end, one=True)

    # e.g., if "."
    except ValueError:
        return None


class FeatureNotFoundError(Exception):
    """
    Error to be raised when an ID is not in the database.
    """
    def __init__(self, feature_id):
        Exception.__init__(self)
        self.feature_id = feature_id

    def __str__(self):
        return self.feature_id


class DuplicateIDError(Exception):
    pass


class AttributeStringError(Exception):
    pass


def _jsonify(x):
    """Use most compact form of JSON"""
    return simplejson.dumps(x, separators=(',', ':'))


def _unjsonify(x):
    """Convert JSON string to an ordered defaultdict."""
    return simplejson.loads(
        x, object_pairs_hook=DefaultListOrderedDict)


def _dict_to_fields(d, jsonify=True):
    """
    Convert dict to tuple, for faster sqlite3 import
    """
    x = []
    for k in constants._keys:
        v = d[k]
        if jsonify and (k in ('attributes', 'extra')):
            x.append(_jsonify(v))
        else:
            x.append(v)
    return tuple(x)


def asinterval(feature):
    """
    Converts a gffutils.Feature to a pybedtools.Interval
    """
    import pybedtools
    return pybedtools.create_interval_from_list(str(feature).split('\t'))


class DefaultOrderedDict(OrderedDict):
    """
    OrderedDict that is also a defaultdict.

    From http://stackoverflow.com/a/6190500.
    (using callable() instead of Callable as per the comments to that answer)
    """
    def __init__(self, default_factory=None, *a, **kw):
        if (
            default_factory is not None and
            not callable(default_factory)
        ):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.items())


class DefaultListOrderedDict(DefaultOrderedDict):
    """
    defaultdict(list), but ordered.
    """
    def __init__(self, *a, **kw):
        super(DefaultListOrderedDict, self).__init__(list, *a, **kw)
        self.default_factory = list


if __name__ == "__main__":
    d = DefaultListOrderedDict([('a', 1), ('b', 2)])
