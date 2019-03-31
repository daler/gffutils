import six
import collections

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections
from gffutils import constants


# collections.MutableMapping is apparently the best way to provide dict-like
# interface (http://stackoverflow.com/a/3387975)
class Attributes(collectionsAbc.MutableMapping):
    def __init__(self, *args, **kwargs):
        """
        An Attributes object acts much like a dictionary.  However, values are
        always stored internally as lists, even if a single value is provided.

        Whether or not you get a list back depends on the
        `constants.always_return_list` setting, which can be set on-the-fly.
        If True, then one-item lists are returned.  This is best shown with an
        example:

        Set up an Attributes object:

            >>> attr = Attributes()

        Set the "Name" attribute with a string:

            >>> attr['Name'] = 'gene1'

        This is stored internally as a list, and by default, we'll get a list
        back:

            >>> assert attr['Name'] == ['gene1']

        The same thing happens if we set it with a list in the first place:

            >>> attr['Name'] = ['gene1']
            >>> assert attr['Name'] == ['gene1']

        Now, change the setting so that upon access, single-value lists are
        returned as the first item.

            >>> constants.always_return_list = False
            >>> assert attr['Name'] == 'gene1'

        Change it back again:

            >>> constants.always_return_list = True
            >>> assert attr['Name'] == ['gene1']

        """
        self._d = dict()
        self.update(*args, **kwargs)

    def __setitem__(self, k, v):
        if not isinstance(v, (list, tuple)):
            v = [v]
        self._d[k] = v

    def __getitem__(self, k):
        v = self._d[k]
        if constants.always_return_list:
            return v
        if isinstance(v, list) and len(v) == 1:
            v = v[0]
        return v

    def __delitem__(self, key):
        del self._d[key]

    def __iter__(self):
        return iter(self.keys())

    def __len__(self):
        return len(self._d)

    def keys(self):
        return self._d.keys()

    def values(self):
        return [self.__getitem__(k) for k in self.keys()]

    def items(self):
        r = []
        for k in self.keys():
            r.append((k, self.__getitem__(k)))
        return r

    def __str__(self):
        s = []
        for i in self.items():
            s.append("%s: %s" % i)
        return '\n'.join(s)

    def update(self, *args, **kwargs):
        for k, v in six.iteritems(dict(*args, **kwargs)):
            self[k] = v


# Useful for profiling: which dictionary-like class to store attributes in.
# This is used in Feature below and in parser.py

dict_class = Attributes
#dict_class = dict
#dict_class = helper_classes.DefaultOrderedDict
#dict_class = collections.defaultdict
#dict_class = collections.OrderedDict
#dict_class = helper_classes.DefaultListOrderedDict
