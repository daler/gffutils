"""
This module provides different kinds of iterators, all wrapped by the
DataIterator class which should generally be the only one used in practice.

The BaseIterator class allows "peeking" `checklines` lines into the data --
even if it's a consumable iterator -- in order to figure out what the dialect
is and therefore decide whether the data is GFF or GTF format, which is
important for figuring out how to construct the database.

"""
import os
import tempfile
import itertools
from gffutils.feature import feature_from_line
from gffutils.interface import FeatureDB
from gffutils import helpers
from textwrap import dedent
import six
from six.moves.urllib.request import urlopen
if six.PY3:
    from urllib import parse as urlparse
else:
    import urlparse


def peek(it, n):
    _peek = []
    for _ in range(n):
        try:
            _peek.append(six.next(it))
        except StopIteration:
            break
    return _peek, itertools.chain(_peek, it)


class Directive(object):
    def __init__(self, line):
        self.info = line


class _BaseIterator(object):
    def __init__(self, data, checklines=10, transform=None,
                 force_dialect_check=False, dialect=None):
        """
        Base class for iterating over features.  In general, you should use
        DataIterator -- so see the docstring of class for argument
        descriptions.


        All subclasses -- _FileIterator, _URLIterator, _FeatureIterator,
        _StringIterator -- gain the following behavior:

            - self.current_item and self.current_item_number are set on every
              iteration.  This is very useful for debugging, or reporting to
              the user exactly what item or line number caused the issue.

            - transform a Feature before it gets yielded, filter out a Feature

            - auto-detect dialect by peeking `checklines` items into the
              iterator, and then re-reading those, applying the detected
              dialect.  If multiple dialects are found, use
              helpers._choose_dialect to figure out the best one.

            - keep track of directives

        """
        self.data = data
        self.checklines = checklines
        self.current_item = None
        self.current_item_number = None
        self.dialect = None
        self._observed_dialects = []
        self.directives = []
        self.transform = transform
        self.warnings = []

        if force_dialect_check and dialect is not None:
            raise ValueError("force_dialect_check is True, but a dialect "
                             "is provided")
        if force_dialect_check:
            # In this case, self.dialect remains None.  When
            # parser._split_keyvals gets None as a dialect, it tries to infer
            # a dialect.
            self._iter = self._custom_iter()
        elif dialect is not None:
            self._observed_dialects = [dialect]
            self.dialect = helpers._choose_dialect(self._observed_dialects)
            self._iter = self._custom_iter()
        else:
            # Otherwise, check some lines to determine what the dialect should
            # be
            self.peek, self._iter = peek(self._custom_iter(), checklines)
            self._observed_dialects = [i.dialect for i in self.peek]
            self.dialect = helpers._choose_dialect(self._observed_dialects)

    def _custom_iter(self):
        raise NotImplementedError("Must define in subclasses")

    def __iter__(self):
        for i in self._iter:
            i.dialect = self.dialect
            if self.transform:
                i = self.transform(i)
                if i:
                    yield i
            else:
                yield i

    def _directive_handler(self, directive):
        self.directives.append(directive[2:])


class _FileIterator(_BaseIterator):
    """
    Subclass for iterating over features provided as a filename
    """
    def open_function(self, data):
        data = os.path.expanduser(data)
        if data.endswith('.gz'):
            import gzip
            return gzip.open(data)
        return open(data)

    def _custom_iter(self):
        valid_lines = 0
        for i, line in enumerate(self.open_function(self.data)):
            if isinstance(line, six.binary_type):
                line = line.decode('utf-8')
            line = line.rstrip('\n\r')
            self.current_item = line
            self.current_item_number = i

            if line == '##FASTA' or line.startswith('>'):
                return

            if line.startswith('##'):
                self._directive_handler(line)
                continue

            if line.startswith(('#')) or len(line) == 0:
                continue

            # (If we got here it should be a valid line)
            valid_lines += 1
            yield feature_from_line(line, dialect=self.dialect)


class _UrlIterator(_FileIterator):
    """
    Subclass for iterating over features provided as a URL
    """
    def open_function(self, data):
        response = urlopen(data)

        # ideas from
        # http://stackoverflow.com/a/17537107
        # https://rationalpie.wordpress.com/2010/06/02/\
        #               python-streaming-gzip-decompression/
        if data.endswith('.gz'):
            import zlib
            d = zlib.decompressobj(16 + zlib.MAX_WBITS)
            READ_BLOCK_SIZE = 1024

            def _iter():
                last_line = ""
                while True:
                    data = response.read(READ_BLOCK_SIZE)
                    if not data:
                        break
                    data = "".join((last_line, d.decompress(data).decode()))
                    lines = data.split('\n')
                    last_line = lines.pop()
                    for line in lines:
                        yield line + '\n'
                yield last_line
            return _iter()

        else:
            return response


class _FeatureIterator(_BaseIterator):
    """
    Subclass for iterating over features that are already in an iterator
    """
    def _custom_iter(self):
        for i, feature in enumerate(self.data):
            self.current_item = feature
            self.current_item_number = i
            yield feature


class _StringIterator(_FileIterator):
    """
    Subclass for iterating over features provided as a string (e.g., from
    file.read())
    """
    def _custom_iter(self):
        self.tmp = tempfile.NamedTemporaryFile(delete=False)
        data = dedent(self.data)
        if isinstance(data, six.text_type):
            data = data.encode('utf-8')
        self.tmp.write(data)
        self.tmp.close()
        self.data = self.tmp.name
        for feature in super(_StringIterator, self)._custom_iter():
            yield feature
        os.unlink(self.tmp.name)


def is_url(url):
    """
    Check to see if a URL has a valid protocol.

    Parameters
    ----------
    url : str or unicode

    Returns
    -------
    True if `url` has a valid protocol False otherwise.
    """
    try:
        return urlparse.urlparse(url).scheme in urlparse.uses_netloc
    except:
        return False


def DataIterator(data, checklines=10, transform=None,
                 force_dialect_check=False, from_string=False, **kwargs):
    """
    Iterate over features, no matter how they are provided.

    Parameters
    ----------
    data : str, iterable of Feature objs, FeatureDB
        `data` can be a string (filename, URL, or contents of a file, if
        from_string=True), any arbitrary iterable of features, or a FeatureDB
        (in which case its all_features() method will be called).

    checklines : int
        Number of lines to check in order to infer a dialect.

    transform : None or callable
        If not None, `transform` should accept a Feature object as its only
        argument and return either a (possibly modified) Feature object or
        a value that evaluates to False.  If the return value is False, the
        feature will be skipped.

    force_dialect_check : bool
        If True, check the dialect of every feature.  Thorough, but can be
        slow.

    from_string : bool
        If True, `data` should be interpreted as the contents of a file rather
        than the filename itself.

    dialect : None or dict
        Provide the dialect, which will override auto-detected dialects.  If
        provided, you should probably also use `force_dialect_check=False` and
        `checklines=0` but this is not enforced.
    """

    _kwargs = dict(data=data, checklines=checklines, transform=transform,
                   force_dialect_check=force_dialect_check, **kwargs)
    if isinstance(data, six.string_types):
        if from_string:
            return _StringIterator(**_kwargs)
        else:
            if os.path.exists(data):
                return _FileIterator(**_kwargs)
            elif is_url(data):
                return _UrlIterator(**_kwargs)
    elif isinstance(data, FeatureDB):
        _kwargs['data'] = data.all_features()
        return _FeatureIterator(**_kwargs)

    else:
        return _FeatureIterator(**_kwargs)
