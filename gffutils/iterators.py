import os
import gzip
import zlib
import tempfile
import itertools
from gffutils.feature import feature_from_line, Feature
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


# String, URL, and File iterators each yield strings.
# DataIterator dispatches to the correct generator function, and provides it to
# BaseIterator.
#
# BaseIterator handles dialect detection, directives, transformation, etc.

def lines_from_file(data):
    """
    Yield lines from a file (gzipped or not), stopping when we see a line start
    with "##FASTA" or ">"
    """
    if data.endswith('.gz'):
        _open = gzip.open
    else:
        _open = open
    with _open(data) as in_:
        for line in in_:
            if isinstance(line, six.binary_type):
                line = line.decode('utf-8')
            line = line.rstrip('\n\r')
            if line.startswith('##'):
                yield line
                continue
            if line.startswith(('##FASTA', '>')):
                raise StopIteration
            if (len(line) == 0) or (line[0] == '#'):
                continue
            yield line


def lines_from_url(data):
    """
    Yield lines from a URL (gzipped or not, detected by ".gz" extension).
    """
    response = urlopen(data)
    if data.endswith('.gz'):
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        READ_BLOCK_SIZE = 1024
        last_line = ""
        while True:
            data = response.read(READ_BLOCK_SIZE)
            if not data:
                break
            data = "".join((last_line, d.decompress(data).decode()))
            lines = data.split('\n')
            last_line = lines.pop()
            for line in lines:
                if line.startswith('##'):
                    yield line
                    continue
                if line.startswith(('##FASTA', '>')):
                    raise StopIteration
                if (len(line) == 0) or (line[0] == '#'):
                    continue
                yield line + '\n'
        yield last_line

    else:
        for line in response:
            line = line.decode()
            if line.startswith('##'):
                yield line
                continue
            if line.startswith(('##FASTA', '>')):
                raise StopIteration
            if (len(line) == 0) or (line[0] == '#'):
                continue
            yield line


def lines_from_string(data):
    """
    Yield lines from a string.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    data = dedent(data)
    if isinstance(data, six.text_type):
        data = data.encode('utf-8')
    with open(tmp, 'wb') as out:
        out.write(data)
    for line in lines_from_file(tmp):
        yield line
    os.unlink(tmp)


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


class DataIterator(object):
    def __init__(self, data, checklines=10, transform=None,
                 force_dialect_check=False, from_string=False, dialect=None):
        """
        Iterate over features, no matter how they are provided.

        Parameters
        ----------
        data : str, iterable of Feature objs, FeatureDB
            `data` can be a string (filename, URL, or contents of a file, if
            from_string=True), any arbitrary iterable of features, or
            a FeatureDB (in which case its all_features() method will be
            called).

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
            If True, `data` should be interpreted as the contents of a file
            rather than the filename itself.

        dialect : None or dict
            Provide the dialect, which will override auto-detected dialects.
            If provided, you should probably also use
            `force_dialect_check=False` and `checklines=0` but this is not
            enforced.
        """

        if isinstance(data, six.string_types):
            if from_string:
                data = lines_from_string(data)
            else:
                if os.path.exists(data):
                    data = lines_from_file(data)
                elif is_url(data):
                    data = lines_from_url(data)
        elif isinstance(data, FeatureDB):
            data = data.all_features()
        elif isinstance(data, (list, tuple)):
            data = iter(data)
        self.data = data
        self.checklines = checklines
        self.current_item = None
        self.current_item_number = 0
        self.dialect = None
        self._observed_dialects = []
        self.transform = transform
        self.directives = []
        self.warnings = []

        if force_dialect_check and dialect is not None:
            raise ValueError("force_dialect_check is True, but a dialect"
                             " is provided")
        if force_dialect_check:
            # In this case, self.dialect remains None.  When
            # parser._split_keyvals gets None as a dialect, it tries to infer
            # a dialect.
            pass
        elif dialect is not None:
            self._observed_dialects = [dialect]
            self.dialect = helpers._choose_dialect(self._observed_dialects)

        else:
            self.peek, self.data = peek(self.data, checklines)
            for i in self.peek:
                if isinstance(i, six.string_types):
                    i = feature_from_line(i)
                self._observed_dialects.append(i.dialect)
            self.dialect = helpers._choose_dialect(self._observed_dialects)

    def _directive_handler(self, directive):
        self.directives.append(directive[2:])

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        while True:
            line = six.next(self.data)

            if not isinstance(line, Feature):

                if line.startswith('##'):
                    self._directive_handler(line)
                    continue

                feature = feature_from_line(line, dialect=self.dialect)
            else:
                feature = line

            self.current_item = line
            self.current_item_number += 1

            if self.transform is not None:
                return self.transform(feature)
            return feature

if __name__ == "__main__":
    fn = '/home/ryan/proj/gffutils/gffutils/test/data/FBgn0031208.gtf'
    d = DataIterator(fn)
