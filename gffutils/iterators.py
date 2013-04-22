import os
import tempfile
import itertools
from feature import Feature, feature_from_line


def peek(it, n):
    _peek = []
    for _ in range(n):
        try:
            _peek.append(it.next())
        except StopIteration:
            break
    return _peek, itertools.chain(_peek, it)


class Directive(object):
    def __init__(self, line):
        self.info = line


class BaseIterator(object):
    def __init__(self, data, checklines=10, transform=None,
                 force_dialect_check=False):
        self.data = data
        self.checklines = checklines
        self.current_item = None
        self.current_item_number = None
        self.dialect = None
        self._observed_dialects = []
        self.directives = []
        self.transform = transform
        self.warnings = []

        if not force_dialect_check:
            self.peek, self._iter = peek(self._custom_iter(), checklines)
            self._observed_dialects = [i.dialect for i in self.peek]
            self.dialect = self._choose_dialect(self._observed_dialects)
        else:
            self._iter = self._custom_iter()

    def _custom_iter(self):
        raise NotImplementedError("Must define in subclasses")

    def __iter__(self):
        for i in self._iter:
            if self.transform:
                yield self.transform(i)
            else:
                yield i

    def _choose_dialect(self, dialects):
        # Do something easy for now; eventually some heuruistics here? Or at
        # least choose the most common dialect observed.  Can use
        # helpers.dialect_compare.
        try:
            return dialects[0]
        except IndexError:
            return None

    def _directive_handler(self, directive):
        self.directives.append(directive)


class FileIterator(BaseIterator):
    def _custom_iter(self):
        valid_lines = 0
        for i, line in enumerate(open(self.data)):

            line = line.rstrip('\n\r')
            self.current_item = line
            self.current_item_number = i

            if line == '##FASTA' or line.startswith('>'):
                raise StopIteration

            if line.startswith('##'):
                self._directive_handler(line)
                continue

            if line.startswith(('#')) or len(line) == 0:
                continue

            # (If we got here it should be a valid line)
            valid_lines += 1
            yield feature_from_line(line, dialect=self.dialect)


class FeatureIterator(BaseIterator):
    def _custom_iter(self):
        for i, feature in enumerate(self.data):
            self.current_item = feature
            self.current_item_number = i
            yield feature


class StringIterator(FileIterator):
    def _custom_iter(self):
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(self.data)
        tmp.close()
        self.data = tmp.name
        for feature in super(StringIterator, self)._custom_iter():
            yield feature
        os.unlink(tmp.name)
