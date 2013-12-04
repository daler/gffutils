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
from feature import Feature, feature_from_line
import helpers
from textwrap import dedent


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
                 force_dialect_check=False, dialect=None):
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
                yield self.transform(i)
            else:
                yield i

    def _directive_handler(self, directive):
        self.directives.append(directive[2:])


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
        self.tmp = tempfile.NamedTemporaryFile(delete=False)
        self.tmp.write(dedent(self.data))
        self.tmp.close()
        self.data = self.tmp.name
        for feature in super(StringIterator, self)._custom_iter():
            yield feature
        os.unlink(self.tmp.name)


def DataIterator(data, checklines=10, transform=None,
                 force_dialect_check=False, from_string=False, **kwargs):
    _kwargs = dict(data=data, checklines=checklines, transform=transform,
                   force_dialect_check=force_dialect_check, **kwargs)
    if isinstance(data, basestring):
        if from_string:
            return StringIterator(**_kwargs)
        else:
            return FileIterator(**_kwargs)
    else:
        return FeatureIterator(**_kwargs)
