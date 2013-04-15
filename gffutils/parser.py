# Portions copied over from BCBio.GFF.GFFParser

import re
import urllib
import copy
import constants
from helpers import DefaultOrderedDict, example_filename
import helpers
import bins

import logging


formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

gff3_kw_pat = re.compile('\w+=')


def _reconstruct(keyvals, dialect):
    """
    Reconstructs the original attributes string according to the dialect.
    """
    if not dialect:
        raise helpers.AttributeStringError
    if not keyvals:
        return ""
    parts = []

    # May need to split multiple values into multiple key/val pairs
    if dialect['repeated keys']:
        items = []
        for key, val in keyvals.items():
            if len(val) > 1:
                for v in val:
                    items.append((key, [v]))
            else:
                items.append((key, val))
    else:
        items = keyvals.items()

    for key, val in items:

        # Multival sep is usually a comma:
        val_str = dialect['multival separator'].join(val)
        if val_str:

            # Surround with quotes if needed
            if dialect['quoted GFF2 values']:
                val_str = '"%s"' % val_str

            # Typically "=" for GFF3 or " " otherwise
            part = dialect['keyval separator'].join([key, val_str])
        else:
            part = key
        parts.append(part)

    # Typically ";" or "; "
    parts_str = dialect['field separator'].join(parts)

    # Sometimes need to add this
    if dialect['trailing semicolon']:
        parts_str += ';'

    return parts_str


# TODO:
# Cythonize -- profiling shows that the bulk of the time is spent on this
# function...
def _split_keyvals(keyval_str, dialect=None):
    """
    Given the string attributes field of a GFF-like line, split it into an
    attributes dictionary and a "dialect" dictionary which contains information
    needed to reconstruct the original string.

    Lots of logic here to handle all the corner cases.

    If `dialect` is None, then do all the logic to infer a dialect from this
    attribute string.

    Otherwise, use the provided dialect (and return it at the end).
    """
    infer_dialect = False
    if dialect is None:
        # Make a copy of default dialect so it can be modified as needed
        dialect = copy.copy(constants.dialect)
        infer_dialect = True

    quals = DefaultOrderedDict(list)
    if keyval_str is None:
        return quals

    # ensembl GTF has trailing semicolon
    if keyval_str[-1] == ';':
        keyval_str = keyval_str[:-1]
        if infer_dialect:
            dialect['trailing semicolon'] = True

    # GFF2/GTF has a semicolon with at least one space after it.
    # Spaces can be on both sides (e.g. wormbase)
    # GFF3 works with no spaces.
    # So split on the first one we can recognize...
    if infer_dialect:
        for sep in (' ; ', '; ', ';'):
            parts = keyval_str.split(sep)
            if len(parts) > 1:
                dialect['field separator'] = sep
                break
    else:
        parts = keyval_str.split(dialect['field separator'])

    # Is it GFF3?  They have key-vals separated by "="
    if gff3_kw_pat.match(parts[0]):
        key_vals = [p.split('=') for p in parts]
        if infer_dialect:
            dialect['fmt'] = 'gff3'
            dialect['keyval separator'] = '='

    # Otherwise, key-vals separated by space.  Key is first item.
    else:
        if infer_dialect:
            dialect['keyval separator'] = " "
        pieces = []
        for p in parts:
            # Fix misplaced semicolons in keys in some GFF2 files
            if p and p[0] == ';':
                p = p[1:]
                if infer_dialect:
                    dialect['leading semicolon'] = True
            pieces.append(p.strip().split(' '))
        key_vals = [(p[0], " ".join(p[1:])) for p in pieces]

    for item in key_vals:

        # Easy if it follows spec
        if len(item) == 2:
            key, val = item

        # Only key provided?
        else:
            assert len(item) == 1, item
            key = item[0]
            val = ''

        # Is the key already in there?
        if infer_dialect and key in quals:
            dialect['repeated keys'] = True

        # Remove quotes in GFF2
        if (len(val) > 0 and val[0] == '"' and val[-1] == '"'):
            val = val[1:-1]
            if infer_dialect:
                dialect['quoted GFF2 values'] = True
        if val:
            # TODO: if there are extra commas for a value, just use empty
            # strings
            # quals[key].extend([v for v in val.split(',') if v])
            vals = val.split(',')
            if (len(vals) > 1) and dialect['repeated keys']:
                raise helpers.AttributeStringError(
                    "Internally inconsistent attributes formatting: "
                    "some have repeated keys, some do not.")
            quals[key].extend(vals)

        else:
            # since it's a defaultdict, asking for a key creates the empty list
            quals[key]

    for key, vals in quals.items():

        # TODO: urllib.unquote breaks round trip invariance for "hybrid1.gff3"
        # test file.  This is because the "Note" field has %xx escape chars,
        # but "Dbxref" has ":" which, if everything were consistent, should
        # have also been escaped.
        #
        # (By the way, GFF3 spec says only literal use of \t, \n, \r, %, and
        # control characters should be encoded)
        #
        # Solution 1: don't unquote
        # Solution 2: store, along with each attribute, whether or not it
        #             should be quoted later upon reconstruction
        # Solution 3: don't care about invariance

        # unquoted = [urllib.unquote(v) for v in vals]

        quals[key] = vals

    if (
        (dialect['keyval separator'] == ' ') and
        (dialect['quoted GFF2 values'])
    ):
        dialect['fmt'] = 'gtf'

    return quals, dialect


class Parser(object):
    _keys = ['seqid', 'source', 'feature', 'start', 'end', 'score', 'strand',
             'frame', 'attributes', 'extra']

    _keys = constants._gffkeys_extra

    def __init__(self, filename, checklines=10, transform=None,
                 from_string=False):
        """
        Parse a GFF or GTF file, yielding dictionaries of each line.


        Parameters
        ----------
        `filename` : str

            Filename to parse

        `checklines` : int

            Number of valid lines to check in order to infer a dialect.

        `transform` : callable

            Function (or other callable object) that accepts a dictionary and
            returns a dictionary

        `from_string`: bool

            If True, then consider `filename` to actually be a string
            containing the data.  Tabs are optional for the first 8 fields --
            space-separated works, too.  Useful for testing.

        Notes
        -----
        * Commented lines are ignored
        * Directives are stored in the `directives` attribute
        * Will raise StopIteration when:
            * end of the file is reached
            * `##FASTA` is the only string in the current line
            * current line starts with `>` (again indicating a FASTA record).

        Examples
        --------
        >>> from gffutils.helpers import example_filename
        >>> from gffutils.parser import Parser
        >>> filename = example_filename('F3-unique-3.v2.gff')
        >>> p = Parser(filename)
        >>> linedict = iter(p).next()

        Fields are accessed by their key, and values are always strings:

        >>> linedict['seqid']
        '3_336_815_F3'

        >>> linedict['start']
        '55409'

        >>> linedict['source']
        'solid'

        The dialect -- including the format, along with other details about
        reconstructing the attributes -- is in the `dialect` attribute of the
        parser:

        >>> p.dialect['fmt']
        'gff3'

        The first valid line sets the dialect for the file; inconsistencies
        with this are recorded in self.warnings and can be displayed with
        self.show_warnings().

        Directives are added to self.directives:

        >>> p.directives[0]
        'solid-gff-version 0.2'
        """
        self.filename = filename
        self.from_string = from_string
        self.dialect = None
        self.directives = []
        self.warnings = []
        self.checklines = checklines
        self.transform = transform
        self.current_line = None
        self.current_line_number = None

    def _iter_data(self):
        """
        Returns an iterable of data, whether a filename or a string of data was
        supplied.
        """
        if self.from_string:
            for line in self.filename.splitlines(True):
                items = line.split(None, 8)
                yield '\t'.join(items)
        else:
            for line in open(self.filename):
                yield line

    def __iter__(self):
        return self._parse()

    def _sniff(self):
        for i, line in enumerate(self):
            if i > self.checklines:
                break
        return self.dialect

    def show_warnings(self):
        """
        If there were warnings generated while parsing the file, print them.
        """
        for i in self.warnings:
            logger.warn(i)

    def _dialect_compare(self, dialect, linenumber):
        """
        Compares an observed dialect with self.dialect and records any
        differences in self.warnings.
        """
        orig = set(self.dialect.items())
        new = set(dialect.items())
        removed = dict(list(orig.difference(new)))
        added = dict(list(new.difference(orig)))
        msg = ["Inconsistent attribute formatting:"]
        msg.append("\tline %s in %s" % (linenumber + 1, self.filename))
        for k in removed.keys():
            msg.append(
                '\t[%s]: "%s" --> "%s"' % (k, removed[k], added[k]))
        self.warnings.append('\n'.join(msg))

    def _valid_line_count(self):
        """
        Iterates through the entire file to get the number of valid GFF/GTF
        lines.  Useful if you want feedback about percentage complete.
        """
        n_lines = 0
        for line in self._iter_data():
            if line.startswith(('##FASTA', '>')):
                return n_lines

            # This also ignores directive lines when getting the valid line
            # count.
            if line.startswith('#'):
                continue

            if len(line.rstrip('\n\r')) == 0:
                continue

            n_lines += 1
        return n_lines

    def _parse(self):
        """
        Iterate through the lines of the file.
        """
        valid_lines = 0
        for i, line in enumerate(self._iter_data()):

            line = line.rstrip('\n\r')
            self.current_line = line
            self.current_line_number = i

            if line == '##FASTA' or line.startswith('>'):
                raise StopIteration

            if line.startswith('##'):
                self._directive_handler(line)
                continue

            if line.startswith(('#')) or len(line) == 0:
                continue

            # (If we got here it should be a valid line)

            fields = line.rstrip('\n\r').split('\t')

            # Only infer dialect if we're below the number of lines
            if valid_lines < self.checklines:
                attrs, dialect = _split_keyvals(fields[8])

                # self.dialect should only be None the first time
                if self.dialect is None:
                    self.dialect = dialect

                # It's possible that previous lines didn't have any keys
                # with multiple values, which means that "repeated keys"
                # never got the chance to be set to True. So here we assume
                # that if it was seen once, it should be applied to all
                # lines.
                repeated_keys = dialect['repeated keys'] \
                    or self.dialect['repeated keys']
                dialect['repeated keys'] \
                    = self.dialect['repeated keys'] \
                    = repeated_keys

                # Do the expensive comparison if they don't match
                if dialect != self.dialect:
                    self._dialect_compare(dialect, i)

            # Otherwise if we already have a dialect then just use that.
            else:
                attrs, dialect = _split_keyvals(fields[8],
                                                dialect=self.dialect)

            d = dict(zip(self._keys[:8], fields[:8]))
            d[self._keys[8]] = attrs
            d[self._keys[9]] = fields[9:]

            if self.transform is not None:
                d = self.transform(d)
            yield d

        # Some feedback -- mimic R's "show_warnings"
        if len(self.warnings) > 0:
            logger.warn(
                "%s warnings reported for file %s; use show_warnings() to see "
                "them" % (len(self.warnings), self.filename))

    def _directive_handler(self, line):
        """
        Simply stores directives (minus the "##").
        """
        self.directives.append(line[2:].rstrip('\n\r'))
