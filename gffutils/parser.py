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
    if not keyval_str:
        return quals, dialect


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
