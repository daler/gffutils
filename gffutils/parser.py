# Portions copied over from BCBio.GFF.GFFParser

import re
import copy
from gffutils import constants
from gffutils.exceptions import AttributeStringError

import logging

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

gff3_kw_pat = re.compile('\w+=')


def _reconstruct(keyvals, dialect, keep_order=False,
                 sort_attribute_values=False):
    """
    Reconstructs the original attributes string according to the dialect.

    Parameters
    ==========
    keyvals : dict
        Attributes from a GFF/GTF feature

    dialect : dict
        Dialect containing info on how to reconstruct a string version of the
        attributes

    keep_order : bool
        If True, then perform sorting of attribute keys to ensure they are in
        the same order as those provided in the original file.  Default is
        False, which saves time especially on large data sets.

    sort_attribute_values : bool
        If True, then sort values to ensure they will always be in the same
        order.  Mostly only useful for testing; default is False.
    """
    if not dialect:
        raise AttributeStringError()
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
        items = list(keyvals.items())

    def sort_key(x):
        # sort keys by their order in the dialect; anything not in there will
        # be in arbitrary order at the end.
        try:
            return dialect['order'].index(x[0])
        except ValueError:
            return 1e6

    if keep_order:
        items.sort(key=sort_key)

    for key, val in items:

        # Multival sep is usually a comma:
        if val:
            if sort_attribute_values:
                val = sorted(val)

            val_str = dialect['multival separator'].join(val)

            if val_str:

                # Surround with quotes if needed
                if dialect['quoted GFF2 values']:
                    val_str = '"%s"' % val_str

                # Typically "=" for GFF3 or " " otherwise
                part = dialect['keyval separator'].join([key, val_str])
        else:
            if dialect['fmt'] == 'gtf':
                part = dialect['keyval separator'].join([key, '""'])
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
    from gffutils import feature
    quals = feature.dict_class()
    if not keyval_str:
        return quals, dialect

    # If a dialect was provided, then use that directly.
    if not infer_dialect:
        if dialect['trailing semicolon']:
            keyval_str = keyval_str.rstrip(';')

        parts = keyval_str.split(dialect['field separator'])

        kvsep = dialect['keyval separator']
        if dialect['leading semicolon']:
            pieces = []
            for p in parts:
                if p and p[0] == ';':
                    p = p[1:]
                pieces.append(p.strip().split(kvsep))
                key_vals = [(p[0], " ".join(p[1:])) for p in pieces]

        if dialect['fmt'] == 'gff3':
            key_vals = [p.split(kvsep) for p in parts]
        else:
            leadingsemicolon = dialect['leading semicolon']
            pieces = []
            for i, p in enumerate(parts):
                if i == 0 and leadingsemicolon:
                    p = p[1:]
                pieces.append(p.strip().split(kvsep))
                key_vals = [(p[0], " ".join(p[1:])) for p in pieces]

        quoted = dialect['quoted GFF2 values']
        for item in key_vals:
            # Easy if it follows spec
            if len(item) == 2:
                key, val = item

            # Only key provided?
            elif len(item) == 1:
                key = item[0]
                val = ''

            else:
                key = item[0]
                val = dialect['keyval separator'].join(item[1:])

            try:
                quals[key]
            except KeyError:
                quals[key] = []

            if quoted:
                if (len(val) > 0 and val[0] == '"' and val[-1] == '"'):
                    val = val[1:-1]

            if val:
                # TODO: if there are extra commas for a value, just use empty
                # strings
                # quals[key].extend([v for v in val.split(',') if v])
                vals = val.split(',')
                quals[key].extend(vals)

        return quals, dialect

    # If we got here, then we need to infer the dialect....
    #
    # Reset the order to an empty list so that it will only be populated with
    # keys that are found in the file.
    dialect['order'] = []

    # ensembl GTF has trailing semicolon
    if keyval_str[-1] == ';':
        keyval_str = keyval_str[:-1]
        dialect['trailing semicolon'] = True

    # GFF2/GTF has a semicolon with at least one space after it.
    # Spaces can be on both sides (e.g. wormbase)
    # GFF3 works with no spaces.
    # So split on the first one we can recognize...
    for sep in (' ; ', '; ', ';'):
        parts = keyval_str.split(sep)
        if len(parts) > 1:
            dialect['field separator'] = sep
            break

    # Is it GFF3?  They have key-vals separated by "="
    if gff3_kw_pat.match(parts[0]):
        key_vals = [p.split('=') for p in parts]
        dialect['fmt'] = 'gff3'
        dialect['keyval separator'] = '='

    # Otherwise, key-vals separated by space.  Key is first item.
    else:
        dialect['keyval separator'] = " "
        pieces = []
        for p in parts:
            # Fix misplaced semicolons in keys in some GFF2 files
            if p and p[0] == ';':
                p = p[1:]
                dialect['leading semicolon'] = True
            pieces.append(p.strip().split(' '))
        key_vals = [(p[0], " ".join(p[1:])) for p in pieces]

    for item in key_vals:

        # Easy if it follows spec
        if len(item) == 2:
            key, val = item

        # Only key provided?
        elif len(item) == 1:
                key = item[0]
                val = ''

        # Pathological cases where values of a key have within them the key-val
        # separator, e.g., 
        #  Alias=SGN-M1347;ID=T0028;Note=marker name(s): T0028 SGN-M1347 |identity=99.58|escore=2e-126
        else:
            key = item[0]
            val = dialect['keyval separator'].join(item[1:])

        # Is the key already in there?
        if key in quals:
            dialect['repeated keys'] = True
        else:
            quals[key] = []

        # Remove quotes in GFF2
        if len(val) > 0 and val[0] == '"' and val[-1] == '"':
            val = val[1:-1]
            dialect['quoted GFF2 values'] = True
        if val:
            # TODO: if there are extra commas for a value, just use empty
            # strings
            # quals[key].extend([v for v in val.split(',') if v])
            vals = val.split(',')
            if (len(vals) > 1) and dialect['repeated keys']:
                raise AttributeStringError(
                    "Internally inconsistent attributes formatting: "
                    "some have repeated keys, some do not.")
            quals[key].extend(vals)

        # keep track of the order of keys
        dialect['order'].append(key)

    #for key, vals in quals.items():
    #
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

        #quals[key] = vals

    if (
        (dialect['keyval separator'] == ' ') and
        (dialect['quoted GFF2 values'])
    ):
        dialect['fmt'] = 'gtf'

    return quals, dialect
