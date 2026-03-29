# Portions copied over from BCBio.GFF.GFFParser

import re
import collections
from urllib import parse
from gffutils import constants
from gffutils.exceptions import AttributeStringError

import logging

formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

# Regex for each separator that will be tested
quoted_semicolon_patterns = dict()

for sep in (" ; ", "; ", ";"):
    quoted_semicolon_patterns[sep] = re.compile(
        rf"""
            {re.escape(sep)}   # The separator we're considering (escaped for VERBOSE mode)
            (?=                # Positive lookahead: does remaining content match?
                (?:            # Start non-capturing group
                    [^"]       # Either: match any character that is NOT a quote
                    |          # OR
                    "[^"]*"    # Match a complete quoted string, specifically:
                               #   - opening quote ", followed by
                               #   - zero or more non-quote characters [^"]*
                               #   - followed by closing quote "
                )*             # Repeat the above pattern zero or more times
                $              # Until we reach the end of the string
            )                  # End of lookahead
        """,
        re.VERBOSE,
    )

# Encoding/decoding notes
# -----------------------
# From
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#description-of-the-format:
#
#       GFF3 files are nine-column, tab-delimited, plain text files.
#       Literal use of tab, newline, carriage return, the percent (%) sign,
#       and control characters must be encoded using RFC 3986
#       Percent-Encoding; no other characters may be encoded. Backslash and
#       other ad-hoc escaping conventions that have been added to the GFF
#       format are not allowed. The file contents may include any character
#       in the set supported by the operating environment, although for
#       portability with other systems, use of Latin-1 or Unicode are
#       recommended.
#
#           tab (%09)
#           newline (%0A)
#           carriage return (%0D)
#           % percent (%25)
#           control characters (%00 through %1F, %7F)
#
#       In addition, the following characters have reserved meanings in
#       column 9 and must be escaped when used in other contexts:
#
#           ; semicolon (%3B)
#           = equals (%3D)
#           & ampersand (%26)
#           , comma (%2C)
#
#
# See also issue #98.
#
# Note that spaces are NOT supposed to be encoded. Yet some GFF files have
# spaces encoded anyway; in these cases round-trip invariance will not hold
# since the %20 will be decoded but not re-encoded.
_to_quote = "\n\t\r%;=&,"
_to_quote += "".join([chr(i) for i in range(32)])
_to_quote += chr(127)


# Caching idea from urllib.parse.Quoter, which uses a defaultdict for
# efficiency. Here we're sort of doing the reverse of the "reserved" idea used
# there.
class Quoter(collections.defaultdict):
    def __missing__(self, b):
        if b != "" and b in _to_quote:
            res = "%{:02X}".format(ord(b))
        else:
            res = b
        self[b] = res
        return res


quoter = Quoter()


def _split_keyvals(keyval_str, dialect=None):
    """
    Dialect detection requires partially parsing the attributes.
    """
    from gffutils import feature

    quals = feature.dict_class()

    if not keyval_str:
        return quals, dialect

    infer_dialect = False
    if dialect is None:
        infer_dialect = True
        dialect = {}

    # No known cases yet of different multival separator
    dialect["multival separator"] = ","

    # Detection for these dialect fields can work on the full attribute
    # string. Other detection needs to wait until we've further parsed the
    # attributes.
    if infer_dialect:
        dialect["trailing semicolon"] = keyval_str[-1] == ";"
        dialect["leading semicolon"] = keyval_str[0] == ";"
        semicolon_in_quotes = False
        sep = None
        for sep in (" ; ", "; ", ";"):
            parts = keyval_str.split(sep)
            if len(parts) > 1:
                # If naive split differs from more expensive regex, we infer there was
                # a semicolon within quoted value and we'll have to use the expensive
                # method later
                parts_regex = re.split(quoted_semicolon_patterns[sep], keyval_str)
                if parts != parts_regex:
                    semicolon_in_quotes = True
                break
        dialect["semicolon in quotes"] = semicolon_in_quotes
        dialect["field separator"] = sep

    if dialect["trailing semicolon"]:
        keyval_str = keyval_str.rstrip(";")

    if dialect["leading semicolon"]:
        keyval_str = keyval_str.lstrip(";")

    if dialect["semicolon in quotes"]:
        parts = re.split(
            quoted_semicolon_patterns[dialect["field separator"]], keyval_str
        )
    else:
        parts = keyval_str.split(dialect["field separator"])

    # The next stage of dialect inference works on the 'parts' -- unsplit
    # keyval pairs -- like:
    #
    #    parts = ["ID=001", "Name=gene1"]
    #
    # or
    #
    #    parts = ["gene_id ENSG001", "gene_biotype protein_coding"]
    #
    if infer_dialect:
        dialect["fmt"] = "gff3"

        # Note: so far, have not found cases where we need to check more than
        # the first item
        if "=" in parts[0]:
            dialect["fmt"] = "gff3"
            dialect["keyval separator"] = "="
        else:
            dialect["fmt"] = "gtf"
            dialect["keyval separator"] = " "

    # Now we split
    #
    #    parts = ["ID=001", "Name=gene1"]
    #
    # into
    #
    #   key_val_tuples = [("ID", "001"), ("Name", "gene1")]
    #
    # in a dialect-dependent manner.
    kvsep = dialect["keyval separator"]
    key_val_tuples = [p.split(kvsep) for p in parts]

    # With the split keys we can detect whether any are repeated
    if infer_dialect:
        keys = [i[0] for i in key_val_tuples]
        dialect["repeated keys"] = len(keys) != len(set(keys))

        # For dialect detection, this will help figure out if there is
        # inconsistent quoting across values. It will only be used in the loop
        # below if infer_dialect is True
        quoted_values = []

    # Now work splitting the keys if needed.
    for i in key_val_tuples:

        if len(i) == 2:
            # Easy, on-spec case
            key, val = i

        elif len(i) == 1:
            # By convention, no value becomes an empty string, e.g. when done
            # parsing,
            #
            #   "ID=001;is_gene;"
            #
            # will end up as:
            #
            #   {"ID": "001", "is_gene": ""}
            key = i[0]
            val = ""

        else:
            # Multiple *spaces* within quoted values are joined back together
            # without requiring a regex, in contrast to when there's *field*
            # separator like a semicolon in the values.
            #
            # That is:
            #
            #   attributes = 'gene_description "an important gene"; gene_id "g001"'
            #
            # when split on spaces, becomes
            #
            #   key_val_tuples = [("gene_description", "an", "important", "gene"), ("gene_id", "g001")]
            #
            # so here when we only keep the first token as a key, that first
            # key/val pair will become:
            #
            #   {
            #     "gene_description": ["an important gene"],
            #     "gene_id": ["g001"],
            #   }
            #
            # Another pathological case, this time for GFF3:
            #
            #   Alias=SGN-M1347;ID=T0028;Note=marker name(s): T0028 SGN-M1347 |identity=99.58|escore=2e-126
            #
            # will become the following:
            #
            #   {
            #     "Alias": ["SGN-M1347"],
            #     "ID": ["T0028"],
            #     "Note": ["marker name(s): T0028 SGN-M1347 |identity=99.58|escore=2e-126"],
            #   }
            #
            key = i[0]
            val = kvsep.join(i[1:])

        # By convention all values are lists, even if there's only one value
        # (or even no values)
        if key not in quals:
            quals[key] = []

        # This will run on every value, accumulating in quoted_values to check
        # later for consistency
        if infer_dialect:
            quoted = len(val) > 0 and val[0] == '"' and val[-1] == '"'
            quoted_values.append(quoted)
            dialect["quoted GFF2 values"] = quoted

        if dialect["quoted GFF2 values"] and val:
            val = val.strip('"')

        if val:
            # For repeated keys dialect, don't split on an internal comma. That is,
            #
            #   attributes = 'db_xref="g01,g02"; db_xref="XYZ"'
            #
            # becomes:
            #
            #   {
            #     "db_xref": ["g01,g02", "XYZ"]
            #    }
            #
            if dialect.get("repeated keys"):
                quals[key].append(val)

            # Otherwise, split but only if it's a comma without a space. So:
            #
            #    attributes = 'db_xref="g01,g02"'
            #
            # becomes
            #    {
            #      "db_xref": ["g01", "g02"]
            #    }
            # but
            #
            #    attributes = 'description="kinase, subunit 1"'
            #                                      ^ note the space here
            # becomes
            #    {
            #      "description": ["kinase, subunit 1"]
            #    }
            #
            else:
                # E.g. the "kinase, subunit 1" example above
                if ", " in val:
                    quals[key].append(val)
                else:
                    quals[key].extend(val.split(","))

    # If there was inconsistent quoting, we fall back to "not quoted" so
    # as to avoid incorrectly stripping off first and last quotes.
    if infer_dialect and len(set(quoted_values)) > 1:
        # Prior behavior was to use whatever the first value used
        dialect["quoted GFF2 values"] = quoted_values[0]

        # Though there could be an argument for considering quotes in mixed
        # cases to be part of the string, though technically they should be
        # %-encoded if so.
        # dialect["quoted GFF2 values"] = False

    # Handle unquoting of %-encoded values
    if not constants.ignore_url_escape_characters and dialect["fmt"] == "gff3":
        for key, vals in quals.items():
            unquoted = [parse.unquote(v) for v in vals]
            quals[key] = unquoted

    # Now that we're not supporting old Python versions we can rely on dict
    # insertion order
    if infer_dialect:
        dialect["order"] = list(quals.keys())

    return quals, dialect


def _reconstruct(keyvals, dialect, keep_order=False, sort_attribute_values=False):
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

    # Re-encode when reconstructing attributes
    if constants.ignore_url_escape_characters or dialect["fmt"] != "gff3":
        attributes = keyvals
    else:
        attributes = {}
        for k, v in keyvals.items():
            attributes[k] = []
            for i in v:
                attributes[k].append("".join([quoter[j] for j in i]))

    # May need to split multiple values into multiple key/val pairs
    if dialect["repeated keys"]:
        items = []
        for key, val in attributes.items():
            if len(val) > 1:
                for v in val:
                    items.append((key, [v]))
            else:
                items.append((key, val))
    else:
        items = list(attributes.items())

    def sort_key(x):
        # sort keys by their order in the dialect; anything not in there will
        # be in arbitrary order at the end.
        try:
            return dialect["order"].index(x[0])
        except ValueError:
            return 1e6

    if keep_order:
        items.sort(key=sort_key)

    for key, val in items:

        # Multival sep is usually a comma:
        if val:
            if sort_attribute_values:
                val = sorted(val)

            val_str = dialect["multival separator"].join(val)

            if val_str:

                # Surround with quotes if needed
                if dialect["quoted GFF2 values"]:
                    val_str = '"%s"' % val_str

                # Typically "=" for GFF3 or " " otherwise
                part = dialect["keyval separator"].join([key, val_str])
            else:
                part = key
        else:
            if dialect["fmt"] == "gtf":
                # By convention, GTF attributes with no value are reconstructed
                # with an empty string. E.g.:
                #   'gene_id "gene1"; is_gene;'
                #
                # becomes
                #
                #   {
                #     "gene_id": "gene1",
                #     "is_gene": ""
                #    }
                #
                # and is printed as:
                #
                #   'gene_id "gene1"; is_gene "";'
                part = dialect["keyval separator"].join([key, '""'])
            else:
                part = key
        parts.append(part)

    # Typically ";" or "; "
    parts_str = dialect["field separator"].join(parts)

    # Sometimes need to add this
    if dialect["trailing semicolon"]:
        parts_str += ";"

    return parts_str
