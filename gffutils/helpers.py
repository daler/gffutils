import gzip
import os
import tempfile
import gffutils

def example_filename(fn):
    here = os.path.dirname(__file)
    return os.path.join(here, 'test', 'data', fn)


def asinterval(feature):
    """
    Convert a gffutils.Feature object to a pybedtools.Interval, enabling full
    use of pybedtools functionality.  Requires pybedtools.
    """
    try:
        import pybedtools
    except ImportError:
        raise ImportError('Please install pybedtools to convert to '
                          'pybedtools.Interval objects')
    return pybedtools.create_interval_from_list(
            str(feature).split('\t'))


class FeatureNotFoundError(Exception):
    """
    Error to be raised when an ID is not in the database.
    """
    def __init__(self, feature_id):
        Exception.__init__(self)
        self.feature_id = feature_id

    def __str__(self):
        return self.feature_id


def clean_gff(fn, newfn=None, addchr=False, featuretypes_to_remove=None,
        chroms_to_ignore=[], sanity_check=True):
    """
    Helps prepare an optionally gzipped (detected by extension) GFF or GTF file
    *fn* for import into a database. The new, cleaned file will be saved as
    *newfn* (by default, *newfn* will be *gfffn* plus a ".cleaned" extension).

    Use the inspect_featuretypes() function in this module to determine what
    featuretypes are contained in the GFF; then you can filter out those that
    are not interesting by including them in the list of
    *featuretypes_to_remove*

    Supply a list of chromosome in `chroms_to_ignore` to exclude features on
    those chromosomes.

    If *addchr* is True, the string "chr" will be prepended to the chromosome
    name.

    If *sanity_check* is True, only features that pass a simple check will be
    output.  Currently the only sanity check is that coordinates are not
    negative and that feature starts are not greater than feature stop coords.

    Also, some GFF files have FASTA sequence at the end.  This function will
    remove the sequence as well by stopping iteration over the input file when
    it hits a ">" as the first character in the line.
    """
    if newfn is None:
        newfn = fn + '.cleaned'
    if featuretypes_to_remove is None:
        featuretypes_to_remove = []

    fout = open(newfn, 'w')
    if fn.endswith('.gz'):
        infile = gzip.open(fn)
    else:
        infile = open(fn)
    for line in infile:
        if line.startswith('>'):
            break
        L = line.split('\t')

        try:
            if L[0] in chroms_to_ignore:
                continue
            if L[2] in featuretypes_to_remove:
                continue
        except IndexError:
            continue

        if sanity_check:
            start = int(L[3])
            stop = int(L[4])
            if (start < 0) or (stop < 0) or (start > stop):
                continue

        if addchr:
            fout.write('chr' + line)
        else:
            fout.write(line)

    fout.close()


def inspect_featuretypes(gfffn):
    """
    Returns a list of unique featuretypes in the GFF file, *gfffn*.  Useful for
    when you're about to clean up a GFF file and you want to know which
    featuretypes to exclude in the cleaned GFF.
    """
    featuretypes = set()
    for line in open(gfffn):
        L = line.split('\t')
        try:
            featuretypes = featuretypes.union([L[2]])
        except IndexError:
            continue
    return list(featuretypes)

##
## Helpers for gffutils-cli
##
## TODO: move clean_gff here?
##
def get_db_fname(gff_fname,
                 ext=".db"):
    """
    Get db fname for GFF. If the database has a .db file,
    return that. Otherwise, create a named temporary file,
    serialize the db to that, and return the name of the file.

    TODO: Add parameter to control whether or not the .db
    is created or not.
    """
    if not os.path.isfile(gff_fname):
        # Not sure how we should deal with errors normally in
        # gffutils -- Ryan?
        raise Exception, "GFF %s does not exist." %(gff_fname)
        return None
    candidate_db_fname = "%s.%s" %(gff_fname, ext)
    if os.path.isfile(candidate_db_fname):
        # Standard .db file found, so return it
        return candidate_db_fname
    # Otherwise, we need to create a temporary but non-deleted
    # file to store the db in. It'll be up to the user
    # of the function the delete the file when done.
    ## NOTE: Ryan must have a good scheme for dealing with this
    ## since pybedtools does something similar under the hood, i.e.
    ## creating temporary files as needed without over proliferation
    db_fname = tempfile.NamedTemporaryFile(delete=False)
    # Create the database for the gff file (suppress output
    # when using function internally)
    gffutils.create_db(gff_fname, db_fname.name,
                       verbose=False)
    return db_fname


def sanitize(db_fname):
    """
    Sanitize given GFF db. Return a generator of sanitized
    records?

    TODO: Do something with negative coordinates
    """
    db = gffutils.FeatureDB(db_fname)
    for rec in db:
        start = int(rec.start)
        stop = int(rec.stop)
        new_rec = rec.copy()
        if start > stop:
            start
        print "Rec: ", rec

    
