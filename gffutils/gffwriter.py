##
## GFF Writer (writer): serializing gffutils records as GFF text files.
##
import six
import tempfile
import shutil
from time import strftime, localtime
from gffutils.version import version


class GFFWriter:
    """
    Simple GFF writer class for serializing gffutils
    records to a file.

    Parameters:
    -----------

    out: string or file-like object
        If a string, parsed as a filename, otherwise, a file-like object to
        write to.

    with_header: bool
        If True, output a header file for the GFF.  The header indicates the
        file was created with gffutils and includes a timestamp and version
        info.

    in_place: bool
        If True and if `out` is a filename, then write the file in place (uses
        named temporary files.)

    TODO: Add make separate GTFWriter class or add support
    for GTF here.
    """
    def __init__(self, out, with_header=True, in_place=False):
        self.out = out
        self.with_header = with_header
        self.in_place = in_place
        # Temporary file to be used (only applies when in_place is True)
        self.temp_file = None
        # Output stream to write to
        self.out_stream = None
        if isinstance(out, six.string_types):
            if self.in_place:
                # Use temporary file
                self.temp_file = tempfile.NamedTemporaryFile(delete=False)
                self.out_stream = open(self.temp_file.name, "w")
            else:
                # Just use the filename given
                self.out_stream = open(self.out, "w")
        else:
            # Assumed to be a write-able stream
            if self.in_place:
                # The in_place parameter is undefined for
                # streams, since no filenames are involved
                raise ValueError("Cannot use 'in_place' when writing to "
                                 "a stream.")
            self.out_stream = out
        # write header if asked
        if self.with_header:
            timestamp = strftime("%Y-%m-%d %H:%M:%S", localtime())
            header = "#GFF3 file (created by gffutils (v%s) on %s)" \
                % (version, timestamp)
            self.out_stream.write("%s\n" % (header))

    def write_rec(self, rec):
        """
        Output record to file.
        """
        self.out_stream.write("%s\n" % rec)

    def write_recs(self, recs):
        """
        Output several records to file.
        """
        for rec in recs:
            self.write_rec(rec)

    def write_gene_recs(self, db, gene_id):
        """
        NOTE: The goal of this function is to have a canonical ordering when
        outputting a gene and all of its records to a file.  The order is
        intended to be:

        gene
          # mRNAs sorted by length, with longest mRNA first
          mRNA_1
            # Exons of mRNA, sorted by start position (ascending)
            exon_1
              # Children of exon, sorted by start position
              exon_child_1
              exon_child_2
            exon_2
              ...
            # Non-exonic children here
            ...
          mRNA_2
            ...
          # Non-mRNA children here
          ...

        Output records of a gene to a file, given a GFF database
        and a gene_id. Outputs records in canonical order: gene record
        first, then longest mRNA, followed by longest mRNA exons,
        followed by rest, followed by next longest mRNA, and so on.

        Includes the gene record itself in the output.

        TODO: This probably doesn't handle deep GFF hierarchies.
        """
        gene_rec = db[gene_id]
        # Output gene record
        self.write_rec(gene_rec)
        # Get each mRNA's lengths
        mRNA_lens = {}
        c = list(db.children(gene_id, featuretype="mRNA"))
        for mRNA in db.children(gene_id, featuretype="mRNA"):
            mRNA_lens[mRNA.id] = \
                sum(len(exon) for exon in db.children(mRNA,
                                                      featuretype="exon"))
        # Sort mRNAs by length
        sorted_mRNAs = \
            sorted(mRNA_lens.items(), key=lambda x: x[1], reverse=True)
        for curr_mRNA in sorted_mRNAs:
            mRNA_id = curr_mRNA[0]
            mRNA_rec = db[mRNA_id]
            # Write mRNA record to file
            self.write_rec(mRNA_rec)
            # Write mRNA's children records to file
            self.write_mRNA_children(db, mRNA_id)
        # Write non-mRNA children of gene (only level1)
        for gene_child in db.children(gene_id, level=1):
            if gene_child.featuretype != "mRNA":
                self.write_rec(gene_child)

    def write_mRNA_children(self, db, mRNA_id):
        """
        Write out the children records of the mRNA given by the ID
        (not including the mRNA record itself) in a canonical
        order, where exons are sorted by start position and given
        first.
        """
        mRNA_children = db.children(mRNA_id, order_by='start')
        nonexonic_children = []
        for child_rec in mRNA_children:
            if child_rec.featuretype == "exon":
                self.write_rec(child_rec)
                self.write_exon_children(db, child_rec)
            else:
                nonexonic_children.append(child_rec)
        self.write_recs(nonexonic_children)

    def write_exon_children(self, db, exon_id):
        """
        Write out the children records of the exon given by
        the ID (not including the exon record itself).
        """
        exon_children = db.children(exon_id, order_by='start')
        for exon_child in exon_children:
            self.write_rec(exon_child)

    def close(self):
        """
        Close the stream. Assumes stream has 'close' method.
        """
        self.out_stream.close()
        # If we're asked to write in place, substitute the named
        # temporary file for the current file
        if self.in_place:
            shutil.move(self.temp_file.name, self.out)
