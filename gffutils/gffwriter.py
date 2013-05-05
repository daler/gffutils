##
## GFF Writer (writer): serializing gffutils records as GFF text files.
##
##   Dear Sir or Madam, will you read my code?
##   It took me years to write, will you take a look?
##   It's based on code by a man named Daler
##   And I need a job, so I want to be a GFF Writer, GFF Writer.
##
import os
import sys
import time
from time import strftime, gmtime
from gfffeature import GFFFile, Feature
from helpers import FeatureNotFoundError, asinterval


class GFFWriter:
    """
    Simple GFF writer class for serializing gffutils
    records to a file.

    TODO: Add make separate GTFWriter class or add support
    for GTF here.
    """
    def __init__(self, out, with_header=True):
        # write-able output stream
        self.out = out
        # whether or not to write header for GFF file
        self.with_header = with_header
        # write header if asked
        if self.with_header:
            timestamp = strftime("%Y-%m-%d %H:%M:%S", gmtime())
            header = "#GFF3 file (created by gffutils on %s)" %(timestamp)
            self.out.write("%s\n" %(header))


    def write_rec(self, rec):
        """
        Output record to file.
        """
        rec_line = rec.tostring()
        self.out.write("%s\n" %(rec_line))


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
              # Children of exon, no particular order
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
        print "GENE ID: %s" %(gene_id)
        print "mRNA CHILDREN OF GENE: "
        c = list(db.children(gene_id, featuretype="mRNA"))
        print "C: ", c
        print [x.featuretype for x in c]
        for mRNA in db.children(gene_id, featuretype="mRNA"):
            mRNA_lens[mRNA.id] = \
                sum(len(exon) for exon in db.children(mRNA,
                                                      featuretype="exon"))
        # Sort mRNAs by length
        sorted_mRNAs = \
            sorted(mRNA_lens.items(), key=lambda x: x[1], reverse=True)
        print "sorted MRNA: ", sorted_mRNAs
        raise Exception
        for curr_mRNA in sorted_mRNAs:
            mRNA_id = curr_mRNA[0]
            mRNA_rec = db[mRNA_id]
            # Write mRNA record to file
            self.write_rec(mRNA_rec)
            print "writing MRNA --> %s" %(mRNA_id)
            # Write mRNA's children records to file
            self.write_mRNA_children(db, mRNA_id)
        # Write non-mRNA children of gene
        for gene_child in db.children(gene_id):
            if gene_child.featuretype != "mRNA":
                self.write_rec(gene_child)


    def write_mRNA_children(self, db, mRNA_id):
        """
        Write out the children records of the mRNA given by the ID
        (not including the mRNA record itself) in a canonical
        order, where exons are sorted by start position and given
        first.
        """
        mRNA_children = db.children(mRNA_id)
        nonexonic_children = []
        # Write out the exons first, sorted by position
        exon_starts = {}
        for child_rec in mRNA_children:
            if child_rec.featuretype == "exon":
                # Record start positions of all exons
                exon_starts[child_rec.id] = \
                    child_rec.start
            else:
                nonexonic_children.append(child_rec)
        sorted_exons = \
            sorted(exon_starts.items(), key=lambda x: x[1])
        print "SORTED EXONS: ", sorted_exons
        for curr_exon in sorted_exons:
            exon_id = curr_exon[0]
            exon_rec = db[exon_id]
            # Write out the exon 
            self.write_rec(exon_rec)
            # Write out exon's chilren
            self.write_exon_children(db, exon_id)
        # Output remaining record types
        self.write_recs(nonexonic_children)


    def write_exon_children(self, db, exon_id):
        """
        Write out the children records of the exon given by
        the ID (not including the exon record itself).
        """
        exon_children = db.children(exon_id)
        for exon_child in exon_children:
            self.write_rec(exon_child)
                

    def close(self):
        """
        Close the stream. Assumes stream has 'close' method.
        """
        self.out.close()

