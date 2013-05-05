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
        self.lines_written = 0
        # write header if asked
        if self.with_header:
            timestamp = strftime("%Y-%m-%d %H:%M:%S", gmtime())
            header = "#gffutils GFF3 file (created %s)" %(timestamp)
        

    def write_rec(self, rec):
        """
        Output record to file.
        """
        rec_line = rec.to_string()
        self.out.write("%s\n" %(rec_line))


    def close(self):
        """
        Close the stream. Assumes stream has 'close' method.
        """
        self.out.close()

