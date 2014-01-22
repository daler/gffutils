#!/usr/bin/env python

import os
import logging
from collections import Counter
import gffutils
from gffutils.create import logger
logger.setLevel(logging.INFO)

usage = """

    Converts a GFF file downloaded from FlyBase into a GTF file suitable for
    use with Cufflinks.  Support for ignoring chromosomes, featuretypes, and
    fixing various issues with the original GFF file.
"""


def read_ignore_list(fn):
    """
    Converts a file into a list
    """
    return [i.strip() for i in open(fn) if not i.startswith('#')]


def clean_gff(gff, cleaned, add_chr=False, chroms_to_ignore=None,
              featuretypes_to_ignore=None):
    """
    Cleans a GFF file by removing features on unwanted chromosomes and of
    unwanted featuretypes.  Optionally adds "chr" to chrom names.
    """
    logger.info("Cleaning GFF")
    chroms_to_ignore = chroms_to_ignore or []
    featuretypes_to_ignore = featuretypes_to_ignore or []
    with open(cleaned, 'w') as fout:
        for i in gffutils.iterators.DataIterator(gff):
            if add_chr:
                i.chrom = "chr" + i.chrom

            if i.chrom in chroms_to_ignore:
                continue

            if i.featuretype in featuretypes_to_ignore:
                continue
            fout.write(str(i) + '\n')
    return cleaned


if __name__ == "__main__":

    import argparse
    ap = argparse.ArgumentParser(usage=usage)
    ap.add_argument('--gff', help='''GFF file downloaded from FlyBase.  Not
                    required if --db already exists''')
    ap.add_argument('--db', required=True, help='''Database that will be
                    created from FlyBase GFF''')
    ap.add_argument('--gtf', required=True, help='''GTF file that will
                    be created''')
    ap.add_argument('--featuretype-ignore', help='''File containing the
                    featuretypes to ignore in the GFF file, one featuretype on
                    each line.  This can speed up processing dramatically''')
    ap.add_argument('--chrom-ignore', help='''File containing chromosomes to
                    ignore in the GFF file.  If --add-chr, then chromosomes in
                    this list should also have "chr"''')
    ap.add_argument('--add-chr', action='store_true', help='''Prepend a "chr"
                    to chromosome names''')
    ap.add_argument('--fix-strand', action='store_true', help='''When exporting
                    the GTF file, check to ensure that all exons of a gene have
                    the same strand.  If not, then issue a warning and fix them
                    so that they all have the most common strand across all the
                    gene's exons.''')
    ap.add_argument('--fix-extent', choices=['gene', 'exon'], help='''When
                    exporting the GTF file, check that all exons fall within
                    the gene's annotated extent.  If not, then issue a warning,
                    and, if --fix-extent=gene, then assume the gene annotation
                    is correct; if --fix-extent=exon then assume the exon
                    annotations are correct and make no changes to the exported
                    exons''')
    args = ap.parse_args()

    if args.gff:
        chroms_to_ignore = []
        if args.chrom_ignore:
            chroms_to_ignore = read_ignore_list(args.chrom_ignore)

        featuretypes_to_ignore = []
        if args.featuretype_ignore:
            featuretypes_to_ignore = read_ignore_list(args.featuretype_ignore)

        cleaned = clean_gff(
            args.gff,
            args.gff + '.cleaned',
            add_chr=args.add_chr,
            chroms_to_ignore=chroms_to_ignore,
            featuretypes_to_ignore=featuretypes_to_ignore)

        db = gffutils.create_db(args.gff + '.cleaned', args.db, verbose=True,
                                id_spec=['ID'], force=True)

    db = gffutils.FeatureDB(args.db)

    logger.info("Creating GTF file")
    with open(args.gtf, 'w') as fout:
        for gene in db.features_of_type('gene'):
            gene_id = gene.id
            transcripts = []
            exons = []
            for child in db.children(gene, level=1):
                transcript_id = child.id
                transcripts.append(child)
                for grandchild in db.children(child, level=1):
                    if grandchild.featuretype not in ['exon', 'CDS']:
                        continue
                    exons.append(grandchild)

            if args.fix_strand:
                c = Counter([i.strand for i in exons])

                if len(c) > 1:
                    # Exons have inconsistent strands.  So assume the most
                    # common strand is the "true" strand to use.
                    most_common = c.most_common()
                    new_strand = most_common[0][0]
                    logger.warning(
                        'Gene %s has inconsistent strands: %s.  '
                        'Changing all exons to "%s" strand'
                        % (gene_id, most_common, new_strand)
                    )
                    new_exons = []
                    for exon in exons:
                        exon.strand = new_strand
                        new_exons.append(exon)
                    exons = new_exons

            if args.fix_extent:
                if len(exons) > 0:
                    exon_extent = [
                        min(i.start for i in exons),
                        max(i.stop for i in exons)
                    ]
                    gene_extent = [gene.start, gene.stop]
                    if exon_extent != gene_extent:
                        logger.warning(
                            "Exons of gene %s do not match gene annotation. "
                            "Gene extent: %s.  Exon extent: %s"
                            % (gene_id, exon_extent, gene_extent)
                        )

                    if args.fix_extent == 'gene':
                        new_exons = sorted(exons, key=lambda x: x.start)
                        new_exons[0].start = gene.start
                        new_exons[-1].stop = gene.stop
                        exons = new_exons

            for exon in exons:
                fields = str(grandchild).split('\t')[:-1]
                attributes = ('gene_id "{0}"; transcript_id "{1}"; '
                              'gene_name "{2}" transcript_type="{3}"'
                              ).format(
                                  gene_id,
                                  transcript_id,
                                  gene.attributes['Name'][0],
                                  child.featuretype)
                fields.append(attributes)
                fout.write('\t'.join(fields) + '\n')
    logger.info("Wrote %s" % args.gtf)
