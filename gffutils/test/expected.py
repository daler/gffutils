# expected data for tests using FBgn0031208.gff and FBgn0031208.gtf files


# list the children and their expected first-order parents for the GFF test file.
GFF_parent_check_level_1 = {'FBtr0300690':['FBgn0031208'],
                            'FBtr0300689':['FBgn0031208'],
                            'CG11023:1':['FBtr0300689','FBtr0300690'],
                            'five_prime_UTR_FBgn0031208:1_737':['FBtr0300689','FBtr0300690'],
                            'CDS_FBgn0031208:1_737':['FBtr0300689','FBtr0300690'],
                            'intron_FBgn0031208:1_FBgn0031208:2':['FBtr0300690'],
                            'intron_FBgn0031208:1_FBgn0031208:3':['FBtr0300689'],
                            'FBgn0031208:3':['FBtr0300689'],
                            'CDS_FBgn0031208:3_737':['FBtr0300689'],
                            'CDS_FBgn0031208:2_737':['FBtr0300690'],
                            'exon:chr2L:8193-8589:+':['FBtr0300690'],
                            'intron_FBgn0031208:2_FBgn0031208:4':['FBtr0300690'],
                            'three_prime_UTR_FBgn0031208:3_737':['FBtr0300689'],
                            'FBgn0031208:4':['FBtr0300690'],
                            'CDS_FBgn0031208:4_737':['FBtr0300690'],
                            'three_prime_UTR_FBgn0031208:4_737':['FBtr0300690'],
                           }

# and second-level . . . they should all be grandparents of the same gene.
GFF_parent_check_level_2 = {
                            'CG11023:1':['FBgn0031208'],
                            'five_prime_UTR_FBgn0031208:1_737':['FBgn0031208'],
                            'CDS_FBgn0031208:1_737':['FBgn0031208'],
                            'intron_FBgn0031208:1_FBgn0031208:2':['FBgn0031208'],
                            'intron_FBgn0031208:1_FBgn0031208:3':['FBgn0031208'],
                            'FBgn0031208:3':['FBgn0031208'],
                            'CDS_FBgn0031208:3_737':['FBgn0031208'],
                            'CDS_FBgn0031208:2_737':['FBgn0031208'],
                            'exon:chr2L:8193-8589:+':['FBgn0031208'],
                            'intron_FBgn0031208:2_FBgn0031208:4':['FBgn0031208'],
                            'three_prime_UTR_FBgn0031208:3_737':['FBgn0031208'],
                            'FBgn0031208:4':['FBgn0031208'],
                            'CDS_FBgn0031208:4_737':['FBgn0031208'],
                            'three_prime_UTR_FBgn0031208:4_737':['FBgn0031208'],
                           }

# Same thing for GTF test file . . .
GTF_parent_check_level_1 = {
                            'exon:chr2L:7529-8116:+':['FBtr0300689'],
                            'exon:chr2L:7529-8116:+_1':['FBtr0300690'],
                            'exon:chr2L:8193-9484:+':['FBtr0300689'],
                            'exon:chr2L:8193-8589:+':['FBtr0300690'],
                            'exon:chr2L:8668-9484:+':['FBtr0300690'],
                            'exon:chr2L:10000-11000:-':['transcript_Fk_gene_1'],
                            'exon:chr2L:11500-12500:-':['transcript_Fk_gene_2'],
                            'CDS:chr2L:7680-8116:+':['FBtr0300689'],
                            'CDS:chr2L:7680-8116:+_1':['FBtr0300690'],
                            'CDS:chr2L:8193-8610:+':['FBtr0300689'],
                            'CDS:chr2L:8193-8589:+':['FBtr0300690'],
                            'CDS:chr2L:8668-9276:+':['FBtr0300690'],
                            'CDS:chr2L:10000-11000:-':['transcript_Fk_gene_1'],
                            'FBtr0300689':['FBgn0031208'],
                            'FBtr0300690':['FBgn0031208'],
                            'transcript_Fk_gene_1':['Fk_gene_1'],
                            'transcript_Fk_gene_2':['Fk_gene_2'],
                            'start_codon:chr2L:7680-7682:+':['FBtr0300689'],
                            'start_codon:chr2L:7680-7682:+_1':['FBtr0300690'],
                            'start_codon:chr2L:10000-11002:-':['transcript_Fk_gene_1'],
                            'stop_codon:chr2L:8611-8613:+':['FBtr0300689'],
                            'stop_codon:chr2L:9277-9279:+':['FBtr0300690'],
                            'stop_codon:chr2L:11001-11003:-':['transcript_Fk_gene_1'],
                            }


GTF_parent_check_level_2 = {
                            'exon:chr2L:7529-8116:+':['FBgn0031208'],
                            'exon:chr2L:8193-9484:+':['FBgn0031208'],
                            'exon:chr2L:8193-8589:+':['FBgn0031208'],
                            'exon:chr2L:8668-9484:+':['FBgn0031208'],
                            'exon:chr2L:10000-11000:-':['Fk_gene_1'],
                            'exon:chr2L:11500-12500:-':['Fk_gene_2'],
                            'CDS:chr2L:7680-8116:+':['FBgn0031208'],
                            'CDS:chr2L:8193-8610:+':['FBgn0031208'],
                            'CDS:chr2L:8193-8589:+':['FBgn0031208'],
                            'CDS:chr2L:8668-9276:+':['FBgn0031208'],
                            'CDS:chr2L:10000-11000:-':['Fk_gene_1'],
                            'FBtr0300689':[],
                            'FBtr0300690':[],
                            'transcript_Fk_gene_1':[],
                            'transcript_Fk_gene_2':[],
                            'start_codon:chr2L:7680-7682:+':['FBgn0031208'],
                            'start_codon:chr2L:10000-11002:-':['Fk_gene_1'],
                            'stop_codon:chr2L:8611-8613:+':['FBgn0031208'],
                            'stop_codon:chr2L:9277-9279:+':['FBgn0031208'],
                            'stop_codon:chr2L:11001-11003:-':['Fk_gene_1'],
                           }

expected_feature_counts = {
            'gff3':{'gene':3,
                   'mRNA':4,
                   'exon':6,
                   'CDS':5,
                   'five_prime_UTR':1,
                   'intron':3,
                   'pcr_product':1,
                   'protein':2,
                   'three_prime_UTR':2},
            'gtf':{
                #'gene':3,
                #   'mRNA':4,
                   'CDS':6,
                   'exon':7,
                   'start_codon':3,
                   'stop_codon':3}
            }

expected_features = {'gff3':['gene',
                            'mRNA',
                            'protein',
                            'five_prime_UTR',
                            'three_prime_UTR',
                            'pcr_product',
                            'CDS',
                            'exon',
                            'intron'],
                    'gtf':['gene',
                           'mRNA',
                           'CDS',
                           'exon',
                           'start_codon',
                           'stop_codon']}

