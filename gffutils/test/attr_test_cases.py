"""
Test cases used for the nosetest generator over in parser_test.py.

Each item is a 3-tuple.  The first two items are the attribute string and the
expected parsed dictionary.  The third, if not None, is a reconstructed string
that is acceptable -- this is needed for cases like this:

    name "fgenesh1_pg.C_chr_1000007"; transcriptId 873

where there is not consistency in the quoting of values.  So, in this case, the
following would be an acceptable recontructed string (with quotes around the
873):

    name "fgenesh1_pg.C_chr_1000007"; transcriptId "873"


"""
attrs = [

    # c_elegans_WS199_shortened_gff.txt
    (
        "count=1;gene=amx-2;sequence=SAGE:ggcagagtcttttggca;"
        "transcript=B0019.1",

        {'count': ["1"],
         'gene': ['amx-2'],
         'sequence': ['SAGE:ggcagagtcttttggca'],
         'transcript': ['B0019.1'],
         },

        None,
    ),

    # ensembl_gtf.txt
    (
        'gene_id "Y74C9A.6"; transcript_id "Y74C9A.6"; exon_number "1"; '
        'gene_name "Y74C9A.6"; transcript_name "NR_001477.2";',

        {'gene_id': ['Y74C9A.6'],
         'transcript_id': ['Y74C9A.6'],
         'exon_number': ['1'],
         'gene_name': ['Y74C9A.6'],
         'transcript_name': ['NR_001477.2']},

        None,
    ),

    # F3-unique-3.v2.gff
    (
        'g=A3233312322232122211;i=1;p=1.000;q=23,12,18,17,10,24,19,14,27,9,23'
        ',9,16,20,11,7,8,4,4,14;u=0,0,0,1',

        {'g': ['A3233312322232122211'],
         'i': ['1'],
         'p': ['1.000'],
         'q': ['23',
               '12', '18', '17', '10', '24', '19', '14', '27', '9', '23', '9',
               '16', '20', '11', '7', '8', '4', '4', '14'],
         'u': ['0', '0', '0', '1']},

        None,
    ),


    # glimmer_nokeyval.gff3
    (
        "ID=GL0000006;Name=GL0000006;Lack 3'-end;",

        {'ID': ["GL0000006"],
         'Name': ['GL0000006'],
         "Lack 3'-end": []},

        None,
    ),

    # hybrid1.gff3
    (
        'ID=A00469;Dbxref=AFFX-U133:205840_x_at,Locuslink:2688,Genbank-mRNA:'
        'A00469,Swissprot:P01241,PFAM:PF00103,AFFX-U95:1332_f_at,Swissprot:'
        'SOMA_HUMAN;Note=growth%20hormone%201;Alias=GH1',

        {'ID': ['A00469'],
         'Dbxref': [
             'AFFX-U133:205840_x_at',
             'Locuslink:2688',
             'Genbank-mRNA:A00469',
             'Swissprot:P01241',
             'PFAM:PF00103',
             'AFFX-U95:1332_f_at',
             'Swissprot:SOMA_HUMAN',
         ],
         'Note': ['growth hormone 1'],
         'Alias': ['GH1']},

        'ID=A00469;Dbxref=AFFX-U133:205840_x_at,Locuslink:2688,Genbank-mRNA:'
        'A00469,Swissprot:P01241,PFAM:PF00103,AFFX-U95:1332_f_at,Swissprot:'
        'SOMA_HUMAN;Note=growth hormone 1;Alias=GH1',
    ),

    # jgi_gff2.txt
    #
    # This file is inconsitent with how it quotes values -- integers are not
    # quoted but string values are.  Only way to make this be invariant is to
    # keep track of the "flavor" of each attribute; not sure it's worth the
    # effort / processing time.
    (
        'name "fgenesh1_pg.C_chr_1000007"; transcriptId 873',

        {'name': ['fgenesh1_pg.C_chr_1000007'],
         'transcriptId': ['873']},

        'name "fgenesh1_pg.C_chr_1000007"; transcriptId "873"',
    ),

    # mouse_extra_comma.gff3: extra comma line
    #
    # Note extra empty string in the dictionary's "Parent" field.
    #
    (
        'Name=CDS:NC_000083.5:LOC100040603;Parent=XM_001475631.1,',

        {'Name': ['CDS:NC_000083.5:LOC100040603'],
         'Parent': ['XM_001475631.1', '']},

        None,
    ),

    # mouse_extra_comma.gff3
    #
    # Note the empty ID field.  Compare with the "Lack 3'-end" attribute of
    # glimmer_nokeyval.gff3 above.  Presumably the "Lack 3'-end" field should
    # be interpreted as "True", but an empty ID should be interpreted as "None"
    # or something.
    #
    # Furthermore, the "Lack 3'-end" has no trailing "=", but the "ID" field
    # here does.
    #
    # In both cases, the dictionary entry is simply an empty list; it's just in
    # the reconstruction where things get tricky.
    (
        'ID=;Parent=XM_001475631.1',

        {'ID': [],
         'Parent': ['XM_001475631.1']},

        'ID;Parent=XM_001475631.1',
    ),

    # ncbi_gff3.txt
    (
        'ID=NC_008596.1:speB:unknown_transcript_1;Parent=NC_008596.1:speB;'
        'locus_tag=MSMEG_1072;EC_number=3.5.3.11;note=identified%20by%20mat'
        'ch%20to%20protein%20family%20HMM%20PF00491%3B%20match%20to%20prote'
        'in%20family%20HMM%20TIGR01230;transl_table=11;product=agmatinase;p'
        'rotein_id=YP_885468.1;db_xref=GI:118469242;db_xref=GeneID:4535378;'
        'exon_number=1',

        {'ID': ['NC_008596.1:speB:unknown_transcript_1'],
         'Parent': ['NC_008596.1:speB'],
         'locus_tag': ['MSMEG_1072'],
         'EC_number': ['3.5.3.11'],
         'note': ['identified by match to protein family HMM P'
                  'F00491; match to protein family HMM TIGR01'
                  '230'],
         'transl_table': ['11'],
         'product': ['agmatinase'],
         'protein_id': ['YP_885468.1'],
         'db_xref': ['GI:118469242', 'GeneID:4535378'],
         'exon_number': ['1'],
         },

        'ID=NC_008596.1:speB:unknown_transcript_1;Parent=NC_008596.1:speB;'
        'locus_tag=MSMEG_1072;EC_number=3.5.3.11;note=identified by mat'
        'ch to protein family HMM PF00491%3B match to prote'
        'in family HMM TIGR01230;transl_table=11;product=agmatinase;p'
        'rotein_id=YP_885468.1;db_xref=GI:118469242;db_xref=GeneID:4535378;'
        'exon_number=1',
    ),

    # wormbase_gff2_alt.txt
    #
    (
        'CDS "cr01.sctg102.wum.2.1"',

        {'CDS': ['cr01.sctg102.wum.2.1']},

        None,
    ),

]
