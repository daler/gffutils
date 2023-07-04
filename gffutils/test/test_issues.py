"""
Tests for specific issues and pull requests
"""


import os
import tempfile
import difflib
from textwrap import dedent
import gffutils
from gffutils import feature
from gffutils import merge_criteria as mc

import pytest

def test_issue_79():
    gtf = gffutils.example_filename("keep-order-test.gtf")
    db = gffutils.create_db(
        gtf,
        "tmp.db",
        disable_infer_genes=False,
        disable_infer_transcripts=False,
        id_spec={"gene": "gene_id", "transcript": "transcript_id"},
        merge_strategy="create_unique",
        keep_order=True,
        force=True,
    )

    exp = open(gtf).read()
    obs = "\n".join([str(i) for i in db.all_features()])
    exp_1 = exp.splitlines(True)[0].strip()
    obs_1 = obs.splitlines(True)[0].strip()
    print("EXP")
    print(exp_1)
    print("OBS")
    print(obs_1)
    print("DIFF")
    print("".join(difflib.ndiff([exp_1], [obs_1])))
    assert obs_1 == exp_1


def test_issue_82():
    # key-val separator is inside an unquoted attribute value
    x = (
        "Spenn-ch12\tsgn_markers\tmatch\t2621812\t2622049\t.\t+\t.\t"
        "Alias=SGN-M1347;ID=T0028;Note=marker name(s): T0028 SGN-M1347 |identity=99.58|escore=2e-126"
    )
    y = feature.feature_from_line(x)
    assert y.attributes["Note"] == [
        "marker name(s): T0028 SGN-M1347 |identity=99.58|escore=2e-126"
    ]

    gffutils.create_db(gffutils.example_filename("keyval_sep_in_attrs.gff"), ":memory:")


def test_issue_85():
    # when start or stop was empty, #85 would fail Should now work with
    # blank fields
    f = feature.feature_from_line("\t".join([""] * 9))

    # or with "." placeholders
    f = feature.feature_from_line("\t".join(["."] * 9))


def test_issue_105():
    fn = gffutils.example_filename("FBgn0031208.gtf")
    home = os.path.expanduser("~")
    newfn = os.path.join(home, ".gffutils.test")
    with open(newfn, "w") as fout:
        fout.write(open(fn).read())
    f = gffutils.iterators.DataIterator(newfn)
    for i in f:
        pass
    os.unlink(newfn)


def test_issue_107():
    s = dedent(
        """
        chr1\t.\tgene\t10\t15\t.\t+\t.\tID=b;
        chr1\t.\tgene\t1\t5\t.\t-\t.\tID=a;
        chr2\t.\tgene\t25\t50\t.\t-\t.\tID=c;
        chr2\t.\tgene\t55\t60\t.\t-\t.\tID=d;
        """
    )
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, "w") as fout:
        fout.write(s + "\n")
    db = gffutils.create_db(tmp, ":memory:")
    interfeatures = list(
        db.interfeatures(db.features_of_type("gene", order_by=("seqid", "start")))
    )
    assert [str(i) for i in interfeatures] == [
        "chr1\tgffutils_derived\tinter_gene_gene\t6\t9\t.\t.\t.\tID=a-b;",
        "chr2\tgffutils_derived\tinter_gene_gene\t51\t54\t.\t-\t.\tID=c-d;",
    ]


def test_issue_119():
    # First file has these two exons with no ID:
    #
    #   chr2L FlyBase exon  8193  8589  .  +  .  Parent=FBtr0300690
    #   chr2L FlyBase exon  7529  8116  .  +  .  Name=CG11023:1;Parent=FBtr0300689,FBtr0300690
    #
    db0 = gffutils.create_db(gffutils.example_filename("FBgn0031208.gff"), ":memory:")

    # And this one, a bunch of reads with no IDs anywhere
    db1 = gffutils.create_db(
        gffutils.example_filename("F3-unique-3.v2.gff"), ":memory:"
    )

    # When db1 is updated by db0
    db2 = db1.update(db0)
    assert (
        db2._autoincrements == db1._autoincrements == {"exon": 2, "read": 112}
    ), db2._autoincrements

    assert len(list(db0.features_of_type("exon"))) == 6

    # Now we update that with db0 again
    db3 = db2.update(db0, merge_strategy="replace")

    # Using the "replace" strategy, we should have only gotten another 2 exons
    assert len(list(db3.features_of_type("exon"))) == 8

    # Make sure that the autoincrements for exons jumped by 2
    assert (
        db2._autoincrements == db3._autoincrements == {"exon": 4, "read": 112}
    ), db2._autoincrements

    # More isolated test, merging two databases each created from the same file
    # which itself contains only a single feature with no ID.
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, "w") as fout:
        fout.write("chr1\t.\tgene\t10\t15\t.\t+\t.\t\n")

    db4 = gffutils.create_db(tmp, tmp + ".db")
    db5 = gffutils.create_db(tmp, ":memory:")

    assert db4._autoincrements == {"gene": 1}
    assert db5._autoincrements == {"gene": 1}

    db6 = db4.update(db5)

    db7 = gffutils.FeatureDB(db4.dbfn)

    # both db4 and db6 should now have the same, updated autoincrements because
    # they both point to the same db.
    assert db6._autoincrements == db4._autoincrements == {"gene": 2}

    # But db5 was created independently and should have unchanged autoincrements
    assert db5._autoincrements == {"gene": 1}

    # db7 was created from the database pointed to by both db4 and db6. This
    # tests that when a FeatureDB is created it should have the
    # correctly-updated autoincrements read from the db
    assert db7._autoincrements == {"gene": 2}


def test_pr_131():
    db = gffutils.create_db(gffutils.example_filename("FBgn0031208.gff"), ":memory:")

    # previously would raise ValueError("No lines parsed -- was an empty
    # file provided?"); now just does nothing
    db2 = db.update([])


def test_pr_133():
    # Previously, merge_attributes would not deep-copy the values from the
    # second dict, and when the values are then modified, the second dict is
    # unintentionally modified.
    d1 = {"a": [1]}
    d2 = {"a": [2]}
    d1a = {"a": [1]}
    d2a = {"a": [2]}
    d3 = gffutils.helpers.merge_attributes(d1, d2)
    assert d1 == d1a, d1
    assert d2 == d2a, d2


def test_pr_139():
    db = gffutils.create_db(gffutils.example_filename("FBgn0031208.gff"), ":memory:")
    exons = list(db.features_of_type("exon"))
    inter = list(db.interfeatures(exons))

    # previously, the first exon's attributes would show up in subsequent merged features
    first_name = exons[0].attributes["Name"][0]
    for i in inter[1:]:
        if "Name" in i.attributes:
            assert first_name not in i.attributes["Name"], str(i)


def test_pr_144():
    # previously this would fail with:
    #   UnboundLocalError: local variable 'part' referenced before assignment
    f = gffutils.Feature(attributes={"a": [""]})

    # Make sure everything got converted correctly
    assert f.attributes["a"] == [""]
    assert str(f) == ".	.	.	.	.	.	.	.	a"
    g = gffutils.feature.feature_from_line(str(f))
    assert g == f


def test_pr_172():
    line = (
        "NC_049222.1\tGnomon\tgene\t209085\t282880\t.\t-\t.\t"
        'gene_id "ENPP1_3"; transcript_id ""; db_xref "GeneID:100856150";'
        'db_xref "VGNC:VGNC:40374"; gbkey "Gene"; gene "ENPP1"; '
        'gene_biotype "protein_coding";\n'
    )
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, "w") as fout:
        fout.write(line)
    db = gffutils.create_db(tmp, ":memory:")


def test_pr_171():
    q = gffutils.parser.Quoter()
    assert q.__missing__("\n") == "%0A"
    assert q.__missing__("a") == "a"

    assert q.__missing__("") == ""


def test_issue_129():

    # thanks @Brunox13 for the detailed notes on #129

    line = 'chr1\tdemo\tstart_codon\t69091\t69093\t.\t+\t.\tgene_id "demo";\n'
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, "w") as fout:
        fout.write(line)
    db = gffutils.create_db(tmp, ":memory:")

    # ASCII art to visualize each test (coords are along the top, from 69087 to
    # 69090). The tests slide a 4-bp region over the original 3-bp start codon.

    # 7 8 9 0 1 2 3 4 5 6 7
    #         | | |         Orig feature
    # | | | |               Test feature
    res = list(db.region(region=("chr1", 69087, 69090), featuretype="start_codon"))
    assert len(res) == 0

    # NOTE: prior to #162, this did not return anything
    # 7 8 9 0 1 2 3 4 5 6 7
    #         | | |         Orig feature
    #   | | | |             Test feature
    res = list(db.region(region=("chr1", 69088, 69091), featuretype="start_codon"))
    assert len(res) == 1

    # 7 8 9 0 1 2 3 4 5 6 7
    #         | | |         Orig feature
    #     | | | |           Test feature
    res = list(db.region(region=("chr1", 69089, 69092), featuretype="start_codon"))
    assert len(res) == 1

    # 7 8 9 0 1 2 3 4 5 6 7
    #         | | |         Orig feature
    #       | | | |         Test feature
    res = list(db.region(region=("chr1", 69090, 69093), featuretype="start_codon"))
    assert len(res) == 1

    # 7 8 9 0 1 2 3 4 5 6 7
    #         | | |         Orig feature
    #         | | | |       Test feature
    res = list(db.region(region=("chr1", 69091, 69094), featuretype="start_codon"))
    assert len(res) == 1

    # 7 8 9 0 1 2 3 4 5 6 7
    #         | | |         Orig feature
    #           | | | |     Test feature
    res = list(db.region(region=("chr1", 69092, 69095), featuretype="start_codon"))
    assert len(res) == 1

    # NOTE: priro to #162, this did not return anything
    # 7 8 9 0 1 2 3 4 5 6 7
    #         | | |         Orig feature
    #             | | | |   Test feature
    res = list(db.region(region=("chr1", 69093, 69096), featuretype="start_codon"))
    assert len(res) == 1

    # 7 8 9 0 1 2 3 4 5 6 7
    #         | | |         Orig feature
    #               | | | | Test feature
    res = list(db.region(region=("chr1", 69094, 69097), featuretype="start_codon"))
    assert len(res) == 0


def test_issue_128():
    # In #128, some lines had separators of "; " and some with ";". The first
    # one in the file would win. Now the detection pays more attention to lines
    # with more attributes to make it work properly
    gff = gffutils.example_filename('gms2_example.gff3')
    db = gffutils.create_db(gff, ":memory:", force=True)
    expected = {
        'ID': ['1'],
        'Parent': ['gene_1'],
        'gene_type': ['native'],
        'partial': ['11'],
        'gc': ['33'], 
        'length': ['363'],
    }
    assert dict(db['1'].attributes) == expected


def test_issue_157():
    # With the merge overhaul, children_bp incorrectly still used ignore_strand.
    db = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'), ":memory:")
    gene = next(db.features_of_type('gene'))
    children = list(db.children(gene, featuretype='exon'))

    # Modify the last one to have a different strand so we can test the
    # ignore_strand argument.
    children[-1].strand = '-'
    db.update(children[-1:], merge_strategy='replace')

    # and, since updating has been problematic in the past, double-check again
    # that the strand is changed in the db.
    assert list(db.children(gene, featuretype='exon'))[-1].strand == '-'
    cbp1 = db.children_bp(gene, child_featuretype='exon')

    # Previously this would give:
    #   TypeError: merge() got an unexpected keyword argument 'ignore_strand'
    #
    # Now changing to ValueError and suggesting a fix. 
    with pytest.raises(ValueError):
        db.children_bp(gene, child_featuretype='exon', merge=True, ignore_strand=True)
    with pytest.raises(ValueError):
        db.children_bp(gene, ignore_strand=True, nonexistent=True)
    with pytest.raises(TypeError):
        db.children_bp(gene, nonexistent=True)

    # The way to do it now is the following (we can omit the mc.feature_type
    # since we're preselecting for exons anyway):
    db.children_bp(gene, child_featuretype='exon', merge=True, merge_criteria=(mc.overlap_end_inclusive))


def test_issue_159():
    db = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'), ":memory:")
    fasta = gffutils.example_filename('dm6-chr2L.fa')
    for transcript, seq in gffutils.helpers.canonical_transcripts(db, fasta):
        pass


def test_issue_164():
    # Something strange with the original annotation, but seems fine at least
    # after pasting in the offending genes from the GitHub comments.
    db = gffutils.create_db(
        gffutils.example_filename('sharr.gtf'),
        ':memory:',
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'},
        merge_strategy='create_unique',
        keep_order=True)


def test_issue_166():
    # Added the new FeatureDB.seqids() method.
    db = gffutils.create_db(gffutils.example_filename('nonascii'), ':memory:')
    seqs = list(db.seqids())
    assert seqs == ['2L', '2R', '3L', '3R', 'X'], seqs


def test_issue_167():
    # Previously was causing sqlite3.InterfaceError
    db = gffutils.create_db(gffutils.example_filename('issue167.gff'), ':memory:')


def test_issue_174():
    db = gffutils.create_db(
        gffutils.example_filename('issue174.gtf'),
        ':memory:',
        merge_strategy='warning',
    )
    introns = [f for f in db.create_introns()]
    observed = [i.attributes['exon_number'] for i in introns]
    assert observed[7] == ['8', '9']
    assert observed[8] == ['10', '9']
    assert observed[9] == ['10', '11']

    # Now do the same thing, but with the new numeric_sort arg
    introns = [f for f in db.create_introns(numeric_sort=True)]
    observed = [i.attributes['exon_number'] for i in introns]
    assert observed[7] == ['8', '9']
    # This should be fixed:
    assert observed[8] == ['9', '10'] 
    assert observed[9] == ['10', '11']

def test_issue_197():

    # Previously this would fail with ValueError due to using the stop position
    # of the last item on the previous chrom as the start position.

    db = gffutils.create_db(gffutils.example_filename('issue_197.gff'), ':memory:', merge_strategy='error')
    genes = list(db.features_of_type('gene'))
    igss = list( db.interfeatures(genes,new_featuretype='intergenic_space') )


    # Prior to PR #219, multiple IDs could be created by interfeatures, which
    # in turn was patched here by providing the transform to db.update. With
    # #219, this ends up being a no-op because ID is a single value by the time
    # it gets to the transform function.
    #
    # However, keeping the test as-is to ensure backward-compatibility.
    def transform(f):
        f['ID'] = [ '-'.join(f.attributes['ID']) ]
        return f

    db = db.update(igss, transform=transform,  merge_strategy='error')

    obs = list(db.features_of_type('intergenic_space'))
    for i in obs:
        print(i)

    assert [str(i) for i in obs] == [
        'tig00000492\tgffutils_derived\tintergenic_space\t47236\t47350\t.\t-\t.\tID=ctg492.gene0001-ctg492.gene0002;Name=gene0001,gene0002',
        'tig00000492\tgffutils_derived\tintergenic_space\t48257\t49999\t.\t-\t.\tID=ctg492.gene0002-gene0;Name=gene0002',
        'tig00000492\tgffutils_derived\tintergenic_space\t50050\t50054\t.\t-\t.\tID=gene3-gene4',
        'tig00000492\tgffutils_derived\tintergenic_space\t50071\t50071\t.\t-\t.\tID=gene4-gene5',
        'tig00000492\tgffutils_derived\tintergenic_space\t50076\t50089\t.\t-\t.\tID=gene5-gene6',
    ]

def test_issue_198():
    line = 'NC_000001.11	BestRefSeq	gene	14362	29370	.	-	.	gene_id "WASH7P"; transcript_id ""; db_xref "GeneID:653635"; db_xref "HGNC:HGNC:38034"; description "WASP family homolog 7, pseudogene"; gbkey "Gene"; gene "WASH7P"; gene_biotype "transcribed_pseudogene"; gene_synonym "FAM39F"; gene_synonym "WASH5P"; pseudo "true";'

    # Original issue #198 is that this fails with:
    #
    #   gffutils.exceptions.AttributeStringError: Internally inconsistent
    #   attributes formatting: some have repeated keys, some do not.
    #
    # This is because the dialect inference sees the two db_xref keys, and
    # correctly assumes the dialect uses repeated keys rather than
    # multiple, comma-separated values -- but there's a comma in the
    # description.
    #
    # So we need to figure out the best way of interpreting a comma in cases
    # like this. It seems like the best solution is to assume that the presence
    # of repeated keys always wins.
    f = feature.feature_from_line(line)

    assert f.attributes['description'] == ['WASP family homolog 7, pseudogene']

    # If we remove one of the db_xref keys, then the parser sees the comma and
    # figures it's a multivalue key.
    line = 'NC_000001.11	BestRefSeq	gene	14362	29370	.	-	.	gene_id "WASH7P"; transcript_id ""; db_xref "GeneID:653635"; description "WASP family homolog 7, pseudogene"; gbkey "Gene"; gene "WASH7P"; gene_biotype "transcribed_pseudogene"; gene_synonym "FAM39F"; gene_synonym "WASH5P"; pseudo "true";'
    f = feature.feature_from_line(line)

    # Previous result, note leading space --------------------------->| |
    # assert f.attributes['description'] == ['WASP family homolog 7', ' pseudogene']
    assert f.attributes['description'] == ['WASP family homolog 7, pseudogene']

    # But removing that space before "pseudogene" means it's interpreted as
    # a multivalue attribute
    line = 'NC_000001.11	BestRefSeq	gene	14362	29370	.	-	.	gene_id "WASH7P"; transcript_id ""; db_xref "GeneID:653635"; description "WASP family homolog 7,pseudogene"; gbkey "Gene"; gene "WASH7P"; gene_biotype "transcribed_pseudogene"; gene_synonym "FAM39F"; gene_synonym "WASH5P"; pseudo "true";'
    f = feature.feature_from_line(line)
    assert f.attributes['description'] == ['WASP family homolog 7', 'pseudogene']

    # Confirm behavior of corner cases like a trailing comma
    line = "chr17	RefSeq	CDS	6806527	6806553	.	+	0	Name=CDS:NC_000083.5:LOC100040603;Parent=XM_001475631.1,"
    f = feature.feature_from_line(line)
    assert f.attributes['Parent'] == ['XM_001475631.1', '']


def test_issue_207():

    def _check(txt, expected_keys, dialect_trailing_semicolon):
        db = gffutils.create_db(txt.replace(' ', '\t'), ':memory:', from_string=True)
        assert [list(f.attributes.keys()) for f in db.all_features()] == expected_keys
        assert db.dialect['trailing semicolon'] == dialect_trailing_semicolon

    # All lines have trailing semicolon
    _check(
        txt=dedent("""\
        chr1 AUGUSTUS gene 68330 73621 1 - . ID=g1903;
        chr1 AUGUSTUS mRNA 68330 73621 1 - . ID=g1903.t1;Parent=g1903;
        chr1 Pfam protein_match 73372 73618 1 - . ID=g1903.t1.d1;Parent=g1903.t1;
        chr1 Pfam protein_hmm_match 73372 73618 1 - . ID=g1903.t1.d1.1;Parent=g1903.t1.d1;
        """),
        expected_keys = [
            ['ID'],
            ['ID', 'Parent'],
            ['ID', 'Parent'],
            ['ID', 'Parent'],
        ],
        dialect_trailing_semicolon=True
    )

    # First two lines have trailing semicolon. However, the heuristics of
    # dialect selection, which favor attributes with more values (assuming more
    # information), decides that this file does NOT have trailing semicolons.
    _check(
        txt=dedent("""\
        chr1 AUGUSTUS gene 68330 73621 1 - . ID=g1903;
        chr1 AUGUSTUS mRNA 68330 73621 1 - . ID=g1903.t1;Parent=g1903;
        chr1 Pfam protein_match 73372 73618 1 - . ID=g1903.t1.d1;Parent=g1903.t1
        chr1 Pfam protein_hmm_match 73372 73618 1 - . ID=g1903.t1.d1.1;Parent=g1903.t1.d1
        """),
        expected_keys = [
            ['ID', ''],
            ['ID', 'Parent', ''],
            ['ID', 'Parent'],
            ['ID', 'Parent'],
        ],
        dialect_trailing_semicolon=False,
    )

    # APPARENTLY INCONSISTENT: The only thing difference here is that the
    # Parent attribute has been removed, otherwise matches above (first two
    # have trailing semicolon). But now there are no empty keys.
    #
    # This is expected behavior, because there are no attributes with more keys
    # as above to give higher weight, and to break the tie between with and
    # without trailing semicolon, falls back to first dialect observed.
    _check(
        txt=dedent("""\
        chr1 AUGUSTUS gene 68330 73621 1 - . ID=g1903;
        chr1 AUGUSTUS mRNA 68330 73621 1 - . ID=g1903.t1;
        chr1 Pfam protein_match 73372 73618 1 - . ID=g1903.t1.d1
        chr1 Pfam protein_hmm_match 73372 73618 1 - . ID=g1903.t1.d1.1
        """),
        expected_keys=[
            ['ID'],
            ['ID'],
            ['ID'],
            ['ID']
        ],
        dialect_trailing_semicolon=True,
    )

    # We can convince the heuristics to think there should be NO trailing
    # semicolon by giving one more line as evidence. Only difference is from
    # above is the last line.
    _check(
        txt=dedent("""\
        chr1 AUGUSTUS gene 68330 73621 1 - . ID=g1903;
        chr1 AUGUSTUS mRNA 68330 73621 1 - . ID=g1903.t1;
        chr1 Pfam protein_match 73372 73618 1 - . ID=g1903.t1.d1
        chr1 Pfam protein_hmm_match 73372 73618 1 - . ID=g1903.t1.d1.1
        chr1 Pfam protein_hmm_match 73372 73618 1 - . ID=g1904.t1.d1.1
        """),
        expected_keys=[
            ['ID', ''],
            ['ID', ''],
            ['ID'],
            ['ID'],
            ['ID'],
        ],
        dialect_trailing_semicolon=False,
    )


    # Again seems inconsistent at first, but heuristics break ties by
    # preferring first dialect, which here is no trailing semicolon.
    _check(
        txt=dedent("""\
        chr1 AUGUSTUS gene 68330 73621 1 - . ID=g1903
        chr1 AUGUSTUS mRNA 68330 73621 1 - . ID=g1903.t1
        chr1 Pfam protein_match 73372 73618 1 - . ID=g1903.t1.d1;
        chr1 Pfam protein_hmm_match 73372 73618 1 - . ID=g1903.t1.d1.1;
        """),
        expected_keys=[
            ['ID'],
            ['ID'],
            ['ID', ''],
            ['ID', '']
        ],
        dialect_trailing_semicolon=False,
    )


def test_issue_213():
    # GFF header directives seem to be not parsed when building a db from
    # a file, even though it seems to work fine from a string.
    data = dedent(
        """
    ##gff-version 3
    .	.	.	.	.	.	.	.
    .	.	.	.	.	.	.	.
    .	.	.	.	.	.	.	.
    .	.	.	.	.	.	.	.
    """
    )

    # Ensure directives are parsed from DataIterator
    it = gffutils.iterators.DataIterator(data, from_string=True)
    assert it.directives == ["gff-version 3"]


    # Ensure they're parsed into the db from a string
    db = gffutils.create_db(data, dbfn=":memory:", from_string=True, verbose=False)
    assert db.directives == ["gff-version 3"], db.directives

    # Ensure they're parsed into the db from a file
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, "w") as fout:
        fout.write(data + "\n")
    db = gffutils.create_db(tmp, ":memory:")
    assert db.directives == ["gff-version 3"], db.directives
    assert len(db.directives) == 1

    # Ensure they're parsed into the db from a file, and going to a file (to
    # exactly replicate example in #213)
    db = gffutils.create_db(tmp, dbfn='issue_213.db', force=True)
    assert db.directives == ["gff-version 3"], db.directives
    assert len(db.directives) == 1
