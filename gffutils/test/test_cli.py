import subprocess as sp
import gffutils
from gffutils import example_filename, create, feature


def test_issue_224():
    fn = gffutils.example_filename("FBgn0031208.gtf")
    sp.run(["gffutils-cli", "create", "--force", fn])
    p = sp.run(
        ["gffutils-cli", "children", fn + ".db", "FBgn0031208"],
        check=True,
        capture_output=True,
        universal_newlines=True,
    )
    assert (
        p.stdout.splitlines()[0]
        == 'chr2L\tgffutils_derived\tgene\t7529\t9484\t.\t+\t.\tgene_id "FBgn0031208";'
    )
