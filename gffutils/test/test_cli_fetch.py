import os
from pathlib import Path
import subprocess
import sys

import gffutils
import pytest


@pytest.mark.parametrize("input_kind", ["database", "gff", "gff_with_database"])
def test_fetch(tmp_path, input_kind):
    line = "synthetic\ttest\tgene\t1\t9\t.\t+\t.\tID=synthetic_gene"
    gff = tmp_path / "synthetic.gff3"
    gff.write_text(line + "\n")
    filename = gff
    if input_kind != "gff":
        filename = Path(str(gff) + "..db")
        db = gffutils.create_db(str(gff), str(filename))
        db.conn.close()
        if input_kind == "gff_with_database":
            filename = gff

    script = Path(gffutils.__file__).parent / "scripts" / "gffutils-cli"
    env = dict(os.environ, TMPDIR=str(tmp_path), TEMP=str(tmp_path), TMP=str(tmp_path))
    result = subprocess.run(
        [sys.executable, str(script), "fetch", str(filename), "synthetic_gene"],
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[-1] == line


def test_fetch_missing_id(tmp_path):
    filename = tmp_path / "synthetic.db"
    db = gffutils.create_db(
        "synthetic\ttest\tgene\t1\t9\t.\t+\t.\tID=synthetic_gene",
        str(filename),
        from_string=True,
    )
    db.conn.close()
    script = Path(gffutils.__file__).parent / "scripts" / "gffutils-cli"
    result = subprocess.run(
        [sys.executable, str(script), "fetch", str(filename), "missing_gene"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == "missing_gene not found\n"
