from pathlib import Path
import subprocess
import sys

import gffutils
import pytest


@pytest.mark.parametrize("command", ["children", "parents"])
@pytest.mark.parametrize("explicit_id", [False, True])
def test_relations_default_ids(tmp_path, command, explicit_id):
    data = (
        "synthetic\ttest\tgene\t1\t9\t.\t+\t.\tID=synthetic_gene\n"
        "synthetic\ttest\tmRNA\t1\t9\t.\t+\t.\tID=synthetic_tx;Parent=synthetic_gene\n"
        "synthetic\ttest\texon\t1\t9\t.\t+\t.\tID=synthetic_exon;Parent=synthetic_tx\n"
        "synthetic\ttest\tgene\t20\t29\t.\t+\t.\tID=unrelated_gene\n"
    )
    filename = tmp_path / "synthetic.db"
    db = gffutils.create_db(data, str(filename), from_string=True)
    db.conn.close()
    script = Path(gffutils.__file__).parent / "scripts" / "gffutils-cli"
    args = [sys.executable, str(script), command, str(filename)]
    if explicit_id:
        args.append("synthetic_gene" if command == "children" else "synthetic_exon")
    result = subprocess.run(args, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    ids = {
        gffutils.feature.feature_from_line(line)["ID"][0]
        for line in result.stdout.splitlines()
    }
    expected = {"synthetic_gene", "synthetic_tx", "synthetic_exon"}
    if command == "children" and not explicit_id:
        expected.add("unrelated_gene")
    assert ids == expected
