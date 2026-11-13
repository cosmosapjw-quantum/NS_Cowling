"""C6C provider-package builder: pinned identity, determinism, C6B read-only."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "build_sharp_maxwell_c6c_provider_package.py"
C6B = ROOT / "artifacts/sharp_maxwell/C6B_COMPOSITION_SEMANTICS_20261110T100954Z"
AUTH = "a09ec30d74b363652d2d1aefe990557c468679dd"


def _have_authority_commit() -> bool:
    r = subprocess.run(["git", "-C", str(ROOT), "cat-file", "-e", AUTH + "^{commit}"],
                       capture_output=True)
    return r.returncode == 0 and C6B.is_dir()


pytestmark = pytest.mark.skipif(not _have_authority_commit(),
                                reason="source authority commit or C6B evidence not available")


def _run(out):
    subprocess.run([sys.executable, str(SCRIPT), "--repo", str(ROOT), "--out-dir", str(out)],
                   check=True, capture_output=True)
    return json.loads((out / "aggregate_statistics.json").read_text())


def test_identity_counts_and_key_statistics(tmp_path):
    rep = _run(tmp_path / "o")
    spg, gp = rep["sources"]["SPG_M1"], rep["sources"]["GPPVA_FSU2"]
    assert spg["git_blob_sha1"] == "ab79da517271bea1f7de37294e87e019de611699" and spg["rows"] == 365
    assert gp["git_blob_sha1"] == "fa1ad2b809eb1c95825e2ee6eb2650cf3f248693" and gp["rows"] == 269
    s7 = spg["aggregates"]["Iphase=7"]
    g7 = gp["aggregates"]["Iphase=7"]
    assert s7["n_rows"] == 148 and s7["n_YI_equals_Ye_plus_Ymu_1e-9"] == 148
    assert s7["n_C3_C4_implied_consistent_1e-3"] == 0
    assert g7["n_rows"] == 52 and g7["n_YI_equals_Ye_plus_Ymu_1e-9"] == 52
    assert g7["n_C3_C4_implied_consistent_1e-3"] == 52
    assert abs(s7["max_abs_r_C3_direct"] - 0.042790704812013594) < 1e-15
    assert abs(g7["max_abs_r_C4_direct"] - 0.98856) < 1e-12


def test_deterministic(tmp_path):
    _run(tmp_path / "a")
    _run(tmp_path / "b")
    for f in sorted((tmp_path / "a").iterdir()):
        assert f.read_bytes() == (tmp_path / "b" / f.name).read_bytes()


def test_refuses_c6b_directory():
    r = subprocess.run([sys.executable, str(SCRIPT), "--repo", str(ROOT), "--out-dir", str(C6B)],
                       capture_output=True)
    assert r.returncode != 0
