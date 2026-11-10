"""Tests for the read-only C6B composition-semantics audit helpers (synthetic rows plus pinned real rows)."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("c6b_audit", ROOT / "scripts" / "audit_sharp_maxwell_composition_semantics.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

OUTER = "0 1 1 0 4 10 0.00000 11 0.00000 0 0.46429 1 0.00000 1 1 56.000 26.000 1.00000"
SPG_INNER = "0 218 1 7 4 10 0.00000 11 0.00000 0 0.06087 1 0.00000 1 1 82.107 24.387 0.06087"
GPPVA_INNER = "1 218 1 7 4 10 0.00000 11 0.00000 0 0.04865 1 0.00000 1 1 740.000 36.000 0.04865"


def row(line):
    return audit.candidate_row(audit.parse_compo_line(line, 1))


def test_git_blob_sha1_known_vector():
    assert audit.git_blob_sha1(b"") == "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391"


def test_parse_grammar_fields():
    r = audit.parse_compo_line(OUTER, 1)
    assert (r["iT"], r["inb"], r["iYq"], r["Iphase"], r["Npairs"], r["Nquad"]) == (0, 1, 1, 0, 4, 1)
    assert [p[0] for p in r["pairs"]] == [10, 11, 0, 1]
    assert r["quads"] == [(1, 56.0, 26.0, 1.0)]


@pytest.mark.parametrize("bad", [OUTER + " 9", "0 1 1 0"])
def test_parse_rejects_trailing_or_short(bad):
    with pytest.raises(ValueError):
        audit.parse_compo_line(bad, 1)


def test_outer_crust_direct_closes_standard_does_not():
    c = row(OUTER)
    assert c["direct_C3_pass"] and c["direct_C4_pass"]
    assert not c["standard_C4_pass"] and c["Yh_standard_Abar_YI"] == 56.0


def test_spg_inner_row_c3_c4_implied_values_incompatible():
    c = row(SPG_INNER)
    assert c["Yh_C4_implied"] == 1.0
    assert abs(c["Yh_C3_implied"] - 0.204939) < 1e-6
    assert abs(c["C3_vs_C4_implied_gap"]) > 0.7
    assert not c["direct_C3_pass"] and not c["direct_C4_pass"]


def test_gppva_inner_row_implied_values_agree_but_raw_field_does_not():
    c = row(GPPVA_INNER)
    assert abs(c["Yh_C3_implied"] - 1.000028) < 1e-6
    assert abs(c["C3_vs_C4_implied_gap"]) < 1e-3
    assert abs(c["YI_raw"] - 1.0) > 0.9


def test_raw_field_equals_charge_fraction_on_inner_rows():
    for line in (SPG_INNER, GPPVA_INNER):
        assert abs(row(line)["YI_minus_Yq"]) < 1e-12


def test_no_phase_branch_in_candidate_row():
    # identical numerics must give identical candidates regardless of Iphase label
    a = row(SPG_INNER)
    b = row(SPG_INNER.replace(" 7 4 10", " 0 4 10", 1))
    assert a == b


@pytest.mark.skipif(not (ROOT / "eos_generator" / "eos_data" / "SPG_M1" / "eos.compo").exists(), reason="raw source absent")
def test_pinned_real_sources_audit():
    manifest = json.loads((ROOT / "data" / "sharp_maxwell_65_cohort_manifest.v1.json").read_text())["source_catalog"]["crust_models"]
    for sid, n, n7 in (("SPG_M1", 365, 148), ("GPPVA_FSU2", 269, 52)):
        expected = {f: v["git_blob_sha"] for f, v in manifest[sid]["files"].items()}
        res = audit.audit_source(sid, ROOT / "eos_generator" / "eos_data" / sid, expected)
        assert res["checks"]["rows"] == n and res["checks"]["compo_keys_equal_thermo_keys_in_order"]
        assert sum(1 for r in res["rows"] if r["Iphase"] == 7) == n7
        assert all(abs(r["YI_minus_Yq"]) < 1e-9 for r in res["rows"] if r["Iphase"] == 7)
