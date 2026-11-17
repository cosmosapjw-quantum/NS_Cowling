"""Delivery contract tests; fixture receipts do not claim physical admission."""
from __future__ import annotations

import csv
import importlib.util
import io
import json
from pathlib import Path
import zipfile

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/package_maxwell_three_families_20261117.py"
SPEC = importlib.util.spec_from_file_location("maxwell_packaging", SCRIPT)
packager = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(packager)


def table(directory, densities, pressure_scale=1.0):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "eos.nb").write_text("1\n" + str(len(densities)) + "\n" + "\n".join(map(str, densities)) + "\n")
    (directory / "eos.t").write_text("1\n1\n0.0\n")
    (directory / "eos.yq").write_text("1\n1\n0.0\n")
    rows = ["939.0 938.0"]
    for i, n in enumerate(densities, 1):
        q1 = pressure_scale * 2 * n
        q3 = (900 + 2 * n + q1) / 939 - 1
        q7 = (900 + 2 * n) / 939 - 1
        rows.append(f"1 {i} 1 {q1:.17g} 0 {q3:.17g} 0 0 {q7:.17g} {q7:.17g} 0")
    (directory / "eos.thermo").write_text("\n".join(rows) + "\n")


def make_fixture(root):
    cores = ["APR_APR", "BBB_BHF2", "GM_GM1", "RG_KDE0v1", "RG_SLY4"]
    crusts = [f"SPG_M{i}" for i in range(1, 6)] + ["GPPVA_" + s for s in ("DD2", "DDME2", "TW", "FSU2", "FSU2R", "FSU2H", "NL3wrL55", "TM1e")]
    records = []
    for family in packager.FAMILIES:
        for model in cores + crusts:
            original = root / "inputs" / family / model
            prepared = root / "prepared" / family / model
            table(original, [.01, .05, .2, 1.0])
            reconstruction = family == "public_reconstructed"
            table(prepared / "compose", [.05, .2] if reconstruction else [.01, .05, .2, 1.0], 1.1 if reconstruction else 1.0)
            record = {"id": family + ":" + model, "family": family, "model": model,
                "role": "core" if model in cores else "crust", "source_path": str(original.relative_to(root)),
                "source_table_sha256": {p.name: packager.sha(p) for p in original.iterdir()},
                "compose_path": str((prepared / "compose").relative_to(root)),
                "compose_sha256": {p.name: packager.sha(p) for p in (prepared / "compose").iterdir()},
                "selection": {"source_start_zero_based": 1 if reconstruction else 0,
                              "source_stop_exclusive": 3 if reconstruction else 4},
                "state": "FINITE_COMPONENTS_AVAILABLE"}
            packager.save(prepared / "SOURCE.json", record)
            records.append(record)
    manifest = {"families": list(packager.FAMILIES), "core_ids": cores, "crust_ids": crusts, "source_records": records}
    manifest_path = root / "prepared/MANIFEST.json"
    packager.save(manifest_path, manifest)
    cells = []
    for family in packager.FAMILIES:
        for i, crust in enumerate(crusts):
            for j, core in enumerate(cores):
                ident = crust + "__" + core
                path = root / "campaign" / family / ident
                result = []
                if i == j == 0:
                    p = path / "candidate_000"
                    table(p / "compose", [.08, .1, .11, .5])
                    candidate = {"status": packager.SUCCESS, "kernel_invoked": True,
                        "admission": {"status": "ADMITTED", "fixture_note": "metadata-contract fixture only"},
                        "serialized_native_verification": {"overall": "PASS"},
                        "output_path": str((p / "compose").relative_to(root)),
                        "output_sha256": {q.name: packager.sha(q) for q in (p / "compose").iterdir()},
                        "density_fm3": [.08, .5], "transition": {"n_low": .1, "n_high": .11},
                        "scope": "synthetic metadata fixture, not actual numerical validation",
                        "crust_component": {"id": "component_000"}, "core_component": {"id": "component_000"}}
                    packager.save(p / "RESULT.json", candidate)
                    result = [candidate]
                elif i == j == 1:
                    p = path / "candidate_000"
                    table(p / "compose", [.08, .1, .11, .5])  # Failed output must never enter accepted payload.
                    candidate = {"status": "REJECTED", "kernel_invoked": True, "category": "TEST_REJECTION"}
                    packager.save(p / "RESULT.json", candidate)
                    result = [candidate]
                cell = {"family": family, "id": ident, "crust": crust, "core": core,
                    "crust_source_id": family + ":" + crust, "core_source_id": family + ":" + core,
                    "status": packager.SUCCESS if i == j == 0 else "NO_MATCH_ON_ADMITTED_COMPONENTS" if result else "SOURCE_HAS_NO_ADMITTED_COMPONENT_PAIR",
                    "generated_count": int(i == j == 0), "kernel_invocations": len(result),
                    "root_search_domain_is_complete_source": False, "candidate_results": result,
                    "source_states": {"crust": "FINITE_COMPONENTS_AVAILABLE", "core": "FINITE_COMPONENTS_AVAILABLE"},
                    "source_interval_failures": {"crust": {}, "core": {}}}
                packager.save(path / "CELL.json", cell)
                cells.append(cell)
    campaign = {"same_family_only": True, "prepared_manifest": "prepared/MANIFEST.json",
        "prepared_manifest_sha256": packager.sha(manifest_path), "cells": cells,
        "summary": {f: {"generated_tables": 1} for f in packager.FAMILIES}}
    path = root / "campaign/CAMPAIGN.json"
    packager.save(path, campaign)
    return path


def test_complete_census_exact_bytes_and_separate_diagnostic_views(tmp_path):
    campaign = make_fixture(tmp_path)
    receipt = packager.package(campaign, tmp_path / "delivery", root=tmp_path, figure=False)
    catalog = packager.load(tmp_path / "delivery/CATALOG.json")
    assert len(catalog["cells"]) == 195
    assert len(catalog["tables"]) == receipt["accepted_table_count"] == 3
    assert not receipt["source_retaining_csv_admitted_as_EOS"]
    for archive in receipt["archives"]:
        assert archive["summary"]["model_pairs"] == 65
        assert archive["summary"]["generated_tables"] == 1
        with zipfile.ZipFile(tmp_path / archive["path"]) as z:
            assert not any(n.startswith("matched_eos/SPG_M2__BBB_BHF2") for n in z.namelist())
            member = "matched_eos/SPG_M1__APR_APR/candidate_000/compose/eos.thermo"
            original = tmp_path / "campaign" / archive["family"] / "SPG_M1__APR_APR/candidate_000/compose/eos.thermo"
            assert z.read(member) == original.read_bytes()
            diagnostic = list(csv.DictReader(io.StringIO(z.read("diagnostics/SPG_M1__APR_APR/candidate_000/source_retaining_thermo.csv").decode())))
            assert [r["view_segment"] for r in diagnostic[:2]] == ["crust_source_prefix"] * 2
            assert diagnostic[-1]["view_segment"] == "core_source_suffix"
            if archive["family"] == "public_reconstructed":
                assert diagnostic[0]["source_representation"] == "original_source"
                assert diagnostic[1]["source_representation"] == "prepared_public_reconstruction"
                assert float(diagnostic[1]["P_MeV_fm3"]) == pytest.approx(.0055)
                assert diagnostic[-1]["source_representation"] == "original_source"


def test_changed_verified_eos_bytes_are_rejected(tmp_path):
    campaign = make_fixture(tmp_path)
    p = tmp_path / "campaign/native_generated/SPG_M1__APR_APR/candidate_000/compose/eos.thermo"
    p.write_text(p.read_text() + "\n")
    with pytest.raises(ValueError, match="verified output hash changed"):
        packager.package(campaign, tmp_path / "delivery", root=tmp_path, figure=False)


def test_cross_family_reference_is_rejected_before_packaging(tmp_path):
    campaign = make_fixture(tmp_path)
    data = packager.load(campaign)
    data["cells"][0]["core_source_id"] = "public_original:APR_APR"
    packager.save(campaign, data)
    with pytest.raises(ValueError, match="cross-family source reference"):
        packager.package(campaign, tmp_path / "delivery", root=tmp_path, figure=False)
    assert not (tmp_path / "delivery").exists()


def test_row_data_uses_density_indices_and_fortran_exponents(tmp_path):
    directory = tmp_path / "table"
    table(directory, [.01, .05, .2, 1.0])
    expected = packager.row_data(directory)
    (directory / "eos.nb").write_text("0\n3\n1d-2\n5D-2\n2d-1\n1D0\n")
    lines = (directory / "eos.thermo").read_text().splitlines()
    rows = []
    for i, line in enumerate(lines[1:]):
        row = line.split()
        row[1] = str(i)
        for j in range(3, 10):
            row[j] = f"{float(row[j]):.17E}".replace("E", "d" if j % 2 else "D")
        rows.append(row)
    rows.reverse()
    (directory / "eos.thermo").write_text("9.39d2 9.38D2\n" + "\n".join(" ".join(r) for r in rows) + "\n")
    actual = packager.row_data(directory)
    for a, e in zip(actual, expected):
        assert a == pytest.approx(e)
    rows[0][1] = rows[1][1]
    (directory / "eos.thermo").write_text("9.39d2 9.38D2\n" + "\n".join(" ".join(r) for r in rows) + "\n")
    with pytest.raises(ValueError, match="duplicate or missing CompOSE density index"):
        packager.row_data(directory)
