#!/usr/bin/env python3
"""Source identity/support inventory, not a physical Maxwell admission test."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from intake_three_family_sources_20261117 import CORE, SPG, TF, ROOT, table_inventory


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path):
    return str(path.relative_to(ROOT))


def phase_inventory(path: Path):
    base = table_inventory(path)
    n = [float(x) for x in (path / "eos.nb").read_text().split()[2:]]
    rows = [x.split() for x in (path / "eos.thermo").read_text().splitlines() if x.strip()][1:]
    if (path / "eos.compo").is_file():
        phases = [int(x.split()[3]) for x in (path / "eos.compo").read_text().splitlines() if x.strip()]
        if len(phases) != len(n):
            raise ValueError("composition row count differs")
        blocks = []
        start = 0
        for end in range(1, len(n) + 1):
            if end == len(n) or phases[end] != phases[start]:
                blocks.append({"phase_id": phases[start], "start_zero_based": start,
                               "stop_exclusive": end, "rows": end - start,
                               "n_min_fm3": n[start], "n_max_fm3": n[end - 1]})
                start = end
    else:
        blocks = []
    nadd = sorted({int(row[10]) for row in rows})
    extra = [float(row[11]) for row in rows if len(row) > 11]
    base.update(phase_blocks=blocks, nadd_values=nadd,
                q_extra_range=[min(extra), max(extra)] if extra else None,
                source_Q5_range=[min(float(r[7]) for r in rows), max(float(r[7]) for r in rows)],
                thermodynamic_changes="NONE", table_sha256={f.name: sha(f) for f in sorted(path.glob("eos.*")) if f.is_file()})
    return base


def native_entry(model: str):
    if model in CORE:
        path = ROOT / "data/research_matching_20261117/integration/real_pairs_v2/core" / model / "generated"
        receipt = path.parent / "DISPATCH.json"
    elif model in SPG:
        path = ROOT / "data/research_matching_20261117/spg/native_129_repaired" / model / "generated"
        receipt = path.parent / "RECEIPT.json"
    else:
        path = ROOT / "data/research_native_20261116/source_grid_cohort" / model / "attempt_001/output/generated"
        receipt = path.parent.parent / "COHORT_RESULT.json"
    result = {"model": model, "family": "code_generated", "path": relative(path),
              "receipt_path": relative(receipt), "receipt_sha256": sha(receipt), **phase_inventory(path)}
    record = json.loads(receipt.read_text())
    evidence = record.get("evidence", record.get("native_evidence", {}))
    expected = evidence.get("output_sha256", {})
    if not expected:
        expected = {k: v for k, v in record.get("output_sha256", {}).items() if k.startswith("eos.")}
    checks = {name: {"expected": digest, "actual": sha(path / name)} for name, digest in expected.items() if (path / name).is_file()}
    if not checks or any(v["expected"] != v["actual"] for v in checks.values()) or any(not (path / name).is_file() for name in expected):
        raise ValueError("native source output identity failed: " + model)
    result.update(retained_hashes_verified=checks, support=json.loads((path / "SUPPORT.json").read_text()),
                  support_sha256=sha(path / "SUPPORT.json"), request=record.get("request"))
    if model in TF:
        selected = list(csv.DictReader((path / "selected.csv").open()))
        attempts = list(csv.DictReader((path / "attempts.csv").open()))
        lookup = {(r["node"], r["mesh"], r["attempt"]): r for r in attempts}
        selected_attempts = [lookup[(r["node"], r["mesh"], r["attempt"])] for r in selected]
        result.update(selected_geometry_counts={g: sum(r["geometry"] == g for r in selected) for g in sorted({r["geometry"] for r in selected})},
                      selected_state_count=len(selected), selected_status_counts={s: sum(r["status"] == s for r in selected_attempts) for s in sorted({r["status"] for r in selected_attempts})},
                      selected_optimizer_modes=sorted({int(r["optimizer_mode"]) for r in selected_attempts}),
                      maximum_selected_KKT_MeV=max(abs(float(r["KKT"])) for r in selected),
                      maximum_selected_radius_stationarity=max(abs(float(r["radius_stationarity"])) for r in selected),
                      requested_densities=record["request"]["densities"],
                      selected_densities=[float(r["n_b"]) for r in selected],
                      held_nodes=evidence["native_result"]["held_nodes"],
                      historical_native_checks=evidence["native_consistency"]["generated"],
                      source_scope="Inner-crust finite N16 samples; no cross-phase interpolation, lower-crust attachment, global minimum or continuum claim",
                      refined_native_alternatives=[
                          f"data/research_completion_20261116/gppva_tf/final/{model}/compose",
                          f"data/research_finish_20261116/geometry/accelerated_final/{model}/N384/compose",
                          f"data/research_finish_20261116/geometry/branches/{model}/N384"],
                      refinement_note="Alternatives are separate higher-mesh finite samples and branch pieces, not interchangeable full-crust replacements. Use their receipts/support and no cross-gap interpolation.")
    elif model in SPG:
        result.update(model_variant="public_composition_inferred_20261116", source_scope="Native spherical CLD inner crust, candidate inferred surface parameters; original author calibration unresolved",
                      public_thermo_reconstruction=False)
    else:
        result.update(source_scope="Native finite core support; full retained phase branches must remain explicit",
                      model_variant=record["request"]["contract"]["representation"])
    return result


def public_entry(model: str, path: Path):
    result = {"model": model, "family": "public_original", "path": relative(path), **phase_inventory(path)}
    if model in CORE:
        starts = {"APR_APR": 69, "BBB_BHF2": 69, "GM_GM1": 48, "RG_KDE0v1": 966, "RG_SLY4": 143}
        n = [float(x) for x in (path / "eos.nb").read_text().split()[2:]]
        start = starts[model]
        result["homogeneous_core_selection"] = {"start_zero_based": start, "stop_exclusive": len(n),
            "n_min_fm3": n[start], "n_max_fm3": n[-1], "attached_lower_rows_excluded": start,
            "evidence": "Fresh provider page/PDF plus registered model-specific component boundary; exact raw rows preserved",
            "provider_beta_equilibrium": True}
        if model == "GM_GM1":
            result["homogeneous_core_selection"]["warning"] = "Fresh file labels all rows Iphase=3, including attached DH crust. Provider cutoff 0.001 fm^-3, first available grid row above cutoff is n=0.048. Phase3 alone is insufficient."
        if model.startswith("RG_"):
            result["homogeneous_core_selection"]["phase_dictionary"] = "Provider PDF: Iphase4 heavy nuclei; Iphase3 homogeneous matter"
            result["chemical_column_status"] = "Provider says beta equilibrium, stored Q5 nonzero; must resolve representation/thermodynamic semantics, not automatically classify source as non-beta-equilibrated."
        if model == "RG_KDE0v1":
            result["extra_column_semantics"] = "UNKNOWN_PROVIDER_AUXILIARY; neither supplied 3-page PDF, README.txt nor eos.init defines the single extra column. It must not be called sound speed based on numerical resemblance."
            result["auxiliary_projection_policy"] = "Lossless raw plus row-indexed auxiliary sidecar; Nadd=0 only in an explicitly named thermodynamic-only working view; never fabricate or interpret the omitted quantity."
    else:
        result["source_scope"] = "Public full raw crust including attached lower segment; source gap between phases0 and7 is not licensed as a smooth interpolating interval."
        if model in TF:
            result["extra_column_semantics"] = "Nadd categorical geometry; distinct from Iphase; no numerical interpolation"
            result["extra_column_primary_evidence"] = {"path": relative(path / "eos.pdf"), "sha256": sha(path / "eos.pdf"), "page": 2,
                "definition": {"0": "outer crust", "1": "droplets", "2": "rods", "3": "slabs", "4": "tubes", "5": "bubbles", "6": "homogeneous core"}}
            result["auxiliary_projection_policy"] = "Retain exact per-row geometry separately from Iphase; a thermo-only working view may omit the payload only with a declared loss receipt and without erasing its branch-boundary information."
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--intake", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    reference_map = json.loads((a.intake / "REFERENCE_MAP.json").read_text())
    native = [native_entry(m) for m in CORE + SPG + TF]
    public = [public_entry(m, ROOT / reference_map[m]) for m in CORE + SPG + TF]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({"schema": "ns_cowling.maxwell_source_inventory.v1",
        "purpose": "Exact source identity and bounded support; no new physical acceptance claim",
        "script_sha256": sha(Path(__file__)), "core_ids": CORE, "crust_ids": SPG + TF,
        "within_family_pair_count": 65, "requested_family_count": 3,
        "native": native, "public_original": public,
        "reconstructed": "Created separately from public originals with per-change thermodynamic diagnosis; never from native generator output"}, indent=2) + "\n")
    print(json.dumps({"native_models": len(native), "public_models": len(public), "all_native_output_hashes_match": True, "output": str(a.out)}))


if __name__ == "__main__":
    main()
