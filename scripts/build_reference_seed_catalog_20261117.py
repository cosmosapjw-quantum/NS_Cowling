#!/usr/bin/env python3
"""Extract bounded DD2 reference diagnostics and composition seed coordinates.

These are public table coordinates, not solved TF spatial profiles.  No public
pressure, energy, phase endpoint, or composition is forced into the generator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from compare_generated_reference import read_cold_table, continuous_segments

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    tables, records = {}, {}
    for model in ("GPPVA_DD2", "GPPVA_DD2_UNIFIED"):
        table = tables[model] = read_cold_table(args.references / model)
        segments = [dict(phase=int(table.phase[s[0]]), rows=len(s),
                         n_min_fm3=float(table.n[s[0]]), n_max_fm3=float(table.n[s[-1]]))
                    for s in continuous_segments(table, set(table.phase.tolist()))]
        boundaries = []
        for i in np.flatnonzero(table.phase[1:] != table.phase[:-1]) + 1:
            boundaries.append(dict(left_row=int(i-1), right_row=int(i),
                                   left_phase=int(table.phase[i-1]), right_phase=int(table.phase[i]),
                                   left_n_fm3=float(table.n[i-1]), right_n_fm3=float(table.n[i]),
                                   left_P_mev_fm3=float(table.fields["pressure"][i-1]),
                                   right_P_mev_fm3=float(table.fields["pressure"][i]),
                                   status="REFERENCE_BOUNDARY_REQUIRES_EXPLICIT_JOIN_LAW"))
        seeds = []
        crust = np.flatnonzero(table.phase == 7)
        for target in (.002, .01, .03, .05, .067):
            j = int(crust[np.argmin(np.abs(table.n[crust]-target))])
            ye, ym = (float(table.fields[k][j]) for k in ("y_e", "y_mu"))
            if not np.isfinite(ye) or not np.isfinite(ym):
                raise ValueError("seed needs explicitly supplied lepton fractions")
            seeds.append(dict(target_n_fm3=target, actual_n_fm3=float(table.n[j]), row=j,
                              P_mev_fm3=float(table.fields["pressure"][j]),
                              epsilon_mev_fm3=float(table.fields["epsilon"][j]),
                              mu_B_mev=float(table.fields["mu_b"][j]), Ye=ye, Ymu=ym,
                              Yp_total_from_neutrality=ye+ym,
                              derivation="total proton fraction = Ye+Ymu for charge-neutral npe(mu)N matter",
                              scope="composition coordinate for a candidate initial ansatz; no TF fields/cell geometry supplied"))
        records[model] = dict(source=table.metadata, segments=segments, boundaries=boundaries, seeds=seeds)
    a, b = tables["GPPVA_DD2"], tables["GPPVA_DD2_UNIFIED"]
    common = len(a.n)
    if len(b.n) < common or not np.array_equal(a.phase, b.phase[:common]):
        raise ValueError("crust/unified phase order differs; cannot align by retained order")
    density_relative = np.max(np.abs(a.n-b.n[:common])/a.n)
    if density_relative > 1e-6:
        raise ValueError("crust/unified density ordering not aligned within printed precision")
    comparison = {"rows": common, "alignment": "same ordered phase rows; density-rounding difference retained",
                  "maximum_relative_density_difference": float(density_relative), "metrics": {}}
    for field in ("pressure", "epsilon", "mu_b", "enthalpy_b", "y_e", "y_mu"):
        aa, bb = a.fields[field], b.fields[field][:common]
        valid = np.isfinite(aa) & np.isfinite(bb)
        diff = aa[valid]-bb[valid]
        nonzero = np.abs(aa[valid]) > 0
        comparison["metrics"][field] = {"max_abs_difference": float(np.max(np.abs(diff))),
            "max_relative_difference_nonzero": float(np.max(np.abs(diff[nonzero]/aa[valid][nonzero]))) if np.any(nonzero) else None}
    result = dict(schema="ns-cowling.reference-seed-catalog.v1", status="REFERENCE_DIAGNOSTICS_PREPARED",
                  production_solver_run=False, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  models=records, crust_unified_common_rows=comparison,
                  limitations=["no complete TF profile is supplied by these tables",
                               "nuclear group A/Z/abundance semantics are not inferred",
                               "no interpolation across a phase boundary or uncovered density gap",
                               "public cutoff does not determine a generated model's phase transition"])
    (args.out / "SEED_CATALOG.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"status": result["status"], "models": {m:r["segments"] for m,r in records.items()},
                      "common_rows": common, "maximum_relative_density_difference": float(density_relative)}, indent=2))


if __name__ == "__main__":
    main()
