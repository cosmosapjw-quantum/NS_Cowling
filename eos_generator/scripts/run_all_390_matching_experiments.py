#!/usr/bin/env python3
"""Generate the 13x5x6 crust-core EOS matching matrix.

The default output is a strict, thermodynamics-first matched EOS archive.
Failed combinations are never written unless --allow-failed is supplied; when
that flag is used every non-PASS output is tagged diagnostic_only in metadata.
TOV target-mass calculations are opt-in and use the bracketed root solver.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_DIR / "compose_eos_pkg"))

from compose_eos import (
    CRUST_NAMES, CORE_NAMES, METHODS,
    read_compose, write_compose, match_crust_core, check_consistency,
    solve_star_at_target_mass,
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eos-root", type=Path, default=PROJECT_DIR / "eos_data_normalized")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-csv", type=Path)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--n-bridge", type=int, default=40)
    parser.add_argument("--core-compo-mode", choices=("strict_thermo_only", "literature_beta_eq"),
                        default="strict_thermo_only")
    parser.add_argument("--allow-failed", action="store_true",
                        help="write failed EOS for diagnosis; they are tagged diagnostic_only")
    parser.add_argument("--tov-target-mass", type=float, default=None,
                        help="optional single target mass in Msun; solved independently for each output EOS")
    args = parser.parse_args(argv)

    if args.n_bridge < 4:
        parser.error("--n-bridge must be >= 4")
    if not args.eos_root.is_dir():
        parser.error(f"EOS input directory not found: {args.eos_root}")
    if args.output.exists() and any(args.output.iterdir()):
        parser.error(f"output directory is not empty: {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    csv_path = args.summary_csv or (args.output / "Matching_390_Summary_v3_2.csv")
    json_path = args.summary_json or (args.output / "Matching_390_Results_v3_2.json")

    crusts = {name: read_compose(str(args.eos_root / name)) for name in CRUST_NAMES}
    cores = {name: read_compose(str(args.eos_root / name)) for name in CORE_NAMES}
    source_hashes = {
        name: {
            fn: sha256(args.eos_root / name / fn)
            for fn in ("eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo")
            if (args.eos_root / name / fn).is_file()
        }
        for name in (*CRUST_NAMES, *CORE_NAMES)
    }

    rows = []
    counts = {"PASS": 0, "FAIL": 0, "REJECTED": 0, "WRITTEN": 0}
    for core_name in CORE_NAMES:
        for crust_name in CRUST_NAMES:
            for method in METHODS:
                result = {
                    "core": core_name,
                    "crust": crust_name,
                    "method": method,
                    "written": False,
                    "tov_status": "NOT_REQUESTED",
                }
                try:
                    eos, info = match_crust_core(
                        crusts[crust_name], cores[core_name],
                        crust_name, core_name,
                        method=method, n_bridge=args.n_bridge,
                        core_compo_mode=args.core_compo_mode,
                    )
                    report = check_consistency(eos, verbose=False)
                    passed = bool(report["overall"]["pass_"])
                    status = "PASS" if passed else "FAIL"
                    counts[status] += 1
                    failed = list(report["overall"]["failed_checks"])

                    result.update({
                        "status": status,
                        "failed_checks": failed,
                        "composition_available": bool(eos.composition_available),
                        "composition_provenance": info.get("composition_provenance", "none"),
                        "core_composition_input_status": info.get("core_composition_input_status"),
                        "core_composition_input_usable": info.get("core_composition_input_usable"),
                        "nb_transition_left_fm3": float(info["n1"]),
                        "nb_transition_right_fm3": float(info["n2"]),
                        "P_transition_left_MeV_fm3": float(info["P1"]),
                        "P_transition_right_MeV_fm3": float(info["P2"]),
                        "mu_transition_left_MeV": float(info["mu1"]),
                        "mu_transition_right_MeV": float(info["mu2"]),
                        "endpoint_cs2_left": float(info["cs2_1"]),
                        "endpoint_cs2_right": float(info["cs2_2"]),
                        "bridge_cs2_min": float(info["cs2_min_bridge"]),
                        "bridge_cs2_max": float(info["cs2_max_bridge"]),
                        "C2_max_relative_residual": float(report["C2"]["max_rel"]),
                        "C3_max_abs_residual": (float(report["C3"]["max_abs"])
                                                  if report["C3"]["status"] == "CHECKED" else None),
                        "C4_max_abs_residual": (float(report["C4"]["max_abs"])
                                                  if report["C4"]["status"] == "CHECKED" else None),
                        "C10_min_cs2": float(report["C10"]["min_cs2"]),
                        "C10_max_cs2": float(report["C10"]["max_cs2"]),
                        "construction_kind": info["construction_kind"],
                        "is_exact_maxwell": False,
                        "density_jump_fm3": 0.0,
                        "seidov_applicable": False,
                        "phase_transition_order": None,
                        "diagnostic_only": not passed,
                    })

                    if passed or args.allow_failed:
                        out_dir = args.output / f"{crust_name}__{core_name}__{method}"
                        out_dir.mkdir(parents=True, exist_ok=False)
                        write_compose(str(out_dir), eos)
                        metadata = {
                            **result,
                            "schema": "compose-eos-core-crust-match-v3.2",
                            "source_crust": crust_name,
                            "source_core": core_name,
                            "source_hashes": {
                                "crust": source_hashes[crust_name],
                                "core": source_hashes[core_name],
                            },
                            "matching_parameters": {
                                "n_bridge": args.n_bridge,
                                "core_compo_mode": args.core_compo_mode,
                                "dm_crust_MeV": float(info["dm_crust"]),
                                "dm_core_MeV": float(info["dm_core"]),
                            },
                            "interpretation": {
                                "maxwell_token_is_compatibility_label_only": method == "maxwell",
                                "physical_construction": info["construction_kind"],
                                "no_density_jump": True,
                                "seidov_criterion_used": False,
                            },
                        }
                        (out_dir / "MATCHING_METADATA.json").write_text(
                            json.dumps(metadata, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8",
                        )
                        result["output_directory"] = out_dir.name
                        result["written"] = True
                        counts["WRITTEN"] += 1

                        if args.tov_target_mass is not None:
                            try:
                                pc_guess = float(eos.p[min(max(5, eos.n // 2), eos.n - 2)])
                                M, R, k2, lam, dR, pc = solve_star_at_target_mass(
                                    eos, float(info["P1"]), float(info["P2"]),
                                    pc_guess, float(args.tov_target_mass),
                                )
                                result["tov_status"] = "ROOT_SOLVED"
                                result["tov"] = {
                                    "target_mass_Msun": float(args.tov_target_mass),
                                    "actual_mass_Msun": float(M),
                                    "mass_residual_Msun": float(M - args.tov_target_mass),
                                    "central_pressure_MeV_fm3": float(pc),
                                    "radius_km": float(R),
                                    "k2": float(k2),
                                    "Lambda": float(lam),
                                    "crust_thickness_km": float(dR),
                                }
                            except Exception as exc:
                                result["tov_status"] = "NOT_BRACKETED_OR_FAILED"
                                result["tov_error"] = f"{type(exc).__name__}: {exc}"
                    elif not passed:
                        result["tov_status"] = "SKIPPED_FAILED_HARD_GATE"
                except Exception as exc:
                    counts["REJECTED"] += 1
                    result.update({
                        "status": "REJECTED",
                        "error": f"{type(exc).__name__}: {exc}",
                        "diagnostic_only": True,
                    })
                rows.append(result)

    fields = sorted({key for row in rows for key in row if not isinstance(row.get(key), dict) and key != "failed_checks"})
    fields.append("failed_checks")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            serial = {key: row.get(key) for key in fields if key != "failed_checks"}
            serial["failed_checks"] = ";".join(row.get("failed_checks", [])) if isinstance(row.get("failed_checks"), list) else ""
            writer.writerow(serial)

    summary = {
        "schema": "compose-eos-matching-sweep-v3.2",
        "input_root_label": "eos_data_normalized",
        "counts": counts,
        "expected_combinations": len(CRUST_NAMES) * len(CORE_NAMES) * len(METHODS),
        "core_compo_mode": args.core_compo_mode,
        "allow_failed": bool(args.allow_failed),
        "tov_target_mass_Msun": args.tov_target_mass,
        "phase_transition_interpretation": "All six methods construct smooth single-valued barotropes. The maxwell method token is a legacy compatibility label for a smooth near-pressure-plateau bridge; none implements an exact Maxwell density jump.",
        "source_hashes": source_hashes,
        "results": rows,
    }
    json_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"counts": counts, "summary_csv": str(csv_path), "summary_json": str(json_path)}, indent=2))
    return 0 if counts["REJECTED"] == 0 and (args.allow_failed or counts["FAIL"] == 0) else 2


if __name__ == "__main__":
    raise SystemExit(main())
