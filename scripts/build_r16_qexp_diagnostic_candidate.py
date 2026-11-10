#!/usr/bin/env python3
"""Build a diagnostic-only native-knot q-exponential candidate over the 234 EOS cohort.

No authoritative source or projected EOS file is modified.  The candidate is a
separate derived evidence table, not a solver input or a physical EOS release.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from r16_pressure_preserving_q import (  # noqa: E402
    qexp_pressure_increment,
    qexp_shape_parameter,
)

SUMMARY = ROOT / "oscillation_release_v1/data/spectrum_summary_234.csv"
SOURCE_DIR = ROOT / "oscillation_release_v1/data/eos_original_normalized"
PROJECTED_DIR = ROOT / "oscillation_release_v1/data/eos_projected"
SUMMARY_SHA256 = "ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965"
EXPECTED_TABLES = 234
EXPECTED_CELLS = 109776
REGULARIZATION_FRACTION = 1.0e-3
RECONSTRUCTION_TOLERANCE = 5.0e-12
METHODS = ("css", "cubic_F", "cubic_P", "maxwell", "quintic_F", "quintic_P")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def numeric(rows: list[dict[str, str]], field: str) -> list[float]:
    values = [float(row[field]) for row in rows]
    if any(not math.isfinite(x) for x in values):
        raise ValueError(f"{field}: nonfinite values")
    return values


def max_relative(actual: list[float], expected: list[float]) -> float:
    return max(
        abs(a-e) / max(abs(e), 1.0e-300)
        for a, e in zip(actual, expected)
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0].keys()) if rows else ["status"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def build_one(reference: dict[str, str]) -> tuple[dict, list[dict], list[dict]]:
    name = reference["name"]
    source_path = SOURCE_DIR / f"{name}.csv"
    projected_path = PROJECTED_DIR / f"{name}.csv"
    if not source_path.is_file() or not projected_path.is_file():
        raise FileNotFoundError(f"{name}: normalized or projected input is missing")

    source_sha_before = sha256(source_path)
    projected_sha_before = sha256(projected_path)
    if source_sha_before != reference.get("raw_normalized_eos_sha256"):
        raise ValueError(f"{name}: normalized source SHA-256 mismatch")
    if projected_sha_before != reference.get("projected_eos_sha256"):
        raise ValueError(f"{name}: projected table SHA-256 mismatch")

    source_rows = read_csv(source_path)
    projected_rows = read_csv(projected_path)
    if len(source_rows) != len(projected_rows) or len(source_rows) < 3:
        raise ValueError(f"{name}: source/projected row count mismatch")
    for i, (s, p) in enumerate(zip(source_rows, projected_rows), start=1):
        if int(s["source_row"]) != i or int(p["source_row"]) != i:
            raise ValueError(f"{name}: source-row order mismatch at {i}")
        for field in ("source_row", "i_t", "i_nb", "i_yq"):
            if field in s and field in p and s[field] != p[field]:
                raise ValueError(f"{name}: identity field {field} changed at row {i}")

    mu = numeric(source_rows, "mu_b_mev")
    n = numeric(source_rows, "n_b_fm3")
    q = [ni/mui for ni,mui in zip(n,mu)]
    p_source = numeric(source_rows, "pressure_mev_fm3")
    eps_source = numeric(source_rows, "energy_mev_fm3")
    p_projected = numeric(projected_rows, "pressure_mev_fm3")
    eps_projected = numeric(projected_rows, "energy_mev_fm3")
    mu_projected = numeric(projected_rows, "mu_b_mev")
    n_projected = numeric(projected_rows, "n_b_fm3")

    if any(b <= a for a,b in zip(mu,mu[1:])) or any(b <= a for a,b in zip(n,n[1:])):
        raise ValueError(f"{name}: source mu or n is not strictly increasing")
    if any(b < a for a,b in zip(q,q[1:])):
        raise ValueError(f"{name}: source q=n/mu endpoints are not monotone")

    candidate_p = [p_source[0]]
    qlinear_p = [p_source[0]]
    cell_rows = []
    exact_count = regularized_count = 0
    k_values = []
    chi_values = []
    exact_dp_errors = []
    for i in range(len(mu)-1):
        mu0,mu1=mu[i],mu[i+1]
        q0,q1=q[i],q[i+1]
        dp_source=p_source[i+1]-p_source[i]
        exact=qexp_shape_parameter(mu0,mu1,q0,q1,dp_source)
        exact_feasible=bool(exact.get("feasible",False))
        regularized=None
        if exact_feasible:
            selected=exact
            target_dp=dp_source
            exact_count+=1
            reconstructed=qexp_pressure_increment(mu0,mu1,q0,q1,float(exact["k"]))
            exact_dp_errors.append(abs(reconstructed-dp_source)/max(abs(dp_source),1e-300))
        else:
            regularized=qexp_shape_parameter(
                mu0,mu1,q0,q1,dp_source,
                regularization_fraction=REGULARIZATION_FRACTION,
            )
            if regularized.get("k") is None:
                raise ValueError(f"{name}: no finite regularized shape at source cell {i+1}: {regularized}")
            selected=regularized
            target_dp=float(regularized["target_pressure_increment"])
            regularized_count+=1

        k=float(selected["k"])
        reconstructed_target=qexp_pressure_increment(mu0,mu1,q0,q1,k)
        target_error=abs(reconstructed_target-target_dp)/max(abs(target_dp),1e-300)
        if target_error>RECONSTRUCTION_TOLERANCE:
            raise ValueError(f"{name}: qexp target reconstruction mismatch at cell {i+1}: {target_error}")
        chi_lo=float(selected["chi_lower"])
        chi_hi=float(selected["chi_upper"])
        if not (0.0<chi_lo<=1.0 and 0.0<chi_hi<=1.0):
            raise ValueError(f"{name}: qexp endpoint chi outside causal domain at cell {i+1}")

        dp_linear=qexp_pressure_increment(mu0,mu1,q0,q1,0.0)
        candidate_p.append(candidate_p[-1]+target_dp)
        qlinear_p.append(qlinear_p[-1]+dp_linear)
        k_values.append(k)
        chi_values.extend([chi_lo,chi_hi])
        cell_rows.append({
            "name":name,
            "crust":reference["crust"],
            "core":reference["core"],
            "matching_method":reference["matching_method"],
            "lower_source_row":i+1,
            "upper_source_row":i+2,
            "mu_lower_mev":mu0,
            "mu_upper_mev":mu1,
            "q_lower_fm3_per_mev":q0,
            "q_upper_fm3_per_mev":q1,
            "source_delta_p_mev_fm3":dp_source,
            "qexp_candidate_delta_p_mev_fm3":target_dp,
            "qlinear_candidate_delta_p_mev_fm3":dp_linear,
            "candidate_delta_p_change_relative":target_dp/dp_source-1.0,
            "pressure_preserving_exact":exact_feasible,
            "regularized_to_interior_envelope":not exact_feasible,
            "shape_k":k,
            "chi_lower":chi_lo,
            "chi_upper":chi_hi,
            "whole_cell_stable_analytic":True,
            "whole_cell_causal_analytic":True,
            "exact_source_increment_reconstruction_relerr":(
                abs(reconstructed_target-dp_source)/max(abs(dp_source),1e-300)
                if exact_feasible else None
            ),
        })

    candidate_eps=[mui*ni-pi for mui,ni,pi in zip(mu,n,candidate_p)]
    qlinear_eps=[mui*ni-pi for mui,ni,pi in zip(mu,n,qlinear_p)]
    # The projector is exactly the k=0 member because q endpoints are already monotone.
    mismatch_p=max_relative(p_projected,qlinear_p)
    mismatch_eps=max_relative(eps_projected,qlinear_eps)
    mismatch_n=max_relative(n_projected,n)
    mismatch_mu=max_relative(mu_projected,mu)
    if max(mismatch_p,mismatch_eps,mismatch_n,mismatch_mu)>RECONSTRUCTION_TOLERANCE:
        raise ValueError(
            f"{name}: q-linear reference differs from frozen projected table "
            f"(P={mismatch_p:.3e}, eps={mismatch_eps:.3e}, n={mismatch_n:.3e}, mu={mismatch_mu:.3e})"
        )

    candidate_p_monotone=all(b>a for a,b in zip(candidate_p,candidate_p[1:]))
    candidate_eps_positive=all(x>0.0 for x in candidate_eps)
    candidate_eps_monotone=all(b>a for a,b in zip(candidate_eps,candidate_eps[1:]))
    if not candidate_p_monotone or not candidate_eps_positive or not candidate_eps_monotone:
        raise ValueError(f"{name}: qexp diagnostic candidate violates positive monotone P/epsilon")

    source_after=sha256(source_path)
    projected_after=sha256(projected_path)
    if source_sha_before!=source_after or projected_sha_before!=projected_after:
        raise RuntimeError(f"{name}: authoritative source/projected file changed during candidate build")

    pmax=max(p_source)
    epsmax=max(eps_source)
    p_rel=[candidate_p[i]/p_source[i]-1.0 for i in range(len(mu))]
    p_projected_rel=[p_projected[i]/p_source[i]-1.0 for i in range(len(mu))]
    eps_rel=[candidate_eps[i]/eps_source[i]-1.0 for i in range(len(mu))]
    qexp_nodes=[]
    for i in range(len(mu)):
        qexp_nodes.append({
            "name":name,
            "crust":reference["crust"],
            "core":reference["core"],
            "matching_method":reference["matching_method"],
            "source_row":i+1,
            "mu_b_mev":mu[i],
            "n_b_fm3":n[i],
            "q_n_over_mu_fm3_per_mev":q[i],
            "p_source_mev_fm3":p_source[i],
            "p_projected_qlinear_mev_fm3":p_projected[i],
            "p_qexp_diagnostic_mev_fm3":candidate_p[i],
            "epsilon_source_mev_fm3":eps_source[i],
            "epsilon_projected_qlinear_mev_fm3":eps_projected[i],
            "epsilon_qexp_diagnostic_mev_fm3":candidate_eps[i],
            "qexp_pressure_delta_over_source_pct":100.0*p_rel[i],
            "qlinear_pressure_delta_over_source_pct":100.0*p_projected_rel[i],
            "qexp_energy_delta_over_source_pct":100.0*eps_rel[i],
        })

    row217=row_index=216
    row218=217
    candidate_first_law=max(
        abs(e+p-mu0*n0)/max(abs(e+p),1e-300)
        for e,p,mu0,n0 in zip(candidate_eps,candidate_p,mu,n)
    )
    table={
        "name":name,
        "crust":reference["crust"],
        "core":reference["core"],
        "matching_method":reference["matching_method"],
        "row_count":len(source_rows),
        "normalized_source_sha256":source_sha_before,
        "projected_sha256":projected_sha_before,
        "source_unchanged":source_after==source_sha_before,
        "projected_unchanged":projected_after==projected_sha_before,
        "exact_pressure_preserving_cell_count":exact_count,
        "regularized_cell_count":regularized_count,
        "cell_count":len(cell_rows),
        "candidate_q_monotone_at_knots":True,
        "candidate_P_positive_monotone":candidate_p_monotone,
        "candidate_epsilon_positive_monotone":candidate_eps_positive and candidate_eps_monotone,
        "candidate_whole_cell_stable_causal_analytic":all(
            x["whole_cell_stable_analytic"] and x["whole_cell_causal_analytic"] for x in cell_rows
        ),
        "candidate_first_law_max_relative_error":candidate_first_law,
        "qlinear_vs_frozen_projected_P_max_rel_mismatch":mismatch_p,
        "qlinear_vs_frozen_projected_epsilon_max_rel_mismatch":mismatch_eps,
        "qexp_max_abs_delta_P_over_Pmax_pct":100.0*max(abs(a-b) for a,b in zip(candidate_p,p_source))/pmax,
        "qlinear_max_abs_delta_P_over_Pmax_pct":100.0*max(abs(a-b) for a,b in zip(p_projected,p_source))/pmax,
        "qexp_max_abs_local_delta_P_over_P_pct":100.0*max(abs(x) for x in p_rel),
        "qlinear_max_abs_local_delta_P_over_P_pct":100.0*max(abs(x) for x in p_projected_rel),
        "qexp_max_abs_delta_epsilon_over_epsilon_pct":100.0*max(abs(x) for x in eps_rel),
        "row217_qexp_pressure_delta_pct":100.0*p_rel[row217],
        "row218_qexp_pressure_delta_pct":100.0*p_rel[row218],
        "row217_qlinear_pressure_delta_pct":100.0*p_projected_rel[row217],
        "row218_qlinear_pressure_delta_pct":100.0*p_projected_rel[row218],
        "row217_218_qexp_k_min":min(k_values[216:217]) if len(k_values)>216 else None,
        "row217_218_qexp_k_max":max(k_values[216:217]) if len(k_values)>216 else None,
        "row217_218_qexp_chi_lower":cell_rows[216]["chi_lower"],
        "row217_218_qexp_chi_upper":cell_rows[216]["chi_upper"],
        "max_cell_qexp_shape_abs_k":max(abs(x) for x in k_values),
        "status":"PASS_DIAGNOSTIC_QEXP_CANDIDATE",
    }
    return table,cell_rows,qexp_nodes


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=ROOT/"dist/r16-qexp-diagnostic-candidate")
    args=parser.parse_args(argv)

    if sha256(SUMMARY)!=SUMMARY_SHA256:
        raise RuntimeError("frozen spectrum summary SHA-256 mismatch")
    refs=read_csv(SUMMARY)
    names={x["name"] for x in refs}
    if len(refs)!=EXPECTED_TABLES or len(names)!=EXPECTED_TABLES:
        raise RuntimeError("frozen cohort is not exactly 234 unique EOS tables")
    if {x["matching_method"] for x in refs}!=set(METHODS):
        raise RuntimeError("frozen matching method inventory changed")
    if {p.stem for p in SOURCE_DIR.glob("*.csv")}!=names or {p.stem for p in PROJECTED_DIR.glob("*.csv")}!=names:
        raise RuntimeError("normalized/projected table inventory differs from summary")

    tables,cells,nodes,failures=[],[],[],[]
    by_name={x["name"]:x for x in refs}
    for name in sorted(names):
        try:
            table,cell_rows,node_rows=build_one(by_name[name])
            tables.append(table);cells.extend(cell_rows);nodes.extend(node_rows)
        except Exception as exc:
            failures.append({"name":name,"error":f"{type(exc).__name__}: {exc}"[:2000]})

    expected_cells=sum(x["row_count"]-1 for x in tables)
    expected_nodes=sum(x["row_count"] for x in tables)
    all_integrity=(len(tables)==EXPECTED_TABLES and all(x["source_unchanged"] and x["projected_unchanged"] for x in tables))
    all_monotone=(len(tables)==EXPECTED_TABLES and all(
        x["candidate_P_positive_monotone"] and x["candidate_epsilon_positive_monotone"]
        and x["candidate_whole_cell_stable_causal_analytic"] for x in tables
    ))
    all_recon=(len(tables)==EXPECTED_TABLES and all(
        max(x["qlinear_vs_frozen_projected_P_max_rel_mismatch"],
            x["qlinear_vs_frozen_projected_epsilon_max_rel_mismatch"])<=RECONSTRUCTION_TOLERANCE
        for x in tables
    ))
    overall_pass=(
        not failures and len(cells)==expected_cells and len(nodes)==expected_nodes
        and all_integrity and all_monotone and all_recon
    )

    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    write_csv(out/"r16_qexp_diagnostic_candidate_cells.csv",cells)
    write_csv(out/"r16_qexp_diagnostic_candidate_nodes.csv",nodes)
    write_csv(out/"r16_qexp_diagnostic_candidate_tables.csv",tables)

    crust_groups=defaultdict(list)
    for row in tables: crust_groups[row["crust"]].append(row)
    crust_summary=[]
    for crust,rows in sorted(crust_groups.items()):
        crust_summary.append({
            "crust":crust,
            "table_count":len(rows),
            "exact_pressure_preserving_cell_count":sum(x["exact_pressure_preserving_cell_count"] for x in rows),
            "regularized_cell_count":sum(x["regularized_cell_count"] for x in rows),
            "qexp_row218_pressure_delta_pct_mean":sum(x["row218_qexp_pressure_delta_pct"] for x in rows)/len(rows),
            "qexp_row218_pressure_delta_pct_min":min(x["row218_qexp_pressure_delta_pct"] for x in rows),
            "qexp_row218_pressure_delta_pct_max":max(x["row218_qexp_pressure_delta_pct"] for x in rows),
            "qlinear_row218_pressure_delta_pct_mean":sum(x["row218_qlinear_pressure_delta_pct"] for x in rows)/len(rows),
            "qexp_qmax_error_over_Pmax_mean_pct":sum(x["qexp_max_abs_delta_P_over_Pmax_pct"] for x in rows)/len(rows),
            "qlinear_qmax_error_over_Pmax_mean_pct":sum(x["qlinear_max_abs_delta_P_over_Pmax_pct"] for x in rows)/len(rows),
        })
    write_csv(out/"r16_qexp_diagnostic_candidate_crust_summary.csv",crust_summary)

    receipt={
        "schema":"ns-cowling-r16-qexp-diagnostic-candidate-v1",
        "status":"PASS_R16_QEXP_DIAGNOSTIC_CANDIDATE" if overall_pass else "FAIL_R16_QEXP_DIAGNOSTIC_CANDIDATE",
        "source_commit":os.environ.get("GITHUB_SHA","UNAVAILABLE_LOCAL_RUN"),
        "source_summary_sha256":SUMMARY_SHA256,
        "diagnostic_scope":"native-knot q-exponential pressure candidate; not an EOS release or solver input",
        "shape_definition":"q(t)=q0+(q1-q0)*(1-exp(-k*t))/(1-exp(-k)); P is its exact integral; epsilon=mu*n-P",
        "regularization_fraction":REGULARIZATION_FRACTION,
        "source_table_count_expected":EXPECTED_TABLES,
        "source_table_count_audited":len(tables),
        "native_cell_count_expected":expected_cells,
        "native_cell_count_audited":len(cells),
        "native_node_count_expected":expected_nodes,
        "native_node_count_audited":len(nodes),
        "exact_pressure_preserving_cell_count":sum(x["exact_pressure_preserving_cell_count"] for x in tables),
        "regularized_cell_count":sum(x["regularized_cell_count"] for x in tables),
        "candidate_all_tables_P_positive_monotone":all_monotone,
        "candidate_all_tables_epsilon_positive_monotone":all_monotone,
        "candidate_all_cells_stable_causal_by_qprime_identity":all_monotone,
        "candidate_first_law_max_relative_error":max((x["candidate_first_law_max_relative_error"] for x in tables),default=None),
        "qexp_max_abs_delta_P_over_Pmax_mean_pct":sum(x["qexp_max_abs_delta_P_over_Pmax_pct"] for x in tables)/len(tables) if tables else None,
        "qexp_max_abs_delta_P_over_Pmax_max_pct":max((x["qexp_max_abs_delta_P_over_Pmax_pct"] for x in tables),default=None),
        "qlinear_max_abs_delta_P_over_Pmax_mean_pct":sum(x["qlinear_max_abs_delta_P_over_Pmax_pct"] for x in tables)/len(tables) if tables else None,
        "qlinear_max_abs_delta_P_over_Pmax_max_pct":max((x["qlinear_max_abs_delta_P_over_Pmax_pct"] for x in tables),default=None),
        "qexp_max_abs_local_delta_P_over_P_mean_pct":sum(x["qexp_max_abs_local_delta_P_over_P_pct"] for x in tables)/len(tables) if tables else None,
        "qexp_max_abs_local_delta_P_over_P_max_pct":max((x["qexp_max_abs_local_delta_P_over_P_pct"] for x in tables),default=None),
        "qexp_max_abs_delta_epsilon_over_epsilon_mean_pct":sum(x["qexp_max_abs_delta_epsilon_over_epsilon_pct"] for x in tables)/len(tables) if tables else None,
        "qexp_max_abs_delta_epsilon_over_epsilon_max_pct":max((x["qexp_max_abs_delta_epsilon_over_epsilon_pct"] for x in tables),default=None),
        "row217_qexp_pressure_delta_pct_mean":sum(x["row217_qexp_pressure_delta_pct"] for x in tables)/len(tables) if tables else None,
        "row218_qexp_pressure_delta_pct_mean":sum(x["row218_qexp_pressure_delta_pct"] for x in tables)/len(tables) if tables else None,
        "row218_qexp_pressure_delta_pct_range":[min(x["row218_qexp_pressure_delta_pct"] for x in tables),max(x["row218_qexp_pressure_delta_pct"] for x in tables)] if tables else None,
        "row218_qlinear_pressure_delta_pct_mean":sum(x["row218_qlinear_pressure_delta_pct"] for x in tables)/len(tables) if tables else None,
        "source_and_projected_tables_unchanged":all_integrity,
        "qlinear_reconstruction_matches_frozen_projected_input":all_recon,
        "authoritative_eos_files_written":0,
        "synthetic_native_rows_written":0,
        "tov_integrations_executed":0,
        "failures":failures,
    }
    receipt_path=out/"r16_qexp_diagnostic_candidate_receipt.json"
    receipt_path.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    manifest={
        "schema":"ns-cowling-r16-qexp-diagnostic-candidate-manifest-v1",
        "receipt_sha256":sha256(receipt_path),
        "cells_csv_sha256":sha256(out/"r16_qexp_diagnostic_candidate_cells.csv"),
        "nodes_csv_sha256":sha256(out/"r16_qexp_diagnostic_candidate_nodes.csv"),
        "tables_csv_sha256":sha256(out/"r16_qexp_diagnostic_candidate_tables.csv"),
        "crust_summary_csv_sha256":sha256(out/"r16_qexp_diagnostic_candidate_crust_summary.csv"),
        "source_summary_sha256":SUMMARY_SHA256,
        "audited_table_count":len(tables),
        "audited_cell_count":len(cells),
        "audited_node_count":len(nodes),
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2),flush=True)
    return 0 if overall_pass else 1


if __name__=="__main__":
    raise SystemExit(main())
