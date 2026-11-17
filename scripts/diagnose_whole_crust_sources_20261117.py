#!/usr/bin/env python3
"""Independent source audit for the explicitly defined whole-crust baseline.

This tool never constructs or rewrites an EOS.  It reads indexed CoMPOSE rows,
checks candidate P(mu_E) Hermite representations with polynomial extrema, and
records finite endpoint bridge candidates.  A bridge candidate is an explicitly
new barotrope, not a recovered microscopic equilibrium or a Maxwell seam.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from numpy.polynomial import polynomial as poly

ROOT = Path(__file__).resolve().parents[1]
EPS = np.finfo(float).eps


def read_table(path: Path) -> dict:
    tokens = (path / "eos.nb").read_text().split()
    first, last = map(int, tokens[:2])
    n = np.array([float(x.replace("D", "E").replace("d", "e")) for x in tokens[2:]])
    if n.size != last-first+1 or np.any(np.diff(n) <= 0):
        raise ValueError(f"invalid indexed density axis: {path}")
    lines = (path / "eos.thermo").read_text().splitlines()
    mn, mp, il = map(float, lines[0].split())
    q = np.full((len(n), 7), np.nan)
    seen = set()
    for line in lines[1:]:
        if not line.strip():
            continue
        s = line.replace("D", "E").replace("d", "e").split()
        i = int(s[1])-first
        if not 0 <= i < len(n) or i in seen:
            raise ValueError(f"duplicate or invalid density index: {path}")
        seen.add(i)
        q[i] = list(map(float, s[3:10]))
    if len(seen) != len(n) or np.any(~np.isfinite(q)):
        raise ValueError(f"missing or invalid thermo row: {path}")
    phase = np.full(len(n), -1, dtype=int)
    isotopes = [None] * len(n)
    if (path / "eos.compo").exists():
        for line in (path / "eos.compo").read_text().splitlines():
            if not line.strip():
                continue
            s = line.split()
            i = int(s[1])-first
            phase[i] = int(s[3])
            k = 5+2*int(s[4])
            groups = int(s[k])
            if groups:
                # CoMPOSE group: identifier, A_av, Z_av, abundance.  Continuously
                # varying CLD averages are not discrete isotope phase labels.
                mass, charge = float(s[k+2]), float(s[k+3])
                if abs(mass-round(mass)) < 1e-9 and abs(charge-round(charge)) < 1e-9:
                    isotopes[i] = (int(round(mass)), int(round(charge)))
    pressure = n*q[:, 0]
    epsilon = n*mn*(1+q[:, 6])
    mu = mn*(1+q[:, 2])
    mu_e = (epsilon+pressure)/n
    return dict(path=str(path.relative_to(ROOT)), n=n, p=pressure, e=epsilon,
                mu=mu, mu_e=mu_e, phase=phase, isotopes=isotopes,
                q=q, mn=mn, mp=mp, il=int(il))


def cut_table(table: dict, stop: int) -> dict:
    return {k: (v[:stop] if isinstance(v, (np.ndarray, list)) else v)
            for k, v in table.items()}


def polynomial_minimum(coefficients: np.ndarray) -> float:
    """Exact extrema of a floating-point polynomial on normalized [0,1]."""
    points = [0., 1.]
    for root in poly.polyroots(poly.polyder(coefficients)):
        if abs(root.imag) < 1e-11 and 0 < root.real < 1:
            points.append(float(root.real))
    return float(min(poly.polyval(t, coefficients) for t in points))


def hermite_check(n0: float, p0: float, mu0: float,
                  n1: float, p1: float, mu1: float) -> dict:
    values = tuple(map(float, (n0, p0, mu0, n1, p1, mu1)))
    n0, p0, mu0, n1, p1, mu1 = values
    base = dict(n_left=n0, n_right=n1, pressure_left=p0, pressure_right=p1,
                mu_left=mu0, mu_right=mu1)
    if not all(np.isfinite(values)) or min(n0, n1, mu0, mu1) <= 0:
        return dict(base, status="INVALID_ENDPOINT")
    if n1 <= n0:
        return dict(base, status="DENSITY_NOT_INCREASING")
    h = mu1-mu0
    if h <= 0:
        return dict(base, status="MU_NOT_INCREASING")
    dp = p1-p0
    # P(t), t=(mu-mu0)/h; endpoint n=dP/dmu fixed exactly.
    a = np.array([p0, h*n0, 3*dp-h*(2*n0+n1), -2*dp+h*(n0+n1)])
    density = poly.polyder(a)/h
    susceptibility = poly.polyder(density)/h
    causal = poly.polysub(poly.polymul([mu0, h], susceptibility), density)
    min_n = polynomial_minimum(density)
    min_d = polynomial_minimum(susceptibility)
    min_causal = polynomial_minimum(causal)
    status = ("PASS" if min_n > 0 and min_d > 0 and min_causal >= 0
              else "INVALID_HERMITE_REPRESENTATION")
    return dict(base, status=status, mean_density=dp/h,
                necessary_convexity=n0 < dp/h < n1,
                min_density=min_n, min_dn_dmu=min_d,
                min_mu_dn_dmu_minus_n=min_causal)


def intervals(table: dict) -> list[dict]:
    result = []
    for i in range(len(table["n"])-1):
        n0, n1 = table["n"][i:i+2]
        p0, p1 = table["p"][i:i+2]
        m0, m1 = table["mu_e"][i:i+2]
        identity_changed = (table["phase"][i] != table["phase"][i+1]
                            or table["isotopes"][i] != table["isotopes"][i+1])
        mu_tol = 64*EPS*max(abs(m0), abs(m1))
        p_tol = 64*EPS*max(abs(p0), abs(p1))
        if abs(m1-m0) <= mu_tol and abs(p1-p0) <= p_tol and n1 > n0:
            record = dict(status="EXPLICIT_COMMON_TANGENT_PAIR",
                          n_left=float(n0), n_right=float(n1),
                          pressure_left=float(p0), pressure_right=float(p1),
                          mu_left=float(m0), mu_right=float(m1),
                          pressure_residual=float(p1-p0), mu_residual=float(m1-m0))
        else:
            record = hermite_check(n0, p0, m0, n1, p1, m1)
        record.update(left_row=i, right_row=i+1,
                      source_identity_changed=bool(identity_changed),
                      interpolation_authorization=("REQUIRES_EXPLICIT_INTERFACE_OR_BRIDGE"
                          if identity_changed else "SAME_SOURCE_IDENTITY"))
        result.append(record)
    return result


def summary(table: dict) -> dict:
    details = intervals(table)
    bad = [r["left_row"] for r in details
           if r["status"] not in ("PASS", "EXPLICIT_COMMON_TANGENT_PAIR")]
    runs = []
    for i in range(len(table["n"])):
        identity = (int(table["phase"][i]), table["isotopes"][i])
        if runs and runs[-1]["identity"] == identity:
            runs[-1]["last_row"] = i
        else:
            runs.append(dict(identity=identity, first_row=i, last_row=i))
    for run in runs:
        a, b = run["first_row"], run["last_row"]
        run.update(n_range=[float(table["n"][a]), float(table["n"][b])],
                   p_range=[float(table["p"][a]), float(table["p"][b])],
                   mu_e_range=[float(table["mu_e"][a]), float(table["mu_e"][b])])
    source_files = {}
    for name in ("eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo"):
        path = ROOT/table["path"]/name
        if path.exists():
            source_files[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    semantic = hashlib.sha256()
    for key in ("n", "p", "e", "mu", "phase"):
        semantic.update(np.ascontiguousarray(table[key]).tobytes())
    return dict(source_directory=table["path"], rows=len(table["n"]),
                selected_prefix_semantic_sha256=semantic.hexdigest(),
                n_range=[float(table["n"][0]), float(table["n"][-1])],
                mn_MeV=table["mn"], mp_MeV=table["mp"],
                max_abs_mu_e_minus_stored_mu_MeV=float(max(abs(table["mu_e"]-table["mu"]))),
                identity_runs=runs, invalid_hermite_intervals=bad,
                explicit_common_tangent_pairs=[r["left_row"] for r in details
                                              if r["status"] == "EXPLICIT_COMMON_TANGENT_PAIR"],
                intervals=details, source_sha256=source_files)


def bridge_candidates(lower: dict, upper: dict, lower_rows, upper_rows) -> list[dict]:
    candidates = []
    for i in lower_rows:
        for j in upper_rows:
            check = hermite_check(lower["n"][i], lower["p"][i], lower["mu_e"][i],
                                  upper["n"][j], upper["p"][j], upper["mu_e"][j])
            if check["status"] == "PASS":
                check.update(lower_row=int(i), upper_row=int(j),
                             density_width=float(upper["n"][j]-lower["n"][i]))
                candidates.append(check)
    return sorted(candidates, key=lambda x: (x["density_width"], x["upper_row"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=ROOT/"data/whole_crust_20261117/source_audit")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    report = dict(schema="ns-cowling.whole-crust-source-audit.v1",
                  input_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  source_tables={}, bridge_candidates={},
                  method="indexed CoMPOSE dimensional recovery; exact polynomial extrema for cubic P(mu_E)",
                  caveat="Representability diagnostics do not certify microscopic equilibrium or authorize interpolation across isotope boundaries.")
    reference = ROOT/"data/matched_eos_reference_20261117"
    public = cut_table(read_table(reference/"references/SPG_M1"), 217)
    old_outer = cut_table(read_table(ROOT/"data/research_finish_20261116/crust_assembly/restored_spg/SPG_M1/compose"), 547)
    report["source_tables"]["independent_outer_reused_prefix"] = summary(old_outer)
    bridge_csv = []
    for model in range(1, 6):
        name = f"SPG_M{model}"
        source = read_table(ROOT/f"data/research_matching_20261117/spg/native_129_repaired/{name}/generated")
        report["source_tables"][name+"_public_lower"] = summary(cut_table(read_table(reference/f"references/{name}"), 217))
        report["source_tables"][name+"_native"] = summary(source)
        for method in ("quintic_P", "css"):
            path = reference/f"campaign/pairs/{name}__RG_SLY4__{method}/compose"
            report["source_tables"][name+"_matched_"+method] = summary(read_table(path))
        for label, lower, lower_rows, upper_rows in (
                ("public_lower", public, range(217), range(8)),
                ("independent_outer", old_outer, [546], range(len(source["n"])) )):
            candidates = bridge_candidates(lower, source, lower_rows, upper_rows)
            key = label+"__"+name
            report["bridge_candidates"][key] = dict(count=len(candidates),
                minimum_width_candidate=candidates[0] if candidates else None,
                all_candidates=candidates,
                scope="Explicit new cubic barotrope; no energy-reference adjustment; knot endpoints fixed.")
            for row in candidates:
                bridge_csv.append(dict(lower_source=label, upper_source=name, **row))
    report["public_lower_prefixes_semantically_identical"] = len({
        report["source_tables"][f"SPG_M{m}_public_lower"]["selected_prefix_semantic_sha256"]
        for m in range(1, 6)}) == 1
    path = args.output/"SOURCE_AUDIT.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    with (args.output/"BRIDGE_CANDIDATES.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(bridge_csv[0]))
        writer.writeheader()
        writer.writerows(bridge_csv)
    print(json.dumps(dict(report=str(path), source_tables=len(report["source_tables"]),
                          bridge_candidates=len(bridge_csv),
                          public_bad_intervals=report["source_tables"]["SPG_M1_public_lower"]["invalid_hermite_intervals"],
                          independent_outer_bad_intervals=report["source_tables"]["independent_outer_reused_prefix"]["invalid_hermite_intervals"])))


if __name__ == "__main__":
    main()
