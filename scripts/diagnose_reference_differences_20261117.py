#!/usr/bin/env python3
"""Bounded, non-production controls for GM1/SLy4 public-source differences.

The source density nodes are evaluated directly; this calculation performs no
interpolation, retuning of the generator, or attribution to an unknown provider
implementation. Densities are fm^-3 and energies are MeV (including rest energy).
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq, lsq_linear

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import gm1_uniform as gm
from compose_eos import skyrme_uniform as sk
from compose_eos.cold_leptons import free_fermion


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_table(path: Path) -> dict:
    """Read just the cold scalar slice needed by these two diagnostics.

    The density axis may have any integer origin. A duplicate thermodynamic
    density index (multiple T/Yq slices) is rejected rather than averaged.
    Composition is joined by the complete (iT,iN,iYq) index tuple.
    """
    nb_tokens = (path / "eos.nb").read_text().split()
    lo, hi = map(int, nb_tokens[:2])
    axis = np.array([float(x.replace("D", "E")) for x in nb_tokens[2:]])
    if len(axis) != hi - lo + 1:
        raise ValueError(f"{path}: density-axis cardinality mismatch")
    lines = (path / "eos.thermo").read_text().splitlines()
    mn, mp = map(float, lines[0].split()[:2])
    rows = [line.split() for line in lines[1:] if line.strip()]
    keys = [tuple(map(int, row[:3])) for row in rows]
    if (len({key[0] for key in keys}) != 1 or
            len({key[2] for key in keys}) != 1 or
            any(not lo <= key[1] <= hi for key in keys)):
        raise ValueError(f"{path}: invalid scalar slice or density-axis index")
    if len({key[1] for key in keys}) != len(keys):
        raise ValueError(f"{path}: diagnostic requires one scalar T/Yq slice")
    composition = {}
    for line in (path / "eos.compo").read_text().splitlines():
        row = line.split()
        if not row:
            continue
        key = tuple(map(int, row[:3]))
        pairs = int(row[4])
        composition[key] = (int(row[3]), {
            int(row[5 + 2*k]): float(row[6 + 2*k]) for k in range(pairs)})
    order = np.argsort([axis[key[1] - lo] for key in keys])
    q = np.array([[float(x.replace("D", "E")) for x in row[3:10]]
                  for row in rows])[order]
    keys = [keys[j] for j in order]
    n = np.array([axis[key[1] - lo] for key in keys])
    if np.any(np.diff(n) <= 0):
        raise ValueError(f"{path}: nonincreasing density support")
    return {
        "path": str(path), "n": n, "p": n*q[:, 0],
        "epsilon": n*mn*(1+q[:, 6]), "mu_b_stored": mn*(1+q[:, 2]),
        "mu_q_stored": mn*q[:, 3], "mu_l_stored": mn*q[:, 4],
        "mn_header": mn, "mp_header": mp,
        "phase": np.array([composition[key][0] for key in keys]),
        "yp": np.array([composition[key][1].get(11, np.nan) for key in keys]),
        "ye": np.array([composition[key][1].get(0, np.nan) for key in keys]),
        "ymu": np.array([composition[key][1].get(1, np.nan) for key in keys]),
        "sha256": {name: sha(path/name) for name in
                   ("eos.nb", "eos.thermo", "eos.compo")},
    }


def metrics(delta, reference=None) -> dict:
    delta = np.asarray(delta, dtype=float)
    if not np.all(np.isfinite(delta)):
        raise ValueError("nonfinite comparison")
    result = {"max_abs": float(np.max(np.abs(delta))),
              "rms": float(np.sqrt(np.mean(delta**2))),
              "min_signed": float(delta.min()), "max_signed": float(delta.max())}
    if reference is not None:
        reference = np.asarray(reference, dtype=float)
        mask = reference != 0
        result["relative_nonzero_reference"] = metrics(delta[mask]/reference[mask])
    return result


def split_kinetic_delta(n, yp):
    """Analytic change when replacing common K by hbar^2/(2m_n,p)."""
    kn = sk.HBARC**2/(2*sk.MN)
    kp = sk.HBARC**2/(2*sk.MP)
    nn, np_ = n*(1-yp), n*yp
    tau_n = 0.6*sk.FERMI_FACTOR*nn**(5/3)
    tau_p = 0.6*sk.FERMI_FACTOR*np_**(5/3)
    de = (kn-sk.KINETIC_COEFFICIENT)*tau_n + (kp-sk.KINETIC_COEFFICIENT)*tau_p
    dmn = (kn-sk.KINETIC_COEFFICIENT)*sk.FERMI_FACTOR*nn**(2/3)
    dmp = (kp-sk.KINETIC_COEFFICIENT)*sk.FERMI_FACTOR*np_**(2/3)
    return de, (2/3)*de, dmn, dmp


def nuclear_control(n, yp, split):
    s = sk.nuclear_state("RG_SLY4", n*(1-yp), n*yp)
    if split:
        de, dp, dmn, dmp = split_kinetic_delta(n, yp)
        s["energy_density_mev_fm3"] += de
        s["pressure_mev_fm3"] += dp
        s["mu_n_mev"] += dmn
        s["mu_p_mev"] += dmp
    return s


def sly_state(n, split):
    if not split:
        return sk.beta_equilibrium("RG_SLY4", n)

    def evaluate(yp):
        state = nuclear_control(n, yp, split)
        mu = state["mu_n_mev"] - state["mu_p_mev"]
        e = sk.lepton_state(mu, sk.ME)
        m = sk.lepton_state(mu, sk.MMU)
        residual = n*yp-e["density_fm3"]-m["density_fm3"]
        return state, e, m, residual

    nodes = np.linspace(0., .5, 65)
    roots = []
    vals = [evaluate(y)[-1] for y in nodes]
    for a, b, fa, fb in zip(nodes[:-1], nodes[1:], vals[:-1], vals[1:]):
        if fa == 0:
            roots.append(a)
        if fa*fb < 0:
            roots.append(brentq(lambda y: evaluate(y)[-1], a, b,
                                xtol=1e-15, rtol=1e-14))
    if not roots:
        raise RuntimeError("split-kinetic control has no composition root")
    candidates = []
    for yp in roots:
        s, e, m, r = evaluate(yp)
        candidates.append({
            "n_b_fm3": float(n), "Y_p": float(yp),
            "energy_density_mev_fm3": s["energy_density_mev_fm3"]+
                e["energy_density_mev_fm3"]+m["energy_density_mev_fm3"],
            "pressure_mev_fm3": s["pressure_mev_fm3"]+
                e["pressure_mev_fm3"]+m["pressure_mev_fm3"],
            "charge_residual_fm3": float(r),
        })
    return min(candidates, key=lambda s: s["energy_density_mev_fm3"])


def diagnose_sly4(source):
    mask = source["phase"] == 3
    for field in ("yp", "ye", "ymu"):
        if not np.all(np.isfinite(source[field][mask])):
            raise ValueError("SLy4 homogeneous diagnostic needs explicit n,p,e,mu fractions")
    n = source["n"][mask]
    rows = []
    for j in np.flatnonzero(mask):
        density, yp, ye = source["n"][j], source["yp"][j], source["ye"][j]
        mue = free_fermion(density*ye, sk.ME)[2]
        states = [sly_state(density, split) for split in (False, True)]
        r = {"n_b_fm3": float(density), "source_Yp": float(yp),
             "source_pressure_mev_fm3": float(source["p"][j]),
             "source_energy_per_baryon_mev": float(source["epsilon"][j]/density)}
        for label, split, state in zip(("common_K", "split_K"), (False, True), states):
            nuclear = nuclear_control(density, yp, split)
            r[label] = {
                "Yp_difference": float(state["Y_p"]-yp),
                "pressure_difference_mev_fm3": float(state["pressure_mev_fm3"]-source["p"][j]),
                "energy_per_baryon_difference_mev": float(state["energy_density_mev_fm3"]/density-source["epsilon"][j]/density),
                "source_composition_beta_residual_mev": float(nuclear["mu_n_mev"]-nuclear["mu_p_mev"]-mue),
            }
        rows.append(r)
    summaries = {}
    for scope, selected in (("full_public_homogeneous", rows),
                            ("campaign_0p08_to_0p5", [r for r in rows if .08 <= r["n_b_fm3"] <= .5])):
        summaries[scope] = {"rows": len(selected), "density_range_fm3":
                            [selected[0]["n_b_fm3"], selected[-1]["n_b_fm3"]]}
        for label in ("common_K", "split_K"):
            summaries[scope][label] = {name: metrics(
                [r[label][name] for r in selected],
                [r["source_pressure_mev_fm3"] for r in selected]
                if name == "pressure_difference_mev_fm3" else None)
                for name in selected[0][label]}
    delta = np.array([r["common_K"]["energy_per_baryon_difference_mev"] for r in rows])
    sign_changes = int(np.count_nonzero(delta[1:]*delta[:-1] < 0))
    return {
        "status": "CONTROL_EXECUTED_PROVIDER_ATTRIBUTION_UNRESOLVED",
        "scope": "exact public density nodes; explicit homogeneous phase=3; no interpolation",
        "kinetic_coefficients_mev_fm2": {
            "common": sk.KINETIC_COEFFICIENT,
            "split_neutron": sk.HBARC**2/(2*sk.MN),
            "split_proton": sk.HBARC**2/(2*sk.MP)},
        "control_definition": "Only nuclear kinetic K is split by species; the published interaction, rest masses and lepton prescription are unchanged.",
        "source_composition_beta_definition": "mu_n(model,n,Yp_source)-mu_p(model,n,Yp_source)-mu_e(Ye_source*n); rounded source fractions limit accuracy",
        "constant_energy_shift_diagnostic": {
            "source_minus_model_sign_changes": sign_changes,
            "range_mev": float(np.ptp(delta)),
            "statement": "A constant energy-per-baryon offset cannot remove a varying residual and changes no pressure. No offset was applied."},
        "summaries": summaries, "rows": rows,
        "remaining_issue": "A kinetic prescription can be tested and its composition effect measured; this does not identify the provider prescription or explain every pressure/energy residual.",
    }


GM_NAMES = ("C_SIGMA", "C_OMEGA", "C_RHO", "B", "C")


@contextmanager
def gm_parameters(values=None, mass=None):
    saved = {name: getattr(gm, name) for name in (*GM_NAMES, "MASS")}
    try:
        if values is not None:
            for name, value in zip(GM_NAMES, values):
                setattr(gm, name, float(value))
        if mass is not None:
            gm.MASS = float(mass)
        yield
    finally:
        for name, value in saved.items():
            setattr(gm, name, value)


def gm_states(densities, parameters=None, mass=None):
    with gm_parameters(parameters, mass):
        return [gm.beta_state(float(n), include_muons=False) for n in densities]


def diagnose_gm1(source):
    mask = (source["n"] >= .08) & (source["n"] <= .5) & (source["phase"] == 3)
    n, pref, eref = source["n"][mask], source["p"][mask], source["epsilon"][mask]/source["n"][mask]
    base = np.array([getattr(gm, name) for name in GM_NAMES])
    halfbins = np.array([.005, .0005, .0005, 5e-7, 5e-7])
    baseline = gm_states(n)
    header_control = gm_states(n, mass=source["mn_header"])
    train = np.arange(0, len(n), 8)
    if train[-1] != len(n)-1:
        train = np.r_[train, len(n)-1]
    holdout = np.setdiff1d(np.arange(len(n)), train)

    def vector(states):
        return np.r_[[s["epsilon"]/s["n_b"] for s in states], [s["pressure"] for s in states]]

    center = vector([baseline[i] for i in train])
    target = np.r_[eref[train], pref[train]]-center
    jacobian = []
    for k in range(5):
        step = np.eye(5)[k]*halfbins[k]
        jacobian.append((vector(gm_states(n[train], base+step))-
                         vector(gm_states(n[train], base-step)))/2)
    jacobian = np.array(jacobian).T
    fit = lsq_linear(jacobian, target, bounds=(-1, 1), tol=1e-12)
    values = base+fit.x*halfbins
    fitted = gm_states(n, parameters=values)
    rows = []
    for i, density in enumerate(n):
        row = {"n_b_fm3": float(density), "rounding_fit_train": bool(i in train),
               "source_pressure_mev_fm3": float(pref[i]),
               "source_energy_per_baryon_mev": float(eref[i])}
        for label, states in (("published_rounded", baseline),
                              ("header_mass_as_kernel_control", header_control),
                              ("bounded_rounding_control", fitted)):
            row[label] = {
                "pressure_difference_mev_fm3": float(states[i]["pressure"]-pref[i]),
                "energy_per_baryon_difference_mev": float(states[i]["epsilon"]/density-eref[i]),
            }
        rows.append(row)
    summary = {}
    for label in ("published_rounded", "header_mass_as_kernel_control", "bounded_rounding_control"):
        summary[label] = {name: metrics([r[label][name] for r in rows], pref
                                        if name == "pressure_difference_mev_fm3" else None)
                          for name in rows[0][label]}
        summary[label]["heldout_energy_per_baryon_difference_mev"] = metrics(
            [rows[i][label]["energy_per_baryon_difference_mev"] for i in holdout])
    return {
        "status": "TESTED_CONVENTION_AND_ROUNDING_CONTROLS_SOURCE_DIFFERENCE_UNRESOLVED",
        "density_range_fm3": [float(n[0]), float(n[-1])], "rows_count": len(n),
        "species": "npe in every control; public provider specifies npe",
        "kernel_masses_mev": {"baseline_common": gm.MASS,
                              "control_common": source["mn_header"]},
        "rounding_control": {
            "definition": "A linearized least-squares search within each printed Table II half-digit interval, followed by a full nonlinear evaluation at its result. This is a diagnostic, not an adopted retuning or a global impossibility proof.",
            "weight_definition": "Residual E/A divided by 1 MeV and P divided by 1 MeV fm^-3 have equal numerical weights; deterministic every-eighth training rows plus final row; remaining rows held out.",
            "parameters": {name: {"published": float(v), "half_last_digit": float(h),
                                  "controlled_value": float(c), "scaled_change": float(u)}
                           for name, v, h, c, u in zip(GM_NAMES, base, halfbins, values, fit.x)},
            "training_rows": len(train), "heldout_rows": len(holdout),
            "linearized_residual_norm_ratio": float(np.linalg.norm(target-jacobian@fit.x)/np.linalg.norm(target)),
            "active_bounds": fit.active_mask.tolist(),
        },
        "summaries": summary, "rows": rows,
        "conclusion_scope": "Header masses do not establish kernel masses. A tested rounding box is not evidence that all source differences are explained by published precision.",
    }


def analytic_control():
    residuals = []
    for n, yp in ((.06, .03), (.2, .06), (.5, .08), (1., .1)):
        h = n*1e-5
        de, dp, _, _ = split_kinetic_delta(n, yp)
        derivative = (split_kinetic_delta(n+h, yp)[0]/(n+h)-
                      split_kinetic_delta(n-h, yp)[0]/(n-h))/(2*h)
        residuals.append(abs(n*n*derivative-dp)/max(abs(dp), 1e-30))
    maximum = max(residuals)
    if maximum > 1e-8:
        raise AssertionError("split kinetic energy/pressure control failed")
    return {"identity": "delta P = n^2*d(delta epsilon/n)/dn = (2/3)*delta epsilon at fixed composition",
            "density_composition_points": 4, "max_relative_residual": maximum,
            "status": "NUMERICALLY_CHECKED"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("use a fresh diagnosis output directory")
    sources = {model: source_table(args.source_root/model) for model in ("RG_SLY4", "GM_GM1")}
    result = {
        "schema_version": 1, "source_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "calculation_role": "independent Python diagnostic; native production EOS are not modified",
        "source_root": str(args.source_root.resolve()),
        "source_sha256": {m: s["sha256"] for m, s in sources.items()},
        "implementation_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
             (Path(__file__).resolve(), ROOT/"python_frontend/compose_eos/skyrme_uniform.py",
              ROOT/"python_frontend/compose_eos/gm1_uniform.py",
              ROOT/"python_frontend/compose_eos/cold_leptons.py")},
        "analytic_control": analytic_control(),
        "RG_SLY4": diagnose_sly4(sources["RG_SLY4"]),
        "GM_GM1": diagnose_gm1(sources["GM_GM1"]),
    }
    args.output.mkdir(parents=True)
    (args.output/"DIAGNOSIS.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"output": str(args.output), "status": "CONTROLS_EXECUTED",
                      "analytic_control": result["analytic_control"],
                      "SLy4": result["RG_SLY4"]["summaries"],
                      "GM1": result["GM_GM1"]["summaries"]}, indent=2))


if __name__ == "__main__":
    main()
