#!/usr/bin/env python3
"""Minimal energy-primary reconstruction of diagnosed public EOS inconsistencies.

Fortran owns the production derivative calculation. Python diagnoses whether a
pressure repair is justified, preserves untouched fields, and records every
change. A reconstructed candidate is not automatically physically admissible.
No energy shift, model refit, composition repair, clipping, or phase bridging is
performed here. A caller must select documented source domains before calling.
"""
from __future__ import annotations

import argparse
import copy
from decimal import Decimal, InvalidOperation
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any

import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos.eos_data import EOSData

VERSION = "1.0.0"
NATIVE_SOURCE = ROOT / "fortran_eos/app/reconstruct_public_thermo.f90"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native_executable() -> Path:
    override = os.environ.get("COMPOSE_PUBLIC_RECON_EXE")
    if override:
        result = Path(override).resolve()
        if not result.is_file():
            raise FileNotFoundError(f"COMPOSE_PUBLIC_RECON_EXE not found: {result}")
        return result
    folder = ROOT / "build/public_thermo_reconstruction"
    result, signature = folder / "reconstruct_public_thermo", folder / "SOURCE_SHA256"
    source_hash = sha(NATIVE_SOURCE)
    if result.is_file() and signature.is_file() and signature.read_text().strip() == source_hash:
        return result
    compiler = shutil.which(os.environ.get("FC", "gfortran"))
    if compiler is None:
        raise RuntimeError("native public reconstruction needs gfortran; no Python production fallback")
    folder.mkdir(parents=True, exist_ok=True)
    command = [compiler, "-std=f2018", "-O2", "-ffree-line-length-none", "-fimplicit-none",
               "-fcheck=all", str(NATIVE_SOURCE), "-o", str(result)]
    run = subprocess.run(command, capture_output=True, text=True)
    (folder / "BUILD.json").write_text(json.dumps({"command": command, "returncode": run.returncode,
                                                 "stdout": run.stdout, "stderr": run.stderr,
                                                 "source_sha256": source_hash}, indent=2) + "\n")
    if run.returncode:
        raise RuntimeError(f"native reconstruction build failed: {run.stderr}")
    signature.write_text(source_hash + "\n")
    return result


def native_derivatives(n: np.ndarray, binding: np.ndarray, binary: Path | None = None) -> np.ndarray:
    """Columns: dF/dln n, d²F/dln n², independent n*dF/dn control."""
    binary = Path(binary) if binary else native_executable()
    with tempfile.TemporaryDirectory(prefix="public-thermo-reconstruct-") as td:
        source, target = Path(td) / "source.txt", Path(td) / "derivatives.txt"
        with source.open("w") as stream:
            stream.write(f"{len(n)}\n")
            np.savetxt(stream, np.column_stack([n, binding]), fmt="%.17e")
        run = subprocess.run([str(binary), str(source), str(target)], capture_output=True, text=True)
        if run.returncode:
            raise RuntimeError(f"native reconstruction failed ({run.returncode}): {run.stderr or run.stdout}")
        answer = np.loadtxt(target, ndmin=2)
        if answer.shape != (len(n), 3) or not np.all(np.isfinite(answer)):
            raise ValueError("native reconstruction returned invalid derivatives")
        return answer


def _quantum(token: str) -> float:
    try:
        return float(Decimal(1).scaleb(Decimal(token.replace("D", "E").replace("d", "e")).as_tuple().exponent)) / 2
    except (InvalidOperation, OverflowError):
        return float("nan")


def source_rounding(data: EOSData, raw_source_dir: Path | None) -> dict[str, np.ndarray]:
    """Estimate half-last-digit uncertainty, matching physical density values.

    This is a rounding sensitivity estimate, not an uncertainty on the nuclear
    model. Missing raw tokens are explicit and never interpreted as zero noise.
    """
    result = {k: np.full(data.n, np.nan) for k in ("n", "q1", "q3", "q7")}
    if raw_source_dir is None or not (Path(raw_source_dir) / "eos.thermo").is_file():
        return result
    folder = Path(raw_source_dir)
    tokens = (folder / "eos.nb").read_text().split()
    origin, stop = int(tokens[0]), int(tokens[1])
    axis = np.array([float(x.replace("D", "E").replace("d", "e")) for x in tokens[2:]])
    if len(axis) != stop-origin+1:
        raise ValueError("source-rounding density axis is malformed")
    by_n: dict[float, dict[str, float]] = {}
    for line in (folder / "eos.thermo").read_text().splitlines()[1:]:
        row = line.split()
        if not row:
            continue
        i = int(row[1])-origin
        by_n[float(axis[i])] = {"n": _quantum(tokens[2+i]), "q1": _quantum(row[3]),
                               "q3": _quantum(row[5]), "q7": _quantum(row[9])}
    for i, n in enumerate(data.nb):
        record = by_n.get(float(n))
        if record is not None:
            for field, value in record.items():
                result[field][i] = value
    return result


def _segments(phase: np.ndarray) -> list[np.ndarray]:
    edges = np.r_[0, 1+np.flatnonzero(np.diff(phase)), len(phase)]
    return [np.arange(a, b) for a, b in zip(edges[:-1], edges[1:])]


def _quadratic_min(c0: np.ndarray, c1: np.ndarray, c2: np.ndarray) -> np.ndarray:
    result = np.minimum(c0, c0+c1+c2)
    nonzero = np.abs(c2) > np.finfo(float).tiny
    t = np.divide(-c1, 2*c2, out=np.full_like(c1, np.nan), where=nonzero)
    inner = nonzero & (t > 0) & (t < 1)
    result[inner] = np.minimum(result[inner], c0[inner]+c1[inner]*t[inner]+c2[inner]*t[inner]**2)
    return result


def hermite_interval_diagnostics(n: np.ndarray, f: np.ndarray, q1: np.ndarray) -> dict[str, Any]:
    """Independent whole-cell checks of the actual native P(mu) cubic family.

    Used only as a diagnostic; the production Maxwell kernel repeats its own
    acceptance calculation. Index j denotes the interval [j,j+1], zero-based.
    """
    p, mu = n*q1, f+q1
    h = np.diff(mu)
    regular = (h > 0) & (np.diff(n) > 0) & (np.diff(p) > 0)
    mechanical = np.zeros(len(h), bool)
    causal = np.zeros(len(h), bool)
    minimum_dn = np.full(len(h), np.nan)
    minimum_h = np.full(len(h), np.nan)
    if np.any(regular):
        i = np.flatnonzero(regular)
        hh = h[i]
        c = 3*np.diff(p)[i]-hh*(2*n[i]+n[i+1])
        d = -2*np.diff(p)[i]+hh*(n[i]+n[i+1])
        dn0, dn1 = 2*c/hh**2, (2*c+6*d)/hh**2
        min_n = _quadratic_min(n[i], 2*c/hh, 3*d/hh)
        mh = _quadratic_min(2*mu[i]*c/hh**2-n[i], 6*mu[i]*d/hh**2, 3*d/hh)
        md = np.minimum(dn0, dn1)
        scale = np.maximum.reduce([np.abs(n[i]), np.abs(n[i+1]),
                                   np.abs(mu[i]*dn0), np.abs(mu[i+1]*dn1)])
        mechanical[i] = (min_n > 0) & (md > 0)
        causal[i] = mh >= -1e-8*scale
        minimum_dn[i], minimum_h[i] = md, mh
    bad_m = np.flatnonzero(~regular | ~mechanical)
    bad_c = np.flatnonzero(regular & mechanical & ~causal)
    good = regular & mechanical & causal
    ranges = []
    first = None
    for i, passed in enumerate(np.r_[good, False]):
        if passed and first is None:
            first = i
        if not passed and first is not None:
            ranges.append({"first_row": first, "last_row": i,
                           "density_fm3": [float(n[first]), float(n[i])], "rows": i-first+1})
            first = None
    return {"all_intervals_pass": bool(np.all(good)), "intervals": len(h),
            "mechanical_failure_intervals": bad_m.tolist(),
            "causality_failure_intervals": bad_c.tolist(),
            "mechanical_failure_density_fm3": n[bad_m].tolist(),
            "causality_failure_density_fm3": n[bad_c].tolist(),
            "contiguous_admissible_ranges": ranges,
            "scope": "independent diagnostic of P(mu) Hermite cells; native kernel remains authoritative"}


def reconstruct_data(data: EOSData, *, native_binary: Path | None = None,
                     raw_source_dir: Path | str | None = None,
                     model: str | None = None) -> tuple[EOSData, dict[str, Any]]:
    """Return the minimally repaired candidate and an auditable receipt.

    All pressure changes require a measured source inconsistency.
    Reconstruction never repairs Q5 or certifies beta equilibrium.
    """
    if not isinstance(data, EOSData):
        raise TypeError("expected EOSData")
    if data.temperature_mev != 0 or np.any(data.Q2 != 0):
        raise ValueError("energy-primary cold reconstruction requires T=0 and Q2=0")
    if np.any(~np.isfinite(data.nb)) or np.any(data.nb <= 0) or np.any(np.diff(data.nb) <= 0):
        raise ValueError("source density grid must be finite, positive and increasing")
    if np.any(~np.isfinite(data.Q7)) or np.max(np.abs(data.Q6-data.Q7)) > 2e-6:
        raise ValueError("source must provide consistent cold energy Q6=Q7")
    binary = Path(native_binary) if native_binary else native_executable()
    out = copy.deepcopy(data)
    rounding = source_rounding(data, Path(raw_source_dir) if raw_source_dir else None)
    parts: list[dict[str, Any]] = []
    warnings: list[str] = []
    beta_supported = bool(data.beta_equilibrated and np.max(np.abs(data.Q5)) <= 1e-8)
    if not beta_supported:
        warnings.append("SOURCE_BETA_CONVENTION_UNRESOLVED: Q3 closure left unchanged; Q5 is not repaired")
    for indices in _segments(np.asarray(data.Iphase)):
        n = data.nb[indices]
        binding = data.mn*data.Q7[indices]
        f = data.mn+binding
        original_q = data.Q1[indices]
        part: dict[str, Any] = {"phase": int(data.Iphase[indices[0]]), "first_row": int(indices[0]),
                                "last_row": int(indices[-1]), "rows": len(indices),
                                "density_fm3": [float(n[0]), float(n[-1])]}
        part["original_intervals"] = hermite_interval_diagnostics(n, f, original_q)
        # These named public crust datasets attach the BSk22 outer crust under
        # one Iphase=0 label. It contains discrete nuclear-species changes; a
        # single smooth energy spline across it is not a thermodynamic repair.
        # Preserve that inherited block until its internal junction laws have
        # a separate source-backed representation. Campaigns select phase 7
        # inner-crust support before calling this helper.
        if model and model.startswith(("SPG_", "GPPVA_")) and part["phase"] == 0:
            part["decision"] = "UNCHANGED_INHERITED_OUTER_CRUST_TOPOLOGY"
            part["reconstructed_intervals"] = part["original_intervals"]
            parts.append(part)
            continue
        if len(indices) < 7:
            part["decision"] = "UNCHANGED_TOO_SHORT_TO_DIAGNOSE"
            part["reconstructed_intervals"] = part["original_intervals"]
            parts.append(part)
            continue
        deriv = native_derivatives(n, binding, binary)
        proposed, alternate = deriv[:, 0], deriv[:, 2]
        spread = np.abs(proposed-alternate)
        raw_sigma_n, raw_sigma_f = rounding["n"][indices], data.mn*rounding["q7"][indices]
        raw_sigma_q = rounding["q1"][indices]
        printing_known = bool(np.all(np.isfinite(raw_sigma_n+raw_sigma_f+raw_sigma_q)))
        if printing_known:
            x = np.log(n/n[0])
            # Derivative interpolation is linear in the energy samples. The
            # absolute coefficient matrix gives the exact fixed-grid bound;
            # density rounding enters by a first-order equivalent F error.
            sensitivity = CubicSpline(x, np.eye(len(n)), axis=0)(x, 1)
            sample_sigma = raw_sigma_f+np.abs(proposed)*raw_sigma_n/n
            noise = np.abs(sensitivity)@sample_sigma+raw_sigma_q
        else:
            # Explicit conservative floating-point floor when token precision
            # is unavailable. Do not advertise a printed-precision diagnosis.
            min_dx = np.min(np.diff(np.log(n)))
            noise = np.full(len(n), 64*np.finfo(float).eps*max(1., np.max(np.abs(f)))/min_dx)
        tolerance = 5*spread+8*noise+1e-7*np.maximum(np.abs(original_q), 1e-12)
        discrepancy = proposed-original_q
        significant = np.abs(discrepancy) > tolerance
        interior = np.zeros(len(n), bool)
        interior[2:-2] = True
        significant &= interior
        # An unverified beta/lepton convention invalidates the one-charge
        # first-law premise. Its residual may be measured but not "repaired".
        diagnose_pressure = bool(beta_supported and np.count_nonzero(significant) >= 3)
        repair_pressure = diagnose_pressure
        if repair_pressure:
            out.Q1[indices] = proposed
        changed_q1 = out.Q1[indices]-original_q
        mu_expected = f+out.Q1[indices]
        if data.Il == 0:
            mu_expected = mu_expected-data.charge_fraction_grid*data.mn*data.Q4[indices]
        euler_delta = mu_expected-data.mn*(1+data.Q3[indices])
        q3_sigma = rounding["q3"][indices]*data.mn
        q3_sigma = np.where(np.isfinite(q3_sigma), q3_sigma, 1e-10)
        euler_bad = np.abs(euler_delta) > np.maximum(8*q3_sigma, 1e-9)
        repair_euler = bool(beta_supported and (repair_pressure or np.any(euler_bad)))
        if repair_euler:
            out.Q3[indices] = mu_expected/data.mn-1
        decision = ("UNCHANGED_BETA_CONVENTION_UNRESOLVED" if not beta_supported else
                    "ENERGY_PRIMARY_PRESSURE_AND_EULER_REPAIR" if repair_pressure and repair_euler else
                    "ENERGY_PRIMARY_PRESSURE_REPAIR_Q3_CONVENTION_UNRESOLVED" if repair_pressure else
                    "ALGEBRAIC_EULER_REPAIR_ONLY" if repair_euler else "UNCHANGED_THERMO")
        part.update(decision=decision, pressure_inconsistency_diagnosed=diagnose_pressure,
                    printing_precision_known=printing_known,
                    significant_first_law_rows=int(np.count_nonzero(significant)),
                    significant_first_law_row_indices=indices[significant].tolist(),
                    max_abs_first_law_q1_discrepancy_MeV=float(np.max(np.abs(discrepancy))),
                    max_abs_interpolation_control_q1_MeV=float(np.max(spread)),
                    max_printing_sensitivity_q1_MeV=float(np.max(noise)),
                    max_relative_pressure_change=float(np.max(np.abs(changed_q1)/np.maximum(np.abs(original_q), 1e-30))),
                    max_abs_euler_source_discrepancy_MeV=float(np.max(np.abs(euler_delta))),
                    changed_pressure_rows=int(np.count_nonzero(changed_q1)),
                    changed_mu_rows=int(np.count_nonzero(out.Q3[indices] != data.Q3[indices])))
        part["reconstructed_intervals"] = hermite_interval_diagnostics(n, f, out.Q1[indices])
        parts.append(part)
    immutable = ("nb", "mn", "mp", "Q2", "Q4", "Q5", "Q6", "Q7", "q_extra", "Iphase", "Yn", "Yp",
                 "Ye", "Ymu", "Aav", "Zav", "Yh", "Il", "iT", "iYq", "inb_start", "beta_equilibrated")
    preservation = {name: bool(np.array_equal(np.asarray(getattr(data, name)), np.asarray(getattr(out, name)), equal_nan=True))
                    for name in immutable}
    if not all(preservation.values()):
        raise AssertionError("minimal reconstruction changed an immutable source field")
    changed = [name for name in ("Q1", "Q3") if np.any(getattr(data, name) != getattr(out, name))]
    all_shape = all(part["reconstructed_intervals"]["all_intervals_pass"] for part in parts)
    report = {"schema": "ns_cowling.public_thermo_reconstruction.v1", "version": VERSION,
              "model": model, "status": "READY" if changed else "UNCHANGED_THERMO",
              "family": "public_thermo_reconstructed", "source_family": "public_original",
              "orchestration_source_sha256": sha(Path(__file__)),
              "native_source_sha256": sha(NATIVE_SOURCE), "native_binary_sha256": sha(binary),
              "raw_source_dir": str(raw_source_dir) if raw_source_dir else None,
              "changed_fields": changed, "immutable_field_exact_preservation": preservation,
              "physical_energy_offset_MeV_per_baryon": 0.0, "density_grid_changed": False,
              "nuclear_parameters_changed": False, "composition_changed": False,
              "interpolation": "phase-segmented C2 not-a-knot cubic F(log(n_B/n_0)); no smoothing",
              "first_law": "Q1=dF/dln n_B; P=n_B Q1; epsilon=n_B F; mu_E=F+Q1",
              "pressure_repair_trigger": "at least three interior discrepancies >5*coordinate-control +8*printing-sensitivity +1e-7*abs(Q1)",
              "readiness_is_not_admission": True,
              "native_maxwell_source_admission": "NOT_TESTED_BY_THIS_HELPER",
              "all_reconstructed_intervals_diagnostic_pass": all_shape,
              "ready_for_matching": bool(all_shape and beta_supported),
              "source_beta_convention_supported": beta_supported,
              "phase_boundaries_interpolated": False, "warnings": warnings, "segments": parts}
    return out, report


def reconstruct_source(source: Path, output: Path, *, raw_source_dir: Path | None = None,
                       model: str | None = None, native_binary: Path | None = None) -> dict:
    from compose_eos import read_compose, write_compose
    data = read_compose(source)
    rebuilt, report = reconstruct_data(data, raw_source_dir=raw_source_dir or source,
                                       model=model, native_binary=native_binary)
    output.mkdir(parents=True, exist_ok=False)
    write_compose(output, rebuilt)
    report["input_file_sha256"] = {p.name: sha(p) for p in source.iterdir() if p.is_file() and p.name.startswith("eos.")}
    report["output_file_sha256"] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
    (output / "RECONSTRUCTION.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--raw-source", type=Path)
    parser.add_argument("--model")
    parser.add_argument("--native-binary", type=Path)
    args = parser.parse_args()
    report = reconstruct_source(args.source, args.output, raw_source_dir=args.raw_source,
                                model=args.model, native_binary=args.native_binary)
    print(json.dumps({"model": report["model"], "status": report["status"],
                      "changed_fields": report["changed_fields"],
                      "ready_for_matching": report["ready_for_matching"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
