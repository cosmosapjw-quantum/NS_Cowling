"""Installed CLI for causality-aware TOV sequences on cold CompOSE EOS tables."""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

import numpy as np

from .io import read_compose
from .tov import (
    BarotropicEOS, TOVConfig, CausalDomain, find_causal_domain,
    first_law_nb_consistency, solve_tov_sequence, solve_tov_star,
)

FULL_CORE_MODELS = {
    "APR_APR", "BBB_BHF2", "GM_GM1_fixed", "RG_SLY4_fixed", "RG_KDE0v1_fixed"
}


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, sort_keys=True, allow_nan=False)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise


def _json_safe(x):
    if isinstance(x, dict):
        return {str(k): _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_safe(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def _find_model_dirs(root: Path) -> dict[str, Path]:
    found = {}
    for thermo in root.rglob("eos.thermo"):
        directory = thermo.parent
        if directory.name in FULL_CORE_MODELS:
            found[directory.name] = directory
    missing = FULL_CORE_MODELS - set(found)
    if missing:
        raise FileNotFoundError(f"full-core EOS directories missing under {root}: {sorted(missing)}")
    return found


def _provider_mr(directory: Path):
    path = directory / "eos.mr"
    if not path.is_file():
        return None
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            values = [float(t) for t in line.split()]
        except ValueError:
            continue
        if len(values) >= 2 and np.isfinite(values[0]) and np.isfinite(values[1]):
            rows.append((values[0], values[1]))
    if not rows:
        return None
    arr = np.asarray(rows)
    radius, mass = arr[:, 0], arr[:, 1]
    imax = int(np.argmax(mass))
    m_branch, r_branch = mass[:imax + 1], radius[:imax + 1]
    r14 = None
    idx = np.flatnonzero((m_branch[:-1] - 1.4) * (m_branch[1:] - 1.4) <= 0.0)
    if idx.size:
        j = int(idx[0])
        if m_branch[j + 1] != m_branch[j]:
            t = (1.4 - m_branch[j]) / (m_branch[j + 1] - m_branch[j])
            r14 = float(r_branch[j] + t * (r_branch[j + 1] - r_branch[j]))
    return {"Mmax_solar": float(np.max(mass)), "R14_km": r14,
            "n_rows": int(len(rows)), "source": "provider eos.mr; parsed, not recomputed"}


def _run_model(name: str, directory: Path, coarse_points: int, fine_points: int,
               rtol: float, central_nb_min: float) -> dict:
    data = read_compose(str(directory))
    eos = BarotropicEOS(data.nb, data.p, data.eps)
    domain = find_causal_domain(eos, tolerance=1e-8)
    config = TOVConfig(rtol=rtol)
    coarse = solve_tov_sequence(eos, n_points=coarse_points,
                                central_nb_min=central_nb_min, config=config,
                                causal_domain=domain)
    fine = solve_tov_sequence(eos, n_points=fine_points,
                              central_nb_min=central_nb_min, config=config,
                              causal_domain=domain)

    fine_stars = [s for s in fine["stars"] if s["status"] == "SUCCESS"]
    if fine_stars:
        max_star = max(fine_stars, key=lambda s: s["mass_solar"])
        nc_ref = float(max_star["central_nb"])
        tol_runs = [solve_tov_star(eos, nc_ref, config=TOVConfig(rtol=tol),
                                   causal_domain=domain).to_dict()
                    for tol in (1e-6, 2e-8, 1e-10)]
        tolerance_study = {"status": "COMPLETED", "central_nb_fixed_fm3": nc_ref,
                           "runs": tol_runs}
    else:
        tolerance_study = {"status": "NO_SUCCESSFUL_FINE_STARS", "runs": []}

    fl_8 = first_law_nb_consistency(eos, oversampling=8)
    fl_16 = first_law_nb_consistency(eos, oversampling=16)
    fl_32 = first_law_nb_consistency(eos, oversampling=32)
    fl = {
        "preferred": fl_32,
        "quadrature_convergence": {
            "oversampling_8_max_abs_log_ratio": fl_8["max_abs_log_ratio"],
            "oversampling_16_max_abs_log_ratio": fl_16["max_abs_log_ratio"],
            "oversampling_32_max_abs_log_ratio": fl_32["max_abs_log_ratio"],
            "abs_change_16_to_32": abs(fl_32["max_abs_log_ratio"] - fl_16["max_abs_log_ratio"]),
        },
    }

    acausal_diagnostic = None
    if domain.status == "ACASUAL_CUTOFF":
        full = CausalDomain(
            status="ACASUAL_DIAGNOSTIC_ONLY", nb_max=float(data.nb[-1]),
            eps_max=float(data.eps[-1]), pressure_max=float(data.p[-1]),
            cs2_at_limit=float(eos.cs2(data.eps[-1])),
            reason="Intentional extension into a known acausal segment; diagnostic only.")
        seq = solve_tov_sequence(eos, n_points=fine_points,
                                 central_nb_min=central_nb_min, config=config,
                                 causal_domain=full)
        acausal_diagnostic = {
            "status": "NONPHYSICAL_DIAGNOSTIC_ONLY",
            "full_table_max_cs2": float(np.max(eos.cs2_on_grid())),
            "sequence": seq,
            "warning": "Includes c_s^2 > 1; do not interpret as a physical causal EOS prediction.",
        }

    c_m = coarse["summary"].get("Mmax_sampled_solar")
    f_m = fine["summary"].get("Mmax_sampled_solar")
    c_r = coarse["summary"].get("R14_km_linear_on_stable_branch")
    f_r = fine["summary"].get("R14_km_linear_on_stable_branch")
    return {
        "model": name,
        "source_directory": directory.name,
        "source_grid_points": int(data.n),
        "nb_range_fm3": [float(data.nb[0]), float(data.nb[-1])],
        "causal_domain": domain.to_dict(),
        "coarse_sequence": coarse,
        "fine_sequence": fine,
        "grid_convergence": {
            "Mmax_abs_change_solar": abs(f_m - c_m) if c_m is not None and f_m is not None else None,
            "Mmax_rel_change": abs(f_m - c_m) / abs(f_m) if c_m is not None and f_m not in (None, 0.0) else None,
            "R14_abs_change_km": abs(f_r - c_r) if c_r is not None and f_r is not None else None,
        },
        "first_law_number_density_consistency": fl,
        "ode_tolerance_study": tolerance_study,
        "acausal_full_table_diagnostic": acausal_diagnostic,
        "provider_eos_mr": _provider_mr(directory),
        "interpretation": {
            "Mmax": "sampled maximum only; upper-domain-limited sequences do not establish the global maximum",
            "R14": "linear interpolation on the sampled stable branch; null when 1.4 Msun is not bracketed",
            "radius": fine["surface_definition"],
            "causality": "higher-density branch excluded at first PCHIP cs2=1 crossing; no acausal extrapolation",
        },
    }


def build_parser():
    p = argparse.ArgumentParser(description="Run causal TOV sequences on the five full-core EOS tables.")
    p.add_argument("input_root", type=Path, help="root containing extracted EOS directories")
    p.add_argument("--output", type=Path, required=True, help="checkpointed JSON result path")
    p.add_argument("--coarse-points", type=int, default=24)
    p.add_argument("--fine-points", type=int, default=36)
    p.add_argument("--rtol", type=float, default=2e-8)
    p.add_argument("--central-nb-min", type=float, default=0.14)
    return p


def main(argv=None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    if args.coarse_points < 4 or args.fine_points < 4:
        ap.error("coarse/fine point counts must be >= 4")
    if not args.input_root.is_dir():
        ap.error(f"input_root is not a directory: {args.input_root}")
    dirs = _find_model_dirs(args.input_root)
    payload = {
        "schema_version": "3.0",
        "method": "log-log PCHIP barotrope; causal cutoff; SI TOV integration",
        "input_root": args.input_root.resolve().name,
        "settings": {"coarse_points": args.coarse_points,
                     "fine_points": args.fine_points, "rtol": args.rtol,
                     "central_nb_min_fm3": args.central_nb_min},
        "models": {}, "status": "RUNNING",
    }
    _atomic_json(args.output, payload)
    for name in sorted(FULL_CORE_MODELS):
        print(f"[TOV] {name}", flush=True)
        try:
            payload["models"][name] = _run_model(
                name, dirs[name], args.coarse_points, args.fine_points,
                args.rtol, args.central_nb_min)
            payload["models"][name]["run_status"] = "COMPLETED"
        except Exception as exc:
            payload["models"][name] = {
                "run_status": "ERROR", "error_type": type(exc).__name__,
                "error": str(exc), "source_directory": dirs[name].name}
        _atomic_json(args.output, _json_safe(payload))
        print(f"[TOV] {name}: {payload['models'][name]['run_status']}", flush=True)
    payload["status"] = ("COMPLETED" if all(m.get("run_status") == "COMPLETED"
                                             for m in payload["models"].values())
                         else "PARTIAL")
    _atomic_json(args.output, _json_safe(payload))
    print(json.dumps({k: v.get("run_status") for k, v in payload["models"].items()}, indent=2))
    return 0 if payload["status"] == "COMPLETED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
