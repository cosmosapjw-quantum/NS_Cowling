#!/usr/bin/env python3
"""Compare physical cold EOS quantities without changing either table.

This is an independent, bounded raw-CompOSE diagnostic reader, not a production
EOS reader or a thermodynamic admission test.  It supports one T=0, Yq slice,
arbitrary axis origins and the standard Q1..Q7 encoding.  Comparisons use
reference nodes inside pairs of contiguous, explicitly selected phase segments.
Public nuclear-group tuples are deliberately not interpreted here.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.interpolate import CubicSpline, PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
THERMO_FIELDS = ("pressure", "epsilon", "mu_b", "enthalpy_b")
PARTICLES = {10: "y_n_free", 11: "y_p_free", 0: "y_e", 1: "y_mu"}
UNITS = {"n_b": "fm^-3", "pressure": "MeV fm^-3", "epsilon": "MeV fm^-3",
         "mu_b": "MeV", "enthalpy_b": "MeV", "energy_per_baryon": "MeV",
         **{v: "1" for v in PARTICLES.values()}}


@dataclass
class ColdTable:
    path: Path
    n: np.ndarray
    phase: np.ndarray
    fields: dict[str, np.ndarray]
    metadata: dict


def _number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def _lines(path: Path) -> list[list[str]]:
    return [line.split("#", 1)[0].split() for line in path.read_text().splitlines()
            if line.split("#", 1)[0].strip()]


def _axis(path: Path) -> tuple[int, np.ndarray]:
    tokens = [v for row in _lines(path) for v in row]
    if len(tokens) < 3:
        raise ValueError(f"incomplete axis: {path}")
    lo, hi = int(tokens[0]), int(tokens[1])
    values = np.asarray([_number(v) for v in tokens[2:]], dtype=float)
    if hi < lo or len(values) != hi - lo + 1 or not np.all(np.isfinite(values)):
        raise ValueError(f"invalid axis bounds/values: {path}")
    return lo, values


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_cold_table(directory: Path | str) -> ColdTable:
    """Preserve absent particles as NaN and never infer physical zero."""
    directory = Path(directory)
    t0, t = _axis(directory / "eos.t")
    y0, y = _axis(directory / "eos.yq")
    n0, n = _axis(directory / "eos.nb")
    if len(t) != 1 or len(y) != 1 or t[0] != 0.0:
        raise ValueError("comparison requires one explicitly cold scalar T/Yq slice")
    if np.any(n <= 0.0) or np.any(np.diff(n) <= 0.0):
        raise ValueError("density must be positive and strictly increasing")
    rows = _lines(directory / "eos.thermo")
    if len(rows[0]) != 3:
        raise ValueError("expected CompOSE mn,mp,Il header")
    mn, mp = map(_number, rows[0][:2])
    il = int(rows[0][2])
    if not np.isfinite(mn) or not np.isfinite(mp) or mn <= 0 or mp <= 0 or il not in (0, 1):
        raise ValueError("invalid reference-mass/lepton header")
    # Cold provider tables sometimes declare the scalar T=0 axis at 0 but
    # consistently label all thermo/composition rows 1 (e.g. RG SLy4).
    # Accept only one consistent declared-or-plus-one scalar label and record
    # the actual label.  This is an index convention, not a changed T slice.
    tlabels = {int(row[0]) for row in rows[1:] if row}
    if len(tlabels) != 1 or not tlabels.issubset({t0, t0 + 1}):
        raise ValueError("inconsistent or unsupported cold scalar T row labels")
    row_t = next(iter(tlabels))
    q = np.full((len(n), 7), np.nan)
    for row in rows[1:]:
        if len(row) < 11:
            raise ValueError("incomplete eos.thermo row")
        it, ib, iy = map(int, row[:3])
        k = ib - n0
        nadd = int(row[10])
        if it != row_t or iy != y0 or not 0 <= k < len(n):
            raise ValueError("thermodynamic row outside its declared axes")
        if nadd < 0 or len(row) != 11 + nadd or np.isfinite(q[k, 0]):
            raise ValueError("invalid additional-column count or duplicate density row")
        q[k] = [_number(v) for v in row[3:10]]
    if not np.all(np.isfinite(q)):
        raise ValueError("missing/nonfinite thermodynamic row")
    if np.max(np.abs(q[:, 1])) > 1e-12:
        raise ValueError("nonzero entropy in requested cold comparison")
    phase = np.full(len(n), -999999, dtype=int)
    fields = {name: np.full(len(n), np.nan) for name in PARTICLES.values()}
    compo = directory / "eos.compo"
    if not compo.exists():
        raise ValueError("phase-aware comparison needs explicit eos.compo phase labels")
    for row in _lines(compo):
        if len(row) < 6:
            raise ValueError("incomplete composition row")
        it, ib, iy, iphase, npart = map(int, row[:5])
        k = ib - n0
        if it != row_t or iy != y0 or not 0 <= k < len(n) or phase[k] != -999999:
            raise ValueError("composition row outside its axes or repeated")
        if npart < 0 or len(row) < 6 + 2 * npart:
            raise ValueError("incomplete composition particles")
        phase[k] = iphase
        seen = set()
        for j in range(npart):
            species = int(row[5 + 2 * j])
            value = _number(row[6 + 2 * j])
            if species in seen or not np.isfinite(value):
                raise ValueError("duplicate/nonfinite composition species")
            seen.add(species)
            if species in PARTICLES:
                fields[PARTICLES[species]][k] = value
        offset = 5 + 2 * npart
        ngroups = int(row[offset])
        if ngroups < 0 or len(row) != offset + 1 + 4 * ngroups:
            raise ValueError("invalid nuclear-group record length")
    if np.any(phase == -999999):
        raise ValueError("missing composition phase row")
    epsilon = n * mn * (1.0 + q[:, 6])
    pressure = n * q[:, 0]
    fields.update(pressure=pressure, epsilon=epsilon, mu_b=mn * (1.0 + q[:, 2]),
                  enthalpy_b=(epsilon + pressure) / n,
                  energy_per_baryon=epsilon / n)
    metadata = {"path": str(directory.resolve()), "rows": len(n),
                "density_fm3": [float(n[0]), float(n[-1])],
                "axis_origins": {"t": t0, "nb": n0, "yq": y0},
                "actual_scalar_row_labels": {"t": row_t, "yq": y0},
                "temperature_MeV": float(t[0]), "charge_fraction_axis": float(y[0]),
                "mn_reference_MeV": mn, "mp_reference_MeV": mp, "Il": il,
                "max_abs_mu_l_MeV": float(np.max(np.abs(mn * q[:, 4]))),
                "max_abs_enthalpy_minus_mu_b_MeV": float(np.max(np.abs(fields["enthalpy_b"] - fields["mu_b"]))),
                "file_sha256": {p.name: _hash(p) for p in
                                [directory / f"eos.{s}" for s in ("nb", "t", "yq", "thermo", "compo")]}}
    return ColdTable(directory, n, phase, fields, metadata)


def continuous_segments(table: ColdTable, phases: set[int]) -> list[np.ndarray]:
    """Keep repeated phase IDs separated if another phase occurs in between."""
    segments = []
    first = 0
    for stop in range(1, len(table.n) + 1):
        if stop == len(table.n) or table.phase[stop] != table.phase[first]:
            if int(table.phase[first]) in phases:
                segments.append(np.arange(first, stop))
            first = stop
    return segments


def _interpolate_segment(n: np.ndarray, values: np.ndarray, targets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate only contiguous finite subsegments; missingness never bridges."""
    pchip, linear = np.full(len(targets), np.nan), np.full(len(targets), np.nan)
    valid = np.isfinite(values)
    starts = np.flatnonzero(valid & ~np.r_[False, valid[:-1]])
    stops = np.flatnonzero(valid & ~np.r_[valid[1:], False]) + 1
    for start, stop in zip(starts, stops):
        if stop - start < 2:
            # An isolated known value is usable only at that exact density.
            hit = targets == n[start]
            pchip[hit] = linear[hit] = values[start]
            continue
        inside = (targets >= n[start]) & (targets <= n[stop - 1])
        pchip[inside] = PchipInterpolator(n[start:stop], values[start:stop], extrapolate=False)(targets[inside])
        linear[inside] = np.interp(targets[inside], n[start:stop], values[start:stop])
    return pchip, linear


def _summary(n: np.ndarray, generated: np.ndarray, reference: np.ndarray, linear: np.ndarray) -> dict:
    valid = np.isfinite(generated) & np.isfinite(reference)
    if not np.any(valid):
        return {"compared_rows": 0, "status": "NO_COMMON_DEFINED_VALUES"}
    n, g, r, lin = n[valid], generated[valid], reference[valid], linear[valid]
    diff = g - r
    scale = max(float(np.max(np.abs(r))), np.finfo(float).tiny)
    floor = scale * 1e-12
    relative = np.abs(r) > floor
    sensitivity = np.abs(g - lin)
    worst = int(np.argmax(np.abs(diff)))
    result = {"compared_rows": len(n), "density_fm3": [float(n[0]), float(n[-1])],
              "max_abs_difference": float(np.max(np.abs(diff))),
              "max_abs_difference_at_n_fm3": float(n[worst]),
              "signed_difference_range": [float(np.min(diff)), float(np.max(diff))],
              "rms_abs_difference": float(np.sqrt(np.mean(diff ** 2))),
              "relative_denominator_floor": floor,
              "relative_defined_rows": int(np.sum(relative)),
              "interpolation_control_max_abs": float(np.max(sensitivity)),
              "status": "MEASURED_DIFFERENCE_NOT_CAUSAL_ATTRIBUTION"}
    if np.any(relative):
        rel = diff[relative] / r[relative]
        worst_rel = int(np.argmax(np.abs(rel)))
        result.update(max_abs_relative_difference=float(np.max(np.abs(rel))),
                      max_abs_relative_difference_at_n_fm3=float(n[relative][worst_rel]),
                      median_signed_relative_difference=float(np.median(rel)),
                      interpolation_control_max_relative=float(np.max(sensitivity[relative] / np.abs(r[relative]))))
    return result


def compare_tables(generated: ColdTable, reference: ColdTable, *, generated_phases: set[int],
                   reference_phases: set[int], density_min: float, density_max: float,
                   composition_fields: tuple[str, ...] = ()) -> tuple[dict, list[dict]]:
    if generated.metadata["Il"] != reference.metadata["Il"]:
        raise ValueError("different Il values require an explicit physical-content comparison")
    if generated.metadata["charge_fraction_axis"] != reference.metadata["charge_fraction_axis"]:
        raise ValueError("different physical charge-fraction slices")
    fields = (*THERMO_FIELDS, "energy_per_baryon", *composition_fields)
    records = []
    overlaps = []
    for gs in continuous_segments(generated, generated_phases):
        if len(gs) < 2:
            continue
        for rs in continuous_segments(reference, reference_phases):
            lo = max(density_min, float(generated.n[gs[0]]), float(reference.n[rs[0]]))
            hi = min(density_max, float(generated.n[gs[-1]]), float(reference.n[rs[-1]]))
            selected = rs[(reference.n[rs] >= lo) & (reference.n[rs] <= hi)]
            if len(selected) == 0:
                continue
            n = reference.n[selected]
            overlaps.append({"generated_phase": int(generated.phase[gs[0]]),
                             "reference_phase": int(reference.phase[rs[0]]),
                             "support_fm3": [lo, hi], "reference_nodes": len(selected)})
            local = [{"n_b": float(x), "generated_phase": int(generated.phase[gs[0]]),
                      "reference_phase": int(reference.phase[rs[0]])} for x in n]
            for field in fields:
                g, lin = _interpolate_segment(generated.n[gs], generated.fields[field][gs], n)
                r = reference.fields[field][selected]
                for row, gv, rv, lv in zip(local, g, r, lin):
                    row.update({f"generated_{field}": float(gv), f"reference_{field}": float(rv),
                                f"difference_{field}": float(gv - rv), f"linear_{field}": float(lv)})
            # A diagnostic derivative of each retained E/n curve. Exclude two
            # knots at each phase-segment edge; never differentiate across a
            # transition or replace the supplied pressure with this value.
            for name, table, segment in (("generated", generated, gs), ("reference", reference, rs)):
                derivative = np.full(len(n), np.nan)
                control = np.full(len(n), np.nan)
                if len(segment) >= 5:
                    nn = table.n[segment]
                    e = table.fields["energy_per_baryon"][segment]
                    interior = (n >= nn[2]) & (n <= nn[-3])
                    derivative[interior] = n[interior] ** 2 * CubicSpline(nn, e)(n[interior], 1)
                    control[interior] = n[interior] ** 2 * PchipInterpolator(nn, e)(n[interior], 1)
                for row, value, sensitivity in zip(local, derivative, control):
                    row[f"{name}_pressure_from_energy_cubic"] = float(value)
                    row[f"{name}_pressure_from_energy_pchip"] = float(sensitivity)
            records.extend(local)
    if not records:
        raise ValueError("no common nodes within selected contiguous phase supports")
    records.sort(key=lambda r: r["n_b"])
    n = np.asarray([r["n_b"] for r in records])
    metrics = {field: _summary(n, *[np.asarray([r[f"{kind}_{field}"] for r in records])
                                   for kind in ("generated", "reference", "linear")]) for field in fields}
    first_law = {}
    for kind in ("generated", "reference"):
        arrays = [np.asarray([r[f"{kind}_{field}"] for r in records]) for field in
                  ("pressure_from_energy_cubic", "pressure", "pressure_from_energy_pchip")]
        first_law[kind] = _summary(n, *arrays)
    # deltaP_measured = (P_g - P_r). If P_r differs from its E/n derivative,
    # compare that data-level contribution separately. No attribution to an
    # inaccessible author kernel and no additive energy fitting is made.
    pg = np.asarray([r["generated_pressure"] for r in records])
    pr = np.asarray([r["reference_pressure"] for r in records])
    pe = np.asarray([r["reference_pressure_from_energy_cubic"] for r in records])
    control = np.asarray([r["reference_pressure_from_energy_pchip"] for r in records])
    first_law["pressure_difference_after_reference_energy_derivative_diagnostic"] = _summary(
        n, pg - pe + pr, pr, pg - control + pr)
    first_law["definition"] = "P_energy=n_B^2*d(epsilon/n_B)/dn_B; CubicSpline and PCHIP, separately within each phase"
    first_law["edge_policy"] = "omit first/last two knots of each full phase segment from derivative comparisons"
    first_law["residual_definition"] = "(P_generated-P_reference)-(P_reference_energy_derivative-P_reference) = P_generated-P_reference_energy_derivative"
    first_law["interpretation"] = "finite-grid diagnostic; interpolation-control spread is not a rigorous uncertainty bound; no source-kernel causal attribution"
    result = {"generated": generated.metadata, "reference": reference.metadata,
              "selected_generated_phases": sorted(generated_phases),
              "selected_reference_phases": sorted(reference_phases),
              "requested_density_fm3": [density_min, density_max], "overlaps": overlaps,
              "interpolation": "generated PCHIP in n_B at exact reference nodes; linear control; no extrapolation",
              "interpolation_control_is_error_bound": False,
              "mu_b_definition": "mn_reference*(1+Q3); compared separately from (epsilon+P)/n_B",
              "energy_convention": "epsilon includes rest energy, recovered using each file's mn header",
              "composition": {"compared": list(composition_fields),
                              "missing_values": "unavailable; never filled with zero",
                              "nuclear_group_tuples": "excluded: public/native semantics not assumed identical"},
              "metrics": metrics, "first_law_diagnostic": first_law,
              "status": "UNRESOLVED_SOURCE_DIFFERENCE"}
    return result, records


def model_spec(model: str) -> dict:
    if model in ("GM_GM1", "RG_SLY4"):
        return {"generated": ROOT / "data/research_matching_20261117/integration/real_pairs_v2/core" / model / "generated",
                "generated_phases": {6}, "reference_phases": {3}, "density_min": .08, "density_max": .5,
                "composition_fields": ("y_p_free",) if model == "GM_GM1" else ("y_p_free", "y_e", "y_mu")}
    if model in {f"SPG_M{i}" for i in range(1, 6)}:
        return {"generated": ROOT / "data/research_matching_20261117/spg/native_129_repaired" / model / "generated",
                "generated_phases": {1}, "reference_phases": {7}, "density_min": .002, "density_max": .2,
                "composition_fields": ("y_e", "y_mu")}
    raise ValueError(f"unsupported campaign model: {model}")


def _write_plot(output: Path, model: str, records: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = np.asarray([r["n_b"] for r in records])
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), constrained_layout=True)
    labels = {"pressure": r"$100(P_g-P_r)/P_r$ [%]",
              "epsilon": r"$100(\epsilon_g-\epsilon_r)/\epsilon_r$ [%]",
              "mu_b": r"$\mu_{B,g}-\mu_{B,r}$ [MeV]",
              "energy_per_baryon": r"$(\epsilon_g-\epsilon_r)/n_B$ [MeV]"}
    for ax, field in zip(axes.ravel(), ("pressure", "epsilon", "mu_b", "energy_per_baryon")):
        g = np.asarray([r[f"generated_{field}"] for r in records])
        r = np.asarray([r[f"reference_{field}"] for r in records])
        if field in ("pressure", "epsilon"):
            ax.plot(n, 100 * (g - r) / np.where(r == 0, np.nan, r))
        else:
            ax.plot(n, g - r)
        ax.set_ylabel(labels[field])
        ax.axhline(0, color="0.5", lw=.7)
        ax.set_xlabel(r"$n_B$ [fm$^{-3}$]")
        ax.grid(alpha=.25)
    fig.suptitle(f"{model}: generated (g) minus public reference (r)")
    fig.savefig(output / f"{model}.png", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, default=ROOT / "eos_generator/eos_data")
    parser.add_argument("--reference-map", type=Path, help="JSON mapping model names to raw CompOSE directories")
    parser.add_argument("--retained-reference-root", type=Path,
                        help="also compare an older retained raw reference against the selected current reference")
    parser.add_argument("--models", nargs="+", default=["GM_GM1", "RG_SLY4", *[f"SPG_M{i}" for i in range(1, 6)]])
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args(argv)
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("output directory must be new or empty")
    args.output.mkdir(parents=True, exist_ok=True)
    reference_map = json.loads(args.reference_map.read_text()) if args.reference_map else {}
    report = {"schema": "generated-public-eos-comparison-v1", "status": "RUNNING",
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "script_sha256": _hash(Path(__file__)), "units": UNITS, "models": {},
              "scope": "finite tabulated same-phase support; neither source equivalence nor thermodynamic admission"}
    summary_path = args.output / "COMPARISON.json"
    summary_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for model in args.models:
        spec = model_spec(model)
        generated = read_cold_table(spec.pop("generated"))
        reference = read_cold_table(reference_map.get(model, args.reference_root / model))
        result, records = compare_tables(generated, reference, **spec)
        result["model"] = model
        report["models"][model] = result
        with (args.output / f"{model}.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            # Blank is unavailable composition; a numeric zero stays a zero.
            writer.writerows({k: (v if not isinstance(v, float) or np.isfinite(v) else "")
                              for k, v in row.items()} for row in records)
        if args.plots:
            _write_plot(args.output, model, records)
        if args.retained_reference_root:
            retained = read_cold_table(args.retained_reference_root / model)
            historical, historical_rows = compare_tables(
                retained, reference, generated_phases=spec["reference_phases"],
                reference_phases=spec["reference_phases"], density_min=spec["density_min"],
                density_max=spec["density_max"], composition_fields=spec["composition_fields"])
            historical["meaning"] = "generated labels denote older retained reference in this subsection only"
            historical["byte_identity"] = {
                name: old == reference.metadata["file_sha256"][name]
                for name, old in retained.metadata["file_sha256"].items()}
            result["retained_vs_current_reference"] = historical
            with (args.output / f"{model}_reference_history.csv").open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(historical_rows[0]))
                writer.writeheader()
                writer.writerows({k: (v if not isinstance(v, float) or np.isfinite(v) else "")
                                  for k, v in row.items()} for row in historical_rows)
        summary_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(f"{model}: {len(records)} common reference nodes; max |deltaP/P|="
              f"{result['metrics']['pressure'].get('max_abs_relative_difference', float('nan')):.8g}", flush=True)
    report["status"] = "COMPLETE_MEASURED_COMPARISON"
    report["scientific_source_equivalence"] = "UNRESOLVED; see model-specific difference ledger"
    summary_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
