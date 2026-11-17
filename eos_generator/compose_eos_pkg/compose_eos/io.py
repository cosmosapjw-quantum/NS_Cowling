"""Strict, loss-aware reader/writer for one-dimensional cold CompOSE tables.

This module deliberately supports the subset used by this project: one cold
(T=0) density slice, at one fixed CompOSE Yq-grid coordinate, with at most one
nuclear-group quadruple and the particle species n, p, e-, mu-.  CompOSE itself
supports substantially more general tables.  Unsupported dimensions/species
are rejected instead of silently flattened or discarded.

Important composition convention:
  * CompOSE eos.compo stores Y_I = n_I / n_b for a nuclear group.
  * EOSData.Yh is an internal *baryon fraction* of that group/cell.
    Therefore Yh = Abar * Y_I on read and Y_I = Yh / Abar on write.

CompOSE format authority: CompOSE Reference Manual v3.01 (2022), Sec. 4.2,
and CompOSE Quick Guide for Providers (2023).
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
import warnings
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from .eos_data import EOSData

_SUPPORTED_PARTICLE_IDS = {0, 1, 10, 11}  # e-, mu-, n, p
_SUPPORTED_NUCLEAR_GROUP_ID = 1  # all 16 delivered composition tables use this ID


def _as_int(token: str, context: str) -> int:
    """Parse an integer-valued CompOSE field without silently truncating it."""
    value = float(token)
    if not np.isfinite(value) or value != np.rint(value):
        raise ValueError(f"{context}: expected an integer-valued field, got {token!r}")
    return int(value)


def _nonempty_lines(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]


def _load_grid(path: str) -> Tuple[int, int, np.ndarray]:
    """Read a CompOSE grid file: first/last grid indices, then grid values."""
    lines = _nonempty_lines(path)
    if len(lines) < 3:
        raise ValueError(f"{path}: expected index_min, index_max and at least one value")
    try:
        i0 = _as_int(lines[0].split()[0], f"{path}: first grid index")
        i1 = _as_int(lines[1].split()[0], f"{path}: last grid index")
        values = np.asarray([float(x.split()[0]) for x in lines[2:]], dtype=float)
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: malformed CompOSE grid file") from exc
    if i1 < i0:
        raise ValueError(f"{path}: last index {i1} precedes first index {i0}")
    expected = i1 - i0 + 1
    if values.size != expected:
        raise ValueError(
            f"{path}: index range [{i0},{i1}] requires {expected} values, got {values.size}"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{path}: non-finite grid coordinate")
    return i0, i1, values


def _load_thermo(path: str):
    """Parse eos.thermo and apply CompOSE's documented last-duplicate-wins rule."""
    lines = _nonempty_lines(path)
    if not lines:
        raise ValueError(f"{path}: empty thermo table")
    h = lines[0].split()
    if len(h) != 3:
        raise ValueError(f"{path}: first line must contain mn, mp, Il")
    mn, mp = float(h[0]), float(h[1])
    Il = _as_int(h[2], f"{path}: lepton flag Il")
    if not (mn > 0 and mp > 0 and np.isfinite(mn + mp)):
        raise ValueError(f"{path}: invalid nucleon masses")
    if Il not in (0, 1):
        raise ValueError(f"{path}: unsupported CompOSE lepton flag Il={Il}; expected 0 or 1")
    rows: Dict[Tuple[int, int, int], tuple] = {}
    duplicates = 0
    nadd_seen = set()
    for lineno, ln in enumerate(lines[1:], 2):
        f = ln.split()
        if len(f) < 11:
            raise ValueError(f"{path}:{lineno}: short thermo row ({len(f)} columns)")
        try:
            iT, inb, iYq = (_as_int(f[j], f"{path}:{lineno}: index") for j in range(3))
            q = np.asarray([float(x) for x in f[3:10]], dtype=float)
            nadd = _as_int(f[10], f"{path}:{lineno}: Nadd")
            extra = np.asarray([float(x) for x in f[11:]], dtype=float)
        except ValueError as exc:
            raise ValueError(f"{path}:{lineno}: malformed numeric thermo row: {exc}") from exc
        if nadd < 0 or len(extra) != nadd:
            raise ValueError(
                f"{path}:{lineno}: Nadd={nadd}, but row carries {len(extra)} extra values"
            )
        if not np.all(np.isfinite(q)) or not np.all(np.isfinite(extra)):
            raise ValueError(f"{path}:{lineno}: non-finite thermodynamic quantity")
        key = (iT, inb, iYq)
        if key in rows:
            duplicates += 1
        rows[key] = (q, nadd, extra)  # CompOSE: the last identical index triple wins.
        nadd_seen.add(nadd)
    if duplicates:
        warnings.warn(
            f"{path}: {duplicates} duplicate index triple(s); applying CompOSE last-row-wins rule",
            RuntimeWarning,
        )
    if len(nadd_seen) != 1:
        raise ValueError(f"{path}: Nadd changes within one EOS table: {sorted(nadd_seen)}")
    nadd = next(iter(nadd_seen))
    if nadd > 1:
        raise NotImplementedError(
            f"{path}: Nadd={nadd}; this EOSData version preserves at most one optional column"
        )
    return mn, mp, Il, rows, nadd


def _load_compo(path: str):
    """Parse eos.compo keyed by the full (iT,inb,iYq) triple."""
    rows: Dict[Tuple[int, int, int], dict] = {}
    duplicates = 0
    for lineno, ln in enumerate(_nonempty_lines(path), 1):
        f = ln.split()
        if len(f) < 6:
            raise ValueError(f"{path}:{lineno}: short composition row")
        try:
            iT, inb, iYq, iphase, npairs = (
                _as_int(f[j], f"{path}:{lineno}: index/count") for j in range(5)
            )
            if npairs < 0:
                raise ValueError("negative Npairs")
            k = 5
            pairs = {}
            for _ in range(npairs):
                pid = _as_int(f[k], f"{path}:{lineno}: particle id")
                y = float(f[k + 1]); k += 2
                if pid in pairs:
                    raise ValueError(f"duplicate particle id {pid}")
                if pid not in _SUPPORTED_PARTICLE_IDS:
                    raise NotImplementedError(
                        f"particle id {pid} is not supported by this four-species EOSData model"
                    )
                if not np.isfinite(y):
                    raise ValueError("non-finite particle fraction")
                pairs[pid] = y
            if k >= len(f):
                raise ValueError("missing Nquad field")
            nquad = _as_int(f[k], f"{path}:{lineno}: Nquad"); k += 1
            if nquad < 0:
                raise ValueError("negative Nquad")
            if nquad > 1:
                raise NotImplementedError(
                    f"Nquad={nquad}; current EOSData stores one averaged nuclear group only"
                )
            quads = []
            for _ in range(nquad):
                qi = _as_int(f[k], f"{path}:{lineno}: nuclear group id")
                A = float(f[k + 1]); Z = float(f[k + 2]); YI = float(f[k + 3]); k += 4
                if qi != _SUPPORTED_NUCLEAR_GROUP_ID:
                    raise NotImplementedError(
                        f"nuclear group id {qi} is not supported; expected {_SUPPORTED_NUCLEAR_GROUP_ID}"
                    )
                if not np.all(np.isfinite([A, Z, YI])):
                    raise ValueError("non-finite nuclear-group value")
                if A == 0.0 and Z == 0.0 and YI == 0.0:
                    continue
                if A <= 0 or Z < 0 or Z > A or YI < 0:
                    raise ValueError(f"invalid nuclear group (A={A}, Z={Z}, YI={YI})")
                quads.append((qi, A, Z, YI))
            if k != len(f):
                raise ValueError(f"{len(f)-k} trailing field(s) after parsed composition row")
        except IndexError as exc:
            raise ValueError(f"{path}:{lineno}: truncated composition row") from exc
        except ValueError as exc:
            raise ValueError(f"{path}:{lineno}: malformed composition row: {exc}") from exc
        key = (iT, inb, iYq)
        if key in rows:
            duplicates += 1
        rows[key] = dict(Iphase=iphase, pairs=pairs, quads=quads)
    if duplicates:
        warnings.warn(
            f"{path}: {duplicates} duplicate composition index triple(s); applying last-row-wins rule",
            RuntimeWarning,
        )
    return rows


def _one_dimensional_slice(nb_grid, t_grid, yq_grid, thermo_rows):
    nb_i0, nb_i1, nb = nb_grid
    t_i0, t_i1, T = t_grid
    y_i0, y_i1, Yq = yq_grid
    keys = list(thermo_rows)
    t_indices = sorted({k[0] for k in keys})
    y_indices = sorted({k[2] for k in keys})
    if len(t_indices) != 1 or len(y_indices) != 1:
        raise NotImplementedError(
            "This package handles a 1D cold slice only; eos.thermo contains "
            f"{len(t_indices)} temperature and {len(y_indices)} charge-fraction indices. "
            "Select a 1D slice explicitly before reconstruction."
        )
    iT, iYq = t_indices[0], y_indices[0]
    if not (t_i0 <= iT <= t_i1 and y_i0 <= iYq <= y_i1):
        raise ValueError("thermo slice index lies outside eos.t/eos.yq grid index ranges")
    if len(T) != 1 or len(Yq) != 1:
        raise NotImplementedError(
            "The selected table has a multi-point eos.t or eos.yq grid; refusing to flatten it"
        )
    if abs(float(T[iT - t_i0])) > 1e-12:
        raise NotImplementedError(
            f"selected temperature is {T[iT-t_i0]:g} MeV, not a cold T=0 slice; "
            "compose_eos v3.3.1 supports only cold T=0 tables"
        )
    expected = {(iT, ib, iYq) for ib in range(nb_i0, nb_i1 + 1)}
    if set(keys) != expected:
        missing = sorted(expected - set(keys))[:8]
        extra = sorted(set(keys) - expected)[:8]
        raise ValueError(
            f"thermo rows do not exactly cover the 1D density grid; missing={missing}, extra={extra}"
        )
    ordered = [thermo_rows[(iT, ib, iYq)] for ib in range(nb_i0, nb_i1 + 1)]
    return iT, iYq, float(T[iT - t_i0]), float(Yq[iYq - y_i0]), nb, ordered


def read_compose(directory: str) -> EOSData:
    """Read a cold, one-dimensional CompOSE table without lossy flattening."""
    nb_path = os.path.join(directory, "eos.nb")
    t_path = os.path.join(directory, "eos.t")
    yq_path = os.path.join(directory, "eos.yq")
    th_path = os.path.join(directory, "eos.thermo")
    for p in (nb_path, t_path, yq_path, th_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(f"missing required CompOSE file: {p}")
    nb_grid = _load_grid(nb_path)
    t_grid = _load_grid(t_path)
    yq_grid = _load_grid(yq_path)
    mn, mp, Il, thermo_rows, nadd = _load_thermo(th_path)
    iT, iYq, T, Yq, nb, ordered = _one_dimensional_slice(nb_grid, t_grid, yq_grid, thermo_rows)
    if nb.size < 3 or np.any(nb <= 0) or np.any(np.diff(nb) <= 0):
        raise ValueError("eos.nb must contain at least three positive, strictly increasing densities")
    Q = np.vstack([r[0] for r in ordered])
    if nadd:
        q_extra = np.asarray([r[2][0] for r in ordered], dtype=float)
    else:
        q_extra = np.zeros(len(nb), dtype=float)

    cp_path = os.path.join(directory, "eos.compo")
    composition_file_present = os.path.isfile(cp_path)
    if composition_file_present:
        compo_rows = _load_compo(cp_path)
        expected_keys = set(thermo_rows)
        if set(compo_rows) != expected_keys:
            missing = sorted(expected_keys - set(compo_rows))[:8]
            extra = sorted(set(compo_rows) - expected_keys)[:8]
            raise ValueError(
                f"eos.compo index triples do not match eos.thermo; missing={missing}, extra={extra}"
            )
        crows = [compo_rows[(iT, ib, iYq)] for ib in range(nb_grid[0], nb_grid[1] + 1)]
        # A syntactically present eos.compo can legally contain zero particle
        # pairs and zero nuclear groups at every grid point. Preserve phase IDs,
        # but do not advertise physical composition as available in that case.
        composition_available = any(bool(c["pairs"] or c["quads"]) for c in crows)
        Iphase = np.asarray([c["Iphase"] for c in crows], dtype=int)
        Yn = np.asarray([c["pairs"].get(10, 0.0) for c in crows])
        Yp = np.asarray([c["pairs"].get(11, 0.0) for c in crows])
        Ye = np.asarray([c["pairs"].get(0, 0.0) for c in crows])
        Ymu = np.asarray([c["pairs"].get(1, 0.0) for c in crows])
        Aav = np.asarray([c["quads"][0][1] if c["quads"] else 0.0 for c in crows])
        Zav = np.asarray([c["quads"][0][2] if c["quads"] else 0.0 for c in crows])
        YI = np.asarray([c["quads"][0][3] if c["quads"] else 0.0 for c in crows])
        Yh = Aav * YI
        frac_arrays = (Yn, Yp, Ye, Ymu, Yh)
        if any(np.any(a < -1e-10) for a in frac_arrays):
            raise ValueError("eos.compo contains a negative particle/nuclear fraction")
        if np.any(Yh > 1.0 + 5e-3):
            raise ValueError("CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics")
    else:
        composition_available = False
        Iphase = np.zeros(len(nb), dtype=int)
        Yn = Yp = Ye = Ymu = Aav = Zav = Yh = np.zeros(len(nb), dtype=float)

    if not np.isfinite(T) or abs(float(T))>1e-12:
        raise NotImplementedError(
            f'{directory}: compose_eos v3.3.1 supports only cold T=0 tables; got T={T:g} MeV'
        )
    comp_charge_neutral = True
    if composition_available:
        comp_charge_neutral = float(np.max(np.abs(Yp + Zav * YI - Ye - Ymu))) <= 2e-4
    is_cold = np.isfinite(T) and abs(float(T)) <= 1e-12
    is_beta_eq = bool(
        is_cold
        and Il == 1
        and float(np.max(np.abs(Q[:, 4]))) < 1e-8
        and comp_charge_neutral
    )
    return EOSData(
        nb=nb, mn=mn, mp=mp,
        Q1=Q[:, 0], Q2=Q[:, 1], Q3=Q[:, 2], Q4=Q[:, 3],
        Q5=Q[:, 4], Q6=Q[:, 5], Q7=Q[:, 6],
        q_extra=q_extra, Iphase=Iphase,
        Yn=Yn, Yp=Yp, Ye=Ye, Ymu=Ymu,
        Aav=Aav, Zav=Zav, Yh=Yh,
        Il=Il, iT=iT, iYq=iYq, inb_start=int(nb_grid[0]), temperature_mev=T,
        charge_fraction_grid=Yq, q_extra_count=nadd,
        composition_available=composition_available,
        beta_equilibrated=is_beta_eq,
    )


def _extra_count(data: EOSData, copy_from: Optional[str]) -> int:
    if data.q_extra_count is not None:
        nadd = int(data.q_extra_count)
    elif copy_from and os.path.isfile(os.path.join(copy_from, "eos.thermo")):
        _, _, _, _, nadd = _load_thermo(os.path.join(copy_from, "eos.thermo"))
    else:
        nadd = int(data.q_extra.size == data.n and np.any(data.q_extra != 0.0))
    if nadd not in (0, 1):
        raise NotImplementedError("write_compose currently supports Nadd=0 or 1 only")
    if nadd == 1 and data.q_extra.size != data.n:
        raise ValueError("Nadd=1 requires one q_extra value per density point")
    return nadd


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _compose_dir_hashes(directory: str) -> Dict[str, str]:
    hashes: Dict[str, str] = {}
    for fn in ("eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo"):
        p = os.path.join(directory, fn)
        if os.path.isfile(p):
            hashes[fn] = _sha256_file(p)
    return hashes


def _array_diff_stats(actual: np.ndarray, expected: np.ndarray) -> Dict[str, float]:
    a = np.asarray(actual, dtype=float)
    b = np.asarray(expected, dtype=float)
    abs_diff = np.abs(a - b)
    scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-30)
    rel_diff = abs_diff / scale
    return {
        "max_abs": float(np.max(abs_diff)) if abs_diff.size else 0.0,
        "max_rel": float(np.max(rel_diff)) if rel_diff.size else 0.0,
    }


def verify_compose_readback(
    directory: str,
    expected: EOSData,
    *,
    expected_nadd: Optional[int] = None,
    source_dir: Optional[str] = None,
    rtol: float = 1e-10,
    atol: float = 1e-12,
    compo_rtol: float = 1e-9,
    compo_atol: float = 1e-12,
    failure_receipt_path: Optional[str] = None,
) -> Dict[str, object]:
    """Read back a written CompOSE directory natively and verify full fidelity.

    Verifies numerical arrays, Q1..Q7, units and mass headers, Iphase/Nadd,
    composition arrays and presence masks, C3/C4 applicability and status, and
    source/output SHA256 checksums. Preserves a failure receipt before raising
    if any check fails.
    """
    from .consistency import check_consistency

    errors: List[str] = []
    source_hashes = _compose_dir_hashes(source_dir) if source_dir and os.path.isdir(source_dir) else {}
    output_hashes = _compose_dir_hashes(directory)
    receipt: Dict[str, object] = {
        "schema": "compose-eos-readback-receipt-v1",
        "directory": os.path.abspath(directory),
        "source_directory": os.path.abspath(source_dir) if source_dir else None,
        "source_sha256": source_hashes,
        "output_sha256": output_hashes,
        "units": {
            "nb": "fm^-3",
            "temperature_mev": "MeV",
            "mn": "MeV",
            "mp": "MeV",
            "Q1_to_Q7": "dimensionless (CompOSE Manual v3.01 Sec. 4.2)",
        },
        "status": "FAIL",
        "errors": errors,
    }
    try:
        actual = read_compose(directory)
    except Exception as exc:
        errors.append(f"native read_compose failed: {exc}")
        if failure_receipt_path:
            os.makedirs(os.path.dirname(os.path.abspath(failure_receipt_path)), exist_ok=True)
            with open(failure_receipt_path, "w", encoding="utf-8") as fh:
                json.dump(receipt, fh, indent=2, sort_keys=True)
        raise ValueError(f"CompOSE readback failed: {exc}") from exc

    n_exp = expected.n
    n_act = actual.n
    receipt["point_count"] = {"expected": int(n_exp), "actual": int(n_act)}
    if n_act != n_exp:
        errors.append(f"point count mismatch: expected {n_exp}, got {n_act}")

    # Mass & grid headers
    exp_nadd = int(expected_nadd) if expected_nadd is not None else (
        int(expected.q_extra_count)
        if expected.q_extra_count is not None
        else int(expected.q_extra.size == expected.n and np.any(expected.q_extra != 0.0))
    )
    act_nadd = int(actual.q_extra_count if actual.q_extra_count is not None else 0)
    receipt["headers"] = {
        "mn": {"expected": float(expected.mn), "actual": float(actual.mn)},
        "mp": {"expected": float(expected.mp), "actual": float(actual.mp)},
        "Il": {"expected": int(expected.Il), "actual": int(actual.Il)},
        "iT": {"expected": int(expected.iT), "actual": int(actual.iT)},
        "iYq": {"expected": int(expected.iYq), "actual": int(actual.iYq)},
        "inb_start": {"expected": int(expected.inb_start), "actual": int(actual.inb_start)},
        "temperature_mev": {"expected": float(expected.temperature_mev), "actual": float(actual.temperature_mev)},
        "charge_fraction_grid": {
            "expected": float(expected.charge_fraction_grid),
            "actual": float(actual.charge_fraction_grid),
        },
        "Nadd": {"expected": exp_nadd, "actual": act_nadd},
        "beta_equilibrated": {
            "expected_flag": bool(expected.beta_equilibrated),
            "actual": bool(actual.beta_equilibrated),
        },
    }
    for name in ("mn", "mp", "temperature_mev", "charge_fraction_grid"):
        ev = float(getattr(expected, name))
        av = float(getattr(actual, name))
        if not np.isclose(av, ev, rtol=rtol, atol=atol):
            errors.append(f"header {name} mismatch: expected {ev!r}, got {av!r}")
    for name in ("Il", "iT", "iYq", "inb_start"):
        if int(getattr(actual, name)) != int(getattr(expected, name)):
            errors.append(
                f"header {name} mismatch: expected {int(getattr(expected, name))}, "
                f"got {int(getattr(actual, name))}"
            )
    if act_nadd != exp_nadd:
        errors.append(f"Nadd mismatch: expected {exp_nadd}, got {act_nadd}")

    # Q1..Q7 & derived numerical arrays
    q_stats: Dict[str, Dict[str, float]] = {}
    num_stats: Dict[str, Dict[str, float]] = {}
    if n_act == n_exp:
        for name in ("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"):
            ea = np.asarray(getattr(expected, name), dtype=float)
            aa = np.asarray(getattr(actual, name), dtype=float)
            q_stats[name] = _array_diff_stats(aa, ea)
            if not np.allclose(aa, ea, rtol=rtol, atol=atol):
                errors.append(
                    f"{name} readback mismatch (max_abs={q_stats[name]['max_abs']:.3e}, "
                    f"max_rel={q_stats[name]['max_rel']:.3e})"
                )
        num_fields = ["nb", "p", "eps", "mub", "muq", "mul", "F"]
        if exp_nadd == 1:
            num_fields.append("q_extra")
        for name in num_fields:
            ea = np.asarray(getattr(expected, name), dtype=float)
            aa = np.asarray(getattr(actual, name), dtype=float)
            num_stats[name] = _array_diff_stats(aa, ea)
            if not np.allclose(aa, ea, rtol=rtol, atol=atol):
                errors.append(
                    f"array {name} readback mismatch (max_abs={num_stats[name]['max_abs']:.3e}, "
                    f"max_rel={num_stats[name]['max_rel']:.3e})"
                )
    receipt["Q1_Q7_residuals"] = q_stats
    receipt["numerical_array_residuals"] = num_stats

    # Composition arrays, presence mask, and Iphase
    compo_file_exists = os.path.isfile(os.path.join(directory, "eos.compo"))
    exp_comp = bool(expected.composition_available)
    act_comp = bool(actual.composition_available)
    compo_stats: Dict[str, Dict[str, float]] = {}
    receipt["composition"] = {
        "expected_available": exp_comp,
        "actual_available": act_comp,
        "eos_compo_present": compo_file_exists,
        "residuals": compo_stats,
    }
    if act_comp != exp_comp:
        errors.append(f"composition_available mismatch: expected {exp_comp}, got {act_comp}")
    if exp_comp != compo_file_exists:
        errors.append(
            f"eos.compo presence mismatch: composition_available={exp_comp}, file_exists={compo_file_exists}"
        )
    if n_act == n_exp:
        if exp_comp:
            if not np.array_equal(np.asarray(actual.Iphase, dtype=int), np.asarray(expected.Iphase, dtype=int)):
                errors.append("Iphase mismatch on eos.compo readback")
            for name in ("Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh"):
                ea = np.asarray(getattr(expected, name), dtype=float)
                aa = np.asarray(getattr(actual, name), dtype=float)
                compo_stats[name] = _array_diff_stats(aa, ea)
                if not np.allclose(aa, ea, rtol=compo_rtol, atol=compo_atol):
                    errors.append(
                        f"composition array {name} mismatch (max_abs={compo_stats[name]['max_abs']:.3e}, "
                        f"max_rel={compo_stats[name]['max_rel']:.3e})"
                    )
        else:
            if not np.all(np.asarray(actual.Iphase, dtype=int) == 0):
                errors.append("expected zero Iphase on composition-absent readback")
    receipt["Iphase"] = {
        "expected_unique": [int(x) for x in np.unique(expected.Iphase)] if np.asarray(expected.Iphase).size else [],
        "actual_unique": [int(x) for x in np.unique(actual.Iphase)] if np.asarray(actual.Iphase).size else [],
    }

    # C3/C4 applicability and agreement
    if n_act == n_exp and n_act >= 3:
        exp_rep = check_consistency(expected, verbose=False)
        act_rep = check_consistency(actual, verbose=False)
        c34 = {}
        for key in ("C3", "C4"):
            ec, ac = exp_rep[key], act_rep[key]
            exp_app = bool(ec["status"] == "CHECKED")
            act_app = bool(ac["status"] == "CHECKED")
            c34[key] = {
                "expected_status": ec["status"],
                "actual_status": ac["status"],
                "expected_applicable": exp_app,
                "actual_applicable": act_app,
                "expected_pass": bool(ec.get("pass_", False)) if exp_app else None,
                "actual_pass": bool(ac.get("pass_", False)) if act_app else None,
            }
            if ec["status"] != ac["status"] or exp_app != act_app:
                errors.append(
                    f"{key} applicability mismatch: expected ({ec['status']}, {exp_app}), "
                    f"got ({ac['status']}, {act_app})"
                )
            elif exp_app and bool(ec["pass_"]) != bool(ac["pass_"]):
                errors.append(
                    f"{key} pass status mismatch: expected {ec['pass_']}, got {ac['pass_']}"
                )
        receipt["C3_C4_applicability"] = c34

    if errors:
        receipt["status"] = "FAIL"
        if failure_receipt_path:
            os.makedirs(os.path.dirname(os.path.abspath(failure_receipt_path)), exist_ok=True)
            with open(failure_receipt_path, "w", encoding="utf-8") as fh:
                json.dump(receipt, fh, indent=2, sort_keys=True)
        raise ValueError("CompOSE staged readback verification failed: " + "; ".join(errors))

    receipt["status"] = "PASS"
    return receipt


def _replace_file(src: str, dst: str) -> None:
    """Replace dst with src atomically; separated so unit tests can inject mid-commit I/O faults."""
    os.replace(src, dst)


def write_compose(directory: str, data: EOSData,
                  copy_from: Optional[str] = None,
                  *, composition_tolerance: float = 2e-3,
                  receipt_path: Optional[str] = None,
                  _post_stage_hook: Optional[Callable[[str], None]] = None) -> Dict[str, object]:
    """Write a standards-compliant 1D CompOSE table with transactional staging and rollback.

    Internal `Yh` is converted to the CompOSE group number abundance `Y_I=Yh/A`.
    If composition was absent on input, it remains absent. ``copy_from`` is used
    only to infer legacy optional-column metadata; model-specific eos.mr/eos.pdf/
    README/eos.init sidecars are never copied to a modified EOS.
    All files are staged in a temporary directory and committed with automatic
    rollback if any file replacement fails mid-commit.
    """
    n = data.n
    nb = np.asarray(data.nb, dtype=float)
    if n < 3 or nb.shape != (n,) or np.any(~np.isfinite(nb)) or np.any(nb <= 0) or np.any(np.diff(nb) <= 0):
        raise ValueError("EOSData.nb must be finite, positive, and strictly increasing")
    for name in ("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"):
        arr = np.asarray(getattr(data, name), dtype=float)
        if arr.shape != (n,) or not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} must be a finite vector of length {n}")
    nadd = _extra_count(data, copy_from)
    if int(data.Il) not in (0, 1):
        raise ValueError("CompOSE Il must be 0 (no leptons) or 1 (leptons included)")
    if not np.isfinite(data.temperature_mev) or abs(float(data.temperature_mev)) > 1e-12:
        raise ValueError("this package only writes cold T=0 tables")
    if nadd == 1 and (data.q_extra.shape != (n,) or not np.all(np.isfinite(data.q_extra))):
        raise ValueError("q_extra must be a finite length-N vector when Nadd=1")
    iT, iYq = int(data.iT), int(data.iYq)
    inb_start = int(data.inb_start)

    yI = None
    if data.composition_available:
        required = ("Iphase", "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh")
        for name in required:
            if np.asarray(getattr(data, name)).shape != (n,):
                raise ValueError(f"composition field {name} must have length {n}")
        for name in required[1:]:
            if np.any(~np.isfinite(getattr(data, name))):
                raise ValueError(f"composition field {name} contains non-finite values")
        iphase = np.asarray(data.Iphase, dtype=float)
        if np.any(~np.isfinite(iphase)) or np.any(iphase != np.rint(iphase)):
            raise ValueError("Iphase values must be finite integers")
        for name in ("Yn", "Yp", "Ye", "Ymu", "Yh"):
            arr = np.asarray(getattr(data, name), dtype=float)
            if np.any(arr < -1e-10) or np.any(arr > 1.0 + composition_tolerance):
                raise ValueError(f"composition fraction {name} must lie in [0,1]")
        if np.any(data.Aav < 0) or np.any(data.Zav < 0) or np.any(data.Zav > data.Aav + 1e-9):
            raise ValueError("invalid Abar/Zbar in nuclear group")
        has_group = np.asarray(data.Aav > 0)
        if np.any(data.Aav[has_group] < 1.0) or np.any(data.Zav[has_group] > data.Aav[has_group]):
            raise ValueError("a nuclear group requires Abar >= 1 and 0 <= Zbar <= Abar")
        if np.any((~has_group) & (data.Yh > 1e-10)):
            raise ValueError("Yh is nonzero although Abar indicates no nuclear group")
        yI = np.divide(data.Yh, data.Aav, out=np.zeros(n), where=has_group)
        baryon_res = data.Yn + data.Yp + data.Aav * yI - 1.0
        charge_res = data.Yp + data.Zav * yI - data.Ye - data.Ymu
        max_baryon = float(np.max(np.abs(baryon_res)))
        max_charge = float(np.max(np.abs(charge_res)))
        if max_baryon > composition_tolerance or max_charge > composition_tolerance:
            raise ValueError(
                "internal composition does not satisfy CompOSE number-abundance closure: "
                f"max baryon residual={max_baryon:.3e}, charge residual={max_charge:.3e}"
            )

    os.makedirs(directory, exist_ok=True)
    stage_dir = tempfile.mkdtemp(prefix=".compose_stage_", dir=directory)
    backup_dir = tempfile.mkdtemp(prefix=".compose_backup_", dir=directory)
    committed: List[str] = []
    existed: Dict[str, bool] = {}

    try:
        with open(os.path.join(stage_dir, "eos.nb"), "w", encoding="utf-8") as fh:
            fh.write(f"{inb_start:4d}\n{inb_start + n - 1:4d}\n")
            for value in nb:
                fh.write(f" {value:.16E}\n")
        with open(os.path.join(stage_dir, "eos.t"), "w", encoding="utf-8") as fh:
            fh.write(f"{iT:4d}\n{iT:4d}\n {float(data.temperature_mev):.16E}\n")
        with open(os.path.join(stage_dir, "eos.yq"), "w", encoding="utf-8") as fh:
            fh.write(f"{iYq:4d}\n{iYq:4d}\n {float(data.charge_fraction_grid):.16E}\n")

        with open(os.path.join(stage_dir, "eos.thermo"), "w", encoding="utf-8") as fh:
            fh.write(f" {data.mn:.16E} {data.mp:.16E} {int(data.Il)}\n")
            for j in range(n):
                q = [getattr(data, f"Q{k}")[j] for k in range(1, 8)]
                row = f" {iT:d} {inb_start+j:d} {iYq:d} " + " ".join(f"{x:.16E}" for x in q)
                if nadd == 0:
                    fh.write(row + " 0\n")
                else:
                    fh.write(row + f" 1 {float(data.q_extra[j]):.16E}\n")

        target_names = ["eos.nb", "eos.t", "eos.yq", "eos.thermo"]
        if data.composition_available:
            with open(os.path.join(stage_dir, "eos.compo"), "w", encoding="utf-8") as fh:
                for j in range(n):
                    has_group = data.Aav[j] > 0
                    nquad = 1 if has_group else 0
                    fh.write(
                        f" {iT:d} {inb_start+j:d} {iYq:d} {int(data.Iphase[j])} 4 "
                        f"10 {data.Yn[j]:.12E} 11 {data.Yp[j]:.12E} "
                        f"0 {data.Ye[j]:.12E} 1 {data.Ymu[j]:.12E} {nquad}"
                    )
                    if has_group:
                        fh.write(f" 1 {data.Aav[j]:.12E} {data.Zav[j]:.12E} {yI[j]:.16E}")
                    fh.write("\n")
            target_names.append("eos.compo")

        # Never copy source-model sidecars to a modified EOS. In particular,
        # eos.mr and eos.pdf would falsely describe the original crust/core.

        if _post_stage_hook is not None:
            _post_stage_hook(stage_dir)

        failure_receipt = os.path.join(directory, "eos.readback_failure_receipt.json")
        receipt = verify_compose_readback(
            stage_dir,
            data,
            expected_nadd=nadd,
            source_dir=copy_from,
            failure_receipt_path=failure_receipt,
        )

        # Backup existing files before committing staged replacements
        for fn in target_names + (["eos.compo"] if not data.composition_available else []):
            dst = os.path.join(directory, fn)
            existed[fn] = os.path.isfile(dst)
            if existed[fn]:
                shutil.copy2(dst, os.path.join(backup_dir, fn))

        for fn in target_names:
            _replace_file(os.path.join(stage_dir, fn), os.path.join(directory, fn))
            committed.append(fn)

        if not data.composition_available:
            cp_dst = os.path.join(directory, "eos.compo")
            if os.path.exists(cp_dst):
                os.remove(cp_dst)
                committed.append("eos.compo")

        receipt["directory"] = os.path.abspath(directory)
        receipt["output_sha256"] = _compose_dir_hashes(directory)
        if receipt_path is not None:
            os.makedirs(os.path.dirname(os.path.abspath(receipt_path)), exist_ok=True)
            with open(receipt_path, "w", encoding="utf-8") as fh:
                json.dump(receipt, fh, indent=2, sort_keys=True)
        return receipt
    except Exception:
        # Roll back any files modified before the failure
        for fn in committed:
            dst = os.path.join(directory, fn)
            if existed.get(fn, False):
                shutil.copy2(os.path.join(backup_dir, fn), dst)
            elif os.path.exists(dst):
                os.remove(dst)
        raise
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)
        shutil.rmtree(backup_dir, ignore_errors=True)


def write_compose_verified(
    directory: str,
    data: EOSData,
    copy_from: Optional[str] = None,
    *,
    composition_tolerance: float = 2e-3,
    receipt_path: Optional[str] = None,
) -> Dict[str, object]:
    """Write a CompOSE directory with staged native readback and write a JSON receipt."""
    rpath = receipt_path or os.path.join(directory, "eos.readback_receipt.json")
    return write_compose(
        directory,
        data,
        copy_from=copy_from,
        composition_tolerance=composition_tolerance,
        receipt_path=rpath,
    )


__all__ = [
    "read_compose",
    "write_compose",
    "write_compose_verified",
    "verify_compose_readback",
    "_replace_file",
]
