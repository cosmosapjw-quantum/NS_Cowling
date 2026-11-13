"""Phase-aware Cowling mode diagnostics and profile-overlap continuation.

No automatic f/p1/interface-g labels are assigned. The positive norm is
consistent with the project displacement convention:
N=sum_s int h exp(lambda-nu) [|W|^2/r^2+l(l+1)|V|^2] dr.
The cross-background overlap is a chosen continuation similarity, not a
canonical orthogonality relation or an observable.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import math
import numpy as np
from scipy.optimize import linear_sum_assignment


@dataclass(frozen=True)
class PhaseProfile:
    r_m: np.ndarray
    W_m3: np.ndarray
    V_m2: np.ndarray
    enthalpy_j_m3: np.ndarray
    nu: np.ndarray
    lambda_: np.ndarray

    def validated(self, name: str) -> "PhaseProfile":
        arrays = [np.asarray(x) for x in (
            self.r_m, self.W_m3, self.V_m2, self.enthalpy_j_m3, self.nu, self.lambda_
        )]
        if any(a.ndim != 1 for a in arrays):
            raise ValueError(f"{name}: all phase arrays must be one-dimensional")
        n = arrays[0].size
        if n < 3 or any(a.size != n for a in arrays):
            raise ValueError(f"{name}: all phase arrays must have the same length >= 3")
        if any(not np.all(np.isfinite(a)) for a in arrays):
            raise ValueError(f"{name}: phase arrays must be finite")
        r, _, _, h, _, _ = arrays
        if np.any(r < 0.0) or np.any(np.diff(r) <= 0.0):
            raise ValueError(f"{name}: radius must be nonnegative and strictly increasing")
        if np.any(h <= 0.0):
            raise ValueError(f"{name}: enthalpy density must be positive")
        return self


@dataclass(frozen=True)
class ModeProfile:
    core: PhaseProfile
    envelope: PhaseProfile
    frequency_hz: float
    ell: int = 2
    conversion: str = "slow"
    model_id: str = "unspecified"

    def validated(self) -> "ModeProfile":
        self.core.validated("core")
        self.envelope.validated("envelope")
        if not (math.isfinite(self.frequency_hz) and self.frequency_hz > 0.0):
            raise ValueError("frequency_hz must be finite and positive")
        if self.ell < 1:
            raise ValueError("ell must be >= 1")
        if self.conversion not in ("slow", "rapid"):
            raise ValueError("conversion must be 'slow' or 'rapid'")
        rt_i, rt_o = float(self.core.r_m[-1]), float(self.envelope.r_m[0])
        tol = 1e-10 * max(abs(rt_i), abs(rt_o), 1.0)
        if abs(rt_i - rt_o) > tol:
            raise ValueError("core and envelope must share the same interface radius")
        if self.envelope.r_m[-1] <= rt_i:
            raise ValueError("envelope outer radius must exceed interface radius")
        return self

    @property
    def interface_radius_m(self) -> float:
        return float(self.core.r_m[-1])

    @property
    def surface_radius_m(self) -> float:
        return float(self.envelope.r_m[-1])


def _trapz(y: np.ndarray, x: np.ndarray) -> complex | float:
    fn = getattr(np, "trapezoid", np.trapz)
    value = fn(y, x)
    return value.item() if hasattr(value, "item") else value


def _norm_density(phase: PhaseProfile, ell: int) -> np.ndarray:
    r = np.asarray(phase.r_m, dtype=float)
    W, V = np.asarray(phase.W_m3), np.asarray(phase.V_m2)
    h = np.asarray(phase.enthalpy_j_m3, dtype=float)
    nu, lam = np.asarray(phase.nu, dtype=float), np.asarray(phase.lambda_, dtype=float)
    W_over_r_sq = np.zeros(r.size, dtype=float)
    positive = r > 0.0
    W_over_r_sq[positive] = np.abs(W[positive] / r[positive]) ** 2
    if np.any(~positive) and np.any(np.abs(W[~positive]) >
            1e-12 * max(float(np.max(np.abs(W))), np.finfo(float).tiny)):
        raise ValueError("W must vanish at r=0 for a regular central profile")
    return h * np.exp(lam - nu) * (W_over_r_sq + ell * (ell + 1) * np.abs(V) ** 2)


def _interface_polarization(W: complex | float, V: complex | float, rt: float,
                            ell: int, relative_floor: float) -> dict[str, float | bool | None]:
    abs_w, abs_v = float(abs(W)), float(abs(V))
    scale = max(abs_w, abs_v * rt, np.finfo(float).tiny)
    reliable = abs_w > relative_floor * scale
    ratio = math.sqrt(ell * (ell + 1)) * rt * abs_v / abs_w if reliable else None
    return {"Pi": ratio, "reliable": reliable, "abs_W_m3": abs_w, "abs_V_m2": abs_v}


def diagnose_mode(profile: ModeProfile, *, relative_zero_tol: float = 1e-8,
                  interface_shell_fraction: float = 0.03) -> dict[str, object]:
    """Compute norm, localization and interface diagnostics; node counts remain external."""
    profile.validated()
    if not (math.isfinite(relative_zero_tol) and 0.0 <= relative_zero_tol < 0.1):
        raise ValueError("relative_zero_tol must satisfy 0 <= tol < 0.1")
    if not (math.isfinite(interface_shell_fraction) and 0.0 < interface_shell_fraction < 0.5):
        raise ValueError("interface_shell_fraction must lie in (0,0.5)")
    core, env = profile.core, profile.envelope
    rt, radius, ell = profile.interface_radius_m, profile.surface_radius_m, profile.ell
    dc, de = _norm_density(core, ell), _norm_density(env, ell)
    nc, ne = float(_trapz(dc, core.r_m)), float(_trapz(de, env.r_m))
    total = nc + ne
    if not (math.isfinite(total) and total > 0.0):
        raise ArithmeticError("mode displacement norm is nonpositive or nonfinite")
    shell = interface_shell_fraction * radius
    lo, hi = max(float(core.r_m[0]), rt - shell), min(float(env.r_m[-1]), rt + shell)
    cm, em = core.r_m >= lo, env.r_m <= hi
    sc = float(_trapz(dc[cm], core.r_m[cm])) if np.count_nonzero(cm) >= 2 else 0.0
    se = float(_trapz(de[em], env.r_m[em])) if np.count_nonzero(em) >= 2 else 0.0
    wi, wo, vi, vo = core.W_m3[-1], env.W_m3[0], core.V_m2[-1], env.V_m2[0]
    scale_w = max(float(np.max(np.abs(core.W_m3))),
                  float(np.max(np.abs(env.W_m3))), np.finfo(float).tiny)
    imag_w = max(float(np.max(np.abs(np.imag(core.W_m3)))),
                 float(np.max(np.abs(np.imag(env.W_m3)))))
    real_profile = imag_w <= 1e-12 * scale_w
    sign_reversal = (bool(float(np.real(wi)) * float(np.real(wo)) < 0.0)
                     if real_profile and min(abs(wi), abs(wo)) > 0 else None)
    return {
        "model_id": profile.model_id, "frequency_hz": float(profile.frequency_hz),
        "conversion": profile.conversion, "ell": int(ell), "norm": total,
        "norm_core": nc, "norm_envelope": ne,
        "core_trapping_fraction": nc / total, "envelope_trapping_fraction": ne / total,
        "interface_shell_fraction": (sc + se) / total,
        "interface_W_sign_reversal": sign_reversal,
        "interface_W_jump_m3": float(np.real(wo - wi)),
        "interface_V_jump_m2": float(np.real(vo - vi)),
        "interface_Pi_inner": _interface_polarization(wi, vi, rt, ell, relative_zero_tol),
        "interface_Pi_outer": _interface_polarization(wo, vo, rt, ell, relative_zero_tol),
        "spectral_branch_id": None, "physical_character": "unresolved",
        "identity_status": "DIAGNOSTICS_ONLY_NO_AUTOMATIC_F_P_G_LABEL",
    }


def _phase_coordinate(phase: PhaseProfile, region: str, rt: float,
                      radius: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = np.asarray(phase.r_m, dtype=float)
    if region == "core":
        x, drdx = r / rt, rt
    elif region == "envelope":
        x, drdx = (r - rt) / (radius - rt), radius - rt
    else:
        raise ValueError("region must be 'core' or 'envelope'")
    v1 = np.zeros(r.size, dtype=complex)
    mask = r > 0.0
    v1[mask] = np.asarray(phase.W_m3, dtype=complex)[mask] / r[mask]
    vec = np.column_stack((v1, np.asarray(phase.V_m2, dtype=complex)))
    weight = (np.asarray(phase.enthalpy_j_m3, dtype=float)
              * np.exp(np.asarray(phase.lambda_, dtype=float) - np.asarray(phase.nu, dtype=float))
              * drdx)
    return x, vec, weight


def phase_profile_overlap(a: ModeProfile, b: ModeProfile, *, samples: int = 512) -> dict[str, float]:
    """Phase-aligned overlap; only a single global eigenfunction phase is quotiented."""
    a.validated()
    b.validated()
    if a.ell != b.ell:
        raise ValueError("profiles must have the same ell")
    if a.conversion != b.conversion:
        raise ValueError("track slow and rapid conversion branches separately")
    if samples < 16:
        raise ValueError("samples must be >= 16")
    lfac = a.ell * (a.ell + 1)
    terms: dict[str, tuple[complex, float, float]] = {}
    for name in ("core", "envelope"):
        pa, pb = getattr(a, name), getattr(b, name)
        xa, qa, wa = _phase_coordinate(pa, name, a.interface_radius_m, a.surface_radius_m)
        xb, qb, wb = _phase_coordinate(pb, name, b.interface_radius_m, b.surface_radius_m)
        qa[:, 1] *= math.sqrt(lfac)
        qb[:, 1] *= math.sqrt(lfac)
        left, right = max(float(xa[0]), float(xb[0])), min(float(xa[-1]), float(xb[-1]))
        if right <= left:
            raise ValueError(f"no common normalized-coordinate support in {name}")
        x = np.linspace(left, right, samples)
        qai = np.column_stack([
            np.interp(x, xa, qa[:, k].real) + 1j * np.interp(x, xa, qa[:, k].imag)
            for k in range(2)
        ])
        qbi = np.column_stack([
            np.interp(x, xb, qb[:, k].real) + 1j * np.interp(x, xb, qb[:, k].imag)
            for k in range(2)
        ])
        wai, wbi = np.interp(x, xa, wa), np.interp(x, xb, wb)
        cross = complex(_trapz(np.sqrt(wai * wbi)
                               * np.sum(np.conjugate(qai) * qbi, axis=1), x))
        na = float(_trapz(wai * np.sum(np.abs(qai) ** 2, axis=1), x))
        nb = float(_trapz(wbi * np.sum(np.abs(qbi) ** 2, axis=1), x))
        terms[name] = (cross, na, nb)
    # Sum complex core+envelope overlaps first, then take one absolute value.
    # This retains relative phase/sign information across the sharp interface.
    cross_total = abs(sum(v[0] for v in terms.values()))
    na_total, nb_total = sum(v[1] for v in terms.values()), sum(v[2] for v in terms.values())
    overall = cross_total / math.sqrt(na_total * nb_total)
    out = {"overall": float(min(max(overall, 0.0), 1.0))}
    for name, (cross, na, nb) in terms.items():
        ov = abs(cross) / math.sqrt(na * nb) if na > 0.0 and nb > 0.0 else 0.0
        out[name] = float(min(max(ov, 0.0), 1.0))
    return out


def track_spectra(previous: Sequence[ModeProfile], current: Sequence[ModeProfile],
                  branch_ids: Sequence[str], *, min_overlap: float = 0.55,
                  ambiguity_margin: float = 0.05, frequency_weight: float = 0.02,
                  overlap_samples: int = 256) -> dict[str, object]:
    """Track spectral branches by overlap; physical modal character stays separate."""
    if len(previous) != len(branch_ids):
        raise ValueError("branch_ids must have one entry per previous mode")
    if not previous or not current:
        return {"matches": [], "unmatched_previous": list(branch_ids),
                "unmatched_current": list(range(len(current))), "overlap_matrix": []}
    if not (0.0 <= min_overlap <= 1.0 and 0.0 <= ambiguity_margin < 1.0):
        raise ValueError("invalid overlap or ambiguity threshold")
    if frequency_weight < 0.0 or not math.isfinite(frequency_weight):
        raise ValueError("frequency_weight must be finite and nonnegative")
    overlaps, cost = np.empty((len(previous), len(current))), np.empty((len(previous), len(current)))
    for i, p in enumerate(previous):
        p.validated()
        for j, c in enumerate(current):
            c.validated()
            sim = phase_profile_overlap(p, c, samples=overlap_samples)["overall"]
            overlaps[i, j] = sim
            cost[i, j] = (1.0 - sim) + frequency_weight * abs(math.log(c.frequency_hz / p.frequency_hz))
    rows, cols = linear_sum_assignment(cost)
    matches, used_rows, used_cols = [], set(), set()
    for i, j in zip(rows.tolist(), cols.tolist()):
        row_alt, col_alt = np.delete(overlaps[i, :], j), np.delete(overlaps[:, j], i)
        row_margin = float(overlaps[i, j] - np.max(row_alt)) if row_alt.size else 1.0
        col_margin = float(overlaps[i, j] - np.max(col_alt)) if col_alt.size else 1.0
        ambiguous = (overlaps[i, j] < min_overlap or row_margin < ambiguity_margin
                     or col_margin < ambiguity_margin)
        matches.append({
            "previous_index": i, "current_index": j,
            "spectral_branch_id": None if ambiguous else str(branch_ids[i]),
            "overlap": float(overlaps[i, j]), "row_margin": row_margin,
            "column_margin": col_margin,
            "assignment_status": "AMBIGUOUS_DO_NOT_FORCE_ID" if ambiguous else "MATCHED_BY_PROFILE_OVERLAP",
            "physical_character": "unresolved",
        })
        used_rows.add(i)
        used_cols.add(j)
    return {
        "matches": matches,
        "unmatched_previous": [str(branch_ids[i]) for i in range(len(previous)) if i not in used_rows],
        "unmatched_current": [j for j in range(len(current)) if j not in used_cols],
        "overlap_matrix": overlaps.tolist(),
        "cost_matrix": cost.tolist(),
        "interpretation": "spectral-branch continuation only; physical f/p1/interface-g character is not auto-assigned",
    }
