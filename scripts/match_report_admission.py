"""Independent, method-specific admission of native core--crust match reports.

The native kernels (``match-core css`` and ``maxwell-transition``) already
enforce their own gates.  This module is a *second, independent* layer: it
never infers success from ``status == "ok"`` and re-derives interface conditions
from serialized rows (the internal CSV produced by reading the written CompOSE
artefact back). Reported root counts and native interpolation diagnostics are
checked against their declared contracts; this module does not repeat the
native root search. Branch secants are necessary checks, not a certification of
the native higher-order interpolant or its differential first law between rows.

Every tolerance below is a gate that the native kernel already declares; none
is looser than its native counterpart.  ``SERIAL_FLOOR`` only absorbs the
round-off of reconstructing P, epsilon and mu from 17-digit serialized
(n_B, Q1, Q3, Q7) columns.

Row quantities (cold, one thermodynamic slice):
    P      = n_B * Q1
    eps    = n_B * m_n * (Q7 + 1)
    mu_E   = (eps + P) / n_B                      (Euler chemical potential)
    mu_Q3  = (Q3 + 1) m_n [+ Yq * Q4 * m_n if Il == 0]

For Il = 0 CompOSE stores Q3 = (mu_E - Yq mu_q)/m_n - 1 (charge-shifted), so
the charge shift is undone before Q3 is compared with mu_E.
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

SCHEMA = "ns-cowling.match-report-admission.v1"

CSS_CONSTRUCTION_KIND = "css_single_speed_with_crust_density_jump"

# Native declared gates (cited source in parentheses).
FIRST_LAW_GATE = 5.0e-9           # compose_eos_matching.f90 match_crust_core
CSS_CORE_ENDPOINT_GATE = 2.0e-10  # compose_eos_matching.f90 construct_bridge('css'), x max(1,|x|)
PRESSURE_REL_GATE = 1.0e-10       # compose_eos_phase_transition.f90 pressure_relative_tolerance
MU_REL_GATE = 1.0e-10             # compose_eos_phase_transition.f90 mu_tol
LATENT_IDENTITY_GATE = 1.0e-10    # compose_eos_phase_transition.f90 identity_residual
EULER_GATE = 5.0e-4               # compose_eos_phase_transition.f90 euler_tolerance
CAUSAL_GATE = 1.0e-8              # compose_eos_phase_transition.f90 causal_tolerance
# Round-off floor for quantities rebuilt from 17-significant-digit columns.
SERIAL_FLOOR = 64.0 * 2.220446049250313e-16

_COLUMNS = ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7")
_SPEED_KEYS = ("css_single_sound_speed", "cs2_min", "cs2_max", "cs2_left", "cs2_right")

_CSS_REQUIRED = (
    "method", "construction_kind", "n", "n_bridge", "n1", "n2", "P1", "P2", "mu1", "mu2",
    "cs2_left", "cs2_right", "cs2_min", "cs2_max", "css_single_sound_speed",
    "css_transition_nb_low", "css_transition_nb_high", "css_transition_delta_nb",
    "css_transition_delta_eps", "css_transition_epsilon_low",
    "css_transition_epsilon_high", "css_transition_pressure", "css_transition_mu",
    "css_transition_low_index", "css_transition_high_index", "css_core_bridge_index",
    "first_law_relative_residual", "source_euler_relative_residual",
)

_MAXWELL_REQUIRED = (
    "n", "pressure_crossings", "low_phase_rows", "high_phase_rows", "mu_transition",
    "pressure_transition", "n_low", "n_high", "epsilon_low", "epsilon_high",
    "density_jump", "latent_energy_density", "maxwell_identity_residual",
    "max_euler_residual_low", "max_euler_residual_high", "min_cs2_low",
    "max_cs2_low", "min_cs2_high", "max_cs2_high",
)


# --------------------------------------------------------------------- rows

def rows_from_table(table: Iterable[Sequence[Any] | Mapping[str, Any]],
                    meta: Mapping[str, Any]) -> list[dict[str, float]]:
    """Convert internal-CSV rows (sequence in column order, or mapping) to physics rows."""
    try:
        mn = float(meta["mn"])
        il = float(meta["il"])
        yq = float(meta["charge_fraction_grid"])
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError("invalid thermodynamic normalization metadata") from exc
    if not _finite(mn, il, yq) or mn <= 0.0 or il not in (0.0, 1.0):
        raise ValueError("require finite mn>0, Il in {0,1}, and finite Yq")
    out = []
    for row_number, raw in enumerate(table, 1):
        try:
            if isinstance(raw, Mapping):
                vals = {c: float(raw[c]) for c in _COLUMNS}
            else:
                vals = {c: float(raw[i]) for i, c in enumerate(_COLUMNS)}
        except (KeyError, IndexError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"invalid thermodynamic columns at row {row_number}") from exc
        if not _finite(*vals.values()):
            raise ValueError(f"nonfinite thermodynamic columns at row {row_number}")
        n = vals["nb"]
        p = n * vals["Q1"]
        eps = n * mn * (vals["Q7"] + 1.0)
        mu_e = (eps + p) / n if n != 0.0 else math.nan
        mu_q3 = (vals["Q3"] + 1.0) * mn
        if il == 0:
            mu_q3 += yq * vals["Q4"] * mn
        out.append({"n": n, "P": p, "eps": eps, "mu_e": mu_e, "mu_q3": mu_q3})
    return out


def _parse_metadata(path: Path) -> dict[str, Any]:
    meta: dict[str, Any] = {}
    text = path.read_text(encoding="utf-8")
    for key in ("mn", "il", "charge_fraction_grid"):
        m = re.search(rf"^\s*{key}\s*=\s*([^\s,/]+)", text, flags=re.MULTILINE | re.IGNORECASE)
        if m is None:
            raise ValueError(f"metadata key {key!r} missing in {path}")
        token = m.group(1).replace("D", "E").replace("d", "e")
        meta[key] = int(token) if key == "il" else float(token)
    return meta


def load_rows(csv_path: Path, nml_path: Path) -> tuple[list[dict[str, float]], dict[str, Any]]:
    meta = _parse_metadata(Path(nml_path))
    with Path(csv_path).open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = [c.strip() for c in next(reader, [])]
        if (not header or len(set(header)) != len(header) or "" in header
                or not set(_COLUMNS).issubset(header)):
            raise ValueError("missing or duplicate internal CSV columns")
        table = []
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f"internal CSV row width mismatch at line {line}")
            table.append(dict(zip(header, row)))
    return rows_from_table(table, meta), meta


# ------------------------------------------------------------------ helpers

def _finite(*xs: Any) -> bool:
    try:
        return all(math.isfinite(float(x)) for x in xs)
    except (TypeError, ValueError, OverflowError):
        return False


def _close(a: float, b: float, rel: float, scale: float | None = None) -> tuple[bool, float]:
    """|a-b| <= rel*max(|a|,|b|) + SERIAL_FLOOR*scale. Returns (ok, relative residual)."""
    if not _finite(a, b):
        return False, math.inf
    mag = max(abs(a), abs(b))
    scale = mag if scale is None else abs(scale)
    diff = abs(a - b)
    ok = diff <= rel * mag + SERIAL_FLOOR * scale
    return ok, diff / mag if mag > 0.0 else diff


def _jsonable(v: Any) -> Any:
    if isinstance(v, float) and not math.isfinite(v):
        return repr(v)
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    return v


class _Checks:
    def __init__(self) -> None:
        self.checks: dict[str, dict[str, Any]] = {}

    def add(self, name: str, passed: bool, **detail: Any) -> None:
        entry = {"pass": bool(passed)}
        entry.update({k: _jsonable(v) for k, v in detail.items()})
        self.checks[name] = entry

    def result(self, method: str) -> dict[str, Any]:
        failed = sorted(k for k, v in self.checks.items() if not v["pass"])
        return {
            "schema": SCHEMA,
            "method": method,
            "status": "ADMITTED" if self.checks and not failed else "REJECTED",
            "failed_checks": failed,
            "checks": self.checks,
            "gates": {
                "first_law_relative_residual": FIRST_LAW_GATE,
                "css_core_endpoint": CSS_CORE_ENDPOINT_GATE,
                "pressure_relative": PRESSURE_REL_GATE,
                "mu_relative": MU_REL_GATE,
                "latent_identity_relative": LATENT_IDENTITY_GATE,
                "euler_relative": EULER_GATE,
                "causal": CAUSAL_GATE,
                "serial_floor": SERIAL_FLOOR,
            },
        }


def _missing(report: Mapping[str, Any], keys: Sequence[str]) -> list[str]:
    missing = []
    for k in keys:
        if k not in report:
            missing.append(k)
        elif k in ("method", "construction_kind") or k in _SPEED_KEYS:
            continue  # strings, or numerics whose range is judged by a named check
        elif not _finite(report[k]):
            missing.append(k)
    return missing


def _int_index(report: Mapping[str, Any], key: str) -> int | None:
    try:
        f = float(report.get(key))
    except (TypeError, ValueError, OverflowError):
        return None
    if not math.isfinite(f) or f != int(f):
        return None
    return int(f)


def _branch_checks(c: _Checks, rows: list[dict[str, float]], name: str) -> None:
    """Necessary row/secant checks, not certification of a higher-order interpolant."""
    worst_cs2 = -math.inf
    stable = causal = len(rows) >= 2
    first_bad = None
    for a, b in zip(rows, rows[1:]):
        dn, dp, de = b["n"] - a["n"], b["P"] - a["P"], b["eps"] - a["eps"]
        if not (_finite(dn, dp, de) and dn > 0.0 and dp > 0.0 and de > 0.0):
            stable = False
            first_bad = first_bad or [a["n"], b["n"]]
            continue
        cs2 = dp / de
        worst_cs2 = max(worst_cs2, cs2)
        if cs2 > 1.0 + CAUSAL_GATE:
            causal = False
    c.add(f"{name}_stable", stable, rows=len(rows), first_violation_nb=first_bad)
    c.add(f"{name}_causal", causal, max_secant_cs2=worst_cs2, gate=CAUSAL_GATE)


def _physical_rows(c: _Checks, rows: list[dict[str, float]]) -> bool:
    """Check all retained rows before divisions, logs or fractional powers."""
    invalid = []
    for i, row in enumerate(rows):
        if (not _finite(*(row.get(k) for k in ("n", "P", "eps", "mu_e", "mu_q3")))
                or any(row[k] <= 0.0 for k in ("n", "P", "eps", "mu_e"))):
            invalid.append(i)
    valid = bool(rows) and not invalid
    c.add("serialized_physical_rows", valid, invalid_rows=invalid)
    if not valid:
        return False
    worst = max(_close(row["mu_q3"], row["mu_e"], EULER_GATE)[1] for row in rows)
    c.add("serialized_euler_mu", worst <= EULER_GATE, max_rel_residual=worst, gate=EULER_GATE)
    return True


def _float_or_nan(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError, OverflowError):
        return math.nan


def _report_vs_rows(c: _Checks, pairs: Mapping[str, tuple[float, float, float | None]]) -> None:
    mism = {}
    for key, (a, b, scale) in pairs.items():
        ok, r = _close(a, b, LATENT_IDENTITY_GATE, scale=scale)
        if not ok:
            mism[key] = r
    c.add("report_matches_rows", not mism, mismatched=mism)


# ---------------------------------------------------------------------- CSS

def admit_css(report: Mapping[str, Any], rows: list[dict[str, float]],
              meta: Mapping[str, Any]) -> dict[str, Any]:
    c = _Checks()
    missing = _missing(report, _CSS_REQUIRED)
    c.add("report_fields_present", not missing, missing=missing)
    if missing:
        return c.result("css")

    c.add("method", report["method"] == "css", value=report["method"])
    c.add("construction_kind", report["construction_kind"] == CSS_CONSTRUCTION_KIND,
          value=report["construction_kind"])

    cs2 = _float_or_nan(report["css_single_sound_speed"])
    cs2_ok = _finite(cs2) and 0.0 < cs2 <= 1.0
    c.add("css_sound_speed_range", cs2_ok, value=cs2)
    speeds = [_float_or_nan(report[k]) for k in _SPEED_KEYS[1:]]
    single = _finite(cs2, *speeds) and all(_close(s, cs2, SERIAL_FLOOR)[0] for s in speeds)
    c.add("css_single_sound_speed", single, reported=speeds)

    flr = float(report["first_law_relative_residual"])
    c.add("first_law_relative_residual", 0.0 <= flr <= FIRST_LAW_GATE,
          value=flr, gate=FIRST_LAW_GATE)
    # Output mu_B is re-derived from F + P/n_B, so the serialized rows cannot
    # reveal an inconsistent source chemical potential; the native report can.
    ser = float(report["source_euler_relative_residual"])
    c.add("source_euler_relative_residual", 0.0 <= ser <= EULER_GATE,
          value=ser, gate=EULER_GATE)

    lo = _int_index(report, "css_transition_low_index")
    hi = _int_index(report, "css_transition_high_index")
    core = _int_index(report, "css_core_bridge_index")
    n_bridge = _int_index(report, "n_bridge")
    n_rows = _int_index(report, "n")
    idx_ok = (
        None not in (lo, hi, core, n_bridge, n_rows)
        and n_rows == len(rows) and 0 <= lo and hi == lo + 1 and n_bridge >= 8
        and core == hi + n_bridge + 1 and core < len(rows)
    )
    c.add("row_indices_consistent", idx_ok, low=lo, high=hi, core=core,
          n_bridge=n_bridge, report_n=n_rows, rows=len(rows))
    if not idx_ok or not cs2_ok:
        return c.result("css")
    if not _physical_rows(c, rows):
        return c.result("css")
    _branch_checks(c, rows[:lo + 1], "crust_branch")
    _branch_checks(c, rows[hi:core + 1], "css_branch")
    _branch_checks(c, rows[core:], "core_branch")

    rl, rh, rc = rows[lo], rows[hi], rows[core]
    p_t, mu_t = rl["P"], rl["mu_e"]

    # Crust -> CSS first-order jump: common P and Euler mu, positive jump,
    # Delta eps = mu_t Delta n_B.
    ok, r = _close(rh["P"], rl["P"], PRESSURE_REL_GATE)
    c.add("crust_transition_pressure", ok, rel_residual=r)
    ok, r = _close(rh["mu_e"], rl["mu_e"], MU_REL_GATE)
    c.add("crust_transition_euler_mu", ok, rel_residual=r)
    dn, de = rh["n"] - rl["n"], rh["eps"] - rl["eps"]
    c.add("crust_transition_positive_jump", _finite(dn, de) and dn > 0.0 and de > 0.0,
          delta_nb=dn, delta_eps=de)
    ok, r = _close(de, mu_t * dn, LATENT_IDENTITY_GATE, scale=rh["eps"])
    c.add("crust_transition_latent_identity", ok and de > 0.0, rel_residual=r)

    # One CSS speed over the bridge interior: n = n_h (mu/mu_t)^(1/cs2) and
    # P - P_t = cs2/(1+cs2) (n mu - n_h mu_t) on every serialized bridge row.
    n_h = rh["n"]
    worst = 0.0
    seg_ok = core - hi - 1 == n_bridge
    for row in rows[hi + 1:core]:
        x, y = row["mu_e"] / mu_t, row["n"] / n_h
        if not (_finite(x, y) and x > 1.0 and y > 1.0):
            seg_ok = False
            continue
        implied = math.log(x) / math.log(y)
        p_pred = p_t + cs2 / (1.0 + cs2) * (row["n"] * row["mu_e"] - n_h * mu_t)
        ok1, r1 = _close(implied, cs2, FIRST_LAW_GATE)
        ok2, r2 = _close(row["P"], p_pred, FIRST_LAW_GATE, scale=row["eps"])
        seg_ok = seg_ok and ok1 and ok2
        worst = max(worst, r1, r2)
    c.add("css_segment_single_speed", seg_ok, bridge_rows=core - hi - 1,
          max_rel_residual=worst, gate=FIRST_LAW_GATE)

    # CSS -> core continuity in n_B, P, epsilon and mu_B (native u=1 gate).
    mu_c = rc["mu_e"]
    if _finite(mu_c) and mu_c > 0.0:
        try:
            n_pred = n_h * (mu_c / mu_t) ** (1.0 / cs2)
            p_pred = p_t + cs2 / (1.0 + cs2) * (n_pred * mu_c - n_h * mu_t)
            e_pred = n_pred * mu_c - p_pred
        except OverflowError:
            n_pred = p_pred = e_pred = math.nan
    else:
        n_pred = p_pred = e_pred = math.nan

    def _endpoint(a: float, b: float) -> tuple[bool, float]:
        if not _finite(a, b):
            return False, math.inf
        d = abs(a - b)
        return d <= CSS_CORE_ENDPOINT_GATE * max(1.0, abs(b)) + SERIAL_FLOOR * abs(rc["eps"]), d

    res = {
        "n": _endpoint(rc["n"], n_pred),
        "P": _endpoint(rc["P"], p_pred),
        "eps": _endpoint(rc["eps"], e_pred),
        "mu_vs_report_mu2": _endpoint(mu_c, float(report["mu2"])),
        "n_vs_report_n2": _endpoint(rc["n"], float(report["n2"])),
        "P_vs_report_P2": _endpoint(rc["P"], float(report["P2"])),
    }
    c.add("css_core_continuity", all(v[0] for v in res.values()),
          abs_residuals={k: v[1] for k, v in res.items()}, gate=CSS_CORE_ENDPOINT_GATE)

    _report_vs_rows(c, {
        "css_transition_pressure": (float(report["css_transition_pressure"]), p_t, None),
        "css_transition_mu": (float(report["css_transition_mu"]), mu_t, None),
        "css_transition_nb_low": (float(report["css_transition_nb_low"]), rl["n"], None),
        "css_transition_nb_high": (float(report["css_transition_nb_high"]), rh["n"], None),
        "css_transition_epsilon_low": (float(report["css_transition_epsilon_low"]), rl["eps"], None),
        "css_transition_epsilon_high": (float(report["css_transition_epsilon_high"]), rh["eps"], None),
        "css_transition_delta_nb": (float(report["css_transition_delta_nb"]), dn, rh["n"]),
        "css_transition_delta_eps": (float(report["css_transition_delta_eps"]), de, rh["eps"]),
        "n1": (float(report["n1"]), rl["n"], None),
        "P1": (float(report["P1"]), p_t, None),
        "mu1": (float(report["mu1"]), mu_t, None),
    })
    return c.result("css")


# ------------------------------------------------------------------ Maxwell

def admit_maxwell(report: Mapping[str, Any], rows: list[dict[str, float]],
                  meta: Mapping[str, Any]) -> dict[str, Any]:
    c = _Checks()
    missing = _missing(report, _MAXWELL_REQUIRED)
    c.add("report_fields_present", not missing, missing=missing)
    if missing:
        return c.result("maxwell")

    crossings = _int_index(report, "pressure_crossings")
    c.add("pressure_crossing", crossings == 1, value=crossings)

    n_low = _int_index(report, "low_phase_rows")
    n_high = _int_index(report, "high_phase_rows")
    n_rows = _int_index(report, "n")
    idx_ok = (None not in (n_low, n_high, n_rows) and n_low >= 2 and n_high >= 2
              and n_low + n_high == len(rows) == n_rows)
    c.add("row_indices_consistent", idx_ok, low_phase_rows=n_low,
          high_phase_rows=n_high, report_n=n_rows, rows=len(rows))
    if not idx_ok:
        return c.result("maxwell")
    if not _physical_rows(c, rows):
        return c.result("maxwell")

    rl, rh = rows[n_low - 1], rows[n_low]
    p_t, mu_t = rl["P"], rl["mu_e"]
    ok, r = _close(rh["P"], rl["P"], PRESSURE_REL_GATE)
    c.add("coexistence_pressure", ok, rel_residual=r)
    ok, r = _close(rh["mu_e"], rl["mu_e"], MU_REL_GATE)
    c.add("coexistence_euler_mu", ok, rel_residual=r)
    dn, de = rh["n"] - rl["n"], rh["eps"] - rl["eps"]
    c.add("coexistence_positive_jump", _finite(dn, de) and dn > 0.0 and de > 0.0,
          delta_nb=dn, delta_eps=de)
    ok, r = _close(de, mu_t * dn, LATENT_IDENTITY_GATE, scale=rh["eps"])
    c.add("coexistence_latent_identity", ok and de > 0.0, rel_residual=r)

    idr = float(report["maxwell_identity_residual"])
    c.add("report_identity_residual", 0.0 <= idr <= LATENT_IDENTITY_GATE, value=idr,
          note="native value is built from mu_t*n-P_t; the independent check is "
               "coexistence_latent_identity")

    eul = [float(report["max_euler_residual_low"]), float(report["max_euler_residual_high"])]
    c.add("report_euler_residuals", all(0.0 <= e <= EULER_GATE for e in eul),
          values=eul, gate=EULER_GATE)
    cs = [float(report[k]) for k in ("min_cs2_low", "max_cs2_low", "min_cs2_high", "max_cs2_high")]
    c.add("report_branch_sound_speeds",
          0.0 < cs[0] <= cs[1] <= 1.0 + CAUSAL_GATE
          and 0.0 < cs[2] <= cs[3] <= 1.0 + CAUSAL_GATE, values=cs)
    _branch_checks(c, rows[:n_low], "low_branch")
    _branch_checks(c, rows[n_low:], "high_branch")

    _report_vs_rows(c, {
        "mu_transition": (float(report["mu_transition"]), mu_t, None),
        "pressure_transition": (float(report["pressure_transition"]), p_t, None),
        "n_low": (float(report["n_low"]), rl["n"], None),
        "n_high": (float(report["n_high"]), rh["n"], None),
        "epsilon_low": (float(report["epsilon_low"]), rl["eps"], None),
        "epsilon_high": (float(report["epsilon_high"]), rh["eps"], None),
        "density_jump": (float(report["density_jump"]), dn, rh["n"]),
        "latent_energy_density": (float(report["latent_energy_density"]), de, rh["eps"]),
    })
    return c.result("maxwell")


# ----------------------------------------------------------------- dispatch

def admit_match_report(method: str, report: Mapping[str, Any],
                       csv_path: Path, nml_path: Path) -> dict[str, Any]:
    """Admit a native match report against serialized readback rows."""
    if method not in ("css", "maxwell"):
        c = _Checks()
        c.add("method_supported", False, value=method)
        return c.result(method)
    try:
        rows, meta = load_rows(Path(csv_path), Path(nml_path))
    except (OSError, ValueError, KeyError, UnicodeError, csv.Error) as exc:
        c = _Checks()
        c.add("serialized_rows_readable", False, error=f"{type(exc).__name__}: {exc}")
        return c.result(method)
    if method == "css":
        return admit_css(report, rows, meta)
    return admit_maxwell(report, rows, meta)


# ----------------------------------------------------------- failure taxonomy

_RULES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("implementation_convention", (
        "rest-mass", "mass conventions", "nucleon masses", "lepton-content",
        "charge-fraction conventions")),
    ("implementation_scope", (
        "composition", "nadd", "not yet integrated", "not yet ported",
        "multiple pressure-equality roots", "exactly one isolated",
        "locally charge neutral", "native matching currently requires cold")),
    ("input_contract", (
        "beta-equilibrat", "beta-equilibrium", "requires t=0", "zero entropy",
        "q6 and q7", "q6=q7", "lepton potential", "iphase", "single smooth phase",
        "provider q3", "each phase must have strictly increasing",
        "at least three source points", "core table has no point", "phase eos",
        "not enough points", "invalid crust eos", "invalid core eos")),
    ("physical", (
        "no causal single-speed css", "no pressure equality", "no overlapping",
        "insufficient common", "does not connect a distinct denser",
        "no positive crust-to-bridge density jump", "nonpositive energy density or jump",
        "endpoints must increase", "no convex matching window",
        "endpoint sound speed is outside", "violates 0", "causal bound",
        "fewer than two retained source knots", "coexistence energy densities are nonpositive")),
    ("unresolved_source_or_interpolant_shape", ("interpolation violates",)),
    ("numerical", (
        "first-law integral tolerance", "did not converge", "residual check", "nonfinite",
        "overflow", "has no integral root", "failed continuous core endpoint")),
)

_NATIVE_CATEGORY = {
    "COMPOSITION_UNSUPPORTED": "implementation_scope",
    "MASS_CONVENTION_FAIL": "implementation_convention",
    "BETA_EQUILIBRIUM_FAIL": "input_contract",
    "CAUSALITY_FAIL": "physical",
    "NO_COMMON_MU_SUPPORT": "physical",
    "NO_COEXISTENCE_CROSSING": "physical",
    "MULTIPLE_COEXISTENCE_CROSSINGS": "implementation_scope",
    "NONPOSITIVE_DENSITY_JUMP": "physical",
    "INTERPOLANT_SHAPE_FAIL": "unresolved_source_or_interpolant_shape",
}


def classify_failure(stage: str | None, message: str | None, native_category: str | None) -> str:
    """Separate data-provenance, contract, implementation, numerical and physical failures.

    Unmatched messages are ``unclassified`` -- never silently counted as physics.
    """
    s = (message or "").lower()
    if stage == "source" or "missing compose source" in s:
        return "data_provenance"
    if stage == "admission":
        return "admission"
    if stage in ("readback", "write"):
        return "implementation"
    if stage == "consistency":
        return "native_consistency_gate"
    for label, needles in _RULES:
        if any(n in s for n in needles):
            return label
    if native_category and native_category in _NATIVE_CATEGORY:
        return _NATIVE_CATEGORY[native_category]
    return "unclassified"
