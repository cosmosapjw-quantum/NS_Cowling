"""Select finite TF candidates without hiding unresolved physical seed families.

A family is (spatial dimension, initial topology), not an arbitrary optimizer
restart. Failed attempts remain evidence. An infeasible trial's low energy is
not a competing physical branch. Units are fm^-3 and MeV per baryon.
"""
from __future__ import annotations
import math

STATIONARY = frozenset(("STATIONARY", "UNIFORM_STATIONARY"))
SUCCESS = "LOWEST_EVALUATED_STATIONARY_CANDIDATE"

def _number(row, key):
    try:
        x=float(row[key])
        return x if math.isfinite(x) else None
    except (KeyError,TypeError,ValueError):
        return None

def _family(row):
    try:
        return (int(row["dimension"]),str(row["seed"]))
    except (KeyError,TypeError,ValueError):
        return None

def feasibility(row, n_b):
    """Check comparable density/charge/field feasibility, independently of status."""
    required=("n_b","baryon_constraint_error","charge_constraint_error",
              "scalar_pde_residual","vector_pde_residual","electron_pde_residual")
    values={k:_number(row,k) for k in required}
    if any(v is None for v in values.values()):
        return "UNKNOWN"
    density_tol=max(1e-12,1e-8*n_b)
    if (abs(values["n_b"]-n_b)>density_tol
        or abs(values["baryon_constraint_error"])>density_tol
        or abs(values["charge_constraint_error"])>density_tol
        or abs(values["scalar_pde_residual"])>1e-8
        or abs(values["vector_pde_residual"])>1e-8
        or abs(values["electron_pde_residual"])>2e-10):
        return "INFEASIBLE"
    return "FEASIBLE"

def select_geometry_candidates(attempts, expected_families, *, n_b=None,
                               energy_tolerance=1e-6):
    """Return selected=None on incomplete coverage or a lower unresolved branch.

    best is a diagnostic candidate even on HOLD; only selected is exportable.
    Unknown feasibility cannot certify a candidate or dismiss a lower energy.
    """
    expected=[(int(d),str(seed)) for d,seed in expected_families]
    if not expected or len(set(expected))!=len(expected):
        raise ValueError("distinct nonempty expected geometry families required")
    if not math.isfinite(energy_tolerance) or energy_tolerance<0:
        raise ValueError("finite nonnegative energy tolerance required")
    if n_b is None:
        n_b=next((_number(a,"n_b") for a in attempts if _number(a,"n_b") is not None),None)
    if n_b is None or not math.isfinite(n_b) or n_b<=0:
        raise ValueError("positive requested baryon density required")
    diagnoses=[];valid=[]
    for i,a in enumerate(attempts):
        energy=_number(a,"energy_per_baryon")
        field_status=feasibility(a,n_b)
        stationary=(a.get("numerical_status") in STATIONARY and
                    not a.get("muon_threshold_crossed",False) and
                    energy is not None and _number(a,"epsilon") is not None and
                    field_status=="FEASIBLE")
        diagnoses.append({"attempt":i,"family":_family(a),"feasibility":field_status,
                          "stationary_supported":stationary})
        if stationary:valid.append(i)
    best_index=min(valid,key=lambda i:attempts[i]["energy_per_baryon"]) if valid else None
    best=attempts[best_index] if best_index is not None else None
    resolved={_family(attempts[i]) for i in valid}
    missing=[k for k in expected if k not in resolved]
    lower=[]
    if best is not None:
        for i,a in enumerate(attempts):
            energy=_number(a,"energy_per_baryon")
            if (i not in valid and energy is not None
                and diagnoses[i]["feasibility"]!="INFEASIBLE"
                and energy<best["energy_per_baryon"]-energy_tolerance):
                lower.append(i)
    status=("NO_STATIONARY_CANDIDATE" if best is None else
            "HOLD_UNRESOLVED_GEOMETRY" if missing else
            "LOWER_UNCONVERGED_ATTEMPT" if lower else SUCCESS)
    return {"status":status,"best":best,"selected":best if status==SUCCESS else None,
            "best_attempt_index":best_index,"unresolved_families":missing,
            "lower_unconverged_attempt_indices":lower,"attempt_diagnostics":diagnoses,
            "energy_tolerance_MeV_per_baryon":energy_tolerance,
            "scope":"Finite physical seed-family coverage; no global-minimum certificate"}
