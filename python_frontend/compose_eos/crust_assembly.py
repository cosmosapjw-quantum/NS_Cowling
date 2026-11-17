"""Explicit crust-piece assembly without changing physical energy references.

CompOSE's neutron mass header is a normalization convention. Re-expressing
dimensionless Q fields below preserves dimensional epsilon, pressure and every
chemical potential; it is NOT epsilon -> epsilon + C n_B. A first-order seam
must be supplied as two actual one-sided states at a common positive P and mu.
No interpolation, fitting, energy anchoring or phase optimization occurs here.
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from .eos_data import EOSData

ROW_FIELDS = (
    "nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra",
    "Iphase", "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh",
    "row_status", "Yn_present", "Yp_present", "Ye_present", "Ymu_present",
    "Abar_present", "Zbar_present", "Yh_present", "C3_applicable", "C4_applicable",
)

def _cold_slice(data: EOSData) -> None:
    if data.temperature_mev != 0 or data.q_full.size or data.q_extra_full.size:
        raise ValueError("crust assembly accepts a single cold slice only")
    if data.n < 1 or np.any(~np.isfinite(data.nb)) or np.any(data.nb <= 0):
        raise ValueError("finite positive baryon density required")
    if np.any(np.diff(data.nb) <= 0):
        raise ValueError("branch density must increase strictly")
    if not np.isfinite(data.mn) or not np.isfinite(data.mp) or min(data.mn,data.mp)<=0:
        raise ValueError("positive finite header masses required")
    for name in ROW_FIELDS:
        value=np.asarray(getattr(data,name))
        if value.shape != (data.n,):
            raise ValueError("only one-dimensional row field supported: "+name)
    for name in ("p","eps","mub","muq","mul"):
        if not np.all(np.isfinite(getattr(data,name))):
            raise ValueError("nonfinite physical field: "+name)

def rebase_header(data: EOSData, *, mn: float, mp: float) -> EOSData:
    """Change only serialization normalization; no physical rest energy shift."""
    _cold_slice(data)
    if not np.isfinite(mn) or not np.isfinite(mp) or min(mn,mp)<=0:
        raise ValueError("positive finite target header masses required")
    copied={name:np.array(getattr(data,name),copy=True) for name in ROW_FIELDS}
    copied.update(
        Q3=data.mub/mn-1., Q4=data.muq/mn, Q5=data.mul/mn,
        Q6=data.F/mn-1., Q7=data.eps/(data.nb*mn)-1.,
    )
    return replace(data, mn=float(mn),mp=float(mp),source_dir=None,**copied)

def select_rows(data: EOSData, indices) -> EOSData:
    """Copy row fields together, including sparse composition masks."""
    _cold_slice(data)
    idx=np.atleast_1d(np.arange(data.n)[indices])
    if idx.size==0:
        raise ValueError("empty branch selection")
    out=replace(data,source_dir=None,**{name:np.array(getattr(data,name)[idx],copy=True)
                                      for name in ROW_FIELDS})
    _cold_slice(out)
    return out

def certify_interface(left: EOSData, right: EOSData, *,
                      pressure_rtol: float=1e-7, mu_atol_MeV: float=1e-6,
                      euler_atol_MeV: float=1e-7) -> dict:
    """Validate thermodynamic identities of last/first endpoint rows.

    Local neutrality and beta equilibrium must be established by the upstream
    physical solvers; this helper cannot infer them from sparse TF masks.
    This certifies endpoint identities only, not global minimization or mesh
    convergence. Positive density jump and equal P, mu establish the common
    tangent. The interval between the endpoints is explicitly a phase gap.
    """
    _cold_slice(left);_cold_slice(right)
    if not all(np.isfinite(v) and v>0 for v in (pressure_rtol,mu_atol_MeV,euler_atol_MeV)):
        raise ValueError("positive interface tolerances required")
    nl,nr=float(left.nb[-1]),float(right.nb[0])
    pl,pr=float(left.p[-1]),float(right.p[0])
    ml,mr=float(left.mub[-1]),float(right.mub[0])
    el,er=float(left.eps[-1]),float(right.eps[0])
    if not 0<nl<nr:
        raise ValueError("first-order interface requires positive density jump")
    if min(pl,pr)<=0 or min(ml,mr)<=0:
        raise ValueError("positive interface pressure and chemical potential required")
    p_error=abs(pl-pr)/max(abs(pl),abs(pr))
    mu_error=abs(ml-mr)
    euler=max(abs((el+pl)/nl-ml),abs((er+pr)/nr-mr))
    if p_error>pressure_rtol:
        raise ValueError("interface pressure does not match")
    if mu_error>mu_atol_MeV:
        raise ValueError("interface chemical potential does not match")
    if euler>euler_atol_MeV:
        raise ValueError("endpoint violates cold Euler identity")
    return dict(status="ENDPOINT_COMMON_TANGENT_CHECKED",n_left=nl,n_right=nr,
                pressure_left=pl,pressure_right=pr,mu_left=ml,mu_right=mr,
                relative_pressure_error=p_error,mu_error_MeV=mu_error,
                euler_error_MeV=euler,
                interpolation_across_density_gap="FORBIDDEN",
                global_minimum="NOT_ESTABLISHED_BY_ENDPOINT_CHECK",
                neutrality_and_beta_equilibrium="UPSTREAM_SOLVER_PRECONDITION_NOT_CHECKED_HERE")

def assemble_crust(branches, *, labels, mn: float, mp: float,
                   pressure_rtol: float=1e-7, mu_atol_MeV: float=1e-6,
                   euler_atol_MeV: float=1e-7):
    """Concatenate only branches with explicitly supplied certified interfaces."""
    branches=list(branches);labels=list(labels)
    if len(branches)<2 or len(labels)!=len(branches) or len(set(labels))!=len(labels):
        raise ValueError("two or more uniquely labelled branches required")
    normalized=[rebase_header(b,mn=mn,mp=mp) for b in branches]
    seams=[]
    offset=0
    for index,(left,right) in enumerate(zip(normalized[:-1],normalized[1:])):
        cert=certify_interface(left,right,pressure_rtol=pressure_rtol,
                mu_atol_MeV=mu_atol_MeV,euler_atol_MeV=euler_atol_MeV)
        offset+=left.n
        cert.update(left_label=labels[index],right_label=labels[index+1],
                    left_row_zero_based=offset-1,right_row_zero_based=offset)
        seams.append(cert)
    first=normalized[0]
    result=replace(first,source_dir=None,
        composition_file_present=all(b.composition_file_present for b in normalized),
        composition_available=all(b.composition_available for b in normalized),
        beta_equilibrated=all(b.beta_equilibrated for b in normalized),
        **{name:np.concatenate([getattr(b,name) for b in normalized]) for name in ROW_FIELDS})
    _cold_slice(result)
    return result,dict(status="ASSEMBLED_EXPLICIT_PHASE_PIECES",labels=labels,
        row_count=result.n,interfaces=seams,
        physical_reference_change="NONE; header normalization only",
        interpolation_policy="interpolate each branch separately; never across interfaces")
