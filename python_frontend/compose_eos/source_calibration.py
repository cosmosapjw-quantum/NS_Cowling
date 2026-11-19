"""Read-only source-calibration diagnostics for cold one-charge NS matter.

Never relabel a sampled nuclear functional as the provider's original
unified core without comparisons at IDENTICAL baryon densities.
Units: n [fm^-3], P and epsilon [MeV fm^-3], mu [MeV].
Metric choice does not enter these thermodynamic algebraic calculations.
"""
from __future__ import annotations

import numpy as np

from .whole_crust import PiecewisePolytrope


def _original_rows(n, P, epsilon, mu):
    arrays=[np.asarray(x,dtype=np.float64).reshape(-1) for x in (n,P,epsilon,mu)]
    ns,ps,es,ms=arrays
    if (len(ns)<2 or len(set(map(len,arrays)))!=1
        or any(not np.all(np.isfinite(a)) for a in arrays)
        or np.any(ns<=0) or np.any(ps<=0) or np.any(es<=0)
        or np.any(np.diff(ns)<=0)):
        raise ValueError("positive finite, strictly increasing density source required")
    return ns,ps,es,ms


def phase_contiguous_slices(phase):
    """Return source label, start, stop for each contiguous phase component.

    Repeated phase IDs separated by another phase stay different components.
    This must precede thermo PCHIP or pressure-primary reconstruction.
    """
    values=np.asarray(phase)
    if values.ndim!=1 or len(values)==0 or np.any(~np.isfinite(values)):
        raise ValueError("finite one-dimensional nonempty phase labels required")
    if not np.all(values==values.astype(int)):
        raise ValueError("integer-valued original phase labels required")
    ids=values.astype(int)
    switches=np.r_[0,np.flatnonzero(ids[1:]!=ids[:-1])+1,len(ids)]
    return [(int(ids[a]),int(a),int(b)) for a,b in zip(switches[:-1],switches[1:])]


def first_law_segment_audit(n,P,epsilon,mu_Q3):
    """First-law reconstruction on exactly ONE monotone source phase.

    F(n)=epsilon/n, dF/dn=P/n²; sample pressure is unmodified,
    F(n0) is *exactly* anchored at the source's first epsilon/n.
    No density shifts, no arbitrary energy offset, no across-phase fit.
    """
    n,P,eps,mu=_original_rows(n,P,epsilon,mu_Q3)
    if np.any(np.diff(P)<=0):
        raise ValueError("strictly increasing pressure required on a phase segment")
    model=PiecewisePolytrope(n,P,eps[0]/n[0])
    validation=model.validate()
    reconstructed_epsilon=n*model.F
    reconstructed_mu=model.mu
    raw_euler=(eps+P)/n
    delta_f=model.F-eps/n
    return dict(n_source=n.copy(),
        p_source=P.copy(),epsilon_source=eps.copy(),
        source_mu_Q3=mu.copy(),source_mu_Euler=raw_euler.copy(),
        corrected_pressure=P.copy(),corrected_epsilon=reconstructed_epsilon.copy(),
        corrected_mu=reconstructed_mu.copy(),
        corrected_f_energy_per_baryon=model.F.copy(),
        source_anchor_energy_per_baryon_MeV=float(eps[0]/n[0]),
        max_abs_delta_epsilon_per_baryon_MeV=float(np.max(abs(delta_f))),
        median_abs_delta_epsilon_per_baryon_MeV=float(np.median(abs(delta_f))),
        max_abs_relative_epsilon_change=float(np.max(abs((reconstructed_epsilon-eps)/eps))),
        max_abs_mu_Q3_minus_Euler_MeV=float(np.max(abs(mu-raw_euler))),
        max_abs_mu_Q3_minus_corrected_MeV=float(np.max(abs(mu-reconstructed_mu))),
        max_abs_raw_Euler_minus_corrected_MeV=float(np.max(abs(raw_euler-reconstructed_mu))),
        max_cs2_over_c2=float(validation["max_cs2_over_c2"]),
        min_cs2_over_c2=float(validation["min_cs2_over_c2"]),
        valid_pressure_primary=validation["status"]=="PASS",
        n_points=len(n),
        note="caloric correction is a changed physical EOS, never represented as immutable source")


def source_model_row_audit(original, calculated):
    """Compare actual model RHS evaluations at EXACT source n, not regridded.

    No extrapolation; caller must select one homogeneous source phase and
    supply independent e,p,mu_B for each identical density.
    """
    n,P,eps,mu=_original_rows(original["n"],original["P"],
                              original["epsilon"],original["mu_Q3"])
    if len(calculated)!=len(n):
        raise ValueError("source and model rows must have identical cardinality")
    rows=[];rp=[];re=[];dmu=[]
    for i,model in enumerate(calculated):
        try:
            pm=float(model["P"]); em=float(model["epsilon"]); mm=float(model["mu_B"])
        except (KeyError,TypeError,ValueError) as ex:
            raise ValueError("model must supply P,epsilon,mu_B") from ex
        if not all(np.isfinite([pm,em,mm])) or pm<=0 or em<=0 or mm<=0:
            raise ValueError("model state must be positive finite at the source n")
        if "n" in model and not np.isclose(model["n"],n[i],rtol=0,atol=1e-14):
            raise ValueError("model evaluated at wrong density")
        dp=(pm-P[i])/P[i]; de=(em-eps[i])/eps[i]; dm=mm-mu[i]
        rp.append(float(dp));re.append(float(de));dmu.append(float(dm))
        rows.append(dict(density_fm3=float(n[i]),
            pressure_provider=float(P[i]),pressure_independent=float(pm),
            epsilon_provider=float(eps[i]),epsilon_independent=float(em),
            mu_provider_Q3_MeV=float(mu[i]),
            mu_provider_Euler_MeV=float((eps[i]+P[i])/n[i]),
            mu_independent_MeV=float(mm),pressure_relative_difference=float(dp),
            epsilon_relative_difference=float(de),
            model_mu_minus_provider_q3_MeV=float(dm)))
    return dict(n_evaluated=len(n),density_support_fm3=[float(n[0]),float(n[-1])],
        pressure_abs_rel_max=float(max(abs(np.asarray(rp)))),
        pressure_abs_rel_median=float(np.median(abs(np.asarray(rp)))),
        epsilon_abs_rel_max=float(max(abs(np.asarray(re)))),
        epsilon_abs_rel_median=float(np.median(abs(np.asarray(re)))),
        mu_model_minus_provider_q3_abs_max_MeV=float(max(abs(np.asarray(dmu)))),
        mu_provider_q3_minus_Euler_abs_max_MeV=float(np.max(abs(mu-(eps+P)/n))),
        source_identical=False,
        source_identical_status="NOT_PROVEN_SAME_FUNCTIONAL_MODEL_MUST_CHECK_NUMERIC_RESIDUALS",
        input_grid_policy="exact original homogeneous phase6 source n_B; no interpolation",
        samples=rows)
