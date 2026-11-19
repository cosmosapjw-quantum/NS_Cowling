"""Original-core RMF reconciliation without fitting or changing source EOS.

Cold one-charge, beta-equilibrated thermodynamics (n fm^-3, P and epsilon
MeV/fm^3, mu in MeV). No phase-boundary smoothing or extrapolation.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import CubicSpline


def _array(values,name):
    a=np.asarray(values,dtype=float)
    if a.ndim!=1 or len(a)<2 or not np.all(np.isfinite(a)):
        raise ValueError(name+" must be one-dimensional finite with >=2 entries")
    return a


def cold_pressure_from_energy(n,epsilon,P=None,mu_Q3=None):
    """Source-cold first-law diagnostic, not a replacement EOS.

    d(epsilon/n)/dn = P/n^2. The cubic spline is a sampled derivative,
    not certified physical pressure; discard two endpoints in error summaries.
    """
    n=_array(n,"density")
    epsilon=_array(epsilon,"epsilon")
    if len(n)<4:
        raise ValueError("energy derivative needs at least four original source knots")
    if len(n)!=len(epsilon) or np.any(n<=0) or np.any(epsilon<=0) or np.any(np.diff(n)<=0):
        raise ValueError("positive same-length epsilon and strictly increasing density required")
    F=epsilon/n
    pressure=n*n*CubicSpline(n,F,bc_type="not-a-knot")(n,nu=1)
    if not np.all(np.isfinite(pressure)):
        raise ValueError("source energy derivative not finite")
    result=dict(n=n.copy(),P_from_energy=pressure,
        source_euler_mu=(epsilon+pressure)/n,
        max_inner_abs_P_error=None,
        median_inner_abs_rel_P_error=None,
        max_provider_Q3_Euler_discrepancy_MeV=None,
        interpolation="CubicSpline d(epsilon/n)/dn, exclude outermost two samples for comparison",
        status="SAMPLED_ENERGY_DERIVATIVE_NOT_A_THERMODYNAMIC_REPAIR")
    if P is not None:
        P=_array(P,"pressure")
        if len(P)!=len(n) or np.any(P<=0):
            raise ValueError("positive source pressure matching density required")
        delta=P[2:-2]-pressure[2:-2]
        result["max_inner_abs_P_error"]=float(np.max(np.abs(delta)))
        result["median_inner_abs_rel_P_error"]=float(np.median(np.abs(delta)/P[2:-2]))
        result["median_signed_inner_PminusPderived_MeV_fm3"]=float(np.median(delta))
    if mu_Q3 is not None:
        mu=_array(mu_Q3,"provider chemical potential")
        if len(mu)!=len(n):raise ValueError("chemical potential length mismatch")
        result["max_provider_Q3_Euler_discrepancy_MeV"]=float(
            np.max(np.abs(mu-(epsilon+(P if P is not None else pressure))/n)))
    return result


def rest_mass_shift_pressure_invariant(n,epsilon,P,*,delta_mass_MeV):
    """Derivation: epsilon' = epsilon+n delta_m implies dF'/dn=dF/dn."""
    n=_array(n,"density");epsilon=_array(epsilon,"epsilon");P=_array(P,"pressure")
    if len(n)!=len(epsilon) or len(n)!=len(P):
        raise ValueError("density and thermodynamic arrays must align")
    mass=float(delta_mass_MeV)
    if not np.isfinite(mass):raise ValueError("mass shift finite required")
    return dict(epsilon_after_shift=epsilon+n*mass,P_after_shift=P.copy(),
         delta_energy_per_baryon_MeV=mass,max_abs_pressure_shift=0.,
         thermodynamic_identity="P=n^2 d(epsilon/n)/dn unchanged under constant mass offset")


def compare_phase_homogeneous_states(original,models,*,rounding_half_width=5.01e-6):
    """Same-density raw CompOSE source vs independently solved RMF β state.

    Charge fractions are percentage-free *number fractions*; provider printed
    5 decimal places in eos.compo. A rounded fraction match is NOT proof that
    the bulk source is identical.
    """
    fields=("n","P","epsilon","mu_Q3","Yp","Ye","Ymu")
    arrays={k:_array(original[k],k) for k in fields}
    n=arrays["n"];P=arrays["P"];eps=arrays["epsilon"]
    if (any(len(arrays[k])!=len(n) for k in fields)
        or np.any(n<=0) or np.any(P<=0) or np.any(eps<=0)
        or np.any(np.diff(n)<=0)):
        raise ValueError("original homogeneous rows must be one monotone physical phase")
    if not np.isfinite(rounding_half_width) or rounding_half_width<0:
        raise ValueError("rounding tolerance must be nonnegative finite")
    if len(models)!=len(n):
        raise ValueError("independently evaluated model rows must align at exact source density")
    result=[]
    for i,row in enumerate(models):
        if "n" not in row or abs(float(row["n"])-n[i])>1e-14:
            raise ValueError("independent RMF sample evaluated at different density")
        x={k:float(row[k]) for k in ("P","epsilon","mu_B","Yp","Ye","Ymu")}
        if not all(np.isfinite(list(x.values()))) or x["P"]<=0 or x["epsilon"]<=0:
            raise ValueError("independent RMF state invalid")
        y={k:float(x[k]-arrays[k][i]) for k in ("Yp","Ye","Ymu")}
        result.append(dict(density_fm3=float(n[i]),
            pressure_provider=float(P[i]),pressure_model=x["P"],
            epsilon_provider=float(eps[i]),epsilon_model=x["epsilon"],
            mu_Q3_provider=float(arrays["mu_Q3"][i]),
            mu_Euler_provider=float((eps[i]+P[i])/n[i]),
            mu_B_model=x["mu_B"],fraction_deltas=y,
            pressure_delta=float(x["P"]-P[i]),
            energy_per_baryon_delta_MeV=float((x["epsilon"]-eps[i])/n[i]),
            pressure_relative=float((x["P"]-P[i])/P[i])))
    gap_p=np.array([x["pressure_relative"] for x in result])
    gap_y=np.array([z for x in result for z in x["fraction_deltas"].values()])
    gap_mu=np.array([x["mu_B_model"]-x["mu_Q3_provider"] for x in result])
    return dict(n_evaluated=len(result),
        max_abs_composition_fraction_difference=float(max(abs(gap_y))),
        composition_compatible_with_printed_rounding=bool(
            max(abs(gap_y))<=rounding_half_width),
        provider_fraction_half_decimal_resolution=float(rounding_half_width),
        pressure_max_abs_rel=float(max(abs(gap_p))),
        pressure_median_abs_rel=float(np.median(abs(gap_p))),
        bulk_pressure_source_identical=bool(max(abs(gap_p))<1e-6),
        energy_per_baryon_model_minus_provider_range_MeV=[
            float(min(x["energy_per_baryon_delta_MeV"] for x in result)),
            float(max(x["energy_per_baryon_delta_MeV"] for x in result))],
        mu_model_minus_Q3_max_abs_MeV=float(max(abs(gap_mu))),
        source_euler_minus_Q3_max_abs_MeV=float(max(abs(
            x["mu_Euler_provider"]-x["mu_Q3_provider"]) for x in result)),
        homogeneous_rows=result,
        conclusion="composition parity and bulk pressure identity are independent tests")



def charge_mu_consistency(n,Y_e,Q4,Q5,*,mn,
                          source_fraction_rounding_half_width=5e-6):
    """CompOSE mu_e=(Q5-Q4)*mn compared with free electron Fermi energy.

    CompOSE's Q4 is a *charge* chemical potential, NOT the electron
    chemical potential itself. For cold beta-equilibrated lepton number
    unconstrained, Q5=0 and mu_e=-mn Q4. Test all positive-electron
    phase rows without changing any original provider metadata.
    """
    from .cold_leptons import ELECTRON_MASS,free_fermion
    n=_array(n,"baryon density")
    ye=_array(Y_e,"electron fraction")
    q4=_array(Q4,"Q4 scaled charge potential")
    q5=_array(Q5,"Q5 scaled lepton potential")
    if len(set(map(len,(n,ye,q4,q5))))!=1:
        raise ValueError("CompOSE chemical potential and composition rows must align")
    mass=float(mn)
    width=float(source_fraction_rounding_half_width)
    if not np.isfinite(mass) or mass<=0 or not np.isfinite(width) or width<0:
        raise ValueError("positive reference neutron mass and nonnegative rounding width required")
    if np.any(n<=0) or np.any(ye<=width) or np.any(ye>1):
        raise ValueError("electron Fermi consistency needs positive populated electron fractions")
    rows=[]
    for nb,y,qa,qb in zip(n,ye,q4,q5):
        f=lambda yv:float(free_fermion(float(nb*yv),ELECTRON_MASS)[2])
        mu_ideal=f(y)
        lbound=f(y-width);ubound=f(y+width)
        reported=(qb-qa)*mass
        # allow <=0.00005 MeV for finite-Q4 printed precision and source
        # baryon-density rounding; this is irrelevant versus 1-2 MeV.
        margin=5e-5
        passes=lbound-margin<=reported<=ubound+margin
        rows.append(dict(n_fm3=float(nb),Y_e=float(y),
            mu_e_from_source_Q4_Q5_MeV=float(reported),
            mu_e_from_source_n_e_ideal_MeV=float(mu_ideal),
            source_minus_fermion_mu_e_MeV=float(reported-mu_ideal),
            electron_fraction_rounding_mu_e_interval_MeV=[float(lbound),float(ubound)],
            rounding_compatible=bool(passes),
            source_lepton_mu_Q5_MeV=float(qb*mass)))
    diffs=np.array([z["source_minus_fermion_mu_e_MeV"] for z in rows])
    return dict(n_rows=len(rows),
        max_abs_mu_e_discrepancy_MeV=float(np.max(abs(diffs))),
        median_abs_mu_e_discrepancy_MeV=float(np.median(abs(diffs))),
        signed_mu_e_discrepancy_range_MeV=[float(np.min(diffs)),float(np.max(diffs))],
        within_rounding_rows=sum(int(row["rounding_compatible"]) for row in rows),
        all_rows_rounding_consistent=all(row["rounding_compatible"] for row in rows),
        CompOSE_thermodynamic_definition="mu_e = mu_le - mu_q = mn*(Q5-Q4); Q5=0 for unconstrained cold beta equilibrium",
        test_assumption="electrons form free zero-T relativistic Fermi gas, source Ye electron number fraction",
        original_source_metadata_not_modified=True,
        rows=rows)


def decompose_pressure_difference(n,provider_P,provider_epsilon,
                                  model_P,model_epsilon):
    """Decompose same-density pressure mismatch into energy and identity terms.

    Let R=P-n²d(epsilon/n)/dn. At the same baryon-density knots,
    P_model-P_provider =
       n²d[(epsilon_model-epsilon_provider)/n]/dn + R_model-R_provider.

    Uses the SAME cubic spline derivative for all terms so its algebraic
    balance is an implementation invariant, not a physics source-equivalence
    theorem. Interior-only RMS figures reduce endpoint derivative artifacts.
    Do not apply across sharp phase transitions or density inversions.
    """
    density=_array(n,"density")
    p0=_array(provider_P,"provider pressure")
    e0=_array(provider_epsilon,"provider energy")
    p1=_array(model_P,"independent pressure")
    e1=_array(model_epsilon,"independent energy")
    if (len(density)<6 or any(len(v)!=len(density) for v in (p0,e0,p1,e1))
        or np.any(density<=0) or np.any(np.diff(density)<=0)):
        raise ValueError("strictly increasing positive density and >=6 matched source rows required")
    if any(np.any(v<=0) for v in (p0,e0,p1,e1)):
        raise ValueError("positive original and modeled pressure and energy required")
    d0=density*density*CubicSpline(density,e0/density)(density,nu=1)
    d1=density*density*CubicSpline(density,e1/density)(density,nu=1)
    delta_eps=(e1-e0)/density
    ddelta=density*density*CubicSpline(density,delta_eps)(density,nu=1)
    observed=p1-p0
    R0=p0-d0
    R1=p1-d1
    prediction=ddelta+R1-R0
    residual=observed-prediction
    inside=slice(2,-2)
    rms=lambda v:float(np.sqrt(np.mean(np.asarray(v[inside])**2)))
    denom=rms(observed)
    if denom==0:denom=np.finfo(float).tiny
    return dict(
        density_source=density.copy(),
        pressure_observed_difference=observed,
        pressure_from_energy_difference=ddelta,
        provider_identity_residual=R0,
        model_identity_residual=R1,
        algebraic_balance_error=residual,
        max_abs_balance_error_MeV_fm3=float(np.max(abs(residual))),
        rms_model_source_pressure_difference_MeV_fm3=rms(observed),
        rms_energy_density_functional_term_MeV_fm3=rms(ddelta),
        rms_provider_identity_residual_MeV_fm3=rms(R0),
        rms_model_identity_residual_MeV_fm3=rms(R1),
        rms_energy_term_over_observed=float(rms(ddelta)/denom),
        rms_provider_identity_over_observed=float(rms(R0)/denom),
        rms_model_identity_over_observed=float(rms(R1)/denom),
        source_eos_unchanged=True,
        caveat="CubicSpline sampled density derivative; algebraic identity is not a unique diagnosis of provider-generation implementation")
