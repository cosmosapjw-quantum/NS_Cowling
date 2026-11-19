"""Accuracy-gated sharp relativistic Cowling with independent direct-EOS oracle.

The Fortran DVODE integrator uses pretabulated PCHIP background coefficients.
A tiny ODE residual alone does not certify agreement with the original EOS.
This opt-in interface compares native resolutions and an independent Python
DOP853 calculation that evaluates each source phase directly.

Metric (-,+,+,+), frequencies Hz, n fm^-3, R km, c retained.
"""
from __future__ import annotations
import math
import numpy as np
from .stellar_cowling import CowlingSettings, _nodes
from .stellar_mode_continuation import eigenfunction_overlap
from .stellar_sharp import (SharpContractEnthalpyEOS,
    solve_sharp_background, solve_mode as python_sharp_mode)
from .stellar_sharp_native import (
    NativeSharpCowlingDVODE, solve_sharp_background_native)


def classify_cowling_nodes(w_nodes, v_nodes):
    if (not isinstance(w_nodes,(int,np.integer)) or
        not isinstance(v_nodes,(int,np.integer)) or
        isinstance(w_nodes,bool) or isinstance(v_nodes,bool) or
        min(w_nodes,v_nodes)<0):
        raise ValueError("nonnegative integer W/V node counts required")
    if w_nodes==v_nodes==0:
        return "STRICT_F"
    if w_nodes==v_nodes==1:
        return "STRICT_P1"
    if (w_nodes,v_nodes)==(1,2):
        return "W1_V2_CLASSIFICATION_UNRESOLVED"
    return "UNCLASSIFIED"


def count_reference_nodes_by_phase(profile,transitions):
    """Partition the source-evaluated DOP853 eigenfunction on each side.

    The two values at a sharp jump are not adjacent in a physical phase
    and their sign difference is not an eigenfunction zero.
    """
    if not isinstance(profile,list) or len(profile)<4:
        raise ValueError("profile must have at least four values")
    heights=sorted((float(t["h"]) for t in transitions),reverse=True)
    if any(not math.isfinite(h) for h in heights):
        raise ValueError("nonfinite interface")
    phases=[[] for _ in range(len(heights)+1)]
    for row in profile:
        h=float(row["h"])
        side=row["side"]
        w,v=float(row["W"]),float(row["V_km_inverse"])
        if (not all(math.isfinite(q) for q in (h,w,v)) or
            side not in ("high","low")):
            raise ValueError("invalid reference profile")
        idx=0
        for t in heights:
            tol=1e-12*max(1.,abs(t))
            if h<t-tol or (abs(h-t)<=tol and side=="low"):
                idx+=1
        phases[idx].append((w,v))
    if any(len(p)<2 for p in phases):
        raise ValueError("insufficient points on a physical phase")
    return dict(W_nodes_per_phase=int(sum(_nodes([v[0] for v in p]) for p in phases)),
        V_nodes_per_phase=int(sum(_nodes([v[1] for v in p]) for p in phases)),
        n_phases=len(phases),one_sided_at_every_sharp_interface=True,
        across_interface_zero_counts_included=False)


def select_native_sample_convergence(attempts,*,absolute_tolerance_hz=.03,
                                     max_join=1e-5,minimum_shape_overlap=.98):
    """Candidate native resolution test; source EOS agreement is separate."""
    if (not math.isfinite(absolute_tolerance_hz) or absolute_tolerance_hz<=0
        or not math.isfinite(max_join) or max_join<=0
        or not 0<minimum_shape_overlap<=1):
        raise ValueError("invalid numerical admission gates")
    passed=[]
    for ix,q in enumerate(attempts):
        if q.get("status")!="ADMITTED":
            continue
        freq=float(q.get("frequency_hz",float("nan")))
        join=float(q.get("relative_join_mismatch",float("nan")))
        if math.isfinite(freq) and freq>0 and math.isfinite(join) and join<=max_join:
            passed.append((ix,q))
    report=dict(status="UNRESOLVED",reason="INSUFFICIENT_ADMITTED_MODES",
                native_converged=False,comparison=None)
    if len(passed)<2:
        return report
    (left_index,left),(right_index,right)=passed[-2:]
    if right_index!=len(attempts)-1 or left_index!=right_index-1:
        report["reason"]="NONCONSECUTIVE_RESOLUTION_VALIDATION"
        return report
    if ((left["W_nodes_per_phase"],left["V_nodes_per_phase"])!=
        (right["W_nodes_per_phase"],right["V_nodes_per_phase"])):
        report["reason"]="NODE_PAIR_CHANGED"
        return report
    difference=abs(float(right["frequency_hz"])-float(left["frequency_hz"]))
    report["comparison"]=dict(samples=[left["samples"],right["samples"]],
                              difference_hz=float(difference))
    if difference>absolute_tolerance_hz:
        report["reason"]="COEFFICIENT_SAMPLING_NOT_CONVERGED"
        return report
    try:
        overlap=eigenfunction_overlap(left,right)["absolute_overlap"]
    except (ValueError,KeyError,TypeError) as exc:
        report["reason"]="NO_VALID_PHASEWISE_EIGENFUNCTION_OVERLAP"
        report["overlap_error"]=str(exc)
        return report
    report["comparison"]["eigenfunction_overlap"]=float(overlap)
    if overlap<minimum_shape_overlap:
        report["reason"]="EIGENFUNCTION_SHAPE_DISAGREES"
        return report
    report.update(status="NATIVE_SAMPLED_CANDIDATE_CONVERGED",
                  reason="ALL_NATIVE_SAMPLE_GATES_PASS",native_converged=True)
    return report


def solve_certified_sharp_mode(eos,central_density_fm3,bracket_hz,conversion,*,
    expected,settings=None,native_samples=(4097,8193,16385),
    absolute_tolerance_hz=.03,profile_points=600,max_join=1e-5,
    reference_join_tolerance=1e-6):
    """Try native samples; independently verify using original EOS DOP853.

    Failed native numerical gates never become PASS by simply lowering the
    frequency root tolerance. The reference-only result preserves the failure
    and scientific node-identity gate rather than silently relabelling p1.
    """
    if not isinstance(eos,SharpContractEnthalpyEOS) or not eos.transitions:
        raise TypeError("actual sharp EOS required")
    if conversion not in ("slow","rapid") or expected not in ("f","p1"):
        raise ValueError("unknown conversion or expected mode")
    if not isinstance(native_samples,(list,tuple)) or not native_samples or any(
            not isinstance(v,int) or isinstance(v,bool) or v<129
            for v in native_samples) or any(a>=b for a,b in
                                     zip(native_samples,native_samples[1:])):
        raise ValueError("strict increasing integer sample sizes >=129 required")
    if (not math.isfinite(absolute_tolerance_hz) or absolute_tolerance_hz<=0 or
        not math.isfinite(max_join) or max_join<=0 or
        not math.isfinite(reference_join_tolerance) or reference_join_tolerance<=0):
        raise ValueError("positive finite accuracy tolerances required")
    if not isinstance(profile_points,int) or profile_points<100:
        raise ValueError(">=100 profile points required")
    if not math.isfinite(central_density_fm3) or central_density_fm3<=0:
        raise ValueError("positive finite central density required")
    if len(bracket_hz)!=2:
        raise ValueError("one frequency bracket required")
    lo,hi=map(float,bracket_hz)
    if not 0<lo<hi or not all(math.isfinite(q) for q in (lo,hi)):
        raise ValueError("positive finite ascending frequency bracket required")
    if settings is None:
        settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,
            atol=1e-12,max_step_h=.002)
    if not isinstance(settings,CowlingSettings):
        raise TypeError("CowlingSettings required")
    native_attempts=[]
    internal_modes=[]
    native_masses=[]
    for samples in native_samples:
        receipt=dict(samples=int(samples),status="REJECTED")
        try:
            bg=solve_sharp_background_native(eos,central_density_fm3,
                 surface_h=settings.surface_h,rtol=3e-10,atol=3e-11,
                 max_step_h=.004,native_samples=samples)
            native_masses.append(bg.mass_msun)
            native=NativeSharpCowlingDVODE(bg,settings,samples=samples)
            mode=native.solve_mode((lo,hi),conversion=conversion,
                                   profile_points=profile_points)
            receipt.update(status="ADMITTED",
               frequency_hz=float(mode["frequency_hz"]),
               relative_join_mismatch=float(mode["relative_join_mismatch"]),
               matching_residual=float(mode["matching_residual"]),
               W_nodes_per_phase=int(mode["W_nodes_per_phase"]),
               V_nodes_per_phase=int(mode["V_nodes_per_phase"]),
               profile=mode["profile"])
            internal_modes.append(receipt)
        except (RuntimeError,ValueError,FloatingPointError) as exc:
            receipt["reason"]=type(exc).__name__+": "+str(exc)[:300]
        native_attempts.append({k:v for k,v in receipt.items() if k!="profile"})
    candidate=select_native_sample_convergence(internal_modes,
                         absolute_tolerance_hz=absolute_tolerance_hz,max_join=max_join)
    best_native=internal_modes[-1] if candidate["native_converged"] else None
    try:
        pybg=solve_sharp_background(eos,central_density_fm3,
                 surface_h=settings.surface_h,rtol=3e-10,atol=3e-12,
                 max_step_h=.003)
        oracle=python_sharp_mode(pybg,(lo,hi),settings,
                         conversion=conversion,profile_points=profile_points)
        residual=float(oracle["relative_join_mismatch"])
        if not math.isfinite(residual) or residual>reference_join_tolerance:
            raise RuntimeError("direct-EOS DOP853 oracle failed join admission")
        nodal=count_reference_nodes_by_phase(oracle["profile"],eos.transitions)
        label=classify_cowling_nodes(nodal["W_nodes_per_phase"],
                                     nodal["V_nodes_per_phase"])
        if native_masses and min(abs(pybg.mass_msun-m) for m in native_masses)>3e-6:
            raise RuntimeError("native/Python TOV mass disagreement >3e-6 Msun")
        pyreceipt=dict(status="DIRECT_EOS_REFERENCE_ACCEPTED",
            frequency_hz=float(oracle["frequency_hz"]),
            relative_join_mismatch=residual,
            matching_residual=float(oracle["matching_residual"]),
            mass_MSun=float(pybg.mass_msun),radius_km=float(pybg.radius_km),
            mass_deviation_from_native_Msun=(min(abs(pybg.mass_msun-m)
                for m in native_masses) if native_masses else None),
            reference_nodes=nodal,
            backend="PYTHON_DOP853_DIRECT_EOS")
    except (RuntimeError,ValueError,FloatingPointError) as exc:
        return dict(status="UNRESOLVED_PYTHON_REFERENCE_FAILED",
            frequency_hz=None,selected_backend=None,
            reference_error=type(exc).__name__+": "+str(exc)[:300],
            native_attempts=native_attempts,native_convergence=candidate,
            source_eos_replaced=False)
    if best_native is not None:
        diff=float(best_native["frequency_hz"]-oracle["frequency_hz"])
        native_label=classify_cowling_nodes(best_native["W_nodes_per_phase"],
                                            best_native["V_nodes_per_phase"])
    else:
        diff=None
        native_label=None
    strict_expected=((expected=="p1" and label=="STRICT_P1") or
                     (expected=="f" and label=="STRICT_F"))
    acceptable=(diff is not None and abs(diff)<=absolute_tolerance_hz
                and native_label==label and strict_expected)
    status=("CERTIFIED_NATIVE_REFERENCE_AGREEMENT" if acceptable
            else "REFERENCE_ONLY_NATIVE_NOT_CERTIFIED")
    selected_backend=("FORTRAN_DVODE_PCHIP_CROSSCHECKED" if acceptable
                      else "PYTHON_DOP853_DIRECT_EOS")
    f=(float(best_native["frequency_hz"]) if acceptable
       else float(oracle["frequency_hz"]))
    if (expected=="p1" and label=="STRICT_P1") or (expected=="f" and label=="STRICT_F"):
        identification="STRICT_CONFIRMED"
    elif expected=="p1" and label=="W1_V2_CLASSIFICATION_UNRESOLVED":
        identification="CONDITIONAL_NUMERIC_FAMILY_NOT_STRICT_P1"
    else:
        identification="UNRESOLVED_INCONSISTENT_EXPECTED_NODES"
    return dict(status=status,frequency_hz=f,
        selected_backend=selected_backend,
        central_density_fm3=float(central_density_fm3),
        conversion=conversion,expected=expected,
        bracket_hz=[lo,hi],node_gate=label,
        expected_mode_identification=identification,
        native_attempts=native_attempts,native_convergence=candidate,
        python_reference=pyreceipt,
        native_minus_python_reference_hz=diff,
        native_reference_difference_hz=(abs(diff) if diff is not None else None),
        absolute_tolerance_hz=float(absolute_tolerance_hz),
        reference_does_not_establish_global_mode_family=True,
        reference_physics="equilibrium barotropic relativistic Cowling",
        source_eos_replaced=False)
