#!/usr/bin/env python3
"""Deterministic mock experiment for fixed-gravitational-mass mode response.

This is a mathematical stress test of the implicit-function decomposition,
not a neutron-star EOS or a prediction for physical mode frequencies.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from numpy.polynomial import Polynomial

CORES=("BBB_BHF2","RG_SLY4","GM_GM1")
CRUSTS=tuple(f"crust_{i:02d}" for i in range(1,14))
LAMBDAS=tuple(float(x) for x in np.linspace(0.0,1.0,6))
MASSES=(0.8,1.0,1.2,1.33,1.4,1.6,1.8,2.0)
METHODS=tuple(f"mock_bridge_{i}" for i in range(6))
M_MAX0={"BBB_BHF2":2.30,"RG_SLY4":2.16,"GM_GM1":2.025}
F_BASE={"BBB_BHF2":2.505,"RG_SLY4":2.370,"GM_GM1":1.920}
P1_BASE={"BBB_BHF2":7.008,"RG_SLY4":6.718,"GM_GM1":6.009}


def pair_parameters(core_index: int, crust_index: int, lam: float) -> tuple[float,...]:
    """A smooth synthetic TOV surface with an adjustable first mass fold."""
    phase=2.0*math.pi*crust_index/13.0
    pair_shift=0.004*math.sin(phase+0.7*core_index)
    mmax_shift=0.050*(lam-0.5)+0.009*math.sin(
        2.0*math.pi*lam+0.4*core_index+phase)
    mmax=M_MAX0[CORES[core_index]]+pair_shift+mmax_shift
    curvature=1.0*(1.0+0.035*core_index+0.02*math.cos(phase)
                   +0.06*lam+0.01*lam*math.sin(phase))
    xturn=1.40+0.015*core_index+0.008*math.sin(phase)+0.04*(lam-0.5)
    dmmax=0.050+0.009*2.0*math.pi*math.cos(
        2.0*math.pi*lam+0.4*core_index+phase)
    dcurvature=1.0*(0.06+0.01*math.sin(phase))
    dxturn=0.04
    return mmax,curvature,xturn,dmmax,dcurvature,dxturn


def stable_logpc_at_mass(mass_msun: float, core_index: int,
                         crust_index: int, lam: float) -> float | None:
    """Analytic inverse on the initial rising-mass branch only."""
    mmax,curvature,xturn,*_=pair_parameters(core_index,crust_index,lam)
    if not math.isfinite(mass_msun) or mass_msun<=0 or mass_msun>=mmax:
        return None
    return xturn-math.sqrt(2.0*(mmax-mass_msun)/curvature)


def tov_partials(x: float, lam: float, core_index: int,
                 crust_index: int) -> tuple[float,float]:
    """Return (M_lambda, M_x) for M=Mmax-A(x-xturn)^2/2."""
    mmax,A,xturn,dmmax,dA,dxturn=pair_parameters(core_index,crust_index,lam)
    dx=x-xturn
    m_lambda=dmmax-0.5*dA*dx*dx+A*dx*dxturn
    m_x=-A*dx
    return m_lambda,m_x


def mode_value(x: float, lam: float, core_index: int, crust_index: int,
               mode: str, profile: str="fold") -> float:
    core=CORES[core_index]
    phase=math.sin(2.0*math.pi*crust_index/13.0+0.3*core_index)
    if mode=="f":
        base=F_BASE[core]
        slope_x=0.22*(1.0+0.10*core_index+0.04*phase)
        direct=0.008+0.004*x+0.001*phase
    elif mode=="p1":
        base=P1_BASE[core]
        slope_x=0.75*(1.0+0.12*core_index+0.05*phase)
        direct=0.045+0.018*x+0.006*phase
    else:
        raise ValueError("mode must be 'f' or 'p1'")
    if profile=="fold":
        return (base+slope_x*(x-1.0)+direct*(lam-0.5)
                +0.002*(lam-0.5)*(x-1.0))
    if profile=="turning_point_stationary":
        xturn=pair_parameters(core_index,crust_index,lam)[2]
        dx=x-xturn
        amplitude=0.10 if mode=="f" else 0.40
        return (base+0.015*(lam-0.5)+amplitude*dx*dx
                +0.003*phase*(lam-0.5))
    raise ValueError("unknown synthetic mode profile")


def mode_partials(x: float, lam: float, core_index: int,
                  crust_index: int, mode: str,
                  profile: str="fold") -> tuple[float,float]:
    """Return (f_lambda, f_x), with the other coordinate held fixed."""
    phase=math.sin(2.0*math.pi*crust_index/13.0+0.3*core_index)
    if profile=="fold":
        if mode=="f":
            slope_x=0.22*(1.0+0.10*core_index+0.04*phase)
            direct=0.008+0.004*x+0.001*phase
            direct_x=0.004
        elif mode=="p1":
            slope_x=0.75*(1.0+0.12*core_index+0.05*phase)
            direct=0.045+0.018*x+0.006*phase
            direct_x=0.018
        else:
            raise ValueError("mode must be 'f' or 'p1'")
        f_x=slope_x+direct_x*(lam-0.5)+0.002*(lam-0.5)
        f_lambda=direct+0.002*(x-1.0)
        return f_lambda,f_x
    if profile=="turning_point_stationary":
        xturn=pair_parameters(core_index,crust_index,lam)[2]
        dx=x-xturn
        amplitude=0.10 if mode=="f" else 0.40
        _,_,_,_,_,dxturn=pair_parameters(core_index,crust_index,lam)
        f_x=2.0*amplitude*dx
        f_lambda=0.015-2.0*amplitude*dx*dxturn+0.003*phase
        return f_lambda,f_x
    raise ValueError("unknown synthetic mode profile")


def fixed_mass_response(mass_msun: float, lam: float, core_index: int,
                        crust_index: int, mode: str,
                        profile: str="fold") -> dict | None:
    """Exact IFT decomposition into fixed-Pc and structure-remapping pieces."""
    x=stable_logpc_at_mass(mass_msun,core_index,crust_index,lam)
    if x is None:
        return None
    f_lambda,f_x=mode_partials(x,lam,core_index,crust_index,mode,profile)
    m_lambda,m_x=tov_partials(x,lam,core_index,crust_index)
    if abs(m_x)<1.0e-14:
        raise ZeroDivisionError("fixed-mass IFT derivative is singular at the fold")
    structure=-f_x*m_lambda/m_x
    return {
        "x_logpc":x,
        "f_lambda_fixed_pc":f_lambda,
        "structure_remapping":structure,
        "f_lambda_fixed_mass":f_lambda+structure,
        "M_lambda":m_lambda,
        "M_x":m_x,
    }


def thermodynamic_q_path(npoints: int=401) -> dict:
    """Build a convex path q_lambda(mu) with common endpoint P,n,epsilon.

    With n=mu*q and epsilon=mu*n-P, dP/dmu=n and
    cs^2=q/(q+mu*q_mu). The perturbation polynomial vanishes at both
    endpoints and has zero weighted integral int mu*delta-q dmu, so all
    six lambda slices share both endpoint pressures.
    """
    if npoints<51:
        raise ValueError("npoints must be at least 51")
    t=np.linspace(0.0,1.0,npoints)
    beta=0.2
    T=Polynomial([0.0,1.0])
    Z=Polynomial([1.0,beta])
    base_shape=T*(1.0-T)
    numerator=(Z*T*T*(1.0-T)).integ()
    denominator=(Z*T*(1.0-T)).integ()
    c=(numerator(1.0)-numerator(0.0))/(denominator(1.0)-denominator(0.0))
    perturb=T*(1.0-T)*(T-c)
    q0=Polynomial([0.8,0.3])
    z=1.0+beta*t
    endpoint_pressure=[]
    rows=[]
    q_min=math.inf
    cs_min=math.inf
    cs_max=-math.inf
    endpoint_eps=[]
    for method_index,lam in enumerate(LAMBDAS):
        qpoly=q0+0.1*lam*perturb
        dqdt=qpoly.deriv()
        pressure_primitive=(beta*Z*qpoly).integ()
        p=pressure_primitive(t)-pressure_primitive(0.0)
        q=qpoly(t)
        q_mu=dqdt(t)/beta
        cs2=q/(q+z*q_mu)
        epsilon=z*z*q-p
        if not (np.all(q>0.0) and np.all(q_mu>=0.0)
                and np.all(cs2>0.0) and np.all(cs2<=1.0)):
            raise AssertionError("synthetic q-path lost thermodynamic stability or causality")
        q_min=min(q_min,float(np.min(q)))
        cs_min=min(cs_min,float(np.min(cs2)))
        cs_max=max(cs_max,float(np.max(cs2)))
        endpoint_pressure.append(float(p[-1]))
        endpoint_eps.append(float(epsilon[-1]))
        for j in range(0,npoints,10):
            rows.append({
                "lambda":lam,"method_index":method_index,"t":float(t[j]),
                "mu_over_mu1":float(z[j]),"q":float(q[j]),"n_over_mu1qref":float(z[j]*q[j]),
                "p_dimensionless":float(p[j]),"epsilon_dimensionless":float(epsilon[j]),
                "cs2":float(cs2[j]),
            })
    pspread=max(endpoint_pressure)-min(endpoint_pressure)
    espread=max(endpoint_eps)-min(endpoint_eps)
    if pspread>2.0e-13 or espread>2.0e-13:
        raise AssertionError("endpoint pressure/energy invariance failed along q path")
    return {
        "lambda_slices":len(LAMBDAS),"npoints_per_slice":npoints,
        "q_min":q_min,"cs2_min":cs_min,"cs2_max":cs_max,
        "endpoint_pressure_spread":pspread,"endpoint_energy_spread":espread,
        "weighted_zero_integral_constant_c":float(c),
        "interpretation":"synthetic admissible q-path; not a physical EOS table",
        "sample_rows":rows,
    }


def build_mock(output_dir: Path) -> dict:
    output_dir.mkdir(parents=True,exist_ok=True)
    pair_rows=[]
    for ci,core in enumerate(CORES):
        for si,crust in enumerate(CRUSTS):
            for mass in MASSES:
                samples={"f":[],"p1":[]}
                valid_methods=[]
                for mi,lam in enumerate(LAMBDAS):
                    x=stable_logpc_at_mass(mass,ci,si,lam)
                    if x is None:
                        continue
                    valid_methods.append(mi)
                    samples["f"].append(mode_value(x,lam,ci,si,"f"))
                    samples["p1"].append(mode_value(x,lam,ci,si,"p1"))
                complete=len(valid_methods)==len(LAMBDAS)
                row={
                    "core":core,"crust":crust,"mass_msun":mass,
                    "n_valid_methods":len(valid_methods),"pair_complete":complete,
                    "missing_method_indices":";".join(str(i) for i in range(6) if i not in valid_methods),
                    "f_range_percent":None,"p1_range_percent":None,
                    "f_mean_khz":None,"p1_mean_khz":None,
                }
                if complete:
                    for key in ("f","p1"):
                        avg=float(np.mean(samples[key]))
                        row[f"{key}_mean_khz"]=avg
                        row[f"{key}_range_percent"]=100.0*(max(samples[key])-min(samples[key]))/avg
                pair_rows.append(row)
    pair_path=output_dir/"paired_mock_ranges.csv"
    with pair_path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(pair_rows[0]))
        writer.writeheader();writer.writerows(pair_rows)

    mass_summary=[]
    for mass in MASSES:
        rows=[r for r in pair_rows if r["mass_msun"]==mass]
        complete=[r for r in rows if r["pair_complete"]]
        result={"mass_msun":mass,"n_pairs":len(rows),"n_complete_pairs":len(complete),
                "coverage_fraction":len(complete)/len(rows)}
        for key in ("f","p1"):
            values=[r[f"{key}_range_percent"] for r in complete]
            result[f"{key}_median_range_percent"]=float(np.median(values)) if values else None
            result[f"{key}_p90_range_percent"]=float(np.quantile(values,0.9)) if values else None
            result[f"{key}_max_range_percent"]=float(max(values)) if values else None
        mass_summary.append(result)
    summary_path=output_dir/"mass_coverage_summary.csv"
    with summary_path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(mass_summary[0]))
        writer.writeheader();writer.writerows(mass_summary)

    derivative_rows=[]
    max_fd_error=0.0
    for ci,core in enumerate(CORES):
        for si,crust in enumerate(CRUSTS):
            for mass in (1.0,1.33,1.8):
                lam=0.4
                for mode in ("f","p1"):
                    result=fixed_mass_response(mass,lam,ci,si,mode)
                    if result is None:
                        continue
                    h=1.0e-6
                    xp=stable_logpc_at_mass(mass,ci,si,lam+h)
                    xm=stable_logpc_at_mass(mass,ci,si,lam-h)
                    fp=mode_value(xp,lam+h,ci,si,mode)
                    fm=mode_value(xm,lam-h,ci,si,mode)
                    fd=(fp-fm)/(2.0*h)
                    error=abs(fd-result["f_lambda_fixed_mass"])/max(abs(fd),1.0e-12)
                    max_fd_error=max(max_fd_error,error)
                    derivative_rows.append({
                        "core":core,"crust":crust,"mass_msun":mass,"lambda":lam,"mode":mode,
                        **result,"finite_difference":fd,"relative_error":error,
                    })
    derivative_path=output_dir/"fixed_mass_derivative_checks.csv"
    with derivative_path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(derivative_rows[0]))
        writer.writeheader();writer.writerows(derivative_rows)

    # Approach the same generic fold in two different mode profiles:
    # nonzero f_x gives delta-M^(-1/2); f_x ~ x-x_turn cancels the divergence.
    gaps=np.logspace(-1.0,-4.0,12)
    ci,si,lam=2,6,0.5
    scaling_rows=[]
    fold_abs=[];stationary_abs=[]
    mmax=pair_parameters(ci,si,lam)[0]
    for gap in gaps:
        mass=mmax-float(gap)
        fold=fixed_mass_response(mass,lam,ci,si,"p1","fold")
        stationary=fixed_mass_response(mass,lam,ci,si,"p1","turning_point_stationary")
        fold_abs.append(abs(fold["f_lambda_fixed_mass"]))
        stationary_abs.append(abs(stationary["f_lambda_fixed_mass"]))
        scaling_rows.append({
            "mass_gap_to_mmax_msun":float(gap),
            "fold_abs_fixed_mass_response_khz_per_lambda":fold_abs[-1],
            "stationary_abs_fixed_mass_response_khz_per_lambda":stationary_abs[-1],
            "M_x":fold["M_x"],
            "x_logpc":fold["x_logpc"],
        })
    fold_slope=float(np.polyfit(np.log(gaps),np.log(fold_abs),1)[0])
    stationary_slope=float(np.polyfit(np.log(gaps),np.log(stationary_abs),1)[0])
    scaling_path=output_dir/"turning_point_scaling.csv"
    with scaling_path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(scaling_rows[0]))
        writer.writeheader();writer.writerows(scaling_rows)

    thermo=thermodynamic_q_path()
    thermo_path=output_dir/"thermodynamic_q_path_summary.json"
    thermo_compact={k:v for k,v in thermo.items() if k!="sample_rows"}
    thermo_compact["lambda_values"]=list(LAMBDAS)
    thermo_path.write_text(json.dumps(thermo_compact,indent=2,sort_keys=True)+"\n")
    thermo_rows=thermo["sample_rows"]
    thermo_csv=output_dir/"thermodynamic_q_path_samples.csv"
    with thermo_csv.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(thermo_rows[0]))
        writer.writeheader();writer.writerows(thermo_rows)

    manifest={
        "schema":"ns-cowling-fixed-mass-response-mock-v1",
        "status":"MOCK_ONLY_NOT_A_PHYSICAL_EOS_OR_SPECTRUM",
        "design":{"cores":list(CORES),"crust_count":len(CRUSTS),"crusts":list(CRUSTS),
                  "lambda_values":list(LAMBDAS),"mass_grid_msun":list(MASSES),
                  "crossed_pair_count":len(CORES)*len(CRUSTS),
                  "maximum_pair_method_cells":len(CORES)*len(CRUSTS)*len(LAMBDAS)*len(MASSES)},
        "results":{"pair_mass_cells":len(pair_rows),"derivative_checks":len(derivative_rows),
                   "max_ift_vs_central_fd_relative_error":max_fd_error,
                   "fold_loglog_slope":fold_slope,"stationary_mode_loglog_slope":stationary_slope,
                   "mass_summary":mass_summary,"thermodynamic_path":thermo_compact},
        "formula":"df/dlambda|M = f_lambda|x - f_x M_lambda/M_x, x=ln(Pc/P*)",
        "fold_model":"M=Mmax(lambda)-A(lambda)*(x-xturn(lambda))^2/2",
        "interpretation":{
            "fold_profile":"generic nonzero f_x produces an inverse-square-root response near Mmax",
            "stationary_profile":"f_x proportional to x-xturn cancels the fold divergence",
            "method_labels":"mock_bridge_0..5 are equally spaced samples of a synthetic admissible path, not the six physical NS_Cowling prescriptions",
            "coverage":"incomplete high-mass pairs are reported and excluded, never imputed",
        },
        "files":["paired_mock_ranges.csv","mass_coverage_summary.csv",
                 "fixed_mass_derivative_checks.csv","turning_point_scaling.csv",
                 "thermodynamic_q_path_summary.json","thermodynamic_q_path_samples.csv"],
    }
    manifest_path=output_dir/"mock_manifest.json"
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    if max_fd_error>1.0e-6:
        raise AssertionError(f"IFT finite-difference check failed: {max_fd_error}")
    if not (-0.65<fold_slope<-0.35):
        raise AssertionError(f"generic fold scaling did not approach -1/2: {fold_slope}")
    if abs(stationary_slope)>0.10:
        raise AssertionError(f"stationary-mode cancellation was not demonstrated: {stationary_slope}")
    return manifest


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=Path("dist/mass-response-mock"))
    args=parser.parse_args()
    result=build_mock(args.output_dir)
    print(json.dumps({
        "status":result["status"],
        "derivative_checks":result["results"]["derivative_checks"],
        "max_ift_vs_central_fd_relative_error":result["results"]["max_ift_vs_central_fd_relative_error"],
        "fold_loglog_slope":result["results"]["fold_loglog_slope"],
        "stationary_mode_loglog_slope":result["results"]["stationary_mode_loglog_slope"],
        "mass_coverage":[{"mass_msun":x["mass_msun"],"complete_pairs":x["n_complete_pairs"],
                          "f_median":x["f_median_range_percent"],"p1_median":x["p1_median_range_percent"]}
                         for x in result["results"]["mass_summary"]],
    },indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
