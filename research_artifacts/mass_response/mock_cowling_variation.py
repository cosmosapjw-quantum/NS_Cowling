#!/usr/bin/env python3
"""Synthetic first-variation check for the relativistic Cowling 2x2 shooting system.

This is an ODE sensitivity regression, not a neutron-star model or an identified
f/p1 eigenfrequency. The prescribed coefficient background is intentionally
synthetic; it verifies the tangent equations and residual derivative only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ELL=2
RADIUS=1.0
CENTER_RADIUS=1.0e-4
ROOT_BRACKET=(7.67,8.18)
CHI0=0.3
RTOL=2.0e-9


def _profile(r: float) -> tuple[float,float,float,float]:
    """Return nu, nu', chi-shape g, and its derivative on the toy background."""
    nu=0.5*r*r
    nu_prime=r
    g=r*r*(1.0-r)**2
    return nu,nu_prime,g,0.0


def _integrate(omega: float, parameter: float, *, sensitivities: bool) -> dict:
    if not np.isfinite([omega,parameter]).all() or omega<=0.0:
        raise ValueError("omega must be positive and inputs finite")
    angular=ELL*(ELL+1)

    def rhs(r: float,z: np.ndarray) -> np.ndarray:
        W,V=z[:2]
        nu,nu_prime,g,_=_profile(r)
        chi=CHI0+parameter*g
        if chi<=0.0:
            raise ValueError("toy adiabatic index chi became nonpositive")
        e_to_Lambda=1.0
        alpha=(omega*r)**2*e_to_Lambda*np.exp(-2.0*nu)
        A=np.array([
            [nu_prime/chi,alpha/chi-angular*e_to_Lambda],
            [-e_to_Lambda/(r*r),2.0*nu_prime],
        ])
        y=np.array([W,V])
        if not sensitivities:
            return A@y
        S_parameter=z[2:4]
        S_omega=z[4:6]
        A_parameter=np.array([
            [-nu_prime*g/(chi*chi),-alpha*g/(chi*chi)],
            [0.0,0.0],
        ])
        A_omega=np.array([
            [0.0,2.0*omega*r*r*e_to_Lambda*np.exp(-2.0*nu)/chi],
            [0.0,0.0],
        ])
        return np.concatenate([
            A@y,
            A@S_parameter+A_parameter@y,
            A@S_omega+A_omega@y,
        ])

    r0=CENTER_RADIUS
    # Regular l=2 centre expansion with W=A r^3 and V=-A r^2/2.
    y0=np.array([r0**3,-0.5*r0**2])
    if sensitivities:
        y0=np.concatenate([y0,np.zeros(4)])
    sol=solve_ivp(
        rhs,(r0,RADIUS),y0,method="DOP853",
        rtol=RTOL,atol=RTOL*1.0e-3,max_step=0.025,
    )
    if not sol.success:
        raise RuntimeError(f"DOP853 tangent integration failed: {sol.message}")

    W,V=sol.y[:2,-1]
    nu,nu_prime,g,_=_profile(RADIUS)
    e_to_Lambda=1.0
    alpha=(omega*RADIUS)**2*e_to_Lambda*np.exp(-2.0*nu)
    residual=alpha*V+nu_prime*W
    result={"residual":float(residual)}
    if not sensitivities:
        return result

    S_parameter=sol.y[2:4,-1]
    S_omega=sol.y[4:6,-1]
    alpha_omega=2.0*omega*RADIUS**2*e_to_Lambda*np.exp(-2.0*nu)
    D_parameter=alpha*S_parameter[1]+nu_prime*S_parameter[0]
    D_omega=alpha*S_omega[1]+nu_prime*S_omega[0]+alpha_omega*V
    if not np.isfinite([D_parameter,D_omega]).all() or abs(D_omega)<1e-12:
        raise RuntimeError("simple-root derivative gate failed")
    result.update({
        "D_parameter":float(D_parameter),
        "D_omega":float(D_omega),
        "implicit_frequency_derivative":float(-D_parameter/D_omega),
    })
    return result


def surface_residual(omega: float, parameter: float=0.0) -> float:
    return _integrate(omega,parameter,sensitivities=False)["residual"]


def build_mock(output_dir: Path) -> dict:
    output_dir.mkdir(parents=True,exist_ok=True)
    lo,hi=ROOT_BRACKET
    root=brentq(lambda w:surface_residual(w,0.0),lo,hi,xtol=1e-11,rtol=1e-12)
    tangent=_integrate(root,0.0,sensitivities=True)
    derivative=tangent["implicit_frequency_derivative"]
    checks=[]
    for h in (1.0e-3,5.0e-4,2.5e-4):
        plus=brentq(lambda w:surface_residual(w,+h),lo,hi,xtol=1e-11,rtol=1e-12)
        minus=brentq(lambda w:surface_residual(w,-h),lo,hi,xtol=1e-11,rtol=1e-12)
        finite_difference=(plus-minus)/(2.0*h)
        relative_error=abs(finite_difference-derivative)/max(abs(derivative),1e-14)
        checks.append({
            "step":h,"omega_plus":float(plus),"omega_minus":float(minus),
            "central_difference":float(finite_difference),
            "relative_error":float(relative_error),
        })
    max_error=max(row["relative_error"] for row in checks)
    if abs(tangent["residual"])>1e-9 or max_error>2e-7:
        raise AssertionError("Cowling tangent residual/finite-difference acceptance failed")
    if not all(checks[i+1]["relative_error"]<0.4*checks[i]["relative_error"] for i in range(2)):
        raise AssertionError("central finite-difference error did not show second-order h refinement")

    result={
        "schema":"ns-cowling-cowling-tangent-mock-v1",
        "status":"MOCK_ONLY_NOT_A_PHYSICAL_NEUTRON_STAR_MODE",
        "background":{
            "nu_r":"r^2/2","Lambda_r":"0","radius":RADIUS,
            "chi_r_lambda":"0.3 + lambda*r^2*(1-r)^2",
            "ell":ELL,"c_dimensionless":1.0,"center_radius":CENTER_RADIUS,
            "root_bracket":list(ROOT_BRACKET),"ode_rtol":RTOL,
        },
        "root_frequency_dimensionless":float(root),
        "root_residual_abs":abs(float(tangent["residual"])),
        "D_parameter":tangent["D_parameter"],
        "D_omega":tangent["D_omega"],
        "implicit_frequency_derivative":derivative,
        "finite_difference_checks":checks,
        "max_relative_error":max_error,
        "interpretation":"Forward tangent integration for the smooth ell=2 Cowling coefficient system agrees with independently re-solved roots. No stellar EOS, TOV solution, mode-family label or sharp-interface law is represented.",
    }
    (output_dir/"cowling_tangent_variation.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return result


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=Path("dist/cowling-variation-mock"))
    args=parser.parse_args(argv)
    result=build_mock(args.output_dir)
    print(json.dumps({
        "status":result["status"],
        "root_frequency_dimensionless":result["root_frequency_dimensionless"],
        "root_residual_abs":result["root_residual_abs"],
        "implicit_frequency_derivative":result["implicit_frequency_derivative"],
        "finite_difference_relative_errors":[x["relative_error"] for x in result["finite_difference_checks"]],
        "max_relative_error":result["max_relative_error"],
    },indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
