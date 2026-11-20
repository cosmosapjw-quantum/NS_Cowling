#!/usr/bin/env python3
"""Thin per-star frontend: native Fortran scan, root polishing and profile.

Python supplies the frozen EOS and marshals results. No Python root solver
is called; native DVODE COMMON remains confined to one OS worker process.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))

from compose_eos.stellar_background import ContractEnthalpyEOS
from compose_eos.stellar_cowling_tov_native import solve_background_native
from compose_eos.stellar_cowling_native import NativeCowlingDVODE
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json"
SOURCE_SHA="591ac27aa99448512e291684c98f8273efc729d17c63be88e331c20e83bd2653"


def run(density:float,out:Path):
    if not np.isfinite(density) or density<=0:
        raise ValueError("central density must be finite and positive")
    if hashlib.sha256(SOURCE.read_bytes()).hexdigest()!=SOURCE_SHA:
        raise RuntimeError("original EOS source drift")
    eos=ContractEnthalpyEOS.from_path(SOURCE)
    star=solve_background_native(eos,density,surface_h=1.e-4,
        rtol=1.e-10,atol=1.e-12,max_step_h=.002,native_samples=4097)
    settings=CowlingSettings(surface_h=1.e-4,rtol=1.e-9,
        atol=1.e-12,max_step_h=.002,root_xtol_hz=1.e-5)
    modes=NativeCowlingDVODE(star,settings,samples=2049,method="adams")
    fs=np.linspace(1800.,7200.,121)
    brackets,scan=modes.scan_brackets(fs)
    candidates=[]
    for lo,hi in brackets:
        item={"bracket_hz":[lo,hi]}
        try:
            m=modes.solve_mode((lo,hi),profile_points=240,root_backend="fortran")
        except (ValueError,RuntimeError) as exc:
            item.update(status="REJECTED_NATIVE_ROOT",reason=str(exc))
        else:
            item.update(status="ROOT_COMPUTED_NOT_MODE_CLASSIFIED",
                        frequency_hz=m["frequency_hz"],
                        determinant_residual=m["matching_residual"],
                        surface_relative=m["surface_BC_relative"],
                        join_relative=m["relative_join_mismatch"],
                        W_nodes=m["W_nodes"],V_nodes=m["V_nodes"],
                        root_backend=m["root_backend"],
                        root_iterations=m["root_iterations"])
        candidates.append(item)
    result=dict(schema="ns_cowling.c10.native_root_independent_star.v1",
        source_sha256=SOURCE_SHA,central_density_fm3=float(density),
        mass_msun=float(star.mass_msun),radius_km=float(star.radius_km),
        root_method="Fortran Brent-Dekker + Fortran DVODE",
        frequency_grid_hz=[1800.,7200.,len(fs)],
        n_brackets=len(brackets),n_root_computed=sum(
            item["status"]=="ROOT_COMPUTED_NOT_MODE_CLASSIFIED" for item in candidates),
        n_rejected=sum(item["status"]=="REJECTED_NATIVE_ROOT" for item in candidates),
        candidates=candidates,mode_identification="NOT_CERTIFIED_FROM_SIGN_ONLY",
        bulk_composition_g_modes_supported=False,
        status="COMPLETED_SCOPED_ROOT_CANDIDATE_CENSUS")
    out.parent.mkdir(exist_ok=True,parents=True)
    partial=out.with_name(out.name+".part")
    partial.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    os.replace(partial,out)
    print("NATIVE_C10_ROOT",density,"mass",round(star.mass_msun,8),
          "brackets",len(brackets),"accepted_roots",result["n_root_computed"],
          "rejected",result["n_rejected"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--central-density",type=float,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    run(args.central_density,args.output)
