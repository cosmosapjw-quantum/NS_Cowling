#!/usr/bin/env python3
"""Reproducible source-locked PCHIP coefficient export for standalone Fortran.

This is a frontend input producer, not a Python shooting/root solver.
Fortran reads the phase-separated coefficient contract and computes
the slow/rapid eigenfrequencies without invoking Python.
"""
from __future__ import annotations
import hashlib,io,json,sys,zipfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings

SRC=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
OUT=ROOT/"data/cowling_sharp_native_root_20261120/ORIGINAL_SOURCE_PCHIP_ONEJUMP.dat"
BENCH=ROOT/"data/cowling_sharp_native_root_20261120/PCHIP_ROOT_BENCHMARK.json"

def run():
    raw=SRC.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:
        raise RuntimeError("source-locked EOS archive hash changed")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        with zipfile.ZipFile(io.BytesIO(z.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inside:
            rep=inside.read("FRESH_5x2_WHOLE_STARS/SPG_M2/maxwell/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(rep)))
    star=solve_sharp_background_native(eos,.30385128926615257,
        surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    settings=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,
        max_step_h=.002,root_xtol_hz=1e-6)
    solver=NativeSharpCowlingDVODE(star,settings,
        samples=8193,coefficient_backend="pchip")
    assert len(solver._phases)==2 and len(star.interfaces)==1
    hi,lo=solver._phases
    ht=star.interfaces[0]["h"]
    slow=star.state_at_enthalpy(ht,"low")
    high=star.state_at_enthalpy(ht,"high")
    iface=np.array([slow["r_km"],slow["Phi"],slow["lambda_metric"],
       slow["dPhi_dr_km"],slow["p"],slow["epsilon"],high["epsilon"]],
        dtype=np.float64)
    cases=json.loads(BENCH.read_text())["cases"]
    ordering=(("slow_f",1),("rapid_f",2),("slow_p1",1),("rapid_p1",2))
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open("w") as stream:
        stream.write(f"{len(hi['s'])} {len(lo['s'])}\n")
        stream.write(f"{0.5*hi['s'][-1]:.17e} {solver._smax:.17e} "+
          f"{settings.rtol:.17e} {settings.atol:.17e} "+
          f"{min(solver._smax/200,settings.max_step_h/(2*solver._smax)):.17e} "+
          f"{settings.root_xtol_hz:.17e}\n")
        stream.write(" ".join(f"{v:.17e}" for v in iface)+"\n")
        for phase in hi,lo:
            matrix=np.column_stack((phase["s"],phase["vals"].T,phase["slopes"].T))
            np.savetxt(stream,matrix,fmt="%.17e")
        stream.write(f"{len(ordering)}\n")
        for name,mode in ordering:
            item=cases[name]
            stream.write(f"{mode:d} {item['bracket_hz'][0]:.17e} "+
                f"{item['bracket_hz'][1]:.17e} "+
                f"{item['baseline_frequency_hz']:.17e}\n")
    print("SOURCE_SHA",SHA)
    print("CONTRACT_SHA",hashlib.sha256(rep).hexdigest())
    print("INPUT_SHA",hashlib.sha256(OUT.read_bytes()).hexdigest())
    print("PHASE_SAMPLES",len(hi["s"]),len(lo["s"]))
    print("EXPORTED",OUT)
if __name__=="__main__":run()
