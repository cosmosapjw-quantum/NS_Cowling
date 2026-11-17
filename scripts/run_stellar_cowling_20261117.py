#!/usr/bin/env python3
"""Run bounded exact B-contract l=2 Cowling observables and convergence."""
import argparse,csv,hashlib,json,time
from pathlib import Path
from dataclasses import asdict,replace
import numpy as np
from scipy.optimize import brentq
from compose_eos.stellar_background import ContractEnthalpyEOS,solve_background
from compose_eos.stellar_cowling import CowlingSettings,scan_brackets,solve_mode

ROOT=Path(__file__).resolve().parents[1]
def write_json(p,x):
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--central-density",type=float)
    ap.add_argument("--output",default="data/stellar_consumers_20261117/cowling")
    args=ap.parse_args()
    out=ROOT/args.output;out.mkdir(parents=True,exist_ok=True)
    source=ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json"
    eos=ContractEnthalpyEOS.from_path(source)
    nc=args.central_density
    if nc is None:
        nc=brentq(lambda n:solve_background(eos,n,rtol=2e-10,atol=2e-12).mass_msun-1.4,
                  .5,.6,xtol=2e-10)
    base_bg=solve_background(eos,nc,rtol=1e-10,atol=1e-12)
    base=CowlingSettings(rtol=1e-8,atol=1e-11,surface_h=1e-9,max_step_h=.002)
    scan_settings=replace(base,rtol=1e-5,atol=1e-8,surface_h=1e-7,max_step_h=.005)
    brackets,scan=scan_brackets(base_bg,np.arange(250.,10501.,250.),scan_settings)
    write_json(out/"SCAN.json",dict(settings=asdict(scan_settings),scan=scan,brackets=brackets))
    print("nc",nc,"scan",brackets,flush=True)
    if len(brackets)<2:raise RuntimeError("Scan did not resolve two modes")
    configurations=[
        ("baseline",base,{}),
        ("ode_loose",replace(base,rtol=1e-6,atol=1e-9),{}),
        ("ode_tight",replace(base,rtol=1e-9,atol=1e-12),{}),
        ("surface_h1e-6",replace(base,surface_h=1e-6),{}),
        ("surface_h1e-7",replace(base,surface_h=1e-7),{}),
        ("surface_h1e-8",replace(base,surface_h=1e-8),{}),
        ("match_h007",replace(base,match_h_fraction=.07),{}),
        ("match_h030",replace(base,match_h_fraction=.30),{}),
        ("step_h1e-3",replace(base,max_step_h=.001),{}),
        ("center_10m",base,dict(r0_km=.01)),
        ("center_01m",base,dict(r0_km=.0001)),
        ("background_loose",base,dict(rtol=1e-8,atol=1e-10)),
    ]
    results=[]
    for label,settings,changes in configurations:
        start=time.monotonic()
        kw=dict(rtol=1e-10,atol=1e-12);kw.update(changes)
        bg=base_bg if not changes else solve_background(eos,nc,**kw)
        modes=[]
        for j,bracket in enumerate(brackets[:2]):
            mode=solve_mode(bg,bracket,settings)
            profile=mode.pop("profile")
            with (out/(label+("_f_profile.csv" if j==0 else "_p1_profile.csv"))).open("w") as f:
                writer=csv.DictWriter(f,fieldnames=profile[0].keys());writer.writeheader();writer.writerows(profile)
            mode["name"]="f" if j==0 else "p1"
            if mode["W_nodes"]!=j:
                raise RuntimeError("Unexpected mode node count: "+str(mode))
            mode["residual_accepted"]=(abs(mode["matching_residual"])<=2e-6 and mode["relative_join_mismatch"]<=2e-6)
            modes.append(mode)
        row=dict(label=label,settings=asdict(settings),background=bg.summary(),modes=modes,
                 wall_seconds=time.monotonic()-start)
        results.append(row);write_json(out/(label+".json"),row)
        print(label,[(m["frequency_hz"],m["W_nodes"],m["matching_residual"]) for m in modes],
              "seconds",row["wall_seconds"],flush=True)
    reference=results[0]["modes"]
    for row in results:
        row["frequency_shifts_hz"]=[m["frequency_hz"]-ref["frequency_hz"] for m,ref in zip(row["modes"],reference)]
    summary=dict(schema="ns-cowling.stellar-cowling.v1",backend="Python exact-contract enthalpy DOP853 reference",
                 closure="equilibrium-barotropic perfect fluid: Gamma1=Gamma_eq",
                 convention="(-+++), geometric km, Phi=half temporal metric, frequency at infinity",
                 interface="continuous EOS knots: W,V continuous; no sharp Maxwell jump",
                 source=str(source.relative_to(ROOT)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                 central_density_fm3=nc,background=base_bg.summary(),configurations=results,
                 scope="nonrotating, unmagnetized, no metric perturbation, no elasticity, no superfluid relative motion",
                 unresolved="Frozen composition and composition buoyancy require missing off-equilibrium input; no full-GR mode claim")
    write_json(out/"SUMMARY.json",summary)
    print("Saved",out/"SUMMARY.json",flush=True)
if __name__=="__main__":main()
