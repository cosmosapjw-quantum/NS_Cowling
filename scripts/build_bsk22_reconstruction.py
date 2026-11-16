#!/usr/bin/env python3
"""Validate the optional boundary-anchored BSk22 pressure reconstruction."""
from __future__ import annotations
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.bsk22_thermodynamic_reconstruction import BoundaryAnchoredBSk22
from compose_eos import bsk22_analytic as source
from compose_eos.api import read_compose,write_compose,check_consistency

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--outer-receipt",type=Path,default=ROOT/"data/research_completion_20261116/outer_crust/generated/RECEIPT.json")
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    anchor=json.loads(args.outer_receipt.read_text())["drip"]
    model=BoundaryAnchoredBSk22(outer_density=anchor["n_b"],outer_pressure=anchor["pressure"],outer_mu_b=anchor["g"],outer_epsilon=anchor["epsilon"])
    glx,glw=np.polynomial.legendre.leggauss(96)
    outputs={};solutions={}
    for points in (129,257):
        grid=np.geomspace(model.n0,source.LOWER_INNER_MAX_DENSITY,points)
        data,rows=model.eos_data(grid);target=args.output/str(points);target.mkdir()
        write_compose(target/"compose",data);restored=read_compose(target/"compose")
        readback={}
        for key in ("nb","Q1","Q3","Q4","Q6","Q7","Yn","Yp","Ye","Ymu","Aav","Zav","Yh"):
            a,b=getattr(data,key),getattr(restored,key);readback[key]=float(np.max(abs(a-b)))
            if not np.allclose(a,b,rtol=5e-12,atol=5e-14):raise RuntimeError("native readback mismatch:"+key)
        with (target/"states.csv").open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        gl_errors=[];firstlaw=[];sound=[];pressure_derivative=[]
        for n in grid[1:-1:4]:
            r=model.state(n);h=1e-4
            lo=model.state(n*np.exp(-h));hi=model.state(n*np.exp(h))
            pp=n*(hi["energy_per_baryon"]-lo["energy_per_baryon"])/(2*h)
            firstlaw.append(abs(pp/r["pressure"]-1))
            cs=(hi["pressure"]-lo["pressure"])/(hi["epsilon"]-lo["epsilon"])
            sound.append(abs(cs/r["sound_speed_squared_over_c_squared"]-1))
            dp=(hi["pressure"]-lo["pressure"])/(2*h*n)
            pressure_derivative.append(abs(dp/r["pressure_density_derivative"]-1))
            upper=np.log(n/model.n0);nodes=.5*upper*(glx+1)
            integral=.5*upper*sum(w*source.state(model.n0*np.exp(x))["pressure"]/ (model.n0*np.exp(x)) for x,w in zip(nodes,glw))
            gl_errors.append(abs(model.e0+integral-r["energy_per_baryon"]))
        outputs[str(points)]=dict(points=points,native_readback_max_abs=readback,native_checks=check_consistency(restored,verbose=False),
            firstlaw_finite_difference_max_relative_error=max(firstlaw),
            sound_speed_finite_difference_max_relative_error=max(sound),
            pressure_derivative_max_relative_error=max(pressure_derivative),
            gauss_legendre96_vs_quad_max_energy_difference_MeV=max(gl_errors),
            quadrature_error_estimate_max_MeV=max(r["first_law_quadrature_error_MeV"] for r in rows),
            energy_shift_MeV_range=[min(r["energy_shift_from_source_fit_MeV"] for r in rows),max(r["energy_shift_from_source_fit_MeV"] for r in rows)],
            source_composition_mu_n_residual_MeV_range=[min(r["source_composition_mu_n_residual_MeV"] for r in rows),max(r["source_composition_mu_n_residual_MeV"] for r in rows)],
            cs2_range=[min(r["sound_speed_squared_over_c_squared"] for r in rows),max(r["sound_speed_squared_over_c_squared"] for r in rows)],
            original_pressure_max_absolute_difference=max(abs(r["pressure"]-source.state(r["n_b"])["pressure"]) for r in rows),
            first=rows[0],last=rows[-1])
        if max(firstlaw)>2e-7 or max(sound)>2e-7 or max(pressure_derivative)>2e-7 or max(gl_errors)>1e-9:
            raise RuntimeError("independent first-law/quadrature/derivative verification failed")
        solutions[points]=rows
    common=max(abs(a["energy_per_baryon"]-b["energy_per_baryon"]) for a,b in zip(solutions[129],solutions[257][::2]))
    receipt=dict(status="BOUNDARY_ANCHORED_ANALYTIC_REPRESENTATION_VERIFIED_SCOPED",
        representation="Preserve published BSk22 P(n) and composition; independently anchor and integrate total energy",
        source_recommendation="https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html",
        anchor_input_receipt=str(args.outer_receipt.relative_to(ROOT)) if args.outer_receipt.is_relative_to(ROOT) else str(args.outer_receipt),
        anchor_input_sha256=hashlib.sha256(args.outer_receipt.read_bytes()).hexdigest(),
        anchor=model.anchor_diagnostics(),outputs=outputs,common_grid_energy_difference_MeV=common,
        not_claimed=["microscopic ETF/Strutinsky minimization","new neutron-drip solution","exact original energy-fit reproduction","source chemical-potential-fit stationarity","global phase minimum","core-crust matching"],
        source_code_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/"python_frontend/compose_eos/bsk22_thermodynamic_reconstruction.py",Path(source.__file__)]})
    (args.output/"RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ("status","anchor","common_grid_energy_difference_MeV")},indent=2))
    print(json.dumps({k:{j:v[j] for j in ("firstlaw_finite_difference_max_relative_error","sound_speed_finite_difference_max_relative_error","gauss_legendre96_vs_quad_max_energy_difference_MeV","energy_shift_MeV_range","source_composition_mu_n_residual_MeV_range")} for k,v in outputs.items()},indent=2))
if __name__=="__main__":main()
