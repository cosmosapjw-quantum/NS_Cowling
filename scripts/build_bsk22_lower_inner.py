#!/usr/bin/env python3
"""Generate both published-fit and energy-derivative BSk22 lower-inner representations."""
from __future__ import annotations
import argparse,csv,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import bsk22_analytic as bsk
from compose_eos.api import read_compose,write_compose,check_consistency

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("output",type=Path)
    ap.add_argument("--points",type=int,default=129)
    ap.add_argument("--author-reference-exe",type=Path)
    args=ap.parse_args()
    if args.points<5:ap.error("--points must be >=5")
    args.output.mkdir(parents=True,exist_ok=False)
    ns=np.geomspace(bsk.NEUTRON_DRIP_DENSITY,bsk.LOWER_INNER_MAX_DENSITY,args.points)
    versions={}
    for mode in ("source_fit","energy_derivative"):
        table,rows=bsk.eos_data(ns,thermodynamic_mode=mode)
        target=args.output/mode;target.mkdir()
        write_compose(target/"compose",table);restored=read_compose(target/"compose")
        errors={}
        for field in ("nb","Q1","Q3","Q4","Q6","Q7","Iphase","Yn","Yp","Ye","Ymu","Aav","Zav","Yh"):
            before,after=getattr(table,field),getattr(restored,field)
            errors[field]=float(np.max(abs(before-after)))
            if not np.allclose(before,after,rtol=5e-12,atol=5e-14):
                raise RuntimeError("native composition/thermodynamic readback mismatch: "+field)
        with (target/"states.csv").open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        derivative_errors=[];sound_errors=[]
        for n in ns[1:-1]:
            s=bsk.state(n,thermodynamic_mode=mode);h=1e-5
            lo=bsk.state(n*np.exp(-h),thermodynamic_mode=mode)
            hi=bsk.state(n*np.exp(h),thermodynamic_mode=mode)
            p_energy=n*(hi["energy_per_baryon"]-lo["energy_per_baryon"])/(2*h)
            derivative_errors.append(abs(p_energy/s["pressure_from_energy_derivative"]-1))
            cs=(hi["pressure"]-lo["pressure"])/(hi["epsilon"]-lo["epsilon"])
            sound_errors.append(abs(cs/s["sound_speed_squared_over_c_squared"]-1))
        versions[mode]=dict(
            points=len(rows),native_readback_max_absolute_errors=errors,
            native_checks=check_consistency(restored,verbose=False),
            analytic_energy_derivative_max_relative_error=max(derivative_errors),
            analytic_sound_speed_max_relative_error=max(sound_errors),
            min_sound_speed_squared=min(r["sound_speed_squared_over_c_squared"] for r in rows),
            max_sound_speed_squared=max(r["sound_speed_squared_over_c_squared"] for r in rows),
            first=rows[0],last=rows[-1])
    rows=[bsk.state(n) for n in ns]
    author=None
    if args.author_reference_exe:
        # Avoid source single-precision AXND rounding at exactly neutron drip.
        comparison_ns=np.geomspace(bsk.NEUTRON_DRIP_DENSITY*(1+1e-7),.002,17)
        result=subprocess.run([str(args.author_reference_exe.resolve())],
            input="\n".join(map(str,comparison_ns))+"\n",text=True,capture_output=True,check=True)
        (args.output/"author_reference_stdout.log").write_text(result.stdout)
        values=np.array([[float(v) for v in line.split()] for line in result.stdout.splitlines() if len(line.split())==12])
        if len(values)!=len(comparison_ns):raise RuntimeError("author reference output shape mismatch")
        cases=[]
        for n,ref in zip(comparison_ns,values):
            s=bsk.state(n)
            cases.append(dict(n_b=float(n),
                energy_minus_rest_absolute_difference=s["energy_minus_neutron_rest_MeV"]-ref[1],
                density_relative_difference=s["rho_g_cm3"]/ref[2]-1,
                pressure_relative_difference=s["pressure"]/ref[3]-1,
                chi_relative_difference=s["dlogP_dlogrho"]/ref[4]-1,
                Acluster_relative_difference=s["A_cluster"]/ref[5]-1,
                Acell_relative_difference=s["A_cell"]/ref[6]-1,
                Zcell_difference=s["Z_cell"]-ref[7],Zcluster_difference=s["Z_cluster"]-ref[8],
                mu_n_absolute_difference=s["mu_n_composition_fit"]-bsk.NEUTRON_REST_MEV-ref[9],
                mu_e_relative_difference=s["mu_e"]/ref[10]-1,
                free_neutron_fraction_absolute_difference=s["Y_n_free"]-ref[11]))
        author=dict(cases=cases,known_difference="Paper exponent7/6 versus author Fortran1.16667; author's default-real literal rounding; no coefficients refitted.",
            source_sha256=bsk.AUTHOR_CODE_SHA256)
        (args.output/"author_reference_comparison.json").write_text(json.dumps(author,indent=2)+"\n")
    outer_path=ROOT/"data/research_completion_20261116/outer_crust/generated/RECEIPT.json"
    seam=None
    if outer_path.exists():
        outer=json.loads(outer_path.read_text())["drip"]
        right=bsk.state(outer["n_b"])
        seam=dict(outer=outer,analytic_at_outer_density=right,
            relative_pressure_difference=right["pressure"]/outer["pressure"]-1,
            mu_b_difference_MeV=right["mu_b_euler"]-outer["g"],
            electron_fraction_difference=right["Y_e"]-outer["Y_e"],
            status="NOT_JOINED; source fit and independent outer endpoint differ")
    receipt=dict(status="PUBLISHED_ANALYTIC_REPRESENTATION_VERIFIED_SCOPED",
        scope="BSk22 lower-inner crust only; two explicit thermodynamic representations, composition included",
        limitation="No microscopic ETF/Strutinsky regeneration or automatic join to outer/upper crust",
        source=bsk.SOURCE_URL,author_code=bsk.AUTHOR_CODE_URL,author_code_sha256=bsk.AUTHOR_CODE_SHA256,
        fit_error_from_source="Energy typical<=0.5%, maximum<3%; pressure typical1%, maximum<=4%, strongest near phase boundaries. Derivatives need separate checks.",
        support_fm3=[float(ns[0]),float(ns[-1])],
        pressure_fit_thermodynamic_residual_min=min(r["pressure_fit_thermodynamic_relative_residual"] for r in rows),
        pressure_fit_thermodynamic_residual_max=max(r["pressure_fit_thermodynamic_relative_residual"] for r in rows),
        fractions_closure_max=max(max(abs(r["composition_baryon_residual"]),abs(r["composition_charge_residual"])) for r in rows),
        versions=versions,author_reference_checked=author is not None,outer_endpoint_diagnostic=seam,
        acceptance="Equation implementation, independently checked derivatives and exact composition/readback; no unconditional CompOSE matching gate",
        source_code_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(bsk.__file__),Path(__file__)]},
        head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip())
    (args.output/"RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ("status","support_fm3","pressure_fit_thermodynamic_residual_min","pressure_fit_thermodynamic_residual_max","fractions_closure_max","author_reference_checked")},indent=2))
    print(json.dumps({k:{j:v[j] for j in ("analytic_energy_derivative_max_relative_error","analytic_sound_speed_max_relative_error","min_sound_speed_squared","max_sound_speed_squared")} for k,v in versions.items()},indent=2))
if __name__=="__main__":main()
