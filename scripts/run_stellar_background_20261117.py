#!/usr/bin/env python3
"""Bounded D exact-contract TOV sequence and independent native comparison."""
import csv, hashlib, json, os, platform, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import scipy
from compose_eos.stellar_background import (
    ContractEnthalpyEOS, solve_background, G_SI,C_SI,M_SUN_KG,
    GM_SUN_NOMINAL,MEV_FM3_TO_PA,MEV_FM3_TO_GEOM_KM2)
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"pipeline"))  # avoid scripts/compose_cowling.py shadowing package
OUT=ROOT/"data/stellar_consumers_20261117/tov"
SOURCE=ROOT/"data/core_support_20261117/primary_v1/REPRESENTATION.json"
def dump(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+"\n")
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write_csv(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def normalized_sample(eos,path,subdivisions,hcut):
    # Retain all source-law joins; this is a named native PCHIP interpolation
    # variant, not a rewrite of the exact contract.
    breaks=np.unique(np.r_[hcut,eos.break_h[eos.break_h>hcut]])
    hs=np.unique(np.concatenate([np.linspace(a,b,subdivisions+1)
                                 for a,b in zip(breaks[:-1],breaks[1:])]))
    header=("source_row,i_t,i_nb,i_yq,n_b_fm3,n_b_m3,pressure_mev_fm3,pressure_pa,"
            "energy_mev_fm3,energy_j_m3,mu_b_mev,mu_b_j,q1_p_over_nb_mev,"
            "q2_entropy_per_baryon,q3_mu_b_over_mn_minus_1,q4_mu_q_over_mn,"
            "q5_mu_l_over_mn,q6_f_over_nb_mn_minus_1,q7_e_over_nb_mn_minus_1")
    with path.open("w",newline="") as f:
        f.write(header+"\n")
        w=csv.writer(f)
        for i,h in enumerate(hs):
            s=eos.state_at_enthalpy(h);n,p,e,mu=(s[k] for k in ("n","p","epsilon","mu"))
            w.writerow([i+1,1,i+1,1,n,n*1e45,p,p*MEV_FM3_TO_PA,
                        e,e*MEV_FM3_TO_PA,mu,mu*1.602176634e-13,p/n,
                        0,mu/939.5654133-1,0,0,e/n/939.5654133-1,e/n/939.5654133-1])
    return len(hs)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    start=time.time()
    eos=ContractEnthalpyEOS.from_path(SOURCE)
    provenance=dict(source=str(SOURCE.relative_to(ROOT)),source_sha256=sha(SOURCE),
        physical_csv_sha256=sha(SOURCE.with_name("EOS_PHYSICAL.csv")),
        git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
        invocation=sys.argv,runner_sha256=sha(__file__),
        implementation_sha256=sha(ROOT/"python_frontend/compose_eos/stellar_background.py"),
        test_sha256=sha(ROOT/"python_frontend/tests/test_stellar_background.py"),
        PYTHONPATH=os.environ.get("PYTHONPATH"),
        constants=dict(c_m_s=C_SI,G_SI=G_SI,GM_sun_nominal=GM_SUN_NOMINAL,
                       M_sun_kg=M_SUN_KG,MeV_fm3_to_geom_km2=MEV_FM3_TO_GEOM_KM2),
        backend="exact contract Python DOP853 enthalpy; native independent table comparison",
        source_support_n_fm3=list(eos.density_domain),
        source_support_p_MeV_fm3=list(eos.pressure_domain),
        composition="cold effective beta-equilibrium barotrope; masks unmodified")
    dump("PROVENANCE.json",provenance)
    config=dict(r0_km=.001,rtol=2e-10,atol=1e-12,max_step_h=.002,surface_h=0.)
    if "--native-only" not in sys.argv:
        seq=[]
        for n in np.linspace(.1,1.,37):
            b=solve_background(eos,n,**config)
            row={k:v for k,v in b.summary().items() if k not in ("config","backend")}
            seq.append(row)
            print("sequence",n,b.mass_msun,b.radius_km,flush=True)
        write_csv(OUT/"sequence.csv",seq)
        nc=brentq(lambda n:solve_background(eos,n,**config).mass_msun-1.4,.5,.6,xtol=2e-10)
        star=solve_background(eos,nc,**config)
        dump("CANONICAL_1P4.json",star.summary())
        print("CANONICAL",json.dumps(star.summary()),flush=True)
        hs=np.unique(np.r_[np.linspace(0,star.h_start,1001),eos.break_h[(eos.break_h>0)&(eos.break_h<star.h_start)]])[::-1]
        profile=[star.state_at_enthalpy(h) for h in hs]
        write_csv(OUT/"canonical_profile.csv",profile)
        configs=[
            ("rtol_1e-7",dict(config,rtol=1e-7,atol=1e-10,max_step_h=.01)),
            ("rtol_1e-8",dict(config,rtol=1e-8,atol=1e-11,max_step_h=.005)),
            ("rtol_1e-9",dict(config,rtol=1e-9,atol=1e-12,max_step_h=.003)),
            ("reference",config),
            ("rtol_3e-11",dict(config,rtol=3e-11,atol=2e-13,max_step_h=.001)),
            ("center_10m",dict(config,r0_km=.01)),
            ("center_0p3m",dict(config,r0_km=.0003)),
            ("surface_h1e-4",dict(config,surface_h=1e-4)),
            ("surface_h1e-6",dict(config,surface_h=1e-6)),
            ("surface_h1e-8",dict(config,surface_h=1e-8))]
        convergence=[]
        for label,cfg in configs:
            b=solve_background(eos,nc,**cfg)
            row=dict(label=label,**cfg,mass_msun=b.mass_msun,radius_km=b.radius_km,
                     delta_mass_msun=b.mass_msun-star.mass_msun,
                     delta_radius_m=b.radius_km*1000-star.radius_km*1000)
            convergence.append(row);print("convergence",row,flush=True)
        write_csv(OUT/"convergence.csv",convergence)
    else:
        with (OUT/"sequence.csv").open() as f:
            seq=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
        canonical=json.loads((OUT/"CANONICAL_1P4.json").read_text())
        nc=canonical["central_density_fm3"]
        star=solve_background(eos,nc,**config)
    # Same h cutoff and physical central pressure isolates native interpolation/ODE.
    from compose_cowling.solver import _compile_solver, _parse_tsv
    case_driver,_ = _compile_solver()
    hcut=1e-4
    cutstar=solve_background(eos,nc,**dict(config,surface_h=hcut))
    pc=eos.state_at_enthalpy(star.h_c)["p"]*MEV_FM3_TO_PA
    pcut=eos.state_at_enthalpy(hcut)["p"]*MEV_FM3_TO_PA
    native=[]
    for subdivisions in (1,4):
        path=OUT/("native_exact_samples_sub%d.csv"%subdivisions)
        count=normalized_sample(eos,path,subdivisions,hcut)
        rec=dict(subdivisions=subdivisions,rows=count,table_sha256=sha(path),
                 reference=cutstar.summary(),pc_pa=pc,p_cut_pa=pcut)
        try:
            nml=OUT/("native_sub%d.nml"%subdivisions)
            nml.write_text("&request request_token='D-exact-contract-comparison' /\n&parameters\n"
                +f"table_path='{path}',pc_pa={pc:.17e},p_cut_pa={pcut:.17e},\n"
                +"r0_m=1.0,r_limit_m=1.0e6,rtol=2.0e-10,atol_r_m=1.0e-6,"
                +"atol_mass_kg=1.0e14,atol_nu=1.0e-12,dq_initial=1.0e-5,dq_min=1.0e-16,"
                +"dq_max=0.2,nu_start=0.0,max_steps=100000,emit_profile=.false.\n/\n")
            command=[str(case_driver),"real_table_star",str(nml)]
            run=subprocess.run(command,capture_output=True,text=True)
            (OUT/("native_sub%d.stdout.log"%subdivisions)).write_text(run.stdout)
            (OUT/("native_sub%d.stderr.log"%subdivisions)).write_text(run.stderr)
            rec.update(command=command,exit_code=run.returncode,
                       executable_sha256=sha(case_driver),namelist_sha256=sha(nml))
            result=_parse_tsv(run.stdout)
            if run.returncode or result.get("status")!="ok":
                raise RuntimeError("native TOV rejected: "+run.stdout+run.stderr)
            result["backend"]="Fortran"
            rec.update(status="ok",result=result,
                delta_radius_m=result["radius_m"]-cutstar.radius_km*1000,
                delta_mass_msun=result["mass_kg"]/M_SUN_KG-cutstar.mass_msun)
        except Exception as exc:
            rec.update(status="FAILED",error=repr(exc))
        native.append(rec);print("native",json.dumps(rec),flush=True)
    dump("NATIVE_COMPARISON.json",native)
    slopes=np.diff([r["mass_msun"] for r in seq])/np.diff([r["central_density_fm3"] for r in seq])
    result=dict(status="SUPPORTED_OBSERVABLES_COMPUTED",canonical=star.summary(),
                sequence_points=len(seq),sequence_all_sampled_mass_slopes_positive=bool(np.all(slopes>0)),
                endpoint=seq[-1],last_mass_slope_Msun_per_fm3=float(slopes[-1]),
                maximum_mass_status="NOT_BRACKETED_WITHIN_SOURCE_SUPPORT",
                upper_extension="required to determine maximum mass; not performed",
                causality="inherited exact piecewise continuous validation in WholeCrustBarotrope",
                elapsed_seconds=time.time()-start,native_status=[r["status"] for r in native],
                validation="numerical tolerances, central radius, surface cutoff and independent native comparison")
    dump("RESULT.json",result)
    print("FINISHED",json.dumps(result),flush=True)
if __name__=="__main__":
    main()
