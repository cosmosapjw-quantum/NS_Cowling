#!/usr/bin/env python3
"""Generate cold APR branches, local-neutral Maxwell and two-charge Gibbs cores.

Source APR_APR has no composition. Generated composition is a prediction of the
published Hamiltonian with the declared mass/lepton/phase conventions; it is
not recovered provider data. No crust or core-crust joining is performed.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import scipy
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import apr_uniform as apr
from compose_eos import apr_gibbs as gibbs
from compose_eos.api import read_compose,write_compose,check_consistency
from compose_eos.eos_data import EOSData
from compose_eos.physical_domain import physical_domain_report


def compose_table(rows):
    if any(row["branch"]=="mixed" for row in rows):
        raise ValueError("mixed local-neutral states lack one global charge chemical potential")
    values=lambda name:np.array([row[name] for row in rows],dtype=float)
    n=values("n_b")
    zero=np.zeros(len(rows))
    return EOSData(nb=n,mn=apr.M_N,mp=apr.M_P,Q1=values("pressure")/n,Q2=zero.copy(),
        Q3=values("mu_n")/apr.M_N-1,Q4=-values("mu_e")/apr.M_N,Q5=zero.copy(),
        Q6=values("epsilon")/(n*apr.M_N)-1,Q7=values("epsilon")/(n*apr.M_N)-1,
        q_extra=zero.copy(),Iphase=np.array([{"low":6,"high":7,"gibbs_mixed":8}[row["branch"]] for row in rows]),
        Yn=values("Y_n"),Yp=values("Y_p"),Ye=values("Y_e"),Ymu=values("Y_mu"),
        Aav=np.ones(len(rows)),Zav=zero.copy(),Yh=zero.copy(),
        composition_file_present=True,composition_available=True,beta_equilibrated=True)


def sound_ratio(n,branch,relative_step):
    h=n*relative_step
    left,right=apr.beta_state(max(.01,n-h),branch),apr.beta_state(min(2.,n+h),branch)
    return (right["pressure"]-left["pressure"])/(right["epsilon"]-left["epsilon"])


def causal_crossings(branch,n_min,n_max):
    records=[]
    for step in (1e-4,3e-5,1e-5):
        grid=np.linspace(n_min,n_max,401)
        values=[sound_ratio(float(n),branch,step)-1 for n in grid]
        roots=[]
        for a,b,fa,fb in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
            if fa*fb<0:
                roots.append(float(brentq(lambda n:sound_ratio(n,branch,step)-1,a,b,xtol=1e-10)))
        records.append({"relative_density_difference_step":step,"cs2_equal_1_roots_fm3":roots})
    return {"method":"central thermodynamic derivative on beta-minimized homogeneous branch; one-sided at computational limits",
            "density_scan_fm3":[n_min,n_max],"refinements":records,
            "meaning":"model-functional crossing; not empirical validity endpoint"}


def save_rows(path,rows):
    with path.open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader()
        writer.writerows({key:(json.dumps(value,sort_keys=True) if isinstance(value,dict) else value)
                          for key,value in row.items()} for row in rows)


def maxwell_native_density_grid(density,phase):
    """Pure-branch rows and in-window coexistence endpoints, without extrapolation."""
    left,right=phase["left"]["n_b"],phase["right"]["n_b"]
    endpoints=[n for n in (left,right) if density[0]<=n<=density[-1]]
    return np.unique(np.r_[density[(density<=left)|(density>=right)],endpoints])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--source",type=Path,default=ROOT/"eos_generator/eos_data/APR_APR")
    parser.add_argument("--n-min",type=float,default=.085,
        help="analysis lower bound fm^-3, not a fitted crust transition")
    args=parser.parse_args()
    if not np.isfinite(args.n_min) or not .01<=args.n_min<1.34:
        parser.error("--n-min must be finite and in [.01,1.34)")
    args.output.mkdir(parents=True,exist_ok=False)
    source=read_compose(args.source)
    indices=np.flatnonzero((source.nb>=args.n_min)&(source.nb<=1.34))
    if len(indices)<3:
        raise ValueError("at least three source rows required in the analysis interval")
    density=source.nb[indices]
    phase=apr.maxwell_coexistence()
    gibbs_phase=gibbs.gibbs_coexistence()
    phase_left,phase_right=phase["left"]["n_b"],phase["right"]["n_b"]
    candidate=[apr.equilibrium_state(float(n),phase) for n in density]
    gibbs_candidate=[gibbs.gibbs_equilibrium_state(float(n),gibbs_phase) for n in density]
    endpoint_grid=maxwell_native_density_grid(density,phase)
    crossings={branch:causal_crossings(branch,float(density[0]),float(density[-1]))
               for branch in ("low","high")}
    models={"low_branch":[apr.beta_state(float(n),"low") for n in density],
            "high_branch":[apr.beta_state(float(n),"high") for n in density],
            "local_neutral_Maxwell":[apr.equilibrium_state(float(n),phase) for n in endpoint_grid]}
    gibbs_extra=[row["n_b"] for row in gibbs_phase["volume_fraction_path"]
                 if density[0]<=row["n_b"]<=density[-1]]
    models["bulk_Gibbs"]=[gibbs.gibbs_equilibrium_state(float(n),gibbs_phase)
                          for n in np.unique(np.r_[density,gibbs_extra])]
    roots=crossings["high"]["refinements"][-1]["cs2_equal_1_roots_fm3"]
    if roots:
        models["sampled_causal_Maxwell"]=[row for row in models["local_neutral_Maxwell"]
            if row["n_b"]<roots[0]]
        models["sampled_causal_Gibbs"]=[row for row in models["bulk_Gibbs"]
            if row["n_b"]<roots[0]]
    summaries={}
    for label,rows in models.items():
        folder=args.output/label
        folder.mkdir()
        save_rows(folder/"states.csv",rows)
        data=compose_table(rows)
        write_compose(folder/"compose",data)
        restored=read_compose(folder/"compose")
        for field in ("nb","Q1","Q3","Q4","Q5","Q6","Q7","Yn","Yp","Ye","Ymu"):
            if not np.allclose(getattr(data,field),getattr(restored,field),rtol=5e-12,atol=5e-14):
                raise RuntimeError(f"native readback changed {label}.{field}")
        for field in ("Iphase","Yn_present","Yp_present","Ye_present","Ymu_present",
                      "C3_applicable","C4_applicable"):
            if not np.array_equal(getattr(data,field),getattr(restored,field)):
                raise RuntimeError(f"native readback changed {label}.{field}")
        if "Gibbs" in label:
            expected_phase=np.where(data.nb<=gibbs_phase["left"]["n_b"],6,
                np.where(data.nb>=gibbs_phase["right"]["n_b"],7,8))
            if not np.array_equal(restored.Iphase,expected_phase):
                raise RuntimeError("Gibbs phase interpretation disagrees with density endpoints")
        if not np.all(restored.C3_applicable) or not np.all(restored.C4_applicable):
            raise RuntimeError("generated APR composition closure was unexpectedly masked")
        report={"representation":label,"species":"npemu","points":len(rows),
                "native_composition_readback":"passed",
                "native_phase_counts":{str(int(code)):int(np.count_nonzero(restored.Iphase==code))
                                       for code in np.unique(restored.Iphase)},
                "phase_code_meaning":{"6":"homogeneous low branch","7":"homogeneous high branch","8":"two-charge Gibbs mixture"},
                "native_phase_presence_and_applicability_readback":"passed",
                "closure_applicable_counts":{"C3":int(np.count_nonzero(restored.C3_applicable)),
                                             "C4":int(np.count_nonzero(restored.C4_applicable))},
                "native_checks":check_consistency(restored,verbose=False),
                "physical_domain":physical_domain_report(restored),
                "max_charge_residual_fm3":max(abs(r["charge_residual_fm3"]) for r in rows),
                "max_euler_residual_mev_fm3":max(abs(r["euler_residual_mev_fm3"]) for r in rows)}
        (folder/"RECEIPT.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
        summaries[label]={key:report[key] for key in ("points","native_composition_readback",
                        "max_charge_residual_fm3","max_euler_residual_mev_fm3")}
    comparison=[]
    for index,state,gibbs_state in zip(indices,candidate,gibbs_candidate):
        n=float(source.nb[index])
        item={"source_row_1based":int(index+1),"n_b_fm3":n,
              "provider_pressure_mev_fm3":float(source.p[index]),
              "provider_epsilon_mev_fm3":float(source.eps[index]),
              "provider_mu_b_mev":float(source.mub[index]),
              "provider_composition_available":False,"generated_phase":state["branch"]}
        for label,generated in (("low",apr.beta_state(n,"low")),
                                ("high",apr.beta_state(n,"high")),("maxwell",state),("gibbs",gibbs_state)):
            item[label]={"P_relative_difference":generated["pressure"]/float(source.p[index])-1,
                         "epsilon_relative_difference":generated["epsilon"]/float(source.eps[index])-1,
                         "mu_b_difference_mev":generated["mu_n"]-float(source.mub[index]),
                         "generated_Yp":generated["Y_p"]}
        comparison.append(item)
    mixture=[apr.equilibrium_state(float(n),phase)
             for n in np.linspace(phase_left,phase_right,33)[1:-1]
             if density[0]<=n<=density[-1]]
    (args.output/"mixed_states.json").write_text(json.dumps(mixture,indent=2,allow_nan=False)+"\n")
    receipt={"schema_version":1,"model_id":"APR_APR","sources":apr.SOURCE,
        "head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "generator_sha256":hashlib.sha256(Path(apr.__file__).read_bytes()).hexdigest(),
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "dependency_sha256":{name:hashlib.sha256((ROOT/"python_frontend/compose_eos"/name).read_bytes()).hexdigest()
                             for name in ("cold_leptons.py","physical_domain.py","apr_gibbs.py")},
        "runtime":{"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__},
        "source_input_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(args.source.glob("eos.*")) if p.is_file()},
        "parameters_p1_to_p21":list(apr.PARAMETERS),
        "nucleon_masses_mev":{"neutron":apr.M_N,"proton":apr.M_P},
        "mass_convention":"explicit split nucleon masses; author SNA proton constant differs by 0.0000600 MeV",
        "source_rows_compared":len(comparison),"source_density_interval_fm3":[float(density[0]),float(density[-1])],
        "lower_bound_status":"uniform-core comparison window chosen for this run; provider phase labels unavailable",
        "phase_construction":phase,"gibbs_phase_construction":gibbs_phase,
        "provider_phase_evidence":gibbs.PROVIDER_PHASE_SOURCE,"runs":summaries,
        "causality_crossing_refinement":crossings,
        "source_comparison":comparison,
        "max_absolute_relative_comparison":{
            label:{quantity:max(abs(row[label][quantity]) for row in comparison)
                for quantity in ("P_relative_difference","epsilon_relative_difference")}
            for label in ("low","high","maxwell","gibbs")},
        "source_equivalence_status":"PARTIAL; provider Gibbs prescription identified and implemented in bulk limit, numerical table mismatch remains unresolved",
        "native_Maxwell_representation":"pure-phase endpoints retained only inside requested window; any density jump marked by phase6-to7; mixed states only in JSON because one global local charge potential is undefined",
        "native_Gibbs_representation":"phase6 low, phase8 two-charge bulk Gibbs, phase7 high; common charge chemical potential and averaged species densities retained",
        "scope":"cold uniform APR functional with composition; no crust, no joining, no stellar solve",
        "caveats":["local-neutral Maxwell is distinct from globally neutral Gibbs",
                   "Gibbs calculation omits finite-size surface and Coulomb energies",
                   "minorants tested on sampled branch states; global phase minimum not proven",
                   "high-density acausality retained and localized, not a model-wide rejection",
                   "source table mismatch is reported without parameter tuning"]}
    (args.output/"GENERATION_RECEIPT.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"runs":summaries,"coexistence":{"n_low":phase_left,"n_high":phase_right,
          "P_mev_fm3":phase["pressure"],"mu_b_mev":phase["mu_b"]},
          "max_relative_source_difference":receipt["max_absolute_relative_comparison"],
          "causality_crossings":receipt["causality_crossing_refinement"]},indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
