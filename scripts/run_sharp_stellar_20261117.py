"""Reproduce G's scoped sharp-star calculation and explicit sensitivity controls."""
from pathlib import Path
import hashlib,json,os,sys,time,traceback
from dataclasses import asdict
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,scan_brackets,solve_mode
from compose_eos.stellar_cowling import CowlingSettings

OUT=ROOT/"data/crust_refinement_20261117/stellar"
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
def reference_star(contract,r0=.001,rtol=1e-10,atol=1e-12,step=.002):
    eos=SharpContractEnthalpyEOS.from_path(contract)
    cfg=dict(r0_km=r0,rtol=rtol,atol=atol,max_step_h=step)
    sequence=[solve_sharp_background(eos,n,**cfg).summary() for n in [.08,.16,.32,.5,.8,1.0]]
    brackets=[(a["central_density_fm3"],b["central_density_fm3"]) for a,b in zip(sequence[:-1],sequence[1:])
              if a["mass_msun"]<1.4<b["mass_msun"]]
    if not brackets:raise ValueError("No sampled ascending 1.4-solar-mass bracket in supported density range")
    lo,hi=brackets[0]
    n=brentq(lambda n:solve_sharp_background(eos,n,**cfg).mass_msun-1.4,lo,hi,xtol=1e-10)
    bg=solve_sharp_background(eos,n,**cfg)
    return bg,dict(sequence=sequence,selected_mass_bracket=[lo,hi],selection="first sampled ascending mass bracket; no radial-stability conclusion",background=bg.summary())

def main():
    started=time.time()
    manifest=json.loads((OUT/"FULL_CONTRACT_ADMISSION.json").read_text())
    entries={x["id"]:x for x in manifest["inputs"]}
    selection=[("G_TM1_CANONICAL_N064","primary"),("G_TM1_CANONICAL_N032","primary"),
               ("G_TM1_CANONICAL_N064","next_admissible_join_control")]
    settings=CowlingSettings(rtol=3e-8,atol=3e-11,max_step_h=.001)
    frequencies=[50,100,150,200,300,400,500,650,800,1000,1250,1500,1750,2000,
                 2250,2500,3000,3500,4000,4500,5000,5500,6000,7000,8000,9000,10000]
    source_paths=[ROOT/"python_frontend/compose_eos"/f for f in
                  ["stellar_sharp.py","stellar_background.py","stellar_cowling.py","whole_crust.py"]]
    source_paths+=[Path(__file__),ROOT/"scripts/prepare_sharp_stellar_20261117.py",
                   OUT/"INPUTS.json",OUT/"FULL_CONTRACT_ADMISSION.json"]
    before={str(p.relative_to(ROOT)):sha(p) for p in source_paths}
    result=dict(scope="Conditional perfect-fluid, equilibrium bulk Cowling; slow and rapid interface limits; detected roots in finite scan only",
                mass_target_msun=1.4,settings=asdict(settings),source_hashes_before=before,cases=[])
    main_roots={}
    for index,(key,label) in enumerate(selection):
        admitted=next(a for a in entries[key]["admitted"] if a["label"]==label)
        path=ROOT/admitted["contract"];case_id=key+"__"+label
        bg,record=reference_star(path)
        record.update(id=case_id,contract=admitted["contract"],contract_sha256=sha(path),modes={})
        save(OUT/"calculations"/case_id/"BACKGROUND.json",record)
        print("BACKGROUND",case_id,bg.mass_msun,bg.radius_km,flush=True)
        for conversion in ("slow","rapid"):
            if index==0:
                brackets,scan=scan_brackets(bg,frequencies,settings,conversion=conversion)
            else:
                scan=[];brackets=[]
                for oldroot in main_roots[conversion]:
                    local=np.linspace(.90*oldroot,1.10*oldroot,9)
                    found,values=scan_brackets(bg,local,settings,conversion=conversion)
                    scan.extend(values)
                    if len(found)!=1:
                        raise RuntimeError(f"Expected one controlled nearby root for {case_id} {conversion}: {found}")
                    brackets.extend(found)
            save(OUT/"calculations"/case_id/(conversion+"_SCAN.json"),
                 dict(scan=scan,brackets=brackets,complete_spectrum_claim=False))
            modes=[]
            for k,bracket in enumerate(brackets):
                mode=solve_mode(bg,bracket,settings,conversion=conversion)
                save(OUT/"calculations"/case_id/(conversion+f"_ROOT_{k:02d}.json"),mode)
                summary={k:v for k,v in mode.items() if k!="profile"}
                modes.append(summary)
                print("ROOT",case_id,conversion,mode["frequency_hz"],
                      mode["W_sign_changes_including_interfaces"],flush=True)
            record["modes"][conversion]=modes
            if index==0:main_roots[conversion]=[m["frequency_hz"] for m in modes]
        result["cases"].append(record)
        save(OUT/"STELLAR_RESULTS.json",result)

    # Tighter controls are a numerical sensitivity, distinct from mesh/lower-model changes.
    main=result["cases"][0]
    fine_eos=SharpContractEnthalpyEOS.from_path(ROOT/main["contract"])
    fine_bg=solve_sharp_background(fine_eos,main["background"]["central_density_fm3"],
               r0_km=.0005,rtol=2e-11,atol=2e-13,max_step_h=.001)
    fine_settings=CowlingSettings(rtol=3e-9,atol=3e-12,max_step_h=.0005,surface_h=3e-8)
    refinement=dict(background=fine_bg.summary(),settings=asdict(fine_settings),modes={},
                    varied="r0/2, tighter background and mode tolerances, step/2, surface_h*0.3 at fixed central density")
    for conversion in ("slow","rapid"):
        values=[]
        for k,frequency in enumerate(main_roots[conversion]):
            mode=solve_mode(fine_bg,[frequency*.985,frequency*1.015],fine_settings,conversion=conversion)
            save(OUT/"calculations"/"numerical_refinement"/(conversion+f"_ROOT_{k:02d}.json"),mode)
            values.append({k:v for k,v in mode.items() if k!="profile"})
            print("REFINED_ROOT",conversion,mode["frequency_hz"],flush=True)
        refinement["modes"][conversion]=values
    result["numerical_refinement"]=refinement
    result["source_hashes_after"]={str(p.relative_to(ROOT)):sha(p) for p in source_paths}
    result["source_bytes_unchanged"]=result["source_hashes_after"]==before
    result["wall_time_seconds"]=time.time()-started
    if not result["source_bytes_unchanged"]:raise RuntimeError("Production source changed during calculation")
    save(OUT/"STELLAR_RESULTS.json",result)
    print("STELLAR_FINISHED",result["wall_time_seconds"],flush=True)

if __name__=="__main__":main()
