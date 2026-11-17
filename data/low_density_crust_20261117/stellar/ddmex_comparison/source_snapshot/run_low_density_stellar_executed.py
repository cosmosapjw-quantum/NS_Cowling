"""H: finite, preselected low-density stellar comparisons with immutable receipts."""
from pathlib import Path
import argparse,hashlib,json,sys,time,traceback
from dataclasses import asdict
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main(plan_path,output):
    plan_path=Path(plan_path);out=Path(output);plan=json.loads(plan_path.read_text())
    physics=Path(plan.get("physics_root",ROOT))
    sys.path.insert(0,str(physics/"python_frontend"))
    from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,scan_brackets,solve_mode
    from compose_eos.stellar_cowling import CowlingSettings
    files=[physics/"python_frontend/compose_eos"/s for s in
        ["stellar_sharp.py","stellar_cowling.py","stellar_background.py","whole_crust.py"]]
    files += [Path(__file__),plan_path]+[ROOT/c["contract"] for c in plan["cases"]]
    before={str(p):sha(p) for p in files}
    settings=CowlingSettings(rtol=3e-9,atol=3e-12,max_step_h=.001,root_xtol_hz=1e-8)
    tight=CowlingSettings(rtol=3e-10,atol=3e-13,max_step_h=.0005,root_xtol_hz=1e-10)
    acceptance=plan["acceptance"]
    def passes(m):
        return bool(m["root_converged"] and abs(m["matching_residual"])<acceptance["determinant_abs"]
          and m["relative_join_mismatch"]<acceptance["join_relative"]
          and m["surface_BC_relative"]<acceptance["surface_relative"]
          and all(x["traction_relative_residual"]<acceptance["traction_relative"] and
                  abs(x["kinematic_residual"])<acceptance["kinematic_abs"] for x in m["interface_junctions"]))
    def solve_checked(bg,bracket,cid,conversion,k,cfg):
        dest=out/"calculations"/cid
        mode=solve_mode(bg,bracket,cfg,conversion=conversion)
        original=dest/(conversion+f"_ROOT_{k:02d}.json");save(original,mode)
        used=original;repair=None
        if not passes(mode):
            repair=dest/(conversion+f"_ROOT_{k:02d}_TIGHT.json")
            mode=solve_mode(bg,bracket,tight,conversion=conversion)
            save(repair,mode);used=repair
        result={key:v for key,v in mode.items() if key!="profile"}
        result.update(acceptance="PASS" if passes(mode) else "FAIL",
                      file=str(used.relative_to(ROOT)),file_sha256=sha(used),
                      bounded_tighter_repair=repair is not None)
        print("ROOT",cid,conversion,mode["frequency_hz"],result["acceptance"],flush=True)
        return result
    def reference_star(contract):
        eos=SharpContractEnthalpyEOS.from_path(contract)
        cfg=dict(r0_km=.001,rtol=1e-10,atol=1e-12,max_step_h=.002)
        seq=[solve_sharp_background(eos,n,**cfg).summary() for n in [.08,.16,.32,.5,.8,1.0]]
        br=[(a["central_density_fm3"],b["central_density_fm3"]) for a,b in zip(seq[:-1],seq[1:])
            if a["mass_msun"]<plan["mass_msun"]<b["mass_msun"]]
        if not br:raise ValueError("No sampled ascending target-mass bracket in source support")
        nc=brentq(lambda n:solve_sharp_background(eos,n,**cfg).mass_msun-plan["mass_msun"],*br[0],xtol=1e-10)
        bg=solve_sharp_background(eos,nc,**cfg)
        bar=eos.barotrope
        hl=float(np.log(bar.bridge.mu_left/eos.a));hr=float(np.log(bar.bridge.mu_right/eos.a))
        rl=float(bg.solution.sol(hl)[0]);rr=float(bg.solution.sol(hr)[0])
        return bg,dict(sequence=seq,selected_mass_bracket=br[0],selection="first sampled ascending bracket; no stability inference",
             background=bg.summary(),lower_bridge_radial_width_km=rl-rr,
             outer_and_surface_radial_width_km=bg.radius_km-rl,
             sharp_crust_radial_width_km=bg.radius_km-bg.interfaces[0]["r_km"] if bg.interfaces else None)
    start=time.time();main_roots={};cases=[]
    result=dict(status="RUNNING",scope=plan["scope"],plan=str(plan_path.relative_to(ROOT)),plan_sha256=sha(plan_path),
      mass_target_msun=plan["mass_msun"],settings=asdict(settings),source_hashes_before=before,cases=cases)
    primary_bg=None
    for index,entry in enumerate(plan["cases"]):
        path=ROOT/entry["contract"]
        if sha(path)!=entry["contract_sha256"]:raise ValueError("contract changed after plan")
        bg,record=reference_star(path);record.update(**entry,modes={})
        if index==0:primary_bg=bg
        cases.append(record);save(out/"STELLAR_RESULTS.json",result)
        save(out/"calculations"/entry["id"]/"BACKGROUND.json",record)
        print("BACKGROUND",entry["id"],bg.mass_msun,bg.radius_km,flush=True)
        for conversion in ("slow","rapid"):
            if index==0:
                brackets,scan=scan_brackets(bg,plan["frequencies_hz"],settings,conversion=conversion)
                detected=brackets[:]
                brackets=brackets[:plan["max_roots"][conversion]]
            else:
                scan=[];brackets=[];detected=[]
                for frequency in main_roots[conversion]:
                    found,rows=scan_brackets(bg,np.linspace(.88*frequency,1.12*frequency,7),settings,conversion=conversion)
                    scan.extend(rows);detected.extend(found)
                    if len(found)!=1:
                        raise RuntimeError(f"Continuation for {entry['id']} {conversion}: require exactly one root near {frequency}; {found}")
                    brackets.extend(found)
            save(out/"calculations"/entry["id"]/(conversion+"_SCAN.json"),
                 dict(scan=scan,all_detected_brackets=detected,selected_brackets=brackets,complete_spectrum_claim=False))
            record["modes"][conversion]=[solve_checked(bg,b,entry["id"],conversion,k,settings) for k,b in enumerate(brackets)]
            if index==0:main_roots[conversion]=[v["frequency_hz"] for v in record["modes"][conversion]]
            save(out/"STELLAR_RESULTS.json",result)
    fine_bg=solve_sharp_background(primary_bg.eos,primary_bg.central_density_fm3,
                r0_km=.0005,rtol=2e-11,atol=2e-13,max_step_h=.001)
    fine_cfg=CowlingSettings(rtol=3e-10,atol=3e-13,max_step_h=.0005,surface_h=3e-8,root_xtol_hz=1e-10)
    refinement=dict(background=fine_bg.summary(),settings=asdict(fine_cfg),fixed_central_density=True,modes={})
    for conv in ("slow","rapid"):
        refinement["modes"][conv]=[solve_checked(fine_bg,[.985*f,1.015*f],"numerical_refinement",conv,k,fine_cfg)
                                   for k,f in enumerate(main_roots[conv])]
    result["numerical_refinement"]=refinement
    result["source_hashes_after"]={str(p):sha(p) for p in files}
    result["source_bytes_unchanged"]=before==result["source_hashes_after"]
    allm=[m for c in cases+[refinement] for modes in c["modes"].values() for m in modes]
    result["status"]="PASS_SCOPED" if result["source_bytes_unchanged"] and all(m["acceptance"]=="PASS" for m in allm) and allm else "REJECTED_OR_EMPTY"
    result["wall_time_seconds"]=time.time()-start
    save(out/"STELLAR_RESULTS.json",result)
    print("FINISHED",result["status"],result["wall_time_seconds"],flush=True)
    return 0 if result["status"]=="PASS_SCOPED" else 1
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("plan");parser.add_argument("output")
    args=parser.parse_args();raise SystemExit(main(args.plan,args.output))
