"""G: deterministic full-contract admission with explicitly modeled lower closure."""
from pathlib import Path
import argparse,copy,hashlib,json,os,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.api import read_compose
from compose_eos.whole_crust import WholeCrustBarotrope,HermiteInterval

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run(inputs):
    reference=ROOT/"data/stellar_extension_20261117/eos/fine/REPRESENTATION.json"
    base=json.loads(reference.read_text());bar=WholeCrustBarotrope.from_contract(base)
    left=dict(n=float(bar.outer.nb[-1]),p=float(bar.outer.p[-1]),
              epsilon=float(bar.outer.nb[-1]*bar.outer.F[-1]))
    output=ROOT/"data/crust_refinement_20261117/stellar"
    manifest=[]
    for entry in inputs:
        path=ROOT/entry["compose_path"]
        data=read_compose(path)
        result=json.loads((path.parent/"RESULT.json").read_text())
        native=result.get("native_report") or result["v2_report"]["native_report"]
        interface=int(native["low_phase_rows"])-1
        record=dict(**entry,source_files={p.name:sha(p) for p in sorted(path.glob("*")) if p.is_file()},
                    original_interface_index=interface,lower_join_trials=[],admitted=[])
        for j in range(interface):
            trial=dict(source_row=j,n=float(data.nb[j]),p=float(data.p[j]),epsilon=float(data.eps[j]))
            try:
                bridge=HermiteInterval(left["n"],left["p"],left["epsilon"],data.nb[j],data.p[j],data.eps[j])
                trial["bridge_validation"]=bridge.validate()
                if trial["bridge_validation"]["status"]!="PASS":
                    trial["status"]="BRIDGE_REJECTED"
                else:
                    full=dict(surface=copy.deepcopy(base["surface"]),outer=copy.deepcopy(base["outer"]),
                        lower_bridge=bridge.contract(),upper=dict(n=data.nb[j:].tolist(),p=data.p[j:].tolist(),
                        epsilon=data.eps[j:].tolist(),sharp_interfaces=[interface-j]))
                    WholeCrustBarotrope.from_contract(full)
                    # This join is chosen from source order before computing any stellar observables.
                    label="primary" if not record["admitted"] else "next_admissible_join_control"
                    dest=output/"contracts"/entry["id"]/label/"REPRESENTATION.json"
                    dest.parent.mkdir(parents=True,exist_ok=True)
                    full["provenance"]=dict(source_compose=entry["compose_path"],source_row_start=j,
                        outer_reference=str(reference.relative_to(ROOT)),outer_sha256=sha(reference),
                        choice="first thermodynamically admissible source knot, then next admissible knot",
                        source_values_preserved=True,lower_bridge_is_explicit_phenomenological_closure=True,
                        geometry_global_minimum_proven=False)
                    dest.write_text(json.dumps(full,indent=2)+"\n")
                    trial["status"]="ADMITTED"
                    record["admitted"].append(dict(label=label,source_row_start=j,
                        retained_crust_rows=interface-j+1,contract=str(dest.relative_to(ROOT)),sha256=sha(dest)))
            except (ValueError,RuntimeError) as error:
                trial.update(status="REJECTED",exception_type=type(error).__name__,error=str(error))
            record["lower_join_trials"].append(trial)
            if len(record["admitted"])==2:break
        manifest.append(record)
        print(entry["id"],[(v["label"],v["source_row_start"]) for v in record["admitted"]],flush=True)
    report=dict(scope="Full piece-contract admission under inherited E outer/surface and new explicit lower bridge; no microscopic bridge claim.",
                reference=str(reference.relative_to(ROOT)),reference_sha256=sha(reference),left_state=left,inputs=manifest)
    (output/"FULL_CONTRACT_ADMISSION.json").write_text(json.dumps(report,indent=2)+"\n")
    return report
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("inputs");args=parser.parse_args()
    run(json.loads(Path(args.inputs).read_text()))
