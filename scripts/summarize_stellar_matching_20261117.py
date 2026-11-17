#!/usr/bin/env python3
"""Apply the exact D APR delta to the immutable C 390-cell census."""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args()
    root=args.root.resolve()
    base=root/"data/tf_pair_expansion_20261117/matching/campaign_v1/MATCHING_390_MAP.json"
    delta=root/"data/stellar_consumers_20261117/apr/MATCHING_390_DELTA.json"
    evidence=root/"data/stellar_consumers_20261117/apr/CELLS.json"
    original=json.loads(base.read_text())
    updates=json.loads(delta.read_text())
    cells=original["cells"]
    assert len(cells)==390 and len({c["id"] for c in cells})==390
    by_id={c["id"]:c for c in cells}
    assert len(updates["cells"])==25 and len({c["cell"] for c in updates["cells"]})==25
    for change in updates["cells"]:
        item=by_id[change["cell"]]
        assert item["status"]==change["from"]=="DEFERRED_APR_GIBBS_SUPPORT"
        assert change["native_matching_exit_code"]==0
        item.update(eligible=True,kernel_attempt=True,constructed_readback=True,
                    method_admitted=False,reused_evidence=False,status=change["to"],
                    reason="D native low-branch match and retained APR Gibbs suffix readback")
        item["evidence_variants"].append({
            "path":str(evidence.relative_to(root)),"sha256":sha(evidence),
            "scope":"Retained 6->8->7 finite core assembly; full assembled method/stellar admission not claimed",
            "stage":"D","unique_cell_id":item["id"]})
        item["D_delta"]=change
    summary={
        "requested":len(cells),
        "finite_component_eligible":sum(bool(c["eligible"]) for c in cells),
        "unique_cell_kernel_attempts":sum(bool(c["kernel_attempt"]) for c in cells),
        "unique_cell_constructed_readback":sum(bool(c["constructed_readback"]) for c in cells),
        "unique_cell_method_admitted":sum(bool(c["method_admitted"]) for c in cells),
        "deferred":sum(c["status"].startswith("DEFERRED") for c in cells),
        "by_status":dict(Counter(c["status"] for c in cells)),
        "D_new_unique_cells":25,"D_native_matching_calls_observed":50,
        "D_matching_call_meaning":"25 initial +25 after consumer/diagnostic repair; not 50 unique cells",
        "D_low_branch_CSS_admissions_separate":5,
        "D_full_assembled_APR_method_admissions_claimed":0,
        "D_full_APR_stellar_models_computed":0}
    assert (summary["unique_cell_kernel_attempts"],summary["unique_cell_constructed_readback"],
            summary["unique_cell_method_admitted"],summary["deferred"])==(190,130,25,200)
    result={
        "schema":"ns-cowling.stage-D-matching-union.v1",
        "authority":"C immutable 390-cell map plus D APR cell-ID delta",
        "base_C_commit":"539e6bd8de73f8654d121c0bc30023fb868dc05b",
        "scope":"Finite native-component matching support. These counts do not mean complete whole-star EOS coverage.",
        "source_policy":"Existing C source variants retained, including SPG public-composition-inferred native variant. Separate public provider families are outside these390 IDs.",
        "historical_field_scope":"Unchanged cell reused_evidence values describe their C provenance. D_delta marks only new D cells.",
        "source_hashes":{str(p.relative_to(root)):sha(p) for p in [base,delta,evidence]},
        "summary":summary,"cells":cells}
    out=root/"data/stellar_consumers_20261117/MATCHING_390_MAP.json"
    out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
if __name__=="__main__":
    main()
