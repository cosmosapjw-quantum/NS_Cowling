from pathlib import Path
import json,hashlib
r=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117");b=r/"data/low_density_crust_20261117/verification"
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
before=load(b/"CODE_BEFORE_FULL_SUITE.json");after=load(b/"CODE_AFTER_FULL_SUITE.json");rec=load(b/"FULL_SUITE_RECEIPT.json")
print(type(after).__name__,list(after)[:3])
if all(isinstance(v,str) for v in after.values()):
 changed=[p for p,h in after.items() if sha(r/p)!=h]
else:changed=["schema requires inspection"]
log=(b/"FULL_SUITE.log").read_text().splitlines()[-1]
ok=before==after and not changed and rec["actual_exit_code"]==0 and "1453 passed, 9 skipped, 4 subtests passed" in log
out=dict(status="PASS" if ok else "FAIL",source_count=len(after),current_changed=changed,source_bytes_before_after_equal=before==after,actual_exit_code=rec["actual_exit_code"],stdout_last_line=log,evidence_sha256={str(p.relative_to(r)):sha(p) for p in [b/"CODE_BEFORE_FULL_SUITE.json",b/"CODE_AFTER_FULL_SUITE.json",b/"FULL_SUITE_RECEIPT.json",b/"FULL_SUITE.log"]})
(r/"data/low_density_crust_20261117/review/FULL_SUITE_REVIEW.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2));raise SystemExit(not ok)
