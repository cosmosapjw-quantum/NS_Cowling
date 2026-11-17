"""Final receipt only; no new physical calculation or numerical audit."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'data/crust_refinement_20261117';R=Path(__file__).parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
before=json.loads((BASE/'verification/CODE_BEFORE_FULL_SUITE.json').read_text())
after=json.loads((BASE/'verification/CODE_AFTER_FULL_SUITE.json').read_text())
assert before==after and all(sha(ROOT/p)==v for p,v in after.items())
docs=[ROOT/'docs/crust_refinement_20261117'/p for p in ['INDEPENDENT_REVIEW_KO.md','RESULTS_KO.md','SHARP_STELLAR_KO.md','TF_RESULTS_KO.md','SUPPORT_RESULTS_KO.md','SOURCE_PHYSICS_KO.md']]
assert all(not any(v<32 and v not in (10,) for v in p.read_bytes()) for p in docs)
evidence=[BASE/'stellar/MODE_ADMISSION.json',BASE/'stellar/STELLAR_RESULTS.json',BASE/'stellar/CENTER_GUARD_REPAIR.json',BASE/'verification/FULL_SUITE_RECEIPT.json',BASE/'verification/NATIVE_BUILD_RECEIPT.json',BASE/'G_SUMMARY.json']
receipt=dict(schema='ns_cowling.G.independent_review.v1',status='ACCEPT_SCOPED',fresh_review_passes=1,
 completed_utc=datetime.now(timezone.utc).isoformat(),base_commit='d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a',
 blocking_findings_open=0,closed_findings=['center Taylor crossing a sharp interface','three original mode residuals above gate','TF report escaped equations'],
 final_code_files_matching_full_suite_snapshot=len(after),
 final_source_sha256={'python_frontend/compose_eos/stellar_sharp.py':sha(ROOT/'python_frontend/compose_eos/stellar_sharp.py'),'python_frontend/compose_eos/tf_candidate.py':sha(ROOT/'python_frontend/compose_eos/tf_candidate.py')},
 document_sha256={str(p.relative_to(ROOT)):sha(p) for p in docs},
 production_evidence_sha256={str(p.relative_to(ROOT)):sha(p) for p in evidence},
 review_evidence_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(R.iterdir()) if p.is_file() and p.name!='REVIEW_RECEIPT.json'},
 publication='F published with exact remote ref receipt; G publication requires separate user approval',
 scope='Finite computed branches, declared full-EOS closure, recorded numerical root gates; no global minimum/continuum/spectral-completeness claim')
(R/'REVIEW_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'report_sha256':sha(docs[0]),'receipt_sha256':sha(R/'REVIEW_RECEIPT.json'),'source_files':len(after),'blocking_findings_open':0},indent=2))
