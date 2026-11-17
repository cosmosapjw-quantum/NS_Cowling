from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'python_frontend/tests')]
from test_stellar_sharp import contract
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background
eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(contract()))
rows=[]
for dn in [1e-6,1e-8,1e-10]:
 try:
  bg=solve_sharp_background(eos,5.+dn)
  rows.append(dict(density_offset=dn,status='RETURNED',h_c_minus_transition=bg.h_c-eos.transition_h[0],h_c_minus_start=bg.h_c-bg.h_start,interfaces=len(bg.interfaces),central_phase_swallowed=bool(bg.h_start<eos.transition_h[0]<bg.h_c)))
 except Exception as exc:rows.append(dict(density_offset=dn,status='REJECTED',error=str(exc)))
result=dict(source_sha256=hashlib.sha256((ROOT/'python_frontend/compose_eos/stellar_sharp.py').read_bytes()).hexdigest(),observations=rows)
Path(__file__).with_name('CENTER_INTERFACE_PROBE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
