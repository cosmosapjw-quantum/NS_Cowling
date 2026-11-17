from pathlib import Path
import sys,json,hashlib,ast,copy
ROOT=Path(__file__).resolve().parents[3];sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'python_frontend/tests')]
from test_stellar_sharp import contract
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background
new=ROOT/'python_frontend/compose_eos/stellar_sharp.py'
old=ROOT/'data/crust_refinement_20261117/stellar/production_source/stellar_sharp.py'
tree=ast.parse(new.read_text())
fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='solve_sharp_background')
guard=next(n for n in fn.body if isinstance(n,ast.If) and 'h0 <= ht < hc' in ast.unparse(n.test))
fn.body.remove(guard)
assert ast.dump(tree,include_attributes=False)==ast.dump(ast.parse(old.read_text()),include_attributes=False)
eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(contract()))
rejected=[]
for dn in [1e-8,1e-10]:
 try:solve_sharp_background(eos,5.+dn,r0_km=.001)
 except ValueError as e:
  assert 'smaller r0' in str(e);rejected.append(dict(density_offset=dn,error=str(e)))
 else:raise AssertionError('unresolved interface accepted')
resolved=solve_sharp_background(eos,5.+1e-8,r0_km=.0001)
assert len(resolved.interfaces)==1 and resolved.interfaces[0]['r_km']>.0001
result=dict(status='PASS',source_sha256=hashlib.sha256(new.read_bytes()).hexdigest(),
 original_source_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
 only_added_guard_AST=True,rejected=rejected,resolved_smaller_center_radius_interfaces=resolved.interfaces,
 scope='One bounded repair recheck; original production arithmetic unchanged outside newly rejected edge inputs')
Path(__file__).with_name('CENTER_GUARD_CLOSEOUT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
