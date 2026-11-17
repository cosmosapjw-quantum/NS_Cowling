from pathlib import Path
import ast,json,subprocess,hashlib
ROOT=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117");OUT=ROOT/"data/low_density_crust_20261117/review"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
p="python_frontend/compose_eos/whole_crust.py";old=subprocess.check_output(["git","show","905b845ad0f14d2f1e7f67f9eb802691144850cd:"+p],cwd=ROOT,text=True);new=(ROOT/p).read_text()
def symbols(src):
 return {x.name:x for x in ast.parse(src).body if isinstance(x,(ast.FunctionDef,ast.ClassDef))}
a,b=symbols(old),symbols(new);added=sorted(set(b)-set(a));removed=sorted(set(a)-set(b));changed=[n for n in a if n in b and ast.dump(a[n])!=ast.dump(b[n])]
members={}
for name in changed:
 x={v.name:v for v in a[name].body if isinstance(v,ast.FunctionDef)}
 y={v.name:v for v in b[name].body if isinstance(v,ast.FunctionDef)}
 members[name]=dict(added=sorted(set(y)-set(x)),removed=sorted(set(x)-set(y)),changed=[n for n in x if n in y and ast.dump(x[n])!=ast.dump(y[n])])
test=ROOT/"data/low_density_crust_20261117/thermo_construction/TEST_01.json";receipt=json.loads(test.read_text())
log=(test.parent/receipt["stdout"]).read_text();errors=[]
if sha(ROOT/p)!=receipt["source_sha256"]:errors.append("core differs from tested bytes")
if sha(ROOT/"python_frontend/tests/test_causal_lower_bridge.py")!=receipt["test_sha256"]:errors.append("test differs from executed bytes")
if receipt["returncode"]!=0 or "71 passed" not in log:errors.append("targeted test receipt/log mismatch")
if added!=["CausalExponentialInterval"] or removed or changed!=["WholeCrustBarotrope"]:errors.append("unexpected old symbol change")
if members.get("WholeCrustBarotrope",{}).get("changed")!=["from_contract","_interval_at_pressure"]:errors.append("unexpected consumer change")
report=dict(status="PASS" if not errors else "FAIL",source_sha256=sha(ROOT/p),test_sha256=receipt["test_sha256"],targeted_test_receipt_sha256=sha(test),targeted_test_stdout_sha256=sha(test.parent/receipt["stdout"]),targeted_tests_actual_exit_code=receipt["returncode"],targeted_tests_count=71,ast_delta=dict(added=added,removed=removed,changed=changed,members=members),direct_mathematical_review=[
"causality iff y'=d(n/mu)/dmu>=0 for positive mu,n; strict integral bounds necessary for unequal y endpoints",
"exponential Phi is increasing in x and decreasing in kappa; integral has unique finite endpoint-determined root",
"stable positive/negative/near-zero primitive branches agree with direct integrals",
"susceptibility and causal-margin extrema: 2+mu*kappa/width=0 plus endpoints",
"sound-speed stationary F=y+mu*kappa*A/width has derivative kappa*y/width, hence at most one interior root",
"equal-y luminal branch and explicit saved-kappa verification; unknown/missing opt-in kinds rejected",
"epsilon positivity follows d epsilon/dmu=mu*dn/dmu>0 from positive left endpoint",
"existing Hermite body unchanged; stellar consumers unchanged"],errors=errors)
(OUT/"CORE_REVIEW.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
raise SystemExit(bool(errors))
