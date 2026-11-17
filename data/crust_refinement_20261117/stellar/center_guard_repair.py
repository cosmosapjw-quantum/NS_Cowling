from pathlib import Path
import ast,hashlib,json,shutil,subprocess,os
r=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117")
o=r/"data/crust_refinement_20261117/stellar"
receipt=json.loads((o/"PRODUCTION_RECEIPT.json").read_text())
assert receipt["actual_exit_code"]==0
results=json.loads((o/"STELLAR_RESULTS.json").read_text())
p=r/"python_frontend/compose_eos/stellar_sharp.py"
old=p.read_text()
oldsha=hashlib.sha256(p.read_bytes()).hexdigest()
assert oldsha==results["source_hashes_before"][str(p.relative_to(r))]
snapshot=o/"production_source";snapshot.mkdir(exist_ok=True)
(snapshot/"stellar_sharp.py").write_text(old)
anchor='    y=np.array([r0_km,4*np.pi*ec*r0_km**3/3])'
guard='    if any(h0 <= ht < hc for ht in eos.transition_h):\n        raise ValueError("center Taylor start crosses an interface; use a smaller r0_km")\n'
assert old.count(anchor)==1 and guard not in old
new=old.replace(anchor,guard+anchor)
new_ast=ast.parse(new)
fn=next(n for n in new_ast.body if isinstance(n,ast.FunctionDef) and n.name=="solve_sharp_background")
def is_guard(node):
    return isinstance(node,ast.If) and len(node.body)==1 and isinstance(node.body[0],ast.Raise) and "center Taylor start crosses" in ast.dump(node)
assert sum(is_guard(n) for n in fn.body)==1
fn.body=[n for n in fn.body if not is_guard(n)]
assert ast.dump(new_ast)==ast.dump(ast.parse(old))
backgrounds=[c["background"] for c in results["cases"]]+[results["numerical_refinement"]["background"]]
guard_checks=[]
for bg in backgrounds:
    hs=[j["h"] for j in bg["interfaces"]]
    guard_checks.append(dict(central_density=bg["central_density_fm3"],
        h_c=bg["h_c"],h_start=bg["h_start"],interfaces_h=hs,
        rejected_by_new_guard=any(bg["h_start"]<=h<bg["h_c"] for h in hs)))
assert not any(x["rejected_by_new_guard"] for x in guard_checks)
p.write_text(new)
repaired=dict(original_production_source_sha256=oldsha,
 final_source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
 original_snapshot=str((snapshot/"stellar_sharp.py").relative_to(r)),
 whole_module_AST_equal_after_removing_only_new_guard=True,
 production_inputs_guard_checks=guard_checks,
 scope="Postproduction rejection of unresolved center-interface starts; all previously computed production inputs remain on identical numerical code paths. No claim of rerunning production.")
(o/"CENTER_GUARD_REPAIR.json").write_text(json.dumps(repaired,indent=2)+"\n")
for name in ["g_low_root_tight_control.py","g_low_root_tight_wrapper.py"]:
    shutil.copyfile(r.parent/"scratch"/name,o/name)
doc=r/"docs/crust_refinement_20261117/SHARP_STELLAR_KO.md"
s=doc.read_text()
s=s.replace(r"(h=\int_0^P dP/(\epsilon+P)=\log(\mu/\mu_{\rm vac}))",r"\(h=\int_0^P dP/(\epsilon+P)=\log(\mu/\mu_{\rm vac})\)")
s=s.replace(r"(n,\epsilon)이",r"\(n,\epsilon\)이")
s+=r"""
### 중심 근접 interface의 경계 입력 수리

독립 검토에서 \(h_0\le h_t<h_c\)이면 중심 Taylor 영역 안에 interface가 숨는 입력을 확인했다.
이제 solver는 그 입력을 거절하고 더 작은 \(r_0\)를 요구한다.
최초 regression은 실제 실패했고, 수리 뒤 같은 상태에서 거절과 작은 \(r_0\)의 정상 interface 검출을 각각 확인한다.
이미 실행한 생산 사례에는 새 guard가 모두 false이므로 수치 경로가 같다.
원 생산 source snapshot과 guard만 제거했을 때 전체 module AST가 일치하는 검산을
CENTER_GUARD_REPAIR.json에 보존했다. 전체 생산을 새 코드로 다시 했다고 주장하지 않는다.
"""
doc.write_text(s)
env=dict(os.environ,PYTHONPATH=str(r/"python_frontend"),OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
argv=["/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python","-m","pytest","-q","python_frontend/tests/test_stellar_sharp.py"]
run=subprocess.run(argv,cwd=r,env=env,capture_output=True,text=True)
(o/"CENTER_GUARD_GREEN.log").write_text(run.stdout+run.stderr)
(o/"CENTER_GUARD_GREEN.json").write_text(json.dumps(dict(argv=argv,actual_exit_code=run.returncode),indent=2)+"\n")
print(run.stdout);print("ACTUAL_EXIT",run.returncode)
raise SystemExit(run.returncode)
