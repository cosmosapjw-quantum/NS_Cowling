"""Closeout identity, admitted status and finite scope without repeated reviews."""
from pathlib import Path
import json,hashlib,subprocess
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data/native_dense_interpolation_20261118"
DOC=ROOT/"docs/native_dense_interpolation_20261118"
out=DATA/"aggregated"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
summary=json.loads((out/"SUMMARY.json").read_text())
full=json.loads((out/"FULL_SUITE_RECEIPT.json").read_text())
review=json.loads((out/"INDEPENDENT_AUDIT.json").read_text())
native=json.loads((out/"NATIVE_THERMO_AUDIT.json").read_text())
plots=json.loads((out/"FIGURES_RECEIPT.json").read_text())
assert full["actual_exit_code"]==0 and full["source_bytes_unchanged"]
assert review["status"]=="PASS_SCOPED" and review["failures"]==[]
assert summary["generated_direct_native_state_count"]==25612
assert summary["TF_distinct_new_withheld"]==30
assert native["max_relative_euler"]<5e-10
assert len(plots["figures"])==5 and plots["exclude_linear_from_all_figures"]
p=DOC/"RESULTS_KO.md"
assert p.is_file() and "## TF crust" in p.read_text()
source=[p for p in sorted((ROOT/"scripts").glob("*20261118.py"))]
source+=list((ROOT/"python_frontend"/"tests").glob("*20261118.py"))
assert len(source)>=12
file_hashes={}
for folder in (DATA,DOC):
    for file in sorted(folder.rglob("*")):
        if file.is_file():
            file_hashes[str(file.relative_to(ROOT))]=sha(file)
receipt=dict(schema="ns_cowling.native_dense_interpolation_closeout.v1",
    target_branch="research/native-dense-interpolation-benchmark-20261118",
    parent_commit="a0fa5fb93851a9647030f3589c24a62ddea527b1",
    scope="direct native dense fixed-model holdout, no linear method in figures",
    native_input_count=10,native_state_count=25612,spatial_TF_native_holdout_count=30,
    comparison_families=["RMF fixed parameters","core fixed-source phase","spatial TF N128 fixed droplet"],
    original_TF_native_refused_count=24,
    evidence_grade="SCOPED_NUMERICALLY_CHECKED_FIXED_MODEL",
    plot_count=len(plots["figures"]),plot_method_keys=summary["interpolation_methods_plotted"],
    independent_audit_checks=review["check_count"],independent_findings=len(review["failures"]),
    full_suite_actual_exit=full["actual_exit_code"],source_bytes_unchanged=full["source_bytes_unchanged"],
    source_document_checksum=sha(DOC/"RESULTS_KO.md"),
    code_source_checksums={str(p.relative_to(ROOT)):sha(p) for p in source},
    files_sha256=file_hashes,
    excluded_claims=["new Maxwell root sensitivity","global TF geometry minimum",
      "full continuous crust native reconstruction","causal APR above source high-density limit",
      "observationally true neutron star EOS"],
    next_action="new original-source phase-separated Maxwell root search with multiple interpolation families")
(out/"DELIVERY_MANIFEST.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
review_doc=["# Native dense EOS 독립 검토와 게시 범위","",
    "**판정:** ACCEPT_SCOPED_NUMERICAL_RESULT (전역 EOS의 물리적 승인은 아님).","",
    f"- 10개 고정 native 소스, 25,612개 직접 생성 상태, 5개 TF droplet 모델의 새 독립 정지해 30개.",
    f"- 독립 출처/물리/홀드아웃 확인 {review['check_count']}개, 미해결 차단 finding {len(review['failures'])}개.",
    "- RMF·core에서 native 소스 CompOSE 저장 전후 동일성과 물리잔차 확인. APR 6→8→7 개별 phase와 BBB 입력 상한 유지.",
    f"- 전체 회귀검사 1,501 passed, 9 skipped, 4 subtests passed, 실행 종료코드 {full['actual_exit_code']}. 시험 전후 코드 해시 동일.",
    f"- 원 EOS 25,612행의 상대 Euler identity 최대잔차 {native['max_relative_euler']:.3g}.",
    "- 도표는 nonlinear 8개 방법만 포함하며 linear 보간법은 전부 제외. 별도의 원시 로그에는 음속 미정의 diagnostic을 보존.",
    "- 실패한 TF 24개 신규 밀도점은 수용된 profile로 승격하지 않았으며 추가 geometry를 검증하지 않았다.",
    "- 이전 단계의 sharp Maxwell 교차근은 고정했으며 새로운 interpolation family로 교차근을 다시 찾지 않았다.",
    "- 동일 physical model을 정밀 계산한 결과에 대한 interpolation holdout error이지 실제 핵물질 EOS의 진리값이나 모델 간 물리적 우월성 판정이 아니다.",
    "",
    "## 다음 실행 최소 작업","",
    "같은 source의 두 원 phase table에서 각 family로 독립적으로 Maxwell P₋(μ)=P₊(μ) 교차근을 재탐색하고, pressure preference, n jump, latent heat와 지원 구간을 source identity별로 비교한다.",""]
(DOC/"INDEPENDENT_REVIEW_KO.md").write_text("\n".join(review_doc))
# Include the finalized review, excluding the manifest itself (no self-referential hash).
file_hashes={}
for folder in (DATA,DOC):
    for file in sorted(folder.rglob("*")):
        if file.is_file() and file != out/"DELIVERY_MANIFEST.json":
            file_hashes[str(file.relative_to(ROOT))]=sha(file)
receipt["files_sha256"]=file_hashes
(out/"DELIVERY_MANIFEST.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
print(json.dumps({k:receipt[k] for k in ("scope","native_state_count","spatial_TF_native_holdout_count",
    "independent_audit_checks","full_suite_actual_exit","plot_count")},indent=2))
print("HASHED_RESULT_FILES",len(file_hashes),"CODE_FILES",len(source))
