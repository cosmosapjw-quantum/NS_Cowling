"""Generate compact machine summaries and the Korean finite-domain results note."""
from pathlib import Path
from collections import Counter
import csv
import json
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from maxwell_three_family_campaign import load,save,sha,rel
BASE=ROOT/"data/new_model_maxwell_20261117"
campaign=load(BASE/"campaign/CAMPAIGN.json")
manifest=load(BASE/"prepared/MANIFEST.json")
conv=load(BASE/"matching_validation/convergence/CONVERGENCE.json")
envelope=load(BASE/"matching_validation/ENVELOPE_AUDIT.json")
reject=load(BASE/"matching_validation/REJECTION_DIAGNOSIS.json")
records={r["id"]:r for r in manifest["source_records"]}
success=[]
for cell in campaign["cells"]:
    for i,r in enumerate(cell["candidate_results"]):
        if r["status"]!="GENERATED_VERIFIED_FINITE_DOMAIN":continue
        success.append({"category":cell["category"],"crust":cell["crust"],"core":cell["core"],
            "crust_source_id":cell["crust_source_id"],"core_source_id":cell["core_source_id"],
            "candidate_index":i,"crust_component_id":r["crust_component"]["id"],
            "core_component_id":r["core_component"]["id"],**r["transition"],
            "compose_path":r["output_path"],"compose_sha256":r["output_sha256"],
            "density_support_fm3":r["density_fm3"]})
save(BASE/"campaign/SUCCESS_INDEX.json",{"status":"9_FINITE_DOMAIN_OUTPUTS_VERIFIED",
    "campaign_sha256":sha(BASE/"campaign/CAMPAIGN.json"),"successful_tables":success})
fields=["category","crust","core","candidate_index","crust_component_id","core_component_id",
        "mu_transition","pressure_transition","n_low","n_high","epsilon_low","epsilon_high",
        "density_jump","latent_energy_density","maxwell_identity_residual","compose_path"]
with (BASE/"campaign/SUCCESS_TRANSITIONS.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader();w.writerows(success)
with (BASE/"campaign/MATRIX.csv").open("w",newline="") as f:
    w=csv.writer(f);w.writerow(["category","crust","core","status","generated_tables","candidate_attempts","actual_kernel_invocations"])
    for c in campaign["cells"]:w.writerow([c["category"],c["crust"],c["core"],c["status"],c["generated_count"],len(c["candidate_results"]),c["kernel_invocations"]])
names={"old_native_crust__new_native_core":"기존 native crust × 신규 native core",
       "new_public_crust__old_public_GM1":"신규 public crust × 기존 public GM1",
       "new_public_crust__new_public_core":"신규 public crust × 신규 public core",
       "reference_assisted_hybrid":"명시적 reference-assisted hybrid"}
lines=["# 신규 모델 sharp Maxwell 결과 — 2026-11-17","",
"## 결과와 근거 상태","",
"요청한 120개 모델·변형 조합을 모두 실행했다. 4,197개 connected-component pair를 native Maxwell kernel에 실제 전달했고, 9개 표가 kernel, 독립 interface admission, CompOSE 실제 재읽기, 재직렬화된 native consistency 검사를 통과했다. 기존 195개 조합은 재실행하지 않았다. 새 native core가 포함된 성공은 모두 PK1_CANONICAL이며, 나머지 세 물리 모델과 DDME-X의 별도 진단 변형에서는 이 유한 입력 범위에서 성공하지 않았다.","",
"- **implementation-verified:** 새 parser의 Q1–Q7 및 Nadd5/7 보존, geometry/phase 검증, 120개 origin 계약, geometry/invalid 구간 분리: 5개 targeted test PASS.",
"- **numerically checked:** 9개 finite-domain 표, 독립 polynomial envelope 검사 9/9 PASS, 신규 native core 성공 7개의 같은 component coarse/fine 비교 7/7 유지.",
"- **unresolved / 제한:** crust의 미세물리 ground state, 기존 SPG의 원 논문 surface calibration, 제외된 source 구간, 전역 EOS 및 별 관측량 검증. 이 보고서는 강한 no-go나 whole-star EOS 성립을 주장하지 않는다.","",
"## 고정한 입력과 origin 계약","",
"기존 native crust 13개는 SPG_M1–M5 및 GPPVA_DD2, DDME2, TW, FSU2, FSU2R, FSU2H, NL3wrL55, TM1e다. 이들 및 GM1은 기존 prepared_v2 manifest의 원본·internal SHA-256을 대조한 후 재사용했고, 현재 native executable로 partner-independent interval 진단을 다시 수행했다. 기존 source의 variant·calibration·support 해석은 prepared_v2/MANIFEST.json을 그대로 상속한다. 특히 SPG native는 public_composition_inferred_20261116 candidate를 포함하므로 원 논문 surface 계수의 완전한 재현과 구분해야 한다.","",
"신규 native 입력은 NL3_PROVIDER_PRECISION, PK1_CANONICAL, TM1_CANONICAL, DDMEX_PROVIDER_ROUNDED_N0, DDMEX_ORIGINAL_TABLE_RHO0이다. 물리 모델은 NL3, PK1, TM1, DDME-X 네 개다. DDME-X 두 입력은 n0 문헌 identity를 구분하는 두 변형이며 별도 물리 모델 성공으로 중복 집계하지 않는다. provider forward comparison에서는 ORIGINAL_TABLE_RHO0가 raw source에 더 잘 맞았다. 모든 입력은 생성자가 고정한 벡터로 평가했고 matching을 위한 에너지 이동·파라미터 조절은 없다.","",
"공개 crust는 XMLSLZ_NL3, XMLSLZ_PK1, XMLSLZ_TM1, XMLSLZ_DDME-X, XMLSLZ_PKDD의 nonzero geometrical phase 전체를 출발점으로 한다. 공개 core는 앞의 네 표의 terminal homogeneous Iphase=0 구간이다. 기존 phase7 selector를 적용하지 않았다. Nadd=7(2025) 및 TM1 Nadd=5는 exact auxiliary token과 원문 thermo line을 SOURCE_AUXILIARY.json에 보존했고, native transport는 명시적 Nadd=0 thermo-only projection이다. Native projection의 Q1–Q7 전 열이 raw parsed double과 bitwise equal임을 실제 확인했다. Q3와 Euler chemical potential 차이는 수정하지 않았다.","",
"이 새 캠페인은 public crust와 native core의 reference-assisted hybrid 30개를 현재 계약으로 명시했다. historical maxwell_three_family_campaign.py의 같은-family 제한은 변경하지 않았다. 출처와 role이 고정 목록 밖이면 새 policy도 거절한다.","",
"| 범주 | 요청 셀 | 성공 셀/표 | 실제 kernel 호출 |","|---|---:|---:|---:|"]
for cat,s in campaign["summary"].items():lines.append(f'| {names[cat]} | {s["model_variant_pairs"]} | {s["pairs_with_output"]} / {s["generated_tables"]} | {s["actual_kernel_invocations"]} |')
lines += ["| 합계 | 120 | 9 / 9 | 4197 |","",
"## 검증된 전이 수치","",
"단위는 μ: MeV, P와 Δε: MeV fm⁻³, n: fm⁻³다. native kernel은 μ_E=(ε+P)/n의 P(μ_E) cubic Hermite 표현과 dP/dμ_E=n을 사용한다. 아래 표시 자리수는 readability용이며 원래 17-digit JSON/CSV를 보존했다.","",
"| crust | core와 origin | μt | Pt | nlow | nhigh | Δε |","|---|---|---:|---:|---:|---:|---:|"]
for r in success:
    origin="public" if r["core_source_id"].startswith("public_original") else "native"
    lines.append(f'| {r["crust"]} | {r["core"]} ({origin}) | {r["mu_transition"]:.9f} | {r["pressure_transition"]:.10g} | {r["n_low"]:.10g} | {r["n_high"]:.10g} | {r["latent_energy_density"]:.10g} |')
lines += ["",
"모든 성공은 Pt의 동일성, μt의 동일성, nhigh>nlow, Δε=μtΔn을 actual serialized rows에서 독립 재검산했다. 서로 다른 원본 mass header는 물리 ε·μ를 보존하는 Q 재인코딩만 거쳤다. native_report.physical_energy_shift_mev는 0이다. composition abundance는 matching 출력에서 이용 불가능하며 phase-only eos.compo로 명시된다; 원래 composition은 immutable source에 남는다.","",
"## Envelope 범위와 숨은 재교차","",
"기존 native find_coexistence는 두 component의 전체 common μ 구간을 union knot로 나누고, 각 cubic pressure difference의 quadratic 도함수 임계점을 사용해 모든 monotone subinterval의 root를 분리한다. 성공 조건은 isolated root 하나와 양의 density jump다.","",
"추가로 SciPy CubicHermiteSpline의 독립 coefficient 표현으로 모든 union interval의 cubic polynomial root를 계산했다. ΔP=Pcore−Pcrust의 공통 구간 왼쪽 끝은 음수, 오른쪽 끝은 양수이고, 전이에서 dΔP/dμ>0이며 real root는 하나다. hidden three-root synthetic control도 통과했다. 9개 전이의 independent root와 native root 차이 최댓값은 "+f'{max(r["root_difference_from_native_MeV"] for r in envelope["results"]):.6g}'+" MeV다.","",
"이 envelope 주장은 **같은 두 connected component의 전체 공통 μ 지원 범위**에 한정된다. 다른 geometry, 불안정·비인과·보간 shape 실패 구간 또는 다른 물리 모델의 경쟁까지 포함하는 전역 ground-state 판정은 아니다. 경계 또는 invalid interval을 가로질러 보간하지 않았고, partner를 보고 domain을 잘라 root를 만들지 않았다.","",
"## 실제 수렴 비교","",
"성공한 새 native core는 PK1 하나이므로 그 core를 포함하는 7개 성공 component pair 모두를 독립 생성된 coarse 385-row 표로 다시 계산했다. fine 표는 769 rows, 두 표의 전체 core domain은 n=0.015–1 fm⁻³이다. 같은 phase와 동일 component endpoints를 확인한 뒤 같은 crust component를 그대로 사용했다. 두 해상도 모두 전이를 유지했다.","",
"| 비교량 | 7개 전이 중 최대 상대 차이 |","|---|---:|"]
for key in ("mu_transition","pressure_transition","n_low","n_high","latent_energy_density"):
    lines.append(f'| {key} | {max(r["relative_difference"][key] for r in conv["results"]):.7g} |')
lines += ["",
"이는 해당 homogeneous core discretization 변화에 대한 실제 transition sensitivity다. 두 해상도만으로 asymptotic convergence order 또는 엄밀한 error bound를 주장하지 않는다. 기존 native crust의 mesh 수렴을 새로 측정한 것도 아니다. Public source에는 원문 정밀도와 mesh만 있으므로 grid refinement를 만들어 주장하지 않았고, hybrid 비교에서도 같은 raw public crust를 고정했다.","",
"## 보존한 실패와 공개 표의 제한","",
"9개 성공 외 4,188개 component pair는 다음과 같이 거절됐다. 원래 native category와 메시지는 candidate RESULT.json에 보존했다.","",
"| 진단된 거절 원인 | component pair 수 | 의미 |","|---|---:|---|"]
meanings={"NO_COMMON_MU_SUPPORT":"두 고정 component의 μ 지원 범위가 겹치지 않음",
"NO_COEXISTENCE_CROSSING":"공통 구간에서 pressure equality 없음",
"NONPOSITIVE_DENSITY_JUMP":"원한 forward crust→denser core 방향이 아님",
"MULTIPLE_COEXISTENCE_CROSSINGS":"단일 isolated sharp transition 계약 밖",
"INSUFFICIENT_RETAINED_SOURCE_KNOTS":"전이 한쪽에 보존할 source knot가 두 개 미만",
"OUTPUT_THERMODYNAMIC_C2_REJECTION":"실제 output strict C2 consistency rejection"}
for key,count in reject["diagnosed_rejection_counts"].items():lines.append(f'| {key} | {count} | {meanings[key]} |')
lines += ["",
"native taxonomy가 UNKNOWN으로 반환한 네 건도 모두 진단했다. 세 건은 GPPVA_NL3wrL55×PK1, GPPVA_TM1e×PK1, public NL3×native GM1의 retained knot 부족이다. 마지막 한 건은 public DDME-X×native GM1이며, 원본 crust component를 읽기 전용으로 재검사했을 때 C2 상대 잔차 4.4023902×10⁻⁶가 strict 2×10⁻⁶ gate를 초과했다. 이 component의 max|μQ3−μE|는 0.0041373185 MeV다. input phase Euler gate 통과가 최종 consistency 통과를 보장하지 않는 사례로 남겼다. Q3 또는 물리 에너지를 바꾸거나 tolerance를 늘리지 않았다.","",
"raw public core는 높은 밀도까지 제공되지만 선택한 thermodynamic Hermite interpolant가 원문 정밀도에서 shape·causality 검사를 자주 통과하지 못한다. 아래는 모든 maximal valid ≥3-row component를 보존한 실제 수다. 이는 이 표·표현의 진단이며 모델 전체가 본질적으로 비인과적이라는 결론이 아니다.","",
"| 공개 source | crust components | homogeneous core components |","|---|---:|---:|"]
for model in ("XMLSLZ_NL3","XMLSLZ_PK1","XMLSLZ_TM1","XMLSLZ_DDME-X","XMLSLZ_PKDD"):
    cr=records[f"new_public:{model}:crust"]
    co=records.get(f"new_public:{model}:core")
    lines.append(f'| {model} | {len(cr["components"])} | {len(co["components"]) if co else "범위 밖"} |')
lines += ["",
"20개 public crust×new public core 셀에서 수행한 3,696회 kernel 호출 중 3,685회는 NO_COMMON_MU_SUPPORT, 11회는 NO_COEXISTENCE_CROSSING이었다. 이 유한 source support 밖으로 extrapolate하지 않았다. causal completion과 CSS bridge는 이 sharp Maxwell 캠페인에 사용하지 않았다.","",
"## 산출물과 재현","",
"- data/new_model_maxwell_20261117/prepared/MANIFEST.json: 34개 source/해상도/role record, source hashes, phase diagnostics, frozen components.",
"- data/new_model_maxwell_20261117/campaign/CAMPAIGN.json: 120개 cell과 4,197개 candidate의 실제 결과·오류.",
"- data/new_model_maxwell_20261117/campaign/SUCCESS_INDEX.json: 9개 compose 디렉터리의 exact path와 SHA-256.",
"- data/new_model_maxwell_20261117/campaign/SUCCESS_TRANSITIONS.csv 및 MATRIX.csv: flat numerical summaries.",
"- data/new_model_maxwell_20261117/matching_validation/ENVELOPE_AUDIT.json: 9개 독립 polynomial envelope 검사.",
"- data/new_model_maxwell_20261117/matching_validation/convergence/CONVERGENCE.json: 7개 같은 component fine/coarse 전이.",
"- data/new_model_maxwell_20261117/matching_validation/REJECTION_DIAGNOSIS.json: UNKNOWN 원문 보존 및 네 실제 원인 분류.",
"- 각 성공 candidate의 compose/, expected_internal/, readback_internal/, CONSISTENCY.json, RESULT.json: 출력과 독립 검증 연결.","",
"실행 command는 아래와 같다. prepare는 새 디렉터리에서 한 번 실행하고 execute는 완료 cell checkpoint를 재사용한다. 기존 산출물을 덮어쓰는 rerun은 별도 output version을 먼저 선택해야 한다.","",
"```bash",
"export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1",
"export COMPOSE_EOS_FORTRAN_ROOT=/home/cosmosapjw/ns_cowling_wt/new_model_maxwell_20261117",
"export COMPOSE_EOS_FORTRAN_EXE=/home/cosmosapjw/ns_cowling_wt/stellar_consumers_20261117/build/compose_eos_fortran/compose_eos_fortran",
"/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python -m unittest discover -s python_frontend/tests -p test_new_model_matching.py -v",
"/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/run_new_model_maxwell_20261117.py prepare",
"/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/run_new_model_maxwell_20261117.py execute",
"/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/run_new_model_maxwell_20261117.py convergence",
"/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python data/new_model_maxwell_20261117/matching_validation/audit_success_envelopes.py",
"```","",
"Native executable SHA-256: "+campaign["matching_executable_sha256"]+".",
"Campaign script SHA-256 at production: "+campaign["campaign_script_sha256"]+".",
"Base git HEAD: "+campaign["git_base_head"]+".",
f'Production wall time: {campaign["wall_seconds"]:.3f} s; OMP/BLAS threads=1. Production, preparation, convergence 및 envelope exit code는 모두 0이다.',
"이 subtask에서는 git commit/push를 수행하지 않았다. Parent가 통합 검토·후속 선택을 소유한다.",""]
lines += ["## 통합 검증과 수치 표기 수정", "",
"Parent 통합 suite는 reporting-only 수정 전 1,415 passed, 9 skipped, 4 subtests passed, actual exit 0으로 완료됐고 실행 전후 408개 source hash가 동일했다. 이후 fresh review의 비차단 finding에 따라 near-zero Maxwell identity residual을 상대오차로 나누던 집계만 수정했다. 잔차는 fine/coarse 원값과 절대 차이로 보존하며 물리 전이 수치나 admission은 바뀌지 않았다.", "",
"수정한 함수는 convergence()뿐이다. prepare, prepare_record, execute, run_candidate, make_pairs의 AST identity를 확인했고, 수정한 실제 aggregation AST를 기존 7개 수치에 실행하여 저장된 재집계와 정확히 일치함을 확인했다. primary campaign 당시 runner의 byte snapshot은 matching_validation/run_new_model_maxwell_20261117.production.py.txt에 남아 CAMPAIGN.json의 원래 script SHA를 계속 만족한다. 현재 runner SHA와 수치 재집계 전후 SHA는 CONVERGENCE_REPORTING_REPAIR.json에 별도로 기록했다. 이 수정으로 production이나 coarse Maxwell 계산을 재실행하지 않았다.", ""]
out=ROOT/"docs/new_model_maxwell_20261117/MAXWELL_RESULTS_KO.md"
out.parent.mkdir(parents=True,exist_ok=True);out.write_text("\n".join(lines))
print(json.dumps({"successful_tables":len(success),"report":rel(out),"report_sha256":sha(out)}))
