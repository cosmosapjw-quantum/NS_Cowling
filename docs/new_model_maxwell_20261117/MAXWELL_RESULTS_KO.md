# 신규 모델 sharp Maxwell 결과 — 2026-11-17

## 결과와 근거 상태

요청한 120개 모델·변형 조합을 모두 실행했다. 4,197개 connected-component pair를 native Maxwell kernel에 실제 전달했고, 9개 표가 kernel, 독립 interface admission, CompOSE 실제 재읽기, 재직렬화된 native consistency 검사를 통과했다. 기존 195개 조합은 재실행하지 않았다. 새 native core가 포함된 성공은 모두 PK1_CANONICAL이며, 나머지 세 물리 모델과 DDME-X의 별도 진단 변형에서는 이 유한 입력 범위에서 성공하지 않았다.

- **implementation-verified:** 새 parser의 Q1–Q7 및 Nadd5/7 보존, geometry/phase 검증, 120개 origin 계약, geometry/invalid 구간 분리: 5개 targeted test PASS.
- **numerically checked:** 9개 finite-domain 표, 독립 polynomial envelope 검사 9/9 PASS, 신규 native core 성공 7개의 같은 component coarse/fine 비교 7/7 유지.
- **unresolved / 제한:** crust의 미세물리 ground state, 기존 SPG의 원 논문 surface calibration, 제외된 source 구간, 전역 EOS 및 별 관측량 검증. 이 보고서는 강한 no-go나 whole-star EOS 성립을 주장하지 않는다.

## 고정한 입력과 origin 계약

기존 native crust 13개는 SPG_M1–M5 및 GPPVA_DD2, DDME2, TW, FSU2, FSU2R, FSU2H, NL3wrL55, TM1e다. 이들 및 GM1은 기존 prepared_v2 manifest의 원본·internal SHA-256을 대조한 후 재사용했고, 현재 native executable로 partner-independent interval 진단을 다시 수행했다. 기존 source의 variant·calibration·support 해석은 prepared_v2/MANIFEST.json을 그대로 상속한다. 특히 SPG native는 public_composition_inferred_20261116 candidate를 포함하므로 원 논문 surface 계수의 완전한 재현과 구분해야 한다.

신규 native 입력은 NL3_PROVIDER_PRECISION, PK1_CANONICAL, TM1_CANONICAL, DDMEX_PROVIDER_ROUNDED_N0, DDMEX_ORIGINAL_TABLE_RHO0이다. 물리 모델은 NL3, PK1, TM1, DDME-X 네 개다. DDME-X 두 입력은 n0 문헌 identity를 구분하는 두 변형이며 별도 물리 모델 성공으로 중복 집계하지 않는다. provider forward comparison에서는 ORIGINAL_TABLE_RHO0가 raw source에 더 잘 맞았다. 모든 입력은 생성자가 고정한 벡터로 평가했고 matching을 위한 에너지 이동·파라미터 조절은 없다.

공개 crust는 XMLSLZ_NL3, XMLSLZ_PK1, XMLSLZ_TM1, XMLSLZ_DDME-X, XMLSLZ_PKDD의 nonzero geometrical phase 전체를 출발점으로 한다. 공개 core는 앞의 네 표의 terminal homogeneous Iphase=0 구간이다. 기존 phase7 selector를 적용하지 않았다. Nadd=7(2025) 및 TM1 Nadd=5는 exact auxiliary token과 원문 thermo line을 SOURCE_AUXILIARY.json에 보존했고, native transport는 명시적 Nadd=0 thermo-only projection이다. Native projection의 Q1–Q7 전 열이 raw parsed double과 bitwise equal임을 실제 확인했다. Q3와 Euler chemical potential 차이는 수정하지 않았다.

이 새 캠페인은 public crust와 native core의 reference-assisted hybrid 30개를 현재 계약으로 명시했다. historical maxwell_three_family_campaign.py의 같은-family 제한은 변경하지 않았다. 출처와 role이 고정 목록 밖이면 새 policy도 거절한다.

| 범주 | 요청 셀 | 성공 셀/표 | 실제 kernel 호출 |
|---|---:|---:|---:|
| 신규 public crust × 신규 public core | 20 | 0 / 0 | 3696 |
| 신규 public crust × 기존 public GM1 | 5 | 1 / 1 | 28 |
| 기존 native crust × 신규 native core | 65 | 6 / 6 | 305 |
| 명시적 reference-assisted hybrid | 30 | 2 / 2 | 168 |
| 합계 | 120 | 9 / 9 | 4197 |

## 검증된 전이 수치

단위는 μ: MeV, P와 Δε: MeV fm⁻³, n: fm⁻³다. native kernel은 μ_E=(ε+P)/n의 P(μ_E) cubic Hermite 표현과 dP/dμ_E=n을 사용한다. 아래 표시 자리수는 readability용이며 원래 17-digit JSON/CSV를 보존했다.

| crust | core와 origin | μt | Pt | nlow | nhigh | Δε |
|---|---|---:|---:|---:|---:|---:|
| SPG_M1 | PK1_CANONICAL (native) | 943.292698146 | 0.01744421157 | 0.0112985466 | 0.01963254668 | 7.861401422 |
| SPG_M2 | PK1_CANONICAL (native) | 943.123169258 | 0.01430739828 | 0.009104210537 | 0.01733874598 | 7.766181166 |
| SPG_M3 | PK1_CANONICAL (native) | 943.348062872 | 0.01855074039 | 0.01183363669 | 0.02033609885 | 8.020781206 |
| SPG_M4 | PK1_CANONICAL (native) | 943.314268907 | 0.0178706915 | 0.01117581495 | 0.01990929715 | 8.238418377 |
| SPG_M5 | PK1_CANONICAL (native) | 943.356569520 | 0.01872418426 | 0.01200832548 | 0.02044224522 | 7.956193598 |
| GPPVA_DDME2 | PK1_CANONICAL (native) | 943.120232669 | 0.01425654248 | 0.009540221759 | 0.01729719612 | 7.315759465 |
| XMLSLZ_TM1 | GM_GM1 (public) | 945.433226661 | 0.1409405271 | 0.04496277449 | 0.0542918258 | 8.819995084 |
| XMLSLZ_TM1 | GM_GM1 (native) | 945.255896835 | 0.1330458839 | 0.04407346409 | 0.05354422061 | 8.952288446 |
| XMLSLZ_DDME-X | PK1_CANONICAL (native) | 943.959505651 | 0.03308501435 | 0.01762632517 | 0.02687036277 | 8.725997169 |

모든 성공은 Pt의 동일성, μt의 동일성, nhigh>nlow, Δε=μtΔn을 actual serialized rows에서 독립 재검산했다. 서로 다른 원본 mass header는 물리 ε·μ를 보존하는 Q 재인코딩만 거쳤다. native_report.physical_energy_shift_mev는 0이다. composition abundance는 matching 출력에서 이용 불가능하며 phase-only eos.compo로 명시된다; 원래 composition은 immutable source에 남는다.

## Envelope 범위와 숨은 재교차

기존 native find_coexistence는 두 component의 전체 common μ 구간을 union knot로 나누고, 각 cubic pressure difference의 quadratic 도함수 임계점을 사용해 모든 monotone subinterval의 root를 분리한다. 성공 조건은 isolated root 하나와 양의 density jump다.

추가로 SciPy CubicHermiteSpline의 독립 coefficient 표현으로 모든 union interval의 cubic polynomial root를 계산했다. ΔP=Pcore−Pcrust의 공통 구간 왼쪽 끝은 음수, 오른쪽 끝은 양수이고, 전이에서 dΔP/dμ>0이며 real root는 하나다. hidden three-root synthetic control도 통과했다. 9개 전이의 independent root와 native root 차이 최댓값은 3.21734e-11 MeV다.

이 envelope 주장은 **같은 두 connected component의 전체 공통 μ 지원 범위**에 한정된다. 다른 geometry, 불안정·비인과·보간 shape 실패 구간 또는 다른 물리 모델의 경쟁까지 포함하는 전역 ground-state 판정은 아니다. 경계 또는 invalid interval을 가로질러 보간하지 않았고, partner를 보고 domain을 잘라 root를 만들지 않았다.

## 실제 수렴 비교

성공한 새 native core는 PK1 하나이므로 그 core를 포함하는 7개 성공 component pair 모두를 독립 생성된 coarse 385-row 표로 다시 계산했다. fine 표는 769 rows, 두 표의 전체 core domain은 n=0.015–1 fm⁻³이다. 같은 phase와 동일 component endpoints를 확인한 뒤 같은 crust component를 그대로 사용했다. 두 해상도 모두 전이를 유지했다.

| 비교량 | 7개 전이 중 최대 상대 차이 |
|---|---:|
| mu_transition | 5.771305e-12 |
| pressure_transition | 3.472985e-09 |
| n_low | 2.393783e-09 |
| n_high | 2.12687e-07 |
| latent_energy_density | 4.484829e-07 |

이는 해당 homogeneous core discretization 변화에 대한 실제 transition sensitivity다. 두 해상도만으로 asymptotic convergence order 또는 엄밀한 error bound를 주장하지 않는다. 기존 native crust의 mesh 수렴을 새로 측정한 것도 아니다. Public source에는 원문 정밀도와 mesh만 있으므로 grid refinement를 만들어 주장하지 않았고, hybrid 비교에서도 같은 raw public crust를 고정했다.

## 보존한 실패와 공개 표의 제한

9개 성공 외 4,188개 component pair는 다음과 같이 거절됐다. 원래 native category와 메시지는 candidate RESULT.json에 보존했다.

| 진단된 거절 원인 | component pair 수 | 의미 |
|---|---:|---|
| NO_COEXISTENCE_CROSSING | 328 | 공통 구간에서 pressure equality 없음 |
| NONPOSITIVE_DENSITY_JUMP | 2 | 원한 forward crust→denser core 방향이 아님 |
| NO_COMMON_MU_SUPPORT | 3852 | 두 고정 component의 μ 지원 범위가 겹치지 않음 |
| INSUFFICIENT_RETAINED_SOURCE_KNOTS | 3 | 전이 한쪽에 보존할 source knot가 두 개 미만 |
| OUTPUT_THERMODYNAMIC_C2_REJECTION | 1 | 실제 output strict C2 consistency rejection |
| MULTIPLE_COEXISTENCE_CROSSINGS | 2 | 단일 isolated sharp transition 계약 밖 |

native taxonomy가 UNKNOWN으로 반환한 네 건도 모두 진단했다. 세 건은 GPPVA_NL3wrL55×PK1, GPPVA_TM1e×PK1, public NL3×native GM1의 retained knot 부족이다. 마지막 한 건은 public DDME-X×native GM1이며, 원본 crust component를 읽기 전용으로 재검사했을 때 C2 상대 잔차 4.4023902×10⁻⁶가 strict 2×10⁻⁶ gate를 초과했다. 이 component의 max|μQ3−μE|는 0.0041373185 MeV다. input phase Euler gate 통과가 최종 consistency 통과를 보장하지 않는 사례로 남겼다. Q3 또는 물리 에너지를 바꾸거나 tolerance를 늘리지 않았다.

raw public core는 높은 밀도까지 제공되지만 선택한 thermodynamic Hermite interpolant가 원문 정밀도에서 shape·causality 검사를 자주 통과하지 못한다. 아래는 모든 maximal valid ≥3-row component를 보존한 실제 수다. 이는 이 표·표현의 진단이며 모델 전체가 본질적으로 비인과적이라는 결론이 아니다.

| 공개 source | crust components | homogeneous core components |
|---|---:|---:|
| XMLSLZ_NL3 | 9 | 44 |
| XMLSLZ_PK1 | 6 | 41 |
| XMLSLZ_TM1 | 2 | 6 |
| XMLSLZ_DDME-X | 7 | 41 |
| XMLSLZ_PKDD | 4 | 범위 밖 |

20개 public crust×new public core 셀에서 수행한 3,696회 kernel 호출 중 3,685회는 NO_COMMON_MU_SUPPORT, 11회는 NO_COEXISTENCE_CROSSING이었다. 이 유한 source support 밖으로 extrapolate하지 않았다. causal completion과 CSS bridge는 이 sharp Maxwell 캠페인에 사용하지 않았다.

## 산출물과 재현

- data/new_model_maxwell_20261117/prepared/MANIFEST.json: 34개 source/해상도/role record, source hashes, phase diagnostics, frozen components.
- data/new_model_maxwell_20261117/campaign/CAMPAIGN.json: 120개 cell과 4,197개 candidate의 실제 결과·오류.
- data/new_model_maxwell_20261117/campaign/SUCCESS_INDEX.json: 9개 compose 디렉터리의 exact path와 SHA-256.
- data/new_model_maxwell_20261117/campaign/SUCCESS_TRANSITIONS.csv 및 MATRIX.csv: flat numerical summaries.
- data/new_model_maxwell_20261117/matching_validation/ENVELOPE_AUDIT.json: 9개 독립 polynomial envelope 검사.
- data/new_model_maxwell_20261117/matching_validation/convergence/CONVERGENCE.json: 7개 같은 component fine/coarse 전이.
- data/new_model_maxwell_20261117/matching_validation/REJECTION_DIAGNOSIS.json: UNKNOWN 원문 보존 및 네 실제 원인 분류.
- 각 성공 candidate의 compose/, expected_internal/, readback_internal/, CONSISTENCY.json, RESULT.json: 출력과 독립 검증 연결.

실행 command는 아래와 같다. prepare는 새 디렉터리에서 한 번 실행하고 execute는 완료 cell checkpoint를 재사용한다. 기존 산출물을 덮어쓰는 rerun은 별도 output version을 먼저 선택해야 한다.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export COMPOSE_EOS_FORTRAN_ROOT=/home/cosmosapjw/ns_cowling_wt/new_model_maxwell_20261117
export COMPOSE_EOS_FORTRAN_EXE=/home/cosmosapjw/ns_cowling_wt/stellar_consumers_20261117/build/compose_eos_fortran/compose_eos_fortran
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python -m unittest discover -s python_frontend/tests -p test_new_model_matching.py -v
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/run_new_model_maxwell_20261117.py prepare
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/run_new_model_maxwell_20261117.py execute
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/run_new_model_maxwell_20261117.py convergence
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python data/new_model_maxwell_20261117/matching_validation/audit_success_envelopes.py
```

Native executable SHA-256: 3893fbe061c2564ed4138211f7ea379a973851eb6e785fad22edcec0875152b1.
Campaign script SHA-256 at production: 3e1589704a6d756b0e6cfa572ab8b626bc8015c79889c387edb523a577b3e291.
Base git HEAD: 00bb5e1adcdbb8dfa0a0b2d5d7dd6f9037a4ddc3.
Production wall time: 35.534 s; OMP/BLAS threads=1. Production, preparation, convergence 및 envelope exit code는 모두 0이다.
이 subtask에서는 git commit/push를 수행하지 않았다. Parent가 통합 검토·후속 선택을 소유한다.

## 통합 검증과 수치 표기 수정

Parent 통합 suite는 reporting-only 수정 전 1,415 passed, 9 skipped, 4 subtests passed, actual exit 0으로 완료됐고 실행 전후 408개 source hash가 동일했다. 이후 fresh review의 비차단 finding에 따라 near-zero Maxwell identity residual을 상대오차로 나누던 집계만 수정했다. 잔차는 fine/coarse 원값과 절대 차이로 보존하며 물리 전이 수치나 admission은 바뀌지 않았다.

수정한 함수는 convergence()뿐이다. prepare, prepare_record, execute, run_candidate, make_pairs의 AST identity를 확인했고, 수정한 실제 aggregation AST를 기존 7개 수치에 실행하여 저장된 재집계와 정확히 일치함을 확인했다. primary campaign 당시 runner의 byte snapshot은 matching_validation/run_new_model_maxwell_20261117.production.py.txt에 남아 CAMPAIGN.json의 원래 script SHA를 계속 만족한다. 현재 runner SHA와 수치 재집계 전후 SHA는 CONVERGENCE_REPORTING_REPAIR.json에 별도로 기록했다. 이 수정으로 production이나 coarse Maxwell 계산을 재실행하지 않았다.
