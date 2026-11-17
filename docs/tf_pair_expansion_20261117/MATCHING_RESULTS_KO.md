# C: 실제 지원되는 나머지 core/method matching

신규 80개 cell을 기존 native Fortran으로 한 번 실행했다. 80개 모두 kernel 호출과 serialized readback/native consistency 검증을 통과했다. 이 중 smooth 70개는 `CONSTRUCTED_VERIFIED`, CSS 10개는 독립 interface admission까지 `ADMITTED`다. Smooth 결과를 CSS/Maxwell의 독립 admission과 합치지 않는다.

## 선택한 물리 범위와 source identity

- SPG M1–M5의 `public_composition_inferred_20261116` sealed 생성 자료를 쓴다. 동일 자료의 기존 quintic_P/CSS GM1·SLy4 20개는 다시 생성하지 않았다.
- 신규 cubic_F/cubic_P/quintic_F는 GM1, SLy4, KDE0v1, BBB_BHF2 전부에 연결한다: 5×4×3=60개.
- 신규 quintic_P/CSS는 KDE0v1, BBB_BHF2에 연결한다: 5×2×2=20개.
- Core 입력은 `data/research_matching_20261117/integration/real_pairs_v2/core/`의 실제 129-row, 0.08–0.5 fm^-3 자료다. 모든 신규 출력의 상단은 실제 0.5 fm^-3다. B의 0.02–1.0 fm^-3 core 자료를 이 표에 붙였다고 주장하지 않는다.
- B의 domain audit는 GM1/SLy4/KDE0v1에 대한 기존 근거로 재사용한다. BBB는 별도의 sealed single-phase 자료와 이번 own-source native verifier PASS를 기록한다. BBB에 B domain certification을 전이하지 않았다.
- Outer/native receipt, request, 모든 sealed artifact와 SUPPORT, 단일 contiguous phase를 검사하고 매 cell 전에 다시 identity를 확인했다. 총 9개 source snapshot이다. `UNAPPROVED` declaration을 변경하지 않았다.
- Composition은 명시적으로 thermo-only로 버린다. n_B 단위는 fm^-3, P와 epsilon은 MeV fm^-3, mu는 MeV다. epsilon은 원래 rest energy를 포함한다. 물리적 epsilon→epsilon+Cn 이동은 없다.

## 실제 실행 결과

| method | 신규 kernel | constructed + serialized 검증 | 독립 method admission |
|---|---:|---:|---:|
| cubic_F | 20 | 20 | 0 |
| cubic_P | 20 | 20 | 0 |
| quintic_F | 20 | 20 | 0 |
| quintic_P | 10 | 10 | 0 |
| css | 10 | 10 | 10 |
| 합계 | 80 | 80 | 10 |

실행 exit code는 0, wall time은 8.93 s였다. cell별 결과는 `data/tf_pair_expansion_20261117/matching/campaign_v1/pairs/*/PAIR.json`, 전체 receipt는 `campaign_v1/RECEIPT.json`에 있다. `matching/CATALOG.csv`는 80개 표의 정확한 cell ID, 출력 경로, 행 수, 실제 밀도 범위를 담는다. 실패 cell이나 unchanged failure의 반복 실행은 없다. 실제 적용한 native executable은 B build의 `compose_eos_fortran`이며 SHA-256은 `9aec926552e5f1077b206bb0067e2120f2c991553a493e8ed831daedfabea461`이다. 실행 직전 C/B의 native source 26개에 BYTE_IDENTITY를 확인했다.

## 390개 requested cell의 중복 없는 합집합

`matching/campaign_v1/MATCHING_390_MAP.json`은 original 13 crust ×5 core ×6 method=390개 ID를 그대로 유지한다. 모든 source/domain variant를 `evidence_variants`에 남기고 한 ID는 한 번만 집계한다.

| 항목 | unique model/method cell |
|---|---:|
| requested universe | 390 |
| 해당 finite component/method의 실제 kernel 근거 | 165 |
| constructed + serialized readback 표가 있는 ID | 105 |
| 독립 CSS/Maxwell method admission | 25 |
| construction-only smooth | 80 |
| 이번 신규 실제 kernel | 80 |
| 기존 source-bound 근거를 재사용한 ID | 85 |
| APR retained Gibbs path가 필요한 deferred SPG cell | 25 |
| TF full-source smooth/CSS 지원이 미확정인 deferred cell | 200 |

`eligible=165`는 기록된 finite component/method kernel의 eligibility다. TF full-source production이나 whole-star EOS의 보증이 아니다. 나머지 225개를 영구적인 unsupported로 분류하지 않는다. 선택한 C cohort에서 supported-but-unattempted cell은 0개다. 25개 SPG→APR의 non-Maxwell cell은 retained 6→8→7 경계 처리를 D로 넘긴다. 200개 TF non-Maxwell cell은 8 TF model×5 core×5 method이며 finite TF component의 존재와 전체 source의 smooth/CSS 승인 여부를 구분한다. TF diagnostic/production에 관한 최종 판정은 별도 TF lane과 상위 C 보고서를 따른다.

기존 Maxwell `campaign_v2`의 `native_generated` 65개 pair에는 모두 실제 candidate와 kernel 근거가 있다. 총 **427 component kernel invocations**이고 unique matrix cell은 65개다. 60개는 `NO_MATCH_ON_ADMITTED_COMPONENTS`, 5개는 `GENERATED_VERIFIED_FINITE_DOMAIN`이며 5개 ID에서 실제 표 6개가 나왔다. 이 admission은 원 source의 일부 finite connected components에 한정된다. public_original/public_reconstructed 130개는 이번 native390 합집합에 더하지 않는다. 초기 30-cell campaign의 Maxwell10과 중복되는 ID는 2번 세지 않았다. 기존의 잘못된 pressure/energy-shift census는 재사용하지 않았다.

105개 output ID는 smooth80 + CSS20 + Maxwell5다. 독립 admission25는 CSS20 + Maxwell5다. raw output table 개수는 Maxwell의 다중 output 때문에 unique ID 개수와 다르다.

## 검증과 물리적 유효범위

신규 source snapshot9개와 출력 table80개의 CompOSE payload400파일을 SHA-256으로 확인했다. 재사용한 기존 20개 표와 Maxwell6개 표의 output seal도 확인했다. Maxwell의 native source/prepared payload/receipt seal, 숫자·metadata readback PASS, serialized native verification PASS, 독립 admission ADMITTED를 확인한 뒤에만 표를 합집합에 넣었다. 재검증 때문에 기존 matching physics를 다시 실행하지 않았다.

| smooth method | 최대 analytic first-law relative residual | 최대 bridge c_s²/c² |
|---|---:|---:|
| cubic_F | 3.8377242e-14 | 0.33371499 |
| cubic_P | 9.0409099e-13 | 0.32259057 |
| quintic_F | 1.7315796e-14 | 0.45572286 |
| quintic_P | 2.7705273e-13 | 0.19045897 |

Native hard gate는 C1–C6/C10/C11이고 C7–C9는 finite-difference diagnostic이다. 신규 cohort에서 C8 최대 diagnostic residual은 quintic_F의 0.03935953이다. 이를 analytic first-law 실패나 mesh-independent 정확도 인증과 혼동하지 않는다. 독립적인 smooth-method별 continuum admission 또는 mesh convergence는 이번 산출물에 포함하지 않았다. TOV/Cowling, composition-dependent buoyancy, whole-crust 연결, global minimum의 증거도 이 inner-crust+core 표만으로 주어지지 않는다.

Regression은 TDD RED4 fail→GREEN4 pass, reuse-seal 추가 RED2 fail→matching-focused GREEN15 pass다. 관련한 기존 generated campaign 회귀를 포함한다. 전체 frontend/pipeline suite는 상위 agent가 TF 작업 후 한 번 실행하도록 맡겼다. 해당 결과는 상위 validation record를 따른다.

- `matching/validation/EXECUTION.json`: 실제 명령·executable/source identity·exit code.
- `matching/validation/ARTIFACT_VERIFICATION.json`: 실제 artifact 확인·물리 diagnostic 집계.
- `matching/validation/BBB_SOURCE_NATIVE_VERIFY.json`: BBB own-source 검증.
- `matching/validation/matching_regression.log`: 15 passed.
- `matching/validation/EXECUTED_RUNNER_SNAPSHOT.py`: 실제 80-cell 실행 당시 runner의 byte-exact 보존본. 이 보존본의 SHA-256은 receipt의 실행 당시 runner hash와 같다.

실행 후 runner에는 historical output-seal 확인과 deferred category 표기만 보강했다. 실행 당시와 현재 runner identity를 모두 verification record에 보존했다. Native matching 방정식, 실행된 cell 경로, source table은 바뀌지 않았다. 게시·commit·push·PR은 수행하지 않았다.

## 최대 C8 diagnostic의 제한된 sampling 확인

최대 C8을 보인 정확한 cell `SPG_M2__GM_GM1__quintic_F` 한 개에 대해 같은 source·endpoint·analytic bridge를 유지하고 n_bridge=80,160을 별도 경로에서 추가 실행했다. 이는 이미 집계한 한 matrix cell의 sampling sensitivity이며 coverage나 method admission을 2개 늘리지 않는다. 결과는 `matching/sampling_refinement_v1/REFINEMENT.json`, 실제 재현 명령은 `matching/validation/REFINE_WORST_C8.py`다. 두 run 모두 exit0, serialized/native verification PASS다.

| n_bridge | C7 FD 최대 | C8 FD 최대 | C9 FD 최대 | analytic first-law residual |
|---|---:|---:|---:|---:|
| 40 | 0.0058714458 | 0.0393595273 | 0.0003169324 | 0 |
| 80 | 0.0058714458 | 0.0116674552 | 0.0000898888 | 0 |
| 160 | 0.0058714458 | 0.0033483503 | 0.0000245299 | 0 |

C8 최대 위치는 bridge 내부 n_B≈0.129 fm^-3에 남고, 세 점 비균일-grid 미분의 stencil 간격이 작아지면서 감소한다. 각 단계 감소율은 약3.37,3.48이다. C9도 감소한다. C7 최대는 동일한 SPG source의 첫 행 n_B=0.00024162193 fm^-3, phase1에 있고, unchanged source sampling 때문에 일정하다. Native C7/C9의 endpoint one-sided stencil까지 직접 복원해 이 위치와 residual을 확인했다.

n1,n2,P1,P2,epsilon1,epsilon2,mu1,mu2 및 spline mode는 정확히 같다. 220개 공통 source/endpoint node는 모든 serialized column이 exact equality이고 P·epsilon 차이는 0이다. Bridge의 서로 다른 grid 사이에는 공통 interior knot가 없으므로, knot를 가로지르지 않는 세 u-window에서 q=P/n의 cubic와 F=epsilon/n의 quartic를 각 sampling으로 독립 복원했다. 공통 27개 평가점에서 최대 |ΔP|=5.56e-16 MeV fm^-3, 최대 |Δepsilon|=2.85e-14 MeV fm^-3, polynomial fit 최대 잔차=2.12e-13이다. 이는 동일한 analytic bridge를 더 촘촘히 sampling했다는 근거다. 제한된 세 sampling으로 continuum convergence나 smooth 독립 method admission을 선언하지 않는다.
