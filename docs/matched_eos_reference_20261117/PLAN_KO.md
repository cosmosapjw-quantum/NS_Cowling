# Reference-assisted matched EOS implementation plan

**Goal:** 최신 생성 EOS로 실제 core–crust matching cohort를 만들고 공개 EOS와의 차이를 물리량·모델 가정·수치 오차로 설명한다.

**Architecture:** 생성·matching 물리는 기존 native Fortran을 사용한다. 새 campaign runner는 실제 생성 receipt와 지원구간을 확인하고 기존 API/pipeline을 호출한다. 공개 CompOSE/Zenodo 자료는 별도 reference intake로 보존하여 비교, seed와 분기 진단에 사용한다.

**Tech stack:** Fortran, Python, NumPy/SciPy, CompOSE, Git. 계산은 이 대화에서 제어하는 별도 worktree에서 실행한다.

**Spec:** 2026-11-17 사용자 후속 요청과 `docs/audit_eos_matching_20261117/HANDOFF_CODEX_KO.md`. 최신 시작점은 `32074d48793e265ab21f7b7db0553f6b3d1283dc` / PR16이다.

## 이번 실행의 결정

이전 local 작업의 단일 smoke allowance는 종료된 노드의 이력이다. 그 CLAIM과 산출물은 보존한다.
사용자의 새 후속 실행 요청에 따라 별도 campaign과 출력 경로를 사용한다.
새 실행은 물리 SUPPORT를 변경하거나 과거 결과를 새 실행으로 집계하지 않는다.

첫 cohort는 SPG M1–M5의 명시적 `public_composition_inferred_20261116` variant와 GM1·SLy4 core를 연결한다.
각 pair에 `quintic_P`, `css`, `maxwell`을 요청하므로 10개 model pair, 30개 method cell이다.
기존 sealed SPG_M1–GM1 CSS는 입력·방법·출력 identity를 확인해 재사용하고 별도로 집계한다.
Core 상단은 현행 source의 0.5 fm^-3다. 산출물은 inner-crust+core이며 whole-star EOS와 구분한다.
Maxwell 실패도 실제 결과로 보존하며 동일 입력의 무의미한 반복 실행은 하지 않는다.

비교는 우선 동일 homogeneous 범위의 GM1·SLy4, 그리고 SPG의 공개 crust trajectory를 대상으로 한다.
GPPVA DD2의 공개 crust/unified archive는 TF branch와 lower-segment 진단용으로 내려받고 checksum을 확인한다.
전체 TF 재최적화나 거대 SPG 저자 eos.txt 다운로드를 이 cohort의 선행조건으로 추가하지 않는다.

## 완료된 기반

- PR16의 source seal intake, 390 requested / 1 attempted / 1 admitted 기록을 확인했다.
- 해당 head의 Fortran frontend CI와 sharp-Maxwell CI 성공을 확인했다.
- 새 worktree에 현재 native executable을 빌드했다.

## 실행 순서와 산출물

1. **공개 reference intake**
   - GM1, cold SLy4, SPG M1 및 DD2 crust/unified의 실제 공개 endpoint와 기존 보존 raw 자료를 연결한다.
   - 원 URL/DOI, version, 다운로드 시각, SHA-256, 게시 checksum, 모델·입자·밀도 범위를 기록한다.
   - `data/matched_eos_reference_20261117/references/`에 원 파일과 provenance를 보존한다.
2. **명시적 campaign 실행기**
   - `scripts/generated_pair_campaign.py`와 해당 회귀를 추가한다.
   - 실제 outer/native receipt와 request/artifact/SUPPORT hash를 검증하고 source를 선택한다.
   - 명시적 모델·방법 목록, 새 출력 경로, n_bridge=40을 받아 한 cell씩 실행하고 결과를 즉시 보존한다.
   - request, eligibility, kernel attempt, constructed/readback, method admission, reuse를 분리한다.
3. **실제 matching cohort**
   - 30개 요청 cell을 실행 또는 identity-bound reuse로 평가한다.
   - 각 실제 출력의 CompOSE table, native report, readback, consistency, source identity를 보존한다.
   - CSS/Maxwell의 독립 interface admission과 smooth native construction 검증의 범위를 구분한다.
4. **공개본과의 정량 비교**
   - `scripts/compare_generated_reference.py`와 analytic control을 추가한다.
   - header로 복원한 물리 n, P, epsilon, mu 및 정의된 composition만 비교한다.
   - 같은 연속 phase와 공통 밀도 지원에서 비교하고 interpolation sensitivity를 수치로 남긴다.
   - difference ledger는 confirmed convention/model difference, controlled numerical explanation, unresolved residual을 구분한다.
   - 유의미한 core 차이가 있으면 한 가지 상수·표현 선택을 바꾸는 제한된 controlled calculation으로 확인한다.
5. **결과·검증·게시**
   - 관련 회귀, 실제 산출물 readback, source identity, manifest를 검증한다.
   - 결과 요약, 비교 CSV/그림, machine-readable campaign 상태와 재개 지점을 같은 branch에 게시한다.

## 후속 과학 DAG

| ID | 산출물 | 선행조건 | 완료 판정 |
|---|---|---|---|
| P01 | 현재 source/reference identity와 차이 ledger | 공개자료·기존 receipt | 설명 가능한 비교 조건, 불명 잔차의 크기·구간 명시 |
| P02 | 10 model pair / 30 method-cell 실제 cohort | P01 및 native support | 출력 또는 정확한 실패, kernel attempt와 reuse 분리 |
| P03 | 나머지 core/방법 및 TF의 지원 가능한 구간으로 확장 | P02, 각 source의 실제 지원 | 390 request coverage와 실제 계산·승인을 분리한 census |
| P04 | surface부터 연결되는 whole-crust+core 후보 | 실제 lower/inner 접합, P02/P03 | 경계의 first law·P/mu 법칙과 origin 보존 |
| P05 | APR의 연속 Gibbs 경계별 native/consumer 지원 | 원 construction·endpoint 법칙 | retained 6→8→7 및 piecewise interpolation을 소비자까지 검증 |
| P06 | 필요한 TF continuation·mesh/geometry 계산 | 공개 seed/branch diagnostic, 기존 유한 해 | 독립 residual·수렴·phase 선택 근거 |
| P07 | 물리 범위가 명시된 TOV/Cowling 관측량 | P04/P05 및 perturbation closure | background와 진동 가정별 결과, 안정성 범위 명시 |

## 검증과 차이 설명의 원칙

공개본과의 일치를 생성의 공통 조건으로 강제하지 않는다. 그러나 불일치는 설명 없이 통과시키지 않는다.
다른 calibration/입자 내용/질량 convention/근사라는 문헌상 가능성과, 그 영향이 계산으로 재현됐다는 주장을 구분한다.
남는 잔차는 `UNRESOLVED_SOURCE_DIFFERENCE`와 크기·밀도·영향을 기록한다.
이 상태와 native 열역학 일관성 및 matching 성립 여부는 별도 항목이다.

epsilon은 rest energy 포함, n은 baryon number density, c_s^2/c^2=dP/dε다.
header 재인코딩은 물리 epsilon을 보존하며, 교차점을 만들기 위한 임의 ε→ε+Cn 보정은 사용하지 않는다.
missing composition은 0으로 채우지 않는다. 공개 segment를 직접 사용하는 hybrid에는 행 origin을 표시한다.
sealed input 변경·삭제, 서로 다른 lepton/charge slice, 지원 밖 구간, 이미 존재하는 출력경로,
kernel 호출 이전의 실패를 실제 matching attempt로 잘못 세는 경우를 회귀로 검사한다.

## 2026-11-17 실행 후 상태와 다음 순서

이번 첫 cohort는 완료했다. 요청30 / eligible30 / 신규 kernel29 / 신규 EOS19 / 재사용 EOS1이다.
Smooth10, CSS10의 실제 표와 Maxwell10의 종료 기록은 `deliverable/`과 `campaign/`에 있다.
공개 archive9개를 확보했고 DD2 Zenodo 두 개는 게시 MD5가 일치했다.
최신 공개본 비교, SLy4 kinetic control, GM1의 두 제한된 control, DD2 seed 좌표 추출을 실행했다.
전체 frontend/pipeline은 1260 passed, 9 skipped다. 상세 근거는 `RESULTS_KO.md`에 있다.

P01은 비교·차이 ledger가 완성됐다는 범위에서 완료다. GM1 source 차이의 완전한 원인 귀속과
SPG 원저자 calibration은 미확정이다. P02는 30칸 모두 결과가 있다는 범위에서 완료이며
Maxwell10의 거절을 성공으로 바꾸지 않았다.

### 다음 실행 A: whole-crust+core의 명시적 baseline

먼저 SPG–SLy4 quintic_P를 smooth baseline 후보로 삼고 CSS는 construction sensitivity 후보로
유지한다. 이것은 관측 적합도나 별 안정성을 이미 판정했다는 뜻이 아니다.
공개 BSk22 provider lower segment와 현재 native SPG의 실제 중첩 구간에서 물리량을 복원한다.
공개 표를 붙이는 경로는 `reference-assisted hybrid`로 명명하고 각 row의 public/generated/bridge
origin과 hash를 기록한다. Native-generated whole-crust라고 부르지 않는다.

이번 DD2 archive에서 확인한 `2.69e-4 → 2.00e-3 fm^-3` 사이의 tabulated gap과 phase 경계를
일반 spline으로 감추지 않는다. 각 모델의 실제 overlap 또는 boundary law를 확인하여
연결 endpoint를 선택하고, first law·P/μ 조건·단조성·causality·readback을 검사한다.
표면까지는 명시적 저밀도 EOS 또는 surface cutoff/extension prescription을 추가하고 그 영향의
수렴을 확인한다. 기존 public threshold를 자동으로 생성모델의 전이 밀도로 고정하지 않는다.

**산출물:** surface/outer/inner/core origin이 표시된 EOS, boundary report, 실제 지원구간,
조건부 실패 목록. **종료:** 쓸 수 있는 최소 하나의 whole-crust+core 후보와 그 물리 범위를 확보.
일치하는 경계가 없으면 필요한 구간만 추가 생성하거나 별도 bridge variant를 정의한다.

### 다음 실행 B: 필요한 core 지원 확장과 source-variant 정책

현재 upper bound 0.5 fm^-3는 matching smoke 범위다. 목표 별 모형의 central density를
실제로 수용하도록 각 native core를 확장하되 기계적 불안정·causality 위반·phase 경계를
검사하고 실제 endpoint를 기록한다. 같은 입력의 무차별 390회 재실행을 선행조건으로 두지 않는다.
직접 Maxwell의 두 μ-support 불일치는 유한 core lower bound의 영향을 분리하는 제한된
추가 homogeneous 계산 대상으로 삼을 수 있다. 기존8개 무교차 결과와 섞지 않는다.

SLy4는 현행 fitted common kinetic coefficient를 기본 정의로 유지한다. Public trajectory와의
작은 차이를 kinetic prescription으로 대부분 설명했으므로, 필요한 경우에만 입자별 kinetic
계수를 **명시적 native variant**로 구현하고 request/receipt에 기록한다. 별도 β평형을 풀고,
공통 K 기본값 및 analytic energy/chemical-potential/pressure relation을 회귀로 유지한다.
공개본과 더 잘 맞는다는 이유만으로 기존 fitted functional을 조용히 바꾸지 않는다.

GM1은 정확한 provider input/implementation이 공개된 범위에서 확보되는지를 한정해서 확인한다.
그 정보가 없으면 independent rounded-parameter reconstruction이라는 모델 정의와 현재
최대0.4864% pressure residual을 함께 보존한다. 완전한 source fidelity를 요구하는 분석은
그 attribution이 미확정이라는 상태를 사용하고, 반복 fitting으로 원저자 입력을 가장하지 않는다.

**산출물:** model/request별 실제 density domain과 difference ledger, 필요한 선택적 variant.
**종료:** 해당 별 모형과 matching에 필요한 지원이 확보되거나 정확한 domain failure를 기록.

### 다음 실행 C: TF와 나머지 core/method의 필요한 확대

DD2의 `reference_support/SEED_CATALOG.json`에 있는 density/Ye/Ymu 좌표를 TF ansatz 생성에
사용한다. 없는 meson fields, cell radius, cluster/group 의미를 임의로 보충하지 않는다.
자체 functional로 density continuation을 수행하고 허용 geometry별 에너지·잔차·mesh 및
cell-size refinement를 비교한다. Public pasta/endpoint는 branch 탐색의 참고이며 강제 정답이 아니다.
전체 TF 최적화를 다시 시작하기보다 기존 유한 해에서 필요한 빈 구간을 확장한다.

나머지 native core와 smooth methods는 source seal과 phase 지원이 실제로 제공되는 조합부터
추가한다. 전체 requested390, eligible, attempted, admitted, diagnostic, unsupported를 계속
분리한다. 실패를 지우거나 SUPPORT의 UNAPPROVED 문자열만 바꿔서 coverage를 늘리지 않는다.

**산출물:** TF continuation 증거, 가능한 추가 pair, 390 request coverage의 실제 census.
**종료:** 선택한 물리 범위의 source와 pair가 완결되거나 필요한 입력/계산의 정확한 실패를 보존.

### 다음 실행 D: APR boundary와 TOV/Cowling 소비자 연결

APR의 retained `6→8→7` 경계를 해당 Gibbs construction에 맞춰 표현하고 소비자 interpolation도
같은 경계를 보존하도록 연결한다. 연속 mixed phase에 sharp-Maxwell 조건을 일괄 적용하지 않는다.
Whole-crust와 충분한 core 지원을 갖춘 EOS에서 TOV를 계산한 후 Cowling perturbation의
equilibrium/frozen composition closure 및 phase interface 조건을 명시한다.
Thermo-only barotrope만으로 composition-dependent buoyancy가 정해졌다고 주장하지 않는다.
CSS의 Seidov flag와 실제 stellar sequence/mode 안정성은 각각 보고한다.

**산출물:** EOS identity에 연결된 M–R/central-density/domain 결과와 명시적 closure의 mode 결과.
**종료:** 선택한 관측량의 수치 수렴과 EOS/경계/closure 민감도를 함께 기록.

## 다음 세션에서 재사용할 상태

진행은 이 대화에서 repository와 연결된 계산 환경을 제어하는 방식으로 이어간다.
별도의 local Codex handoff를 선행조건으로 두지 않는다.
새 작업은 branch의 실제 HEAD와 `RESULTS_KO.md`, `campaign/RECEIPT.json`,
`comparison/COMPARISON.json`, `diagnosis/DIAGNOSIS.json` 및 `deliverable/CATALOG.json`부터 읽는다.
이번 완료 cohort를 다시 감사하는 대신 새로 바뀐 입력·지원구간·모델 선택만 계산한다.
