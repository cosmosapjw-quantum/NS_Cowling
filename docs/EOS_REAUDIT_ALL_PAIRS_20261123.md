# NS_Cowling EOS generator·core–crust 전 조합 재감사 — 2026-11-23

> **상태: v2 전수 계산·별도 실패 판정·전체 CI·전달 검증 완료.**
> **Core source 32개 × crust source 73개 × 7방법 = 16,352개 cell**을 새로 계산하여 **9,377개 source-pair/method cell**이 현재 admission을 통과했다. 이 중 **8,045개 endpoint chain / 1,332개 finite component**다. 전체 CI는 **1,907 passed / 35 skipped / 4 subtests passed**로 완료했다. 과거 v1 15,190개 / 8,712개 결과는 별도 기준 실행으로 보존한다. 이번 턴의 공개 data·runtime evidence·private input 백업은 로컬 Dropbox 경로에 기록하고 hash를 확인했다. 게시 버전은 이 문서를 포함하는 research branch와 최종 PR head로 식별한다(§9).

## 1. 무엇을 검증했고 무엇을 고쳤는가

EOS 생성기의 입력·모델 유효범위, raw table의 물리량·질량 기준 해석, core–crust matching의 유한 후보 탐색, bridge 내부의 안정성·인과성, 전체 출력과 serialization의 일관성을 다시 감사했다. 발견한 구현 결함을 수정한 뒤, 과거 matched output을 재사용하지 않고 각 source identity의 Cartesian product를 새 native Fortran 계산에 넣었다. Python은 CLI·입출력·작업 분배·검증 보고·무손실 패키징을 담당한다. 생산 EOS 계산, 변환, bridge 구성, root search와 native physical verification은 Fortran에서 수행한다.

최종 v2에서는 **9,377개 matched artifact가 정해진 admission을 통과**했다. 이 가운데 **8,045개는 선택된 source domain의 양 끝을 연결하는 결과**, **1,332개는 유한 phase component에 한정된 결과**다. 두 집합은 실제 사용 목록에서도 구분한다. 이 숫자는 **source identity × method cell의 수**이며, 수학적으로 동일한 EOS 함수를 semantic deduplication한 개수가 아니다. 어느 쪽도 이 결과만으로 표면부터 중심까지의 전역 EOS, TOV 안정성 또는 Cowling 고유진동 검증을 의미하지 않는다.

| 항목 | 상태 | 근거 |
|---|---|---|
| Generator·raw importer·matching·verifier의 이번 확인 결함 수리 | 지정 회귀에서 검증됨 | §3 및 검증 evidence |
| v1 전 조합 실행 | 완료: 15,190 / 15,190 cell 기록 | `campaign_v1/SUMMARY.json` |
| v1 생존 분류 | 7,434 endpoint chain + 1,278 finite component | v1 method별 coverage 합계 |
| 1 worker / 8 workers 비교 | 112 cell에서 6.9845903배, 56개 생존 EOS byte identity | §7 |
| v1 공개 payload 무손실 압축 | shard payload 기준 10.8381배 축소 | §8 |
| **v2 최종 실행** | **완료 — 16,352개 cell, 9,377개 생존** | `campaign_v2/SUMMARY.json`; §5 |
| **v2 생존 범위 / 공개 전달** | **8,045 endpoint chain + 1,332 finite component; 공개 8,972개 보존** | v2 CSV 목록 / `delivery_v2/`; §8 |
| **전체 CI** | **PASS — 1,907 passed / 35 skipped / 4 subtests passed** | 최종 native/importer freeze 후 실행; §3.4 |
| **이번 턴 백업** | **로컬 기록·hash 확인 완료** | 공개 v1/v2 tar와 runtime/private-input archive; §9 |
| **게시 버전 식별** | 이 문서를 포함한 research branch / 최종 PR head | 실제 commit·PR·remote readback은 외부 `SHIPPING_RECEIPT` 및 최종 전달 기록 |

수학적 유도, 원전, source column·mass convention, 구체적인 실패 재현은 [물리·source 상세 감사](EOS_PHYSICS_SOURCE_AUDIT_20261123.md)에 기록되어 있다. Native matching의 authoritative contract는 [Fortran matching 문서](FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md)다.

## 2. 범위와 source identity

### 2.1 원본을 유지한 새 Cartesian 실행

개별 cell은 **crust source identity × core source identity × method**다. 같은 핵모형 이름이라도 provider original, 보정본, 과거 재구성본, secondary-literature table은 서로 다른 입력 identity다. 따라서 “core 32개”를 “독립된 핵물리 이론 32개”로 해석하지 않는다. 반대로 단순 grid-resolution variant를 새로운 핵모형으로 부풀리지 않는다.

원 raw file과 기존 EOS snapshot은 수정하지 않는다. SHA256, 출처, 선택한 원본 행 범위, 단위, energy convention, source-declared phase 또는 문서화된 geometry 해석을 보존한다. Raw import에서 압력을 spline gradient로 대체하거나, 가공의 zero-density anchor를 삽입하거나, 실패 구간을 제거해 통과 가능한 source를 만들지 않는다.

Source preparation은 미리 정한 domain과 phase view를 materialize하는 과정이다. Physical admission은 이후 native 계산에서 수행한다. `available=true` 또는 import의 `NATIVE_RAW_PROJECTION_READBACK_VERIFIED`는 데이터가 준비되고 읽힌다는 뜻이며, 물리적으로 허용됐다는 뜻이 아니다.

입력 snapshot의 과거 생성 provenance도 그대로 유지한다. 예를 들어 일부 역사적 stationary-RMF snapshot의 Python oracle provenance를 이번에 새로 구현한 Fortran generator 결과로 바꾸어 기술하지 않는다. 이번 계산에서 사용하는 kernel·새 raw conversion과 입력의 역사적 생성 경로는 별도로 추적한다.

### 2.2 v1와 확대된 v2

| 실행 | Core source | Crust source | Source pair | 방법 수 | 요청 cell | 결과 상태 |
|---|---:|---:|---:|---:|---:|---|
| v1 기준 실행 | 31 | 70 | 2,170 | 7 | 15,190 | 완료 |
| v2 최종 범위 | 32 | 73 | 2,336 | 7 | 16,352 | **완료** |

v1 catalog는 [SOURCE_CATALOG.json](../data/eos_matching_reaudit_20261123/SOURCE_CATALOG.json)에 보존한다. 확대 catalog는 별도 `data/eos_matching_reaudit_20261123/SOURCE_CATALOG_EXPANDED.json`을 사용하며, v1 결과·입력 catalog를 덮어쓰지 않는다.

추가 네 identity의 preparation과 source freeze는 완료됐다. Expanded catalog는 **105개 기록 = core 32 + crust 73**, 이 가운데 **102개 available / 3개 blocked**, source phase component는 192개다. 기존 101개 source record는 JSON-value identity를 유지했고, v1 catalog·summary·delivery index의 SHA256도 변하지 않았다. Preparation의 완료는 물리적 생존을 뜻하지 않는다.

| 추가 identity | 역할 | 원 행 → adopted 행 | Adopted n, fm⁻³ | provenance 경계 |
|---|---|---:|---|---|
| `secondary_literature:BBP_Diener2008_raw:crust` | Crust | 35 → 26 | 0.000267–0.0747 | Secondary-literature 전재본 |
| `historical_reconstructed:BBP_pressure_primary:crust` | Crust | 26 → 26 | 0.000266583888518–0.0747 | GM1 representation에 보존된 pressure-based 재구성 prefix |
| `provider_attached:APR_outer_prefix:crust` | Crust | 172 → 49 | 7.92405959867×10⁻¹⁵–0.000257211369325 | Provider EOS에 붙어 있는 outer prefix |
| `external_nuda:2020-MVCD-D1MS-full:core` | Core | 65 → 52 | 0.0838–1.41 | Full table에서 문서화된 core marker 이후 |

D1MS-full의 앞 13개 crust 행은 기존 D1MS crust와 정확한 alias이므로 새 crust record로 중복 등록하지 않았다. 추가 identity·row selection·alias 및 v1 보존 증거는 [SOURCE_IDENTITY_PROOF.json](../data/eos_matching_reaudit_20261123/additional_sources/SOURCE_IDENTITY_PROOF.json)에 있다. Expanded catalog SHA256은 `59a30065b31df341b0b0f03c32784fdd64318f1e855f98af78414f0cce9ee085`다.

이 추가는 기존의 primary-source availability 문제를 없애지 않는다. 다음 세 기록은 계속 blocked 상태로 분리한다.

| 기존 identity | 상태 | 이유 |
|---|---|---|
| `external_nuda:1973-Negele-Vautherin:crust` | `BLOCKED_MISSING_PRESSURE` | 필요한 pressure가 없는 원 source |
| `historical:BPS_outer_only:crust` | `BLOCKED_NO_SELF_CONTAINED_THERMO_TABLE` | 독립적으로 재현할 thermo table 미확보 |
| `historical:BBP_inner_only:crust` | `BLOCKED_NO_SELF_CONTAINED_THERMO_TABLE` | 독립적으로 재현할 primary thermo table 미확보 |

Secondary BBP 자료나 명시적 재구성본이 존재한다는 이유로 마지막 primary-source 기록을 성공으로 바꾸지 않는다.

## 3. 확인한 결함과 수리

### 3.1 생성기·실행 파일·raw import

| 확인된 문제 | 변경과 검증 의미 |
|---|---|
| 기본 native executable이 존재하면 source 변경 여부를 확인하지 않고 재사용했고, 실행 불가능한 파일도 받아들임 | Frontend가 기존 `build.sh --ensure`의 source/compiler/flags signature와 process lock을 이용한다. 반복 호출의 metadata/environment cache는 관련 변경 시 무효화한다. Explicit executable override도 실제 실행 가능한 파일을 요구한다. |
| Public native core API가 자체 선언한 density domain 밖에서도 성공 | Single-state와 bulk entry에서 finite positive density 및 `core_model_bounds`를 강제한다. APR [0.085,1.34], BBB [0.08,1.5], GM1 [0.02,1.0], KDE0v1 [0.02,1.8], SLy4 [0.02,1.4] fm\(^{-3}\)의 기존 선언 범위를 사용한다. 이는 matching을 위한 새 관측 제약이 아니다. |
| Geometric-grid test fixture의 요청 endpoint 1.4가 1.4000000000000001로 생성 | 요청한 양 끝 literal을 정확히 사용한다. 모델 범위나 assertion을 완화하지 않는다. |
| Fortran scalar list-directed parser가 `939 extra`, 복수 token, repetition syntax를 조용히 받아들임 | 질량·선택 domain 인자를 정확히 하나의 finite numeric token으로 검사한다. 잘못된 인자는 output 생성 전에 거부한다. |
| 외부 raw table의 energy·pressure convention 혼동 위험 | 별도 native importer가 exact schema, 모든 원 행의 finite/order 조건, 명시적 domain selection, 정규화와 native readback을 처리한다. 열 추측·정렬·보간·압력 수리·자동 rest-energy shift는 없다. |

Generator/domain 감사는 관측된 **RED 23 failed / 5 passed → GREEN 28 passed**를 보존한다. 관련 기존 회귀를 합친 지정 실행은 **49 passed**였다. Raw importer는 NUDA 세 layout과 사용자 total-energy 두 layout을 포함한 지정 실행에서 **37 passed**를 기록했고, 네 추가 source의 layout·selection contract를 반영한 최종 지정 실행은 **RED 8 failed / 40 passed → GREEN 48 passed**였다. 최종 importer executable SHA256은 `90f24634e20ea73e9629eb6117c85068fc6f3a45575a2af66280780a444a37ee`이며 main matching executable은 v1과 동일한 `d1600879…`다. 이 수치는 해당 테스트 실행의 증거이며 전체 CI 결과와 구별한다.

BSk22/24/25/26에서 upstream loader가 energy/baryon으로 이름 붙인 **column 11은 g/cm³ mass density**였다. 올바른 ETFSI energy offset은 **column 10**, pressure는 column 8이다. Native conversion은 source pressure를 유지하고, 원전·provider header로 확정한 source별 rest-energy convention을 적용한다.

| 원 source 계열 | \(n\) / energy offset / \(P\), 1-based | 더하는 neutron-rest energy |
|---|---|---:|
| GMRS 15개 | 1 / 8 / 9 | 939.565346 MeV |
| Gogny 3개 | 1 / 5 / 6 | 939.0 MeV |
| BSk 4개 | 1 / 10 / 8 | 939.5651828 MeV |
| 사용자 total-energy-density 자료 | layout별 명시 | **0 MeV**; 이미 전체 \(\epsilon\)인 source |

Gogny의 939 MeV는 [Mondal et al. Table II](https://arxiv.org/html/2003.03338v3)의 저자 convention으로 확인했다. BSk는 [수정된 Pearson et al. 문헌](https://arxiv.org/html/1903.04981)과 실제 CompOSE header로 확인했다. Pinned NUDA BSk24/25와 현재 provider의 일부 행은 서로 다른 revision이므로 원본을 교체하지 않는다.

### 3.2 Matching과 verifier

| 확인된 문제 | 변경과 검증 의미 |
|---|---|
| Historical heuristic window 또는 첫 bridge 후보 실패로 source pair를 너무 일찍 거부 | 기존 선호 조건은 탐색 순서로만 사용하고, 최소 retained-knot 수를 충족하는 **유한 supplied endpoint grid 전체**를 fallback 탐색한다. Examined/constructed/rejected/exhausted counter를 남긴다. |
| Near-plateau first-law bridge를 401개 sampling만으로 판정 | \(P(\ln n)\) cubic에 대한 정확한 first-law 해와 exponential-minus-quadratic causal margin을 사용하여 endpoint와 모든 내부 critical point를 검사한다. 표본에서는 통과하나 내부에서 인과성을 위반하는 반례를 native certificate가 거부한다. |
| 동일 endpoint의 repeated derivative 계산과 불필요한 bridge 구성 | Retained physical side의 same-phase endpoint derivative를 cache하고, source prefix/suffix의 기존 positivity·monotonicity·causal-secant 조건과 필요한 causal integral bound를 먼저 검사한다. |
| 긴 method identity의 truncation과 nonfinite constructor 입력 | Full method name을 보존하는 길이와 finite-domain entry 검사를 사용한다. |
| \(d\epsilon/dn\le0\)인 source segment를 mask에서 지워 나머지 점만으로 PASS | Pointwise physical domain과 regular-phase invalid derivative를 별도 실패로 기록한다. `DERIVATIVE_DOMAIN`이 전체 PASS를 막는다. |
| `Il=0`, 고정 \(Y_q\)에서 raw \(\mu_B\)를 density-path potential로 사용 | \(\mu_{\rm eff}=\mu_B+Y_q\mu_q\)를 differential first law·sound speed·potential monotonicity에 사용한다. CompOSE의 raw Q3 algebraic convention 검사는 유지한다. |
| Phase/gap 경계 가까운 endpoint의 3-point one-sided stencil이 경계를 가로질러도 checked로 남음 | 경계를 가로지르는 endpoint stencil까지 mask한다. 검증할 derivative point가 0개인 경우를 PASS로 만들지 않는다. |

Native verifier 전용 **16 assertion**은 bounds/FPE trap debug build와 optimized build 모두 통과했고, 자동 pytest wrapper는 **2 passed**였다. 새 matching native 두 프로그램은 같은 두 build mode에서 **4 pytest case 통과**를 기록했고 기존 adversarial native controls도 실행했다.

기존 Maxwell 회귀의 coarse 6-knot fixture는 이 경계 mask 아래에서 유효 derivative stencil이 없어졌다. 동일한 analytic \(P(\mu)\)와 transition을 유지하는 31-knot fixture로 positive case를 충분히 해상했고, coarse zero-coverage case는 별도 negative test로 보존했다. Production threshold 또는 물리식을 이 테스트를 위해 바꾸지 않았다. 관련 focused tests는 **45 passed**이며, 이후 전체 frontend/pipeline CI도 §3.4의 조건에서 통과했다.

### 3.3 Campaign admission·provenance·개인 자료

각 cell은 native matching receipt, internal consistency, CompOSE write/read, 모든 serialized column·metadata의 readback, written CompOSE 재검증을 통과해야 한다. CSS/Maxwell에는 별도 serialized interface admission이 추가된다. 실패 source, 입력 계약, 물리 조건, interpolation shape, 구현 지원 범위, runtime/protocol 오류를 같은 “물리 불가능”으로 합치지 않는다.

**현재 `check_consistency`의 전체 PASS는 C1·C2·해당 composition 조건 C3/C4·C5·C6·C10·C11 및 jump/domain 조건을 요구한다. C7·C8·C9의 finite-difference first-law residual은 non-veto diagnostic이다.** 즉 source 전체의 연속 first law가 작은 오차로 성립한다는 증명은 admission 정의에 포함되지 않는다. Bridge construction에는 별도의 first-law endpoint closure residual `<=5e-9` 및 continuous shape/stability/causality 조건을 사용한다. 새 bridge의 analytic 조건과 보존한 source table의 finite verification을 구별한다.

Runner는 모든 요청 ID를 기록하고 완료 시 Cartesian 집합과 대조한다. 실행 전후와 resume 시 source·catalog·native binary·runner·admission helper identity를 확인한다. 바뀐 identity가 있으면 canonical survivor 결과로 게시하지 않는다. Private artifact뿐 아니라 private failure trace도 공개 output과 분리한다. Readback 오류에 원 수치가 포함될 수 있으므로 공개 상태표에는 제한된 metadata만 남긴다.

Source별 outcome 표의 `contains_private_pair`는 그 source가 private input과 짝지어진 결과를 포함한다는 뜻이다. 해당 public source 자체가 private이라는 뜻으로 해석하지 않는다.

### 3.4 최종 frontend·pipeline CI

최종 native executable과 importer 및 source preparation을 freeze한 뒤 다음 명령을 실제로 실행했다.

```sh
PYTHONPATH=python_frontend \
COMPOSE_EOS_FORTRAN_EXE=/tmp/ns_eos_reaudit_20261123/build/compose_eos_fortran \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m pytest -q python_frontend/tests pipeline/tests
```

결과는 **1,907 passed / 35 skipped / 4 subtests passed**, return code 0이다. Pytest 보고 시간은 **277.03 s**, wrapper wall time은 **277.638122586 s**이며, 실제 실행 구간은 **2026-11-23 06:54:31–06:59:09 UTC**다. 건너뛴 35개 test를 통과한 test로 계산하지 않는다. Native executable과 importer의 실행 전후 SHA256은 각각 `d1600879…`, `90f24634…`로 동일했다.

Receipt는 `/tmp/ns_eos_reaudit_20261123/CI_FINAL_FRONTEND_PIPELINE_RECEIPT.json`이다. 실제 log SHA256은 `88092ea3532f5fa8c9c9235d3523ec6540db95fe043f783fd8a4372b5eb396f0`다. 이 실행은 frontend/pipeline suite의 구현 회귀 검증이며, 위 source-domain 및 물리적 유효범위의 한계를 확장하지 않는다.

## 4. 일곱 방법과 “생존”의 의미

\(n\)은 baryon density, \(\epsilon\)은 rest energy를 포함한 total energy density, \(F=\epsilon/n\)은 cold energy per baryon이다. 단위는 fm\(^{-3}\), MeV fm\(^{-3}\), MeV이고, native `cs2`는 \(c_s^2/c^2=dP/d\epsilon\)이다.

| Method ID | 구성 |
|---|---|
| `first_law_constrained_hermite_pressure_bridge` | \(P(\ln n)\) cubic와 정확한 first-law energy 적분 |
| `cubic_F` | Free-energy 표현의 cubic bridge |
| `quintic_F` | Endpoint derivative 조건을 포함하는 quintic free-energy bridge |
| `cubic_P` | Pressure–chemical-potential 표현의 cubic bridge |
| `quintic_P` | Pressure–chemical-potential 표현의 quintic bridge |
| `css` | 하나의 constant-sound-speed phase; crust 쪽 한 번의 density jump와 core 쪽 연속 접합 |
| `maxwell` | 두 supplied phase의 pressure·chemical-potential 공존점과 양의 density/energy jump |

Cold one-charge 조건에서
\[
d\epsilon=\mu\,dn,\qquad \epsilon+P=n\mu,\qquad dP=n\,d\mu
\]
를 사용한다. Smooth bridge에는 positive thermodynamic domain, stability, causality와 physical endpoint consistency를 요구한다. Genuine Maxwell interface에는 공통 \(P_t,\mu_t\) 및
\[
\Delta n>0,\quad \Delta\epsilon>0,\quad
\Delta\epsilon=\mu_t\Delta n
\]
가 필요하다. Seidov criterion은 stellar small-core stability의 별도 진단이며 matching의 EOS veto로 사용하지 않는다. 관련 causal endpoint constraints는 [Fortin et al.](https://arxiv.org/pdf/1604.01944)과 일치한다.

**생존은 source-pair/method별로 정해진 admission을 처음 통과한 하나의 deterministic output을 보존했다는 뜻**이다. 가능한 모든 transition density, 모든 연속 interpolation, 모든 phase assignment를 열거했다는 뜻이 아니다. Bridge의 analytic certificate와 retained source의 finite-table verification도 별도 주장이다. 특히 C7–C9가 diagnostic이라는 현재 계약 때문에 `SURVIVED`를 retained source 전체의 continuous first-law 또는 global interpolant certificate로 읽으면 안 된다.

최종 사용자 목록은 다음과 같이 분리한다.

| 파일 | 포함하는 결과 | 해석 |
|---|---|---|
| `campaign_v2/ENDPOINT_CHAIN_SURVIVORS.csv` | `SELECTED_SOURCE_ENDPOINT_CHAIN` | Adopted crust의 시작과 adopted core의 최고 density를 유지한 연결체 |
| `campaign_v2/FINITE_COMPONENT_SURVIVORS.csv` | `FINITE_PHASE_COMPONENT_ONLY` | 지정 phase component만으로 성립하는 local matched artifact |
| `campaign_v2/ALL_CELLS.csv` | 모든 요청 cell | 원 실행의 성공·탈락·blocked accounting |
| `campaign_v2/FINAL_CELL_CLASSIFICATION.csv` | 별도 판정을 반영한 모든 cell | Fixed-prefix rejection과 미지원 multiple-root case를 분리 |
| `campaign_v2/SURVIVORS.json` | Admission 통과 전체 | 상세 method·view·coverage·artifact provenance |

위 v2 목록의 최종 수치는 endpoint chain **8,045개**, finite component **1,332개**다. Endpoint-chain 결과도 **선택된 입력 domain의 양 끝**이라는 범위를 갖는다. 선택한 crust가 inner crust부터 시작하면 stellar surface까지 포함하는 EOS가 되지 않는다. Stellar integration, maximum mass, radial/Cowling mode stability, 관측 제약 또는 global free-energy ground state는 이 campaign으로 증명하지 않는다.

첫 certified bridge 이후 whole-table gate가 실패하면 같은 phase pair의 다른 endpoint 후보가 남아 있을 수 있으므로 `POST_BRIDGE_GATE_FAILURE_SEARCH_INCOMPLETE`를 구별한다. V1에는 해당 outcome이 0개였다. V2에는 secondary BBP raw source에서 **145 trials / 135 cells**가 기록되었고, 이후 모든 대체 endpoint에서도 변하지 않는 retained first-point C11 위반을 별도 Fortran certificate로 확인했다(§6.3). 원 journal의 미완결 분류를 삭제하지 않고 별도 adjudication receipt가 해소 근거를 기록한다. Multiple coexistence roots는 현재 구현 지원 범위로 분류하며, v1과 v2에 각각 남은 **5건**을 보편적 물리 불가능으로 해석하지 않는다.

## 5. 완료된 v1 기준 결과

v1은 8 workers, native child별 BLAS/OpenMP thread 1, bridge output 80개 interior point, operation timeout 180 s로 실행했다. **284.906397796 s**에 15,190개 cell을 기록했고, source·binary identity가 실행 전후 동일했다.

| 방법 | 요청 | 생존 | Endpoint chain | Finite component | No admitted match | Source blocked |
|---|---:|---:|---:|---:|---:|---:|
| First-law pressure bridge | 2,170 | 1,408 | 1,197 | 211 | 669 | 93 |
| `cubic_F` | 2,170 | 1,452 | 1,239 | 213 | 625 | 93 |
| `quintic_F` | 2,170 | 1,452 | 1,239 | 213 | 625 | 93 |
| `cubic_P` | 2,170 | 1,452 | 1,239 | 213 | 625 | 93 |
| `quintic_P` | 2,170 | 1,452 | 1,239 | 213 | 625 | 93 |
| `css` | 2,170 | 1,452 | 1,239 | 213 | 625 | 93 |
| `maxwell` | 2,170 | 44 | 42 | 2 | 2,033 | 93 |
| **합계** | **15,190** | **8,712** | **7,434** | **1,278** | **5,827** | **651** |

이 수치는 완료된 v1 baseline으로 보존한다. V2에 추가된 source를 v1 생존 개수에 단순히 더하지 않고, 확대된 catalog로 전체 Cartesian matrix를 다시 계산한다.

### v2 최종 결과

| 방법 | 요청 | 생존 | Endpoint chain | Finite component | No admitted match | Source blocked |
|---|---:|---:|---:|---:|---:|---:|
| First-law pressure bridge | 2,336 | 1,512 | 1,292 | 220 | 728 | 96 |
| `cubic_F` | 2,336 | 1,564 | 1,342 | 222 | 676 | 96 |
| `quintic_F` | 2,336 | 1,564 | 1,342 | 222 | 676 | 96 |
| `cubic_P` | 2,336 | 1,564 | 1,342 | 222 | 676 | 96 |
| `quintic_P` | 2,336 | 1,564 | 1,342 | 222 | 676 | 96 |
| `css` | 2,336 | 1,564 | 1,342 | 222 | 676 | 96 |
| `maxwell` | 2,336 | 45 | 43 | 2 | 2,195 | 96 |
| **합계** | **16,352** | **9,377** | **8,045** | **1,332** | **6,303** | **672** |

적어도 한 방법이 생존한 source pair는 **1,564개**, 여섯 bridge 방법이 모두 생존한 pair는 **1,512개**, 일곱 방법이 모두 생존한 pair는 **45개**다. Method별 생존 개수가 같더라도 output EOS 함수가 같다는 뜻은 아니다.

| 항목 | v2 값 |
|---|---|
| Catalog | `SOURCE_CATALOG_EXPANDED.json`, SHA256 `59a30065b31df341b0b0f03c32784fdd64318f1e855f98af78414f0cce9ee085` |
| 요청 / 기록 | **16,352 / 16,352**, Cartesian accounting 완료 |
| Public / private 생존 | **8,972 / 405**; private numerical payload는 공개 package에서 제외 |
| Wall time | **297.7550749309885 s**; CI 동시 실행으로 통제된 speed comparison이 아님 |
| Runtime / protocol failure | **0** |
| Native executable SHA256 | `d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588` |
| Runner SHA256 | `063177362d8fdc19bec6777723f5203362e4de8295512bf87cf2492a6b62062c` |
| Source·binary·runner identity | 실행 전후 변경 없음 |
| Post-bridge 재분류 | 원 145 trials / 135 cells를 별도 certificate로 해소; **남은 post-bridge search-incomplete 0**; 원 journal 유지 |
| Delivery index SHA256 | `d9695d8b08bc4ea82c0ba6e815d49bf8d4b53bc972ddbeca5e27ab79a86a06e2` |

[FINAL_CENSUS.json](../data/eos_matching_reaudit_20261123/campaign_v2/FINAL_CENSUS.json)의 최종 분류는 다음과 같다.

| 최종 분류 | Cell 수 |
|---|---:|
| Selected-source endpoint chain 생존 | 8,045 |
| Finite phase component 생존 | 1,332 |
| 그 외 현재 admission 미통과 | 6,163 |
| Fixed retained-prefix C11 rejection | 135 |
| 미지원 multiple coexistence roots | 5 |
| Source 자료 미확보 / blocked | 672 |
| **합계** | **16,352** |

[FINAL_ADJUDICATION.json](../data/eos_matching_reaudit_20261123/campaign_v2/FINAL_ADJUDICATION.json)은 해당 135 cells/145 trials와 원 trace SHA를 certificate에 연결한다. 원 `RESULTS.jsonl`, `SUMMARY.json`, survivor set은 수정하지 않았다. **Endpoint-search 미해소 trial은 0개**, runtime/protocol failure도 0개다. Multiple-root 5개는 여전히 구현 지원 범위의 한계로 남긴다.

Trial failure counts는 한 cell 안의 여러 prescribed phase pair를 포함하므로 cell 수와 합산하지 않는다. `NO_ADMITTED_MATCH`에는 source/input contract와 해당 interpolation 표현의 실패가 포함되며 일반적인 물리적 비존재 정리가 아니다.

## 6. 실제 탈락 사례를 원인별로 확인한 결과

### 6.1 MYN13 → OMEG1

사용자 제공 두 source의 adopted domain은 MYN13 crust 182행, \(n=[6.0222\times10^{-11},0.06757]\) fm\(^{-3}\), OMEG1 core 6,243행, \(n=[0.079,3.2]\) fm\(^{-3}\)다.

원 MYN13 table에는 \(\mu_E=(\epsilon+P)/n\) 감소 구간이 **24개** 있고, 첫 감소가 처음 두 행에서 발생한다. 첫 \(\Delta\mu_E=-0.0028688916867167791\) MeV다. 모든 eligible retained crust prefix가 이 첫 구간을 포함하므로 통과 가능한 prefix endpoint가 0개다. OMEG1은 채택 domain의 모든 6,242 native interval을 통과하고, causal secant 최대값도 약 0.434022다.

여섯 non-Maxwell 방법은 각각 **1,123,380개 endpoint 쌍**을 검사했지만 source gate를 통과하지 못해 bridge를 만들지 않았다. Source gate를 진단 목적으로 제외했을 때 endpoint order와 necessary causal integral bound는 그 모든 쌍에서 통과했다. Maxwell도 low-density source의 monotone chemical potential 조건에서 실패했다. 따라서 이번 7방법 탈락의 직접 원인은 **adopted MYN13 원표의 비단조성**이며, OMEG1 tail 또는 특정 bridge family에 대한 실패 정리가 아니다.

원 decimal 출력 자릿수로 계산한 sensitivity 검사에서 ±0.5 last-place envelope에는 monotone feasibility 위반 2개가 남고 ±1 last-place envelope에서는 0개가 남는다. 이는 source precision에 민감함을 보여준다. 원 생성기의 오차와 rounding convention을 확정한 것은 아니며, 더 넓은 envelope를 admission tolerance로 사용하지 않았다. 원 숫자, source origin, phase label과 threshold를 바꾸지 않고 rejection을 유지했다.

### 6.2 GPPVA_NL3wrL55 → PK1_CANONICAL, Maxwell

V1 `cell_012222`는 native matching에 성공했지만 final C2만 실패했다. 원 provider의 density index 218이 그대로 남은 source/output 첫 행에서
\[
\mu_{\rm Q3}-\mu_{\rm Euler}=-0.0030235134198619562\ {\rm MeV},
\]
relative residual \(3.2139851048496636\times10^{-6}\)을 확인했다. 이 값은 provider input gate \(5\times10^{-4}\)보다 작지만 최종 C2 gate \(2\times10^{-6}\)보다 크다. C10·C11은 통과했고 invalid derivative point는 0개다.

이는 source algebraic consistency rejection이다. 출력 serialization, Maxwell 경계 stencil 또는 finite-difference false failure로 처리하지 않았다. 원 Q3를 Euler 값으로 대체하거나 최종 tolerance를 완화하지 않았다.

두 사례의 source SHA, adopted source indices, native executable identity와 Fortran 진단 증거는 [상세 감사 §5](EOS_PHYSICS_SOURCE_AUDIT_20261123.md)에 있다. Private table과 전체 private interval 수치는 공개 문서·Git에 포함하지 않는다.

### 6.3 Secondary BBP raw의 145개 post-bridge 실패

추가한 `secondary_literature:BBP_Diener2008_raw:crust`에서 **135 cells의 145 native trials**가 bridge 구성 후 C11만 실패했다. 모든 경우가 최소 crust 3행을 보존했고, 다른 endpoint가 남아 있을 가능성을 정직하게 기록하기 위해 원 campaign에는 search-incomplete 분류를 남겼다.

추가 독립 Fortran 검산은 처음 세 행이 모든 eligible retained prefix에 남는다는 사실을 이용했다. 첫 세 \(\mu_E\) 값 자체는 약 **938.8768045, 938.9861022, 940.1564153 MeV로 증가**하지만, 현재 verifier의 3-point one-sided quadratic derivative는 첫 점에서 약 **−1899.11413469 MeV fm³**다. 이것은 “처음 두 tabulated \(\mu\)가 감소한다”는 주장과 다르다.

첫 두 \(\Delta\ln n\)은 0.1589545와 0.2326826으로 고정 gap threshold 0.25보다 작고, phase ID는 모두 0이다. Pressure relative differences 0.05509와 0.08242도 equal-pressure jump trigger \(10^{-10}\)와 명백히 다르다. 그러므로 첫 점은 gap·phase·Maxwell boundary mask에서 빠지지 않는다. Bridge가 시작되는 지점은 최소 세 번째 source 행 뒤이므로 이 첫 stencil을 바꾸지 않는다.

실제 실패에 사용된 **25개 target core view의 mass header·Il convention**으로 다시 encoding해도 첫 \(d\mu/dn\)은
[−1899.1141346935183, −1899.1141346897930] 범위에 남았다. 첫 \(dP/dn,dF/dn,d\epsilon/dn,d\mu/dn\)은 모두 finite이고 \(d\epsilon/dn\approx937.1792\) MeV \(>0\)이므로 invalid-derivative mask도 이 점을 제거하지 않는다. 24개 eligible crust prefix를 native verifier로 직접 확인한 결과도 모두 C11 실패였다. 증명의 핵심은 각 prefix의 마지막 stencil이 아니라 **변하지 않는 첫 stencil**이다.

이에 따라 **현재 immutable source·prescribed view·native C11 operator 아래에서 나머지 endpoint 후보도 통과할 수 없음**을 인증했다. 해당 135 cells는 계속 `NO_ADMITTED_MATCH`이고, 145 trials의 미완결 분류만 별도 `FIXED_RETAINED_PREFIX_C11_REJECTION` certificate로 해소한다. 새로운 physical gate, source 수정, derivative operator 변경 또는 tolerance 완화는 하지 않았다.

별도의 cold first-law 검산에서는 안정한 one-charge branch에 필요한
\(n_1\Delta\mu\le\Delta P\le n_2\Delta\mu\)를 첫 두 source interval이 모두 위반했다. 첫 interval의 \(\Delta P=2.87109417\times10^{-5}\)는 하한 \(2.91824901\times10^{-5}\)보다 작고, 두 번째의 \(4.68113181\times10^{-5}\)는 하한 \(3.66307997\times10^{-4}\)보다 작다. 이 계산은 원 source의 약 3자리 출력 정밀도에 대한 별도 진단이며, 생산 admission에 새 제약을 추가한 것이 아니다. 다른 표현에서 underlying BBP model이 불가능하다는 결론으로 확장하지 않는다.

Certificate는 [audit_evidence/verifier/bbp_raw_prefix_repro/CERTIFICATE.json](../data/eos_matching_reaudit_20261123/audit_evidence/verifier/bbp_raw_prefix_repro/CERTIFICATE.json)이며 SHA256은 `7172f204393c63e6f70e2ab577e950a51d86ed61c29a57518bdab2cf1f73820e`다. `PROOF.md`, `inspect_bbp_prefix.f90`, `NATIVE_PREFIX_CHECK.txt`와 정확한 135 cell/145 trace identity mapping을 함께 보존한다.

## 7. 병렬 계산 성능과 재현성

같은 source·native binary를 사용하는 16 pair × 7 methods = **112개 cell**을 1 worker와 8 workers로 각각 실행했다. 표본은 source pair 14개를 균일 간격으로 고르고 실제 Maxwell 생존 pair의 처음·마지막을 포함했다.

| 측정 | 값 |
|---|---:|
| 1 worker | 9.626607915997738 s |
| 8 workers | 1.3782637919939589 s |
| 실측 speedup | **6.984590302608728배** |
| 생존 cell | **56개** |
| Status·selected endpoint | 완전 일치 |
| 생존 EOS CSV와 NML | **byte-for-byte 동일** |
| Input·binary hash | 변경 없음 |

이는 이 host와 표본에서 한 번의 paired experiment로 얻은 실측이다. 모든 규모에서 같은 가속률이 나온다거나, timing·temporary path를 포함하는 전체 archive의 compressed bytes까지 같다는 주장은 하지 않는다. **V2 전체 실행은 최종 CI와 동시에 진행하므로, v2 wall time을 v1 또는 이 통제된 benchmark와 직접 비교하는 성능 실험으로 해석하지 않는다.**

계산 중에는 8개 독립 native process를 사용하고 각 process의 BLAS/OpenMP thread를 1개로 제한한다. Endpoint derivative cache와 prefix/suffix 검사로 반복 계산을 줄인다. 물리 cutoff를 느슨하게 하거나 Python physical prefilter로 후보를 먼저 제거하여 얻은 속도가 아니다.

Benchmark SSOT:
`/tmp/ns_eos_reaudit_20261123/parallel_benchmark/PARALLEL_EQUIVALENCE.json`.
SHA256:
`2c1401616584609aaeed96e0c3c312b72c3afb7be00e934e774e47066e17d60e`.

### 확대 실행과 기준 실행의 재현성

확대 catalog로 전부 새로 실행한 뒤 v1과 겹치는 **15,190 cells**를 source identity와 method로 대조했다. Status·selected endpoint·native report·verification·failure categories가 일치했고, 공통 생존 **8,712개**의 internal/CompOSE numerical file **60,984개가 byte-for-byte 동일**했다. 차이 0건이다. 이 비교는 archive provenance path와 elapsed time을 제외하며, 기존 output을 재사용했다는 뜻이 아니다.

근거는 [expansion_equivalence/EQUIVALENCE.json](../data/eos_matching_reaudit_20261123/expansion_equivalence/EQUIVALENCE.json)이며 SHA256은 `f4e5dcfb849be9ddf5316a9bcfd89d25ca805c4a5e97a072e0f50863c64af281`이다. 확대로 추가된 1,162 cells는 생존 665개 / 미통과 476개 / blocked 21개로 accounting되었다.

## 8. 무손실 경량화와 결과 사용

### V1 기준 패키지

V1의 공개 생존 artifact **8,318개**를 33개 solid-XZ shard로 묶었고, private 생존 artifact **394개**는 공개 payload에서 제외했다. 개별 candidate의 원래 gzip을 풀었을 때의 **tar bytes 자체**를 보존한다. EOS 재계산·재보간·sampling 감소로 크기를 줄인 것이 아니다.

| V1 패키지 항목 | Bytes |
|---|---:|
| 원 public per-cell gzip 합계 | 717,997,877 |
| 원 decompressed candidate tar 합계 | 5,083,412,480 |
| Solid-XZ shard 33개 합계 | 66,247,512 |
| Delivery index | 6,455,803 |
| Shard + index 합계 | 72,703,315 |

Shard payload만 비교하면 **10.838110825958264배** 축소, 90.7733% 절감이다. Index를 포함한 전체 delivery는 **9.8757배** 축소다. 256 cell/shard, XZ preset 6, packaging workers 4를 사용했다.

Independent verification은 8,318개 원 gzip hash, 33개 shard hash, 모든 raw tar member hash·bytes를 다시 검사했다. Maxwell·CSS·first-law candidate를 실제로 추출해 내부 파일별 bytes와 SHA256도 원 archive와 비교했다. Private candidate 추출 요청은 거부하며 기존 output을 덮어쓰지 않는다. 이 검증은 물리 admission을 다시 계산하는 작업이 아니라 **완료된 candidate 내용의 무손실 전달 검증**이다.

V1 index는 [DELIVERY_INDEX.json](../data/eos_matching_reaudit_20261123/delivery_v1/DELIVERY_INDEX.json)에 보존한다. 그 SHA256은
`e2fbcc6c42433039fde00684d45ba5b3b8268e099507b2882c58b19656c5ea31`이다. V1 index와 기존 33 shard가 변하지 않았음을 v2 packaging 후 다시 확인했다. V1 개별 gzip 원본도 유지한다.

### V2 최종 패키지

최신 공개 numerical delivery는 [delivery_v2/DELIVERY_INDEX.json](../data/eos_matching_reaudit_20261123/delivery_v2/DELIVERY_INDEX.json)과 36개 shard다. 공개 생존 **8,972개**의 complete original tar를 포함하며, private 생존 **405개**의 공개 index에는 metadata만 들어 있다. 전체 private cells는 728개다.

| V2 패키지 항목 | Bytes |
|---|---:|
| 원 public per-cell gzip 합계 | 753,066,466 |
| 원 decompressed candidate tar 합계 | 5,327,052,800 |
| Solid-XZ shard 36개 합계 | 69,840,608 |
| Delivery index | 6,957,074 |
| **Shard + index 합계** | **76,797,682** |

Payload 기준 **10.78264476배** 축소(90.725837% 절감), index를 포함하면 **9.80584891배** 축소(89.802005% 절감)다. V2의 모든 8,972개 원 gzip과 raw tar, 36개 shard를 독립적으로 다시 읽어 hash와 크기를 확인했다. Maxwell `cell_001617`, CSS `cell_000230`, first-law `cell_000225`를 실제 extractor로 복원해 모든 내부 파일 bytes·SHA를 대조했고, private `cell_000442`는 output을 만들지 않고 거부함을 확인했다.

최종 index SHA256은 `d9695d8b08bc4ea82c0ba6e815d49bf8d4b53bc972ddbeca5e27ab79a86a06e2`다. Independent package verification SHA256은 `2d39299261cd08b9710131fceeb95a8864b84a089e27c1efce9ca2cdb6a4e2b3`다. 이 검증은 scientific admission과 구별되는 archive byte-preservation 검증이다.

### V2 공개 candidate 사용 절차

먼저 `ENDPOINT_CHAIN_SURVIVORS.csv`와 `FINITE_COMPONENT_SURVIVORS.csv` 중 필요한 범위를 선택하고 source ID·method·adopted domain을 확인한다. 다음은 실제 복원 검증을 마친 public Maxwell candidate의 사용 예시다. 다른 candidate는 v2 index에 존재하는 public cell ID로 지정한다.

```sh
python3 scripts/extract_eos_reaudit_candidate_20261123.py \
  --index data/eos_matching_reaudit_20261123/delivery_v2/DELIVERY_INDEX.json \
  --cell cell_001617 \
  --out /tmp/restored_eos_candidate
```

Extractor는 selected shard와 tar identity를 확인하고, 새 directory에 해당 candidate의 CompOSE·internal data·verification·provenance 파일을 복원한다. 선택한 output을 stellar calculation에 넣을 때에는 source-domain coverage와 phase-transition 지원 조건을 별도로 확인해야 한다.

## 9. 증거와 최종 마감 항목

### 고정된 v1 실행 identity

| 항목 | SHA256 |
|---|---|
| v1 source catalog | `5e162fb6ec998f06bdf0e88130d22c8391b0635608c56af929d0911e45d66919` |
| v1 native executable | `d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588` |
| v1 runner | `063177362d8fdc19bec6777723f5203362e4de8295512bf87cf2492a6b62062c` |
| v1 SUMMARY | `11e2a8c27efbff3e9f45095745d67e13ff8c1ebc1ba85562677e1952ba0133f0` |
| v1 RESULTS journal | `95556a53f2d3591cccbf65b6bcbd77d487c12316dc0ff680f965879f6304a9cc` |
| v1 SURVIVORS JSON | `e24dea8348f3e99943048d491399697f0088c91460c93030a7e6d34e56f52b6b` |

기준 Git commit은 `e51dadd132125f0c5389da8c92c50b84eeb524e8`, 작업 branch는 `research/eos-core-crust-all-methods-20261123`다. 최종 게시 버전은 **이 문서를 포함하는 해당 research branch와 최종 PR의 head**다. Commit 자기참조를 피하기 위해 최종 commit SHA·PR URL·remote readback은 문서 내부에 미리 쓰지 않고, 외부 `SHIPPING_RECEIPT`와 최종 전달 메시지에 기록한다.

V2 원 journal SHA256은 `c7b30f71ad74df5f9cfe6882fa84f3ee66c2f24fd7d764ad4a651c88aabbf27c`, 원 SUMMARY는 `8329318571475221d29b9c570cee95be114e6e643b969a9f1a77edc4bb49ab28`, 별도 FINAL_CENSUS는 `26ca8b9225a06ca8d6967f26424f65abb97ca4dc15c23508b4c1f90f916e0ff4`다.

### 검증 evidence 위치

| 범위 | 실제 evidence |
|---|---|
| Generator/domain | `/tmp/ns_eos_reaudit_20261123/generator/GENERATOR_AUDIT.md`, `RED.log`, `GREEN.log`, `FOCUSED_REGRESSION_FINAL.log` |
| Native raw import | `generator/IMPORTER_USER_ENERGY_GREEN.log`, `generator/NUDA_NATIVE_IMPORT_LAYOUTS.json` |
| Matching interval·search | `matching/AUDIT_MATCHING.md`, `matching/PYTEST_NATIVE_REGRESSIONS_FINAL.log`, `matching/MATCHING_FINAL_RECEIPT.json` |
| Verifier | `verifier/{red,green-debug,green-opt}/`, `verifier/pytest.log` |
| Maxwell fixture compatibility | `matching/CI_MAXWELL_FIXTURE_REPAIR.json`, `matching/CI_MAXWELL_FIXTURE_ALL_GREEN.log` |
| 지정 integration | `integration/INTEGRATION_AUDIT.md` — 200 passed, 1 skipped |
| 최종 전체 CI | `CI_FINAL_FRONTEND_PIPELINE_RECEIPT.json`, `CI_FINAL_FRONTEND_PIPELINE.log` — 1,907 passed / 35 skipped / 4 subtests passed |
| v1 성능 | `parallel_benchmark/PARALLEL_EQUIVALENCE.json` |
| v1 무손실 전달 | `delivery_evidence/DELIVERY_REPORT.md`, `VERIFICATION.json` 및 18-test delivery 회귀 |
| 구체 source 실패 | `verifier/cell_012222_repro/`, `verifier/user_pair_repro/PUBLIC_ASSESSMENT.json` |
| Fixed-prefix 전 후보 판정 | `verifier/bbp_raw_prefix_repro/CERTIFICATE.json`, `PROOF.md`, `NATIVE_PREFIX_CHECK.txt` |
| V2 무손실 전달 | `delivery_v2_evidence/DELIVERY_REPORT.md`, `VERIFICATION.json` |

위 표는 원 실행의 `/tmp/ns_eos_reaudit_20261123/` 기준 경로다. 공개 가능한 증거의 보존 사본은 [audit_evidence/](../data/eos_matching_reaudit_20261123/audit_evidence/)에 있고, [EVIDENCE_COPY_MANIFEST.json](../data/eos_matching_reaudit_20261123/audit_evidence/EVIDENCE_COPY_MANIFEST.json)이 원 경로와 보존 사본을 연결한다. 일부 log는 `.gz`로 무손실 압축했다. Private numerical payload와 private failure trace는 공개 evidence에 복제하지 않는다.

### 이번 턴 백업

백업 위치는 `/home/cosmosapjw/Dropbox/NS_Cowling_EOS_Reaudit_Backups/20261123/`다. 다음 archive는 실제 local write 후 전체 SHA256을 확인했다.

| Archive | Bytes | SHA256 |
|---|---:|---|
| `public_primary_campaign_v1.tar` | 866,816,000 | `346ab15672d5b94c6b306e476ab04e85a632c913724d9f97a0257c9f91e622e6` |
| `public_expanded_campaign_v2.tar` | 981,012,480 | `c6550721b221ab67df47a5b77f84857b381d57c20a1d3fa734db1d62d23a5313` |
| `runtime_evidence_with_private_inputs.tar.gz` | 20,505,783 | `757d9bc0f3f07d2d969188991401fe1498f00a13cbf0eb99cf6a8e9bd71e45b3` |

Private matched candidates는 같은 backup root의 `private_matched_v2/`에 직접 보관했다. [PRIMARY_BACKUP_RECEIPT.json](../data/eos_matching_reaudit_20261123/audit_evidence/PRIMARY_BACKUP_RECEIPT.json)과 [EXPANDED_BACKUP_RECEIPT.json](../data/eos_matching_reaudit_20261123/audit_evidence/EXPANDED_BACKUP_RECEIPT.json)이 실제 bytes·SHA·경로를 기록한다. **Dropbox cloud sync 완료는 독립적으로 확인하지 않았으므로 local backup과 구별한다.**

### 최종 마감표

| 최종 항목 | 상태 / 확정 값 |
|---|---|
| V2 expanded catalog 준비·identity | **완료** — 105 records, 102 available / 3 blocked; 기존 101 records 보존 |
| V2 16,352개 accounting 및 survivor coverage | **완료** — 9,377 = 8,045 endpoint chain + 1,332 finite component |
| V2 post-bridge incomplete / runtime·protocol | **별도 판정 후 각각 0**; 원 145 trials의 기록과 해소 certificate 보존 |
| V2 delivery packaging·전 member hash 검증 | **완료** — public 8,972 candidates / 36 shards / 76,797,682 bytes |
| 전체 CI 결과 | **PASS** — 1,907 passed / 35 skipped / 4 subtests passed; native/importer 전후 identity 동일 |
| 최종 code/data/document 게시 identity | 이 문서를 포함하는 `research/eos-core-crust-all-methods-20261123` 및 최종 PR head; 외부 `SHIPPING_RECEIPT` 참조 |
| 이번 턴 backup | **로컬 기록·전체 hash 확인 완료**; cloud sync는 미확인 |
| 남겨 둔 source·구현 한계 | **미확보 source 3개 / 미지원 multiple-root cell 5개**; 다른 source·interpolant 실패도 지정 representation 범위의 결과 |

이번 감사는 명시된 source catalog·유한 endpoint grid·현재 admission contract에 대해 계산과 별도 실패 판정을 완료했다. Source 가용성 실패, source precision 한계, 미지원 다중 공존 root, finite interpolation 범위를 남겨 둔 채 “모든 물리적 EOS가 증명되었다”거나 “어떤 matching도 존재할 수 없다”고 일반화하지 않는다.
