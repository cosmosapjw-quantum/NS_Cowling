# NS_Cowling 개선 C — TF continuation과 지원 가능한 matching 확대

2026-11-17 KST. 기준 B commit: `894cbf1c8315f07ee926c29b83766ffdab2ac7b7`.
작업 branch: `research/tf-pair-expansion-20261117`.
B는 사용자 승인에 따라 같은 SHA로 GitHub에 push했고 원격 readback을 확인했다. C의 게시에는 별도 승인이 필요하다.

## 실행 결과

기존 native DD2 유한 해에서 두 저밀도 노드의 TF continuation을 계산하고, N512 droplet 정지해와 mesh/cell-size/열역학 검산을 확보했다. 출처와 단일-phase 지원이 있는 새로운 core/method 조합80개를 native Fortran으로 생성했다. 기존 matching 근거와 합쳐390개 요청 ID의 중복 없는 현재 상태를 기록했다.

기존 A/B primary EOS와 source physics는 보존했다. C의 추가 matching 표는 실제 상단0.5 fm^-3인 inner-crust+core 표이며, B primary의 표면–1.0 fm^-3 지원과 구분한다. 이번 TF 결과는 두 노드의 branch-specific finite-mesh evidence다. 전역 geometry 최소와 연속 TF EOS는 미확정으로 유지했다.

### C 작업별 결과

| 작업 | 산출물과 상태 |
|---|---|
| C1 입력/범위 선택 | 기존 DD2 n=.003 native profile6개와 실제 DD2 seed 좌표, 미계산 지원 조합80개 선택 |
| C2 TF continuation | n=.0025 및 .0020000001 fm^-3에서 native N512 droplet 정지해 확보 |
| C3 수치/물리 검산 | geometry 비교, N64/128/256/512, radius 시작값 ±10%, 독립64-profile 평가, 재최적화 압력 미분 |
| C4 나머지 core/method | 새80개 kernel/output; smooth70 construction, CSS10 independent admission |
| C5 보존/인계 | 390-ID map, 실패·source·실행 기록, 재현 명령, 다음 D 계약 |

A–D 중 C는 이번에 선택한 유한 범위의 실행 단위다. 390개 요청 전체의 해결률과 동일한 진행률로 사용하지 않는다.

## TF: 정의와 실제 계산

T=0, neutrino-free cold nucleonic npe DD2 RMF를 사용했다. Density-dependent vertices와 rearrangement는 local n_n+n_p에서 평가하고 electrons는 screened local solution이다. 공간 미지수는 n_n(r), n_p(r), cell radius R이며 meson/electron fields는 각 functional 평가에서 푼다. 중심 regularity, cell edge의 Neumann 조건, 전역 전하 중성, Coulomb mean-zero gauge를 유지했다.

에너지는 정지 에너지를 포함한다. 단위는 n_B:fm^-3, R:fm, epsilon와 P:MeV fm^-3, chemical potentials:MeV다. ħc=197.3269804 MeV fm과 DD2의 기존 neutron/proton rest energies939.56536/938.27203 MeV를 유지했다. 새로운 에너지 이동이나 parameter fitting은 없다.

공개 DD2 row217은 n_B=.0020000001 fm^-3, Ye=.06988, Ymu=0이라는 초기 평균조성을 제공한다. Native 최소화에서 Ye는 자유롭게 변한다. 공개 phase index7은 inner crust이고 geometry는 별도의 thermo column이다. 공개 자료는 local fields나 Wigner–Seitz radius를 제공하지 않는다. 이 의미는 [공식 DD2 제출 문서](https://compose.obspm.fr/download/1D/Crust/DD2/eos.pdf)와 실제 seed catalog에 결부했다.

연속 밀도 경로의 시작은 기존 native N16, n_B=.003의 실제 보존 profile이다. .0025와 .0020000001로 이동하고 반환된 유한 해를 재시작·refine했다. 고해상도 endpoint에서 .0025로 되돌아간 실제 continuation 경로도 실행 기록에 남겼다. Public pressure/energy를 native 결과로 대입하지 않았다.

### N512 droplet 결과

| n_B [fm^-3] | epsilon [MeV fm^-3] | P [MeV fm^-3] | R [fm] | Ye |
|---:|---:|---:|---:|---:|
|0.0020000001|1.88066653850851|0.00195576004780|42.58456400|0.07005257970|
|0.0025|2.35138626364943|0.00251873421999|41.25776996|0.06164334093|

Endpoint droplet의 mesh 비교:

| N | epsilon [MeV fm^-3] | P [MeV fm^-3] | R [fm] |
|---:|---:|---:|---:|
|64|1.88065861941301|0.00197098930149|42.75740466|
|128|1.88066508195265|0.00195517707387|42.53582878|
|256|1.88066624226361|0.00195579800230|42.57803009|
|512|1.88066653850851|0.00195576004780|42.58456400|

N128→256→512에서 마지막 두 에너지 차이는 .00058015545, .00014812244 MeV/baryon이며 관측 차수는1.96965다. Endpoint 압력의 N256→512 상대 변화는−1.94065e-5(−.001941%)다. .0025 노드에서는 마지막 에너지 차이 .00011589080 MeV/baryon, 상대 압력 변화−1.81018e-4(−.01810%)다. 이 값은 실제 discretization sensitivity이며 엄밀한 continuum error enclosure가 아니다. 압력 차이는 부호가 바뀌므로 에너지와 같은 수렴차수를 가정하지 않았다.

같은 N256 endpoint의 rod−droplet 에너지 차이는+.1414414910 MeV/baryon이다. 계산된 두 branch에서는 droplet이 낮으며 그 차이는 관측 mesh 민감도보다 충분히 크다. 두 branch의 R 초기값을0.9배와1.1배로 바꾸어 실제 재최적화한 네 결과는 같은 basin으로 돌아왔고 |delta epsilon|<=1.18e-13 MeV fm^-3였다. 초기 family명과 최종 geometry는 구분했으며, uniform으로 수렴한 시작값을 독립적인 tube/bubble 최소해로 세지 않았다.

### Stationarity, 실패와 독립 검산

Native 최소화는64회, stationary45회, optimizer mode8로 보존한 원 실행19회다. Process returncode0과 수치 admission은 별개다. 원 mode8 기록을 보존하고 해당 반환 profile에서 baryon 평균을 복원한 새 warm solve를 실제 실행했다. 기존 ftol1e-14, maxiter1800, chemical KKT<.02 MeV, radius residual<2e-4 MeV 및 field 조건을 바꾸지 않았다.

초기 coarse N16에서 보이던 과도하게 낮은 에너지는 mesh 증가에 따라 변하는 underresolution이었다. 다른 해상도의 에너지를 geometry 경쟁값으로 비교하지 않았다. 실행기에는 두 번째 mesh refinement의 출력 label 충돌이 있었고, overwrite guard가 이를 막았다. 해당 구현 결함은 재현 테스트 RED→GREEN으로 수정했으며 원 로그와 과학적 실패 기록을 구분해 보존했다.

64개 모든 반환 profile을 독립 Python functional로 평가했다.

| 독립 검사 | 최대 관측값 |
|---|---:|
| native–oracle epsilon 차이 |7.33e-15 MeV fm^-3|
| local chemical potential 차이 |1.69e-8 MeV|
| log-radius derivative 차이 |7.80e-14 MeV fm^-3|
| electron Poisson residual |1.30e-12 MeV fm^-2|
| baryon constraint residual |2.17e-15 fm^-3|
| charge constraint residual |5.42e-20 fm^-3|

실패 profile도 이 비교에 포함된다. 같은 방정식의 독립 구현 일치와 optimizer 성공 여부, 원저자 source fidelity는 각각 다른 판정이다.

## 열역학과 source 차이

Eliminated Coulomb energy를 E_C, positive meson-gradient magnitudes를 G_sigma,G_omega,G_rho라 쓰면 기존 diffuse TF functional은

```text
d epsilon/d log R = 2(-G_sigma + G_omega + G_rho + E_C)
P = n_B mu_B - epsilon
P_path = n_B^2 d(epsilon/n_B)/dn_B
```

를 준다. 첫 식의 interior stationary 조건은 G_sigma−G_omega−G_rho=E_C다. Wolfram으로 radius scaling, virial identity, fixed-density objective conditioning, nonstationary correction을 포함한 pressure envelope의 잔차4개가 정확히0임을 확인했다. 기호 항등식 검산과 native 수치 검산은 분리했다.

N256 endpoint에서 h=1e-6 fm^-3의 n0+h,n0+2h,n0+4h를 같은 branch에서 실제 재최적화했다. e=epsilon/n_B에 대한 전진3점 미분을 두 step으로 평가했다.

| 압력 [MeV fm^-3] | 값 |
|---|---:|
| Euler |.001955797986885388|
| 재최적화 전진3점, h |.001955797681641538|
| 재최적화 전진3점, 2h |.001955796529652695|

h 결과의 Euler 대비 상대 차이는−1.5607e-7, 두 step 결과의 차이는1.152e-9 MeV fm^-3다. Frozen-profile derivative가 아니다. 이 국소 검사는 노드 사이 전체 밀도 보간의 정확도를 대신하지 않는다.

같은 public seed row와의 endpoint 비교도 별도로 기록했다. N512 native 압력은 public 복원값보다 .216325% 낮다. 그러나 public 복원행 자체의 (epsilon+P−n mu_B)/n은 .08244845 MeV/baryon이다. Native−public 에너지 차이−.08140774 MeV/baryon 전체를 모델 차이로 귀속하지 않았다. Public source convention·tabulation·implementation 차이의 원인은 unresolved로 남기고, reference fitting이나 에너지 이동을 수행하지 않았다.

## Matching 확대와 실제 390-ID 상태

SPG M1–M5의 기존 sealed inferred-composition variant를 사용했다. Source와 단일 phase가 확인된 실제 core table을 native kernel에 전달했다.

| method | 새 kernel/output | 독립 method admission |
|---|---:|---:|
|cubic_F|20|0|
|cubic_P|20|0|
|quintic_F|20|0|
|quintic_P|10|0|
|CSS|10|10|
|합계|80|10|

첫3개 방법은 GM1/SLy4/KDE0v1/BBB 네 core에, 추가 quintic_P/CSS는 KDE0v1/BBB에 적용했다. Smooth70은 native construction+serialized consistency, CSS10은 별도 interface admission이다. Core 입력 .08–.5 fm^-3와 출력 상단 .5 fm^-3를 보존했다. BBB에는 자체 source verifier를 적용했고 B의 세-core 검증을 전이하지 않았다. 모든 matching 결과는 명시적 thermo-only이며 조성에 의존하는 진동 closure를 제공하지 않는다.

| 390개 요청의 현재 상태 | unique model/method ID |
|---|---:|
|전체 requested universe|390|
|해당 finite component/method의 kernel 근거가 있는 ID|165|
|constructed/readback 표가 있는 ID|105|
|독립 CSS/Maxwell admission|25|
|construction-only smooth|80|
|이번 신규 kernel ID|80|
|재사용한 기존 evidence ID|85|
|APR retained Gibbs 지원을 기다리는 non-Maxwell ID|25|
|TF full-source smooth/CSS 지원이 미확정인 ID|200|

105 output ID는 smooth80+CSS20+Maxwell5다. 기존 native Maxwell65개 ID에는427 component invocations가 있고, 5개 output ID에서 실제 표6개가 나온다. Raw table 수, matrix ID 수, component 호출 수를 구분했다. Public 두 family130개와 초기 campaign의 중복 Maxwell ID는 native390 집계에 더하지 않았다. 나머지225는 이유가 명시된 deferred 상태이며 영구적인 물리 no-go 판정이 아니다.

### 큰 유한차분 잔차의 제한된 sampling 검산

최대 C8 잔차를 보인 SPG_M2–GM1 quintic_F 한 cell에서 n_bridge를40→80→160으로 늘렸다. Native 원 연속식/endpoint는 같고 sampling만 바뀐다.

| n_bridge | C8 상대 진단 | C9 상대 진단 | analytic first-law residual |
|---:|---:|---:|---:|
|40|.0393595273|.0003169324|0|
|80|.0116674552|.0000898888|0|
|160|.0033483503|.0000245299|0|

공통 source/endpoint220행은 exact equality이고, bridge 공통27점에서 P 차이는5.56e-16, epsilon 차이는2.85e-14 MeV fm^-3 이하이다. 변하지 않는 C7 잔차 .0058714458은 첫 retained SPG source node의 one-sided stencil에 있음을 재계산했다. 이 추가2회는 새 coverage나 admission으로 세지 않았다. 한 cell의 sampling 검산을 모든 smooth 방법의 독립 continuum admission으로 확장하지 않았다.

## 구현·회귀·독립 검토

Native production physics와 기존 A/B 파일은 변경하지 않았다. 새 코드의 역할은 TF profile continuation adapter/orchestration, supported-pair campaign, 출처/집계/재현 및 그 동작 회귀다.

최종 전체 frontend/pipeline 회귀: **1,378 passed, 9 skipped, 0 failed**, 실제 exit code0, pytest110.89초. 검사 전후 최종 implementation/test source hash가 일치했다. 9개 skip은 준비되지 않은 기존 raw CompOSE fixtures 때문이며 정확한 목록을 FRONTEND_PIPELINE_REPAIRED.log에 보존했다. Matching-focused15개, TF runner2개 회귀도 통과했고 이들을 전체 검사와 중복 합산하지 않았다.

최초 전체 회귀는 새 C 작업트리에 일반 EOS native executable이 없어서 실패했다. Matching 과학 계산은 검증된 B executable, TF는 별도 built driver를 사용했다. Repository CI의 native build 명령으로 C 실행 파일을 만든 뒤 재검사했다. 최초 실패 로그와 환경 진단은 보존했고, 이 실패를 물리적 반례나 solver 수리로 분류하지 않았다.

독립 검토는 원 profile의 가중 제약, mesh/geometry 차이, 압력 미분,390-ID 재집계,506 source/output seals 및 새 sampling 진단을 확인했다. 한 번의 검토에서 현재 범위의 blocking defect는 발견되지 않았다. 검토 종료 때 상위 전체 회귀가 진행 중이었으므로 그 결과는 별도 최종 실행 기록에 따른다.

## 재현과 파일

필요 환경은 기존 project의 NumPy/SciPy/pytest/mpmath, GNU Fortran과 BLAS/LAPACK이다. 실행 환경에서는 기존 mpmath1.3.0을 scoped PYTHONPATH로 사용했고 host 설정이나 package를 바꾸지 않았다.

```bash
bash fortran_eos/build.sh build/compose_eos_fortran
python scripts/expand_supported_pairs_20261117.py \
  --out data/tf_pair_expansion_20261117/matching/reproduction_NEW \
  --native-exe build/compose_eos_fortran/compose_eos_fortran
bash data/tf_pair_expansion_20261117/tf_continuation/REPRODUCE.sh \
  /tmp/tf_C_replay_NEW
PYTHONPATH=python_frontend \
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
python -m pytest -q -rs python_frontend/tests pipeline/tests
```

실행은 repo root에서 하고 새 output 경로를 지정한다. TF REPRODUCE.sh는 이번에 사용한 Python 환경 경로를 보존하므로 다른 환경에서는 해당 PYTHON/PYTHONPATH를 자신의 동일 의존성 환경으로 지정해야 한다. 모든 실제 native 입력과 출력은 별도 보존되어 있다. 전체 재현 스크립트의 일괄 replay를 추가로 다시 수행했다고 주장하지 않는다.

핵심 파일은 아래와 같다.

- `tf_continuation/SUMMARY.json, STATES.json, THERMODYNAMIC_CHECK.json, INDEPENDENT_RESIDUALS.json`: 실제 TF 결과.
- `matching/campaign_v1/RECEIPT.json, MATCHING_390_MAP.json`, `matching/CATALOG.csv`: matching 결과·390-ID 상태.
- `matching/sampling_refinement_v1/REFINEMENT.json`: sampling 진단.
- `B_PUBLICATION_RECEIPT.json, SYMBOLIC_CHECK.json, REFERENCE_ENDPOINT_COMPARISON.json`: 승인된 B 게시·기호 검산·source 비교.
- `CURRENT_OUTCOME.json, NEXT_STEP_D.json, STAGE_MANIFEST.json`: C 최종 상태·D 인계·payload identity.

상대 data 경로의 공통 root는 `data/tf_pair_expansion_20261117/`다. 별도 상세 보고서는 같은 docs directory의 TF_RESULTS_KO, MATCHING_RESULTS_KO, TF_SOURCE_NOTES_KO, INDEPENDENT_REVIEW_KO에 있다.

## 다음 단계 D

다음은 APR의 retained6→8→7 Gibbs 경계와 phase-aware consumer를 연결하고, 충분한 실제 지원을 가진 A/B primary에서 TOV 및 명시한 Cowling closure의 관측량을 계산하는 단계다. C의 TF sparse-node 결과를 연속 stellar EOS로 자동 사용하지 않는다. Equilibrium thermo-only closure와 fully frozen-species buoyancy는 필요한 입력과 해석을 구분한다.

D에서 중심밀도와 EOS 상단의 충분성을 실제 sequence로 확인하고, 필요할 때만 core를 더 확장한다. 수치 수렴과 경계/closure 민감도를 확인한 관측량을 EOS identity에 결부한다. C는 별도 게시 승인을 받은 뒤 게시하며 D는 다음 선택된 실행으로 남아 있다.
