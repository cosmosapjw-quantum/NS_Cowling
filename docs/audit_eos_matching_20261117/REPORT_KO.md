# EOS generator와 core–crust matching 후속 감사

기준일: 2026-11-17 KST. 기준 저장소: [cosmosapjw-quantum/NS_Cowling](https://github.com/cosmosapjw-quantum/NS_Cowling).

## 판정과 범위

감사 기준은 PR15의 실제 head `80b923b7126492544643a46007b45915e7e3f3aa`다.
`main`의 `546435b5bb2687283fea165357f0a01b3b662532`만 조사하면 최신 생성기·matching 수리를 놓친다.
PR15는 `remediation/2026-11-15-m5a-maxwell-eos-contracts`에서
`research/eos390-fortran-unification-v1`로 향하는 열린 draft이며, 기준 head의
[CI #632](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/44000153338)는 성공했다.
이는 이 감사 패치의 CI 결과가 아니라 **감사 시작점의 결과**다.

현재 판단은 **명시된 native 생성 경로와 제한된 matching 구현은 작동하지만, 18개 원저자 full EOS의 재현 및 모든 crust/core 경계의 물리적 연결은 부분 완료**다.
이번에는 실제 반례가 있는 convention·I/O·조성 보존·출력 근거 전달 결함을 수리했다.
원저자 calibration이나 미표본 밀도에서의 평형을 메타데이터 수정으로 완성했다고 주장하지 않는다.

범위는 공개 생성기, native matching, CompOSE 입출력, 직접 연결되는 frontend/pipeline 계약이다.
기존 과학 산출물은 기준 commit의 보고서·receipt를 확인해 재사용했다.
18모델 전체 최적화, 새로운 EOS390 전수 물리 캠페인, stellar spectrum 캠페인은 수행하지 않았다.
생산 물리 계산은 Fortran, Python은 입력·호출·전달·검증이라는 기존 분리를 유지했다.

## 새 결함과 수정

| ID | 재현한 실패 | 수정과 검증 범위 |
|---|---|---|
| A01 | Smooth matching이 서로 다른 `Il` 또는 실제 `Yq` 입력을 허용. 동일 `Il=0,Yq=.1` 입력도 출력 `Yq=0`으로 바뀌어 retained `Q3`가 .001, 즉 939 MeV 기준 .939 MeV만큼 달라짐 | API와 native에서 동일한 물리 slice를 요구. 출력 `Yq` 보존, `Il=0`의 baryon potential 재인코딩 및 source Euler residual 교정. 네 smooth 방법의 독립 analytic source와 실제 readback으로 검증 |
| A02 | 동일한 물리 `Yq`인데 provider-local `iYq`가 3과 9라는 이유로 sharp Maxwell이 거부 | 실제 `Il,Yq` 비교는 유지하고 정수 index equality만 제거. 출력 axis는 low source의 label을 유지 |
| A03 | 비단위 축 label로 matching한 결과가 저장 후 `cold Yq slice could not be resolved`로 읽히지 않음 | CompOSE의 선언된 min/max index와 배열 offset을 구분. cold reader의 잘못된 상한 검사와 native/internal/Python 사이의 선택 label·axis bounds 누락을 수리. declared range 검사는 유지 |
| A04 | 조성이 있는 crust piece와 thermo-only piece를 조립하면 전체 `composition_*` flag가 false가 되어 유효한 sidecar를 저장할 수 없음 | container-level availability에 `any` 사용. 각 행의 presence/applicability mask는 그대로 보존. 양방향 혼합 및 native write/read 검증 |
| A05 | `states.csv`, `first_law.csv`, `attempts.csv`, `profiles.csv`, `selected.csv`가 generation hash inventory에서 누락 | 존재하는 진단 파일을 수집하고 native RESULT 및 실제 outer DISPATCH의 seal을 각각 검사. 이미 seal된 파일의 삭제·변경은 거부. 과거 미봉인 optional 파일의 호환성 유지 |
| A06 | 실제 공개 생성기는 `DISPATCH.json.evidence.output_sha256`에 SUPPORT hash를 저장하지만 matching은 `generated/RESULT.json`만 읽어 seal을 놓침 | `generated/`와 바깥 `DISPATCH.json`의 구조적 관계를 따라 양쪽 seal을 검사. archived absolute path에 의존하지 않음. main/outer/lower_source_fit에서 변경·삭제 반례 검증 |
| A07 | 공개 CLI 문서가 legacy runner 경로·출력명·SPG endpoint를 현재 native CLI와 혼합 | 실제 native entrypoint, flat `generated/eos.*`, profile 진단, `fixed_density_source` 종료점으로 교정. 과거 실행 이력은 그 당시 범위로 보존 |

A06은 실제 저장 산출물에서도 확인했다.
`data/research_matching_20261117/integration/public_spg_m5/DISPATCH.json`은
`SUPPORT.json` hash `01dfff782864a3e26950c08fc5ae1ecef17a891039eba07ec3657a5adf684421`을
가지며 실제 파일과 일치한다. 같은 출력의 `generated/RESULT.json`에는 hash inventory가 없다.
기준 checkout의 `data/` 아래 65개 DISPATCH 중 12개에서 이 바깥 SUPPORT seal을 확인했다.
이는 새로 고안한 가상 receipt에만 해당하는 문제가 아니다.

`composition_file_present`는 파일 전체에 sidecar 정보가 있다는 의미이며, 모든 행에 모든 abundance가 있다는 의미가 아니다.
조성이 없는 행의 false mask는 그대로 남는다. phase-only readback 행은 `PARTIAL`로 표현될 수 있으며,
이를 free neutron/proton 또는 nucleus abundance가 실제 0이라는 물리 정보로 해석하지 않는다.

## 물리·수학적 기준

### 에너지와 chemical potential convention

가정은 $T=0$, baryon number density $n=n_B>0$, rest energy를 포함한
에너지밀도 ε, 하나의 명시된 thermodynamic path다. 단위는

$$
[n]=\mathrm{fm}^{-3},\quad [\epsilon]=[P]=\mathrm{MeV\,fm}^{-3},\quad
[F]=[\bar\mu]=\mathrm{MeV},\quad s=\frac{c_s^2}{c^2}=\frac{dP}{d\epsilon}.
$$

여기서 $F=\epsilon/n$, $\bar\mu=d\epsilon/dn$이며

$$
P=n^2F'(n),\qquad \bar\mu=F+\frac Pn,\qquad dP=n\,d\bar\mu.
$$

중성 total EOS의 beta-equilibrium path에서는 $\bar\mu=\mu_B$다.
반면 lepton이 없는 fixed-$Y_q$ slice의 CompOSE convention은

$$
\bar\mu=\mu_B+Y_q\mu_q,\qquad
Q_3=\frac{\mu_B}{M}-1,\quad Q_4=\frac{\mu_q}{M},\quad M=m_nc^2.
$$

코드의 `mn` 필드는 MeV 단위의 $M$을 담는다. 그러므로 `Il=0`에서 Euler 조합을 그대로
`Q3`에 넣으면 안 된다. A01은 이 차이를 수정한다. 서로 다른 `Il`의 EOS에 lepton energy를
자동으로 더하는 작업은 수행하지 않는다. 정의와 축 인덱스 형식의 원전은
[CompOSE Quick Guides §2.2](https://arxiv.org/html/2311.04715v1)다.

Header 기준을 $M\to M'$로 바꿀 때

$$
Q_7'=\frac{M}{M'}(1+Q_7)-1
$$

로 재인코딩하면 ε가 보존된다. 반면 실제 $\epsilon\to\epsilon+Cn$은

$$
P\to P,\qquad \bar\mu\to\bar\mu+C,\qquad
s\to s\frac{\bar\mu}{\bar\mu+C}
$$

가 되어 gravitating energy와 sound speed를 바꾼다. 교차점을 만들려고 $C$를 맞추지 않는다.
접합 후 chemical potential 보정과 EOS 변경의 관계는
[Fortin et al., PRC 94, 035804 (2016), §III](https://arxiv.org/html/1604.01944)과 대조했다.
위 변환식은 직접 유도하고 기호 검산했다.

### CSS와 Maxwell의 구별

현재 public smooth 방법은 `cubic_F`, `cubic_P`, `quintic_F`, `quintic_P` 네 가지다.
`match_crust_core(method="css")`는 crust 쪽에 density jump가 있고 core 쪽에 연속인
단일 constant-sound-speed 구간을 만든다. `construct_maxwell_transition()`은 별도 API로,
제공된 두 source의 $P(\mu_B)$ 교차와 local charge neutrality 가정 아래 sharp coexistence를 찾는다.
이 두 construction과 과거의 `first_law_constrained_hermite_pressure_bridge`를 같은 Maxwell로 부르지 않는다.

CSS anchor $(n_1,\epsilon_1,P_1)$와 $0<s\le1$에 대해 직접 적분하면

$$
\bar\mu(n)=\bar\mu_1(n/n_1)^s,
$$
$$
\epsilon(n)=\frac{(\epsilon_1+P_1)(n/n_1)^{1+s}+s\epsilon_1-P_1}{1+s},
\qquad P(n)=s\epsilon(n)+P_1-s\epsilon_1.
$$

두 개의 완전 고정된 continuous anchor를 연결하려면

$$
s=\frac{P_2-P_1}{\epsilon_2-\epsilon_1},\qquad
\frac{\bar\mu_2}{\bar\mu_1}=(n_2/n_1)^s
$$

를 동시에 만족해야 한다. 단순히 ΔP/Δε를 정하는 것만으로 충분하지 않다.
현재 crust-side jump construction에서는 $R=\bar\mu_2/\bar\mu_t>1$에 대해

$$
n_+=n_2R^{-1/s},\qquad
P_2-P_t=\frac{s}{1+s}n_2\bar\mu_t\left(R-R^{-1/s}\right),\qquad
\Delta\epsilon=\bar\mu_t(n_+-n_c).
$$

이는 native `single_css_pressure_span`과 일치한다. 일반 matching formalism은
[Canullán-Pascual et al. (2025), §2](https://arxiv.org/html/2502.02373)와 비교했으며,
현재 코드가 선택한 crust-side jump가 논문의 모든 endpoint construction을 구현한 것은 아니다.

Wolfram Language evaluator에서 CSS 기본 identity 다섯 개와 에너지 이동 identity 세 개의
symbolic residual이 모두 0이었다. $s\to0$에서 $P=P_1$,
$\bar\mu=\bar\mu_1$, $\epsilon=\bar\mu_1 n-P_1$도 확인했다.
실행 가능한 입력과 raw 출력은 `data/audit_eos_matching_20261117/native_conventions/WOLFRAM_THERMODYNAMIC_CHECK.json`이다.
WolframContext는 internal error였고 evaluator는 성공했음을 구분해 기록했다.

연속 구간의 first law는 실제 interpolant에 대해

$$
F(n_b)-F(n_a)=\int_{n_a}^{n_b}\frac{P(n)}{n^2}\,dn
$$

등으로 확인한다. density jump에는 $[P]=[\bar\mu]=0$,
$[\epsilon]=\bar\mu_t[n]$를 적용한다. phase boundary를 가로지르는 centered finite difference로
미분형 first law를 검사하지 않는다. 기존 knot별 적분, 다항식 내부 stationary-point 검사,
감소·상수 (P/n) 허용 수리는 유지했다. 부동소수점 극값 검사를 엄밀한 interval proof로 부르지 않는다.

## 최신 생성기의 실제 완료 범위

기준 commit의 `CURRENT_OUTCOME.json`, `RETAINED_COHORT_MATRIX.json`,
`PHASE_RETENTION_REPAIR_KO.md`와 각 원 receipt를 함께 읽어야 한다.
이하 수치는 **이전 계산 근거를 재사용한 결과**이며 이번 감사에서 물리 계산을 다시 수행한 수치가 아니다.

| 계열 | 확보된 표현·실행 근거 | 아직 다른 주장인 것 |
|---|---|---|
| Core 5종 | APR functional+two-charge Gibbs, BBB 에너지 표본의 parabolic/PCHIP 재구성, GM1 rounded RMF, KDE0v1·SLy4 functional | BBB microscopic BHF kernel 재계산·저자 smooth-fit 미분 동일성. APR/KDE0v1/SLy4의 높은 밀도에서 알려진 비인과성 제거 |
| SPG 5종 | 명시적 candidate surface의 CLD; `public_composition_inferred_20261116` 5×129=645행 재최적화 근거 | 원저자 calibration 및 variational/export 동일성. 후보 variant를 원저자 처방으로 승격할 근거 |
| TF 8종 | 요청한 밀도·geometry·seed·mesh에서 실제 spatial RMF-TF 최소화. FSU2R N256의 유한 slab/tube 공존 후보 후속 수리 | 미표본 밀도에서 연속 지원, 모든 geometry의 전역 최소, mesh 연속 극한, 정의되지 않은 free/cluster abundance 분해 |

현재 public SPG 종료점은 동일 밀도의 clustered/uniform energy crossing이다.
이전 실제 진단의 ΔP 약 0.00445–0.00807 MeV fm⁻³ 및
Δμ 약 0.0606–0.0817 MeV는 common-tangent residual이 아니다.
이를 물리적 coexistence endpoint라고 문서에 쓰던 부분을 A07에서 바로잡았다.

`--include-lower-crust`는 outer와 lower-source-fit을 별도 생성한다.
SPG 시작점과 lower 끝 .002 fm⁻³는 겹치고, 공개 TF smoke의 시작 .003과 lower 끝 .002 사이에는 gap이 있다.
따라서 단순 append는 whole-crust assembly가 아니다.
기존에는 제한된 9개 assembly 후보가 확보됐지만, 나머지 구간을 자동으로 승인하지 않는다.
TF의 비정형 `Iphase=9` 및 정의되지 않은 composition mask를 그대로 보존한다.

원전 충실도, 독립 물리 검증, 수치 지원구간, CompOSE reference 비교는 서로 다른 축이다.
CompOSE 파일 readback 성공은 CompOSE reference 수치 동일성의 증거가 아니다.
반대로 특정 명시적 reconstruction을 생성하려면 항상 저자의 microscopic solver 전체를 재현해야 한다는
새로운 공통 gate도 도입하지 않는다. whole-crust reconstruction의 가정 차이는
[CUTER v2, Davis et al. (2025)](https://arxiv.org/html/2506.08658)와 비교했다.

## Matching 잔여 과제와 우선순위

**가장 중요한 다음 구현은 source boundary별 법칙을 native export와 stellar consumer까지 전달하는 것이다.**
최신 소스는 이미 retained `Iphase`를 보존한다. 따라서 다음 작업을 “phase label 보존을 처음 구현”으로 쓰면 중복이다.
전체 APR의 6→8→7 label 보존은 구현됐지만, 각 경계를 넘는 허용 interpolant와 stellar interface law는 아직 선언되지 않았다.
Public full-APR guard와 `UNSUPPORTED_RETAINED_SOURCE_PHASE_BOUNDARY`는 유지된다.
이 상태는 물리적 matching 불가능의 증명이 아니다.

1. **입력 선택의 명시화.** `scripts/match_core_crust.py` 기본값은 8 crust×5 core×2 방법=80이며,
   crust root는 legacy `download/reconstruction_v6_proxy`다. 최신 18개 generator 출력이나 EOS390 전수 계산이라고 보고하면 안 된다.
   실제 source receipt·variant·지원구간을 담은 manifest를 통해 입력을 지정하고, requested/eligible/attempted/admitted를 따로 집계한다.
2. **경계별 실제 계산.** 연속 Gibbs boundary, sharp first-order jump, 단순 source segment 끝을 구분하고,
   source construction으로 정당화한 법칙만 적용한다. lower/inner gap·overlap에서는 실제 branch 공통접선 또는 정의된 reconstruction을 계산한다.
   부족한 입력이나 no crossing 범위를 명시하며 임의 energy shift나 label 삭제로 해결하지 않는다.
3. **Stellar 소비 범위.** background barotrope와 perturbation closure를 구분한다.
   sharp interface의 slow/rapid conversion은 서로 다른 변위 접합조건을 요구할 수 있다
   ([Tonetto & Lugones 2020](https://arxiv.org/abs/2003.01259)).
   현재 자료만으로 frozen-composition Γ₁, buoyancy 또는 elasticity가 주어진다고 주장하지 않는다.
4. **실행 identity 보충.** public dispatcher는 source·request hash를 기록하지만 generator와 readback executable의 SHA-256을 직접 기록하지 않는다.
   이후 일반 provenance 변경 때 두 binary identity를 따로 추가하면 된다. 이번 감사의 실제 binary hash는 별도 실행 receipt에 보존한다.

기존 실제 SPG5×core5×6방법 캠페인은 150시도였고,
smooth/CSS 100개와 APR 단일 branch 5개를 포함한 105개 native readback 근거가 있다.
수리 후 Maxwell 25개는 13개가 제공된 공통 μ 범위 안의 no crossing,
6개가 공통 μ 지원 부재, 1개가 multiple-crossing 미지원, 5개가 full-APR 경계 지원 제한이었다.
admitted 0이라는 결과를 보편적 no-transition으로 확장하지 않는다.

## 검증과 재현

최종 전체 실행은 **1230 passed, 9 skipped, 0 failed**였고,
legacy generator **142 passed**, theory **6 passed**, native matching controls **4개 프로그램 exit 0**이다.
최종 실행 전후 수정 source SHA-256은 일치했으며 generator와 readback 실행 파일의 SHA-256을 따로 기록했다.
최종 source·산출물 manifest 10,102항목의 SHA-256도 확인했다.
이 결과는 수정 worktree의 실행 결과이며 기준 CI #632와 구분한다. 감사 branch의 원격 CI는 별도로 실행하지 않았다.

최종 실행 수치·명령·source/binary identity·실패 원인은
`data/audit_eos_matching_20261117/RESULT.json` 및 `logs/`에 기록한다.
최초 전체 실행에서 나온 11 failures 중 5개는 A03의 실제 I/O 결함,
6개는 기존 test extra `mpmath` 미설치라는 실행환경 문제였다.
이 두 종류를 합쳐 물리 실패로 분류하지 않는다. 실패 로그도 보존한다.
기존 데이터가 staging되지 않아 생긴 9 skips는 성공으로 계산하지 않는다.

실행 환경의 test extra는 `pyproject.toml`에 이미 선언되어 있으므로,
새 환경에서는 `python -m pip install -e '.[test,plot]'`를 사용한다.
실제 수정·명령·종료 조건을 포함한 Codex 인계는 [HANDOFF_CODEX_KO.md](HANDOFF_CODEX_KO.md)를 따른다.

### 변경 종료 기준

확인한 결함별 red→green과 실제 native roundtrip, 관련 전체 회귀,
기존 native 열역학·구간 반례, source manifest 검증을 완료하면 이 감사 수리는 종료한다.
새 calibration 추정이나 whole-star 연구 캠페인을 이번 코드 수리의 완료 조건으로 추가하지 않는다.
