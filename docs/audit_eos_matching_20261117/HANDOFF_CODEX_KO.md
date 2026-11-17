# Codex 인계: EOS generator / core–crust matching

이 문서는 `80b923b7126492544643a46007b45915e7e3f3aa` 이후의 감사 수리를 이어가기 위한 prompt다.
최신 원격 SHA는 실행 시 다시 확인하되, 이 문서와 감사 패치가 포함된 commit을 잃지 않는다.
이미 확인·수리한 항목은 `REPORT_KO.md`와 `data/audit_eos_matching_20261117/RESULT.json`에 있다.

## 그대로 전달할 prompt

~~~text
NS_Cowling의 EOS generator와 core–crust matching 후속 구현을 수행하라.
저장소: https://github.com/cosmosapjw-quantum/NS_Cowling
감사 기준 PR15 head: 80b923b7126492544643a46007b45915e7e3f3aa
감사 수정 branch: audit/eos-matching-20261117

목표는 감사 보고서를 다시 만드는 것이 아니라, 최신 generator 출력이 물리적으로
정당화된 source segment와 interface를 거쳐 matching 및 소비자로 전달되도록
아래의 남은 구현을 실제로 진행하는 것이다. 이미 완료된 수리를 보존하라.

1. 시작점과 SSOT

git fetch 후 PR15, 감사 branch, 관련 CI의 실제 head를 확인하고 정확한 SHA를 기록하라.
기존 dirty checkout은 보존하고 분리된 worktree에서 작업하라.
감사 branch가 PR15보다 앞서면 그 변경을 포함한 상태에서 시작하라.
원격이 더 진행됐다면 이 prompt의 SHA를 최신이라고 가정하지 말고 변경된 의존성만 비교하라.

다음을 순서대로 읽어라.
- docs/audit_eos_matching_20261117/REPORT_KO.md
- data/audit_eos_matching_20261117/RESULT.json
- docs/research_finish_20261116/PUBLIC_GENERATION_CLI.md
- data/research_matching_20261117/CURRENT_OUTCOME.json
- data/research_matching_20261117/generator/RETAINED_COHORT_MATRIX.json
- docs/research_matching_20261117/PHASE_RETENTION_REPAIR_KO.md
- data/research_matching_20261117/phase_retention_repair/assembly_003/RECEIPT.json

기존 18모델 물리 계산은 source·request·artifact 의존성이 같은 범위에서 재사용하라.
실제 물리 입력이나 solver를 바꾼 항목만 필요한 밀도·모델 범위로 재계산하라.

2. 보존할 수정과 물리 계약

감사 패치에는 다음이 실제 반례와 함께 포함돼 있다.
- matching 입력의 동일 Il 및 실제 Yq 계약, Il=0에서 mu_Euler와 mu_B 구분
- output Yq/Q3 보존 및 물리적으로 같은 slice의 서로 다른 iYq 허용
- CompOSE 단일 축의 선언된 min/max label과 배열 offset을 구분하는 reader 및 native/internal/Python metadata 전달 수리
- 혼합 composition coverage의 container availability와 행별 mask 보존
- generator의 states/first_law/attempts/profiles/selected 파일 inventory 및 실제 outer DISPATCH seal 재검증
- 실제 public DISPATCH.evidence.output_sha256의 SUPPORT seal을 matching에 전달
- 현재 native public CLI와 SPG fixed_density_source endpoint에 맞춘 문서

생산 EOS·미분·최소화·matching은 Fortran에 둔다. Python은 frontend/transport와
독립 검증 역할을 유지한다. legacy Python oracle을 생산 물리 backend로 바꾸지 않는다.
T=0, epsilon은 rest energy 포함, n_B는 baryon number density다.
F=epsilon/n, mu_Euler=F+P/n, P=n^2 dF/dn, dP=n dmu_Euler를 사용한다.
Il=0 fixed-Yq에서는 mu_Euler=mu_B+Yq*mu_q이며 Q3는 mu_B를 인코딩한다.
코드의 mn은 MeV 단위 rest-energy reference다. Header reencoding은 물리 epsilon을 보존한다.
epsilon -> epsilon+C*n으로 Maxwell crossing을 인위적으로 만들지 않는다.

source fidelity / independent physics / numerical support / CompOSE comparison을 구분하라.
18개 native route의 실행과 readback은 18개 원저자 full EOS 동일성과 다르다.
SPG의 명시적 candidate variant와 공개 조성 역추정 variant를 원저자 calibration으로 승격하지 않는다.
현재 public SPG endpoint는 fixed-density energy crossing이다. 비영 DeltaP/DeltaMu를
common-tangent residual이라고 이름만 바꾸지 않는다.
TF의 미표본 밀도, finite mesh/seed/geometry 지원 및 비정형 phase9를 보존하라.
누락된 abundance mask를 0 abundance로 채우지 않는다.

3. 첫 실행 목표: 최신 generator 입력을 명시적으로 연결

scripts/match_core_crust.py의 기본 census는 8 crust*5 core*2 methods=80이며
download/reconstruction_v6_proxy를 읽는다. 이것을 최신 18-route 또는 EOS390 결과라고 보고하지 마라.

실제 generator 출력의 model, family, variant, source/request identity, CompOSE 경로,
SUPPORT, density/phase segment, composition policy를 담은 명시적 입력 manifest를 작성하고
pipeline이 이 입력을 받도록 최소 확장하라. 기존 legacy census는 호환성을 유지하라.

현행 13 crust*5 core*(4 smooth+CSS+sharp Maxwell)의 coverage를 다룬다면 390개 요청 cell을
명시하되 requested/eligible/attempted/admitted를 분리하라. 지원구간이 부족한 cell은 정확히
NOT_EVALUATED 또는 typed HOLD로 남기고, 그 수를 실제 실행 성공에 포함하지 마라.
이미 확보된 105개 실제 readback 및 제한된 9개 crust assembly는 dependency identity가
같을 때 재사용하라. 새로운 matrix를 만들었다는 이유만으로 TF 전체 최적화를 반복하지 않는다.

4. 다음 구현 목표: source boundary 법칙과 소비 경로

retained Iphase 보존은 이미 native에 구현됐다. 전체 APR의 6->8->7 label이 보존되는
상태에서 공통 interpolator/stellar intake가 선언되지 않은 경계를 거부하고 있다.
다음 작업을 단순 label 복구로 재정의하지 마라.

source construction과 실제 표본에서 각 경계가 continuous Gibbs endpoint인지,
sharp first-order jump인지, 지식이 없는 segment 끝인지 확인하라.
명시적 boundary metadata를 정의하고 허용된 경계에서만 piecewise interpolation을 하라.
각 연속 구간은 그 실제 interpolant의 first law, stability, causality를 검증하고,
sharp jump는 [P]=[mu]=0 및 [epsilon]=mu_t[n]를 확인하라.
label 값만으로 boundary law를 추정하거나 경계를 가로지르는 FD stencil을 사용하지 않는다.

APR의 연속 Gibbs 경계가 source와 계산으로 확인되면 그 경계 유형에 한해 native export,
CompOSE transport, CSS/TOV consumer를 끝까지 연결하라. 미승인 경계의 기존 typed HOLD는 유지하라.
lower/inner gap 또는 overlap은 실제 branch 접합 계산으로 다뤄라. 근거 없는 append,
energy offset 또는 source label 병합으로 통과시키지 마라.
기존의 명시적 thermo-only 변환은 원본·mask provenance·변환 receipt를 보존하여 사용하고,
조성을 보존한 물리 transition이라고 주장하지 않는다.

sharp Maxwell은 실제 두 source P(mu)의 교차이며, 현재 여러 crossing은 미지원이다.
CSS는 crust 측 density jump와 core 측 연속 연결을 갖는 별도 construction이다.
네 smooth 방법, CSS, true Maxwell, 역사적 FLCHPB의 의미를 혼합하지 않는다.
제공된 common-mu 구간 내 no crossing, common support 부재, 다중 root 미지원,
입력 convention 오류, 모델의 비인과성, 수치 실패, 실행환경 실패를 구분하라.

Cowling까지 진행할 때 background EOS만으로 Gamma1 또는 interface conversion law를
정했다고 주장하지 마라. barotropic f/p 범위를 명시하고, slow/rapid conversion 또는
frozen-composition/elasticity가 필요한 경우 필요한 추가 closure를 구체적으로 남겨라.
이 요청을 g-mode/elasticity 전체 프로젝트로 확장하지 않는다.

5. 검증·종료·인계

먼저 가장 작은 실제 입력 하나를 end-to-end로 완료한 뒤 같은 계약의 지원 대상에 확장하라.
실패를 재현하는 독립 analytic case 또는 실제 입력 case를 먼저 보존하고 최소 수정하라.
수정하지 않은 기능의 reviewer tree나 whole-campaign을 반복하지 않는다.
정의·단위·부호·path·approximation을 명시하고, 구현/수치/모델/환경 실패를 구분하라.

아래 명령과 변경에 직접 관련된 native regression을 실행하라.
비단위 iT/iYq, 동일 Yq/상이한 index, Il0 nonzero Yq, 혼합 composition mask,
실제 DISPATCH seal, 존재/삭제/변경된 diagnostic의 red->green을 유지하라.
합계 PASS만 보고하지 말고 skips와 실제 native binary/source identity를 기록하라.
public generator의 다음 provenance 변경에는 generation executable와 readback executable의
SHA-256을 별도로 남겨라. source hash와 실행 파일 identity를 동일시하지 않는다.

종료 산출물은 실제 코드, 필요한 테스트, 입력 manifest, 실행/실패 receipt,
남은 물리·수치·입력 범위 및 한 개의 업데이트된 handoff다.
해결하지 못한 항목에는 부족한 입력 또는 다음 최소 계산을 명시하라.
로컬 전체 회귀와 원격 CI를 구분하고, 실제 확인하지 않은 통과나 배포를 주장하지 않는다.
기존 역사적 raw evidence는 덮어쓰지 않는다. 새 산출물의 hash를 manifest에 반영하라.
~~~

## 재현 명령

프로젝트 root, test extra와 Fortran/LAPACK 의존성이 준비된 환경에서 실행한다.
`mpmath>=1.3`은 root `pyproject.toml`의 test extra에 이미 포함돼 있다.

~~~bash
python -m pip install -e '.[test,plot]'
bash fortran_eos/build.sh build/compose_eos_fortran
export PYTHONPATH="$PWD/python_frontend:$PWD/scripts"
export COMPOSE_EOS_FORTRAN_ROOT="$PWD"
export COMPOSE_EOS_FORTRAN_EXE="$PWD/build/compose_eos_fortran/compose_eos_fortran"
export COMPOSE_EOS_NATIVE_BUILD_DIR="$PWD/build/compose_eos_fortran"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python -m pytest -q -ra python_frontend/tests pipeline/tests
python fortran_eos/tests/run_matching_native_adversarial.py \
  --build-dir build/matching_controls --log-dir build/matching_controls/logs
python -m pytest -q fortran_eos/tests/matching_theory_soft_barotrope.py \
  fortran_eos/tests/matching_theory_maxwell_topology.py
PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg \
  python -m pytest -q eos_generator/tests
sha256sum -c SHA256SUMS_v4.txt
~~~

검증 범위를 추가할 때는 변경한 물리 의존성이나 남아 있는 구체적 위험을 이유로 기록한다.
