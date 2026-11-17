# NS_Cowling local handoff — 2026-11-17

후속 phase 보존 수리는 [PHASE_RETENTION_REPAIR_KO.md](PHASE_RETENTION_REPAIR_KO.md)와
`data/research_matching_20261117/phase_retention_repair/assembly_003/RECEIPT.json`을 따른다.
native retained label/CompOSE 전달은 구현했고, 미승인 source 경계의 stellar 적용은 typed HOLD다.
public 전체 APR SUPPORT guard는 유지한다. 아래 내용은 앞선 캠페인의 인계 기록이다.

중앙 저장소는 cosmosapjw-quantum/NS_Cowling, 작업 branch는
remediation/2026-11-15-m5a-maxwell-eos-contracts, PR15다.
이번 작업의 시작 commit은 478c7b5de1514ac512f91156f36551478f36f688이다.
이 문서를 포함한 게시 commit의 정확한 SHA와 CI는 최종 인계 메시지에서 확인한다.
다른 dirty checkout을 reset하거나 덮어쓰지 않는다.

## 그대로 사용할 인계 prompt

~~~text
NS_Cowling의 remediation/2026-11-15-m5a-maxwell-eos-contracts를 이어라.
먼저 git fetch와 원격 HEAD/PR15/CI를 확인하고, dirty worktree를 보존한 별도 checkout에서 진행하라.
docs/research_matching_20261117/{REPORT_KO.md,HANDOFF_KO.md},
data/research_matching_20261117/CURRENT_OUTCOME.json,
data/research_matching_20261117/generator/RETAINED_COHORT_MATRIX.json부터 읽어라.
이전 whole-campaign을 무조건 재실행하지 말고 물리 의존성과 산출물 hash를 대조하여 재사용하라.

생산 물리·미분·최소화·matching은 Fortran, Python은 frontend/transport다.
독립 연구 oracle과 legacy reference를 생산 backend로 바꾸지 않는다.
원전 충실도, 독립 물리 검증, 수치 지원구간, CompOSE 비교를 분리한다.
CompOSE 숫자 완전일치, 전체 microscopic BHF/ETFSI 재계산, global minimum의 일반 증명을
모든 명시적 representation 생성의 공통 gate로 추가하지 않는다.

18개 native 생성 경로는 실제 실행·readback 근거가 있지만 18개 원저자 full EOS의
동일성과 모든 연속 crust 구간은 완료된 것이 아니다.
SPG public_composition_inferred_20261116은 공개 조성에서 역추정한 조건부 surface input이며,
기본 AME2408_candidate를 대체하지 않는다. 원저자 calibration/표면 미분 동일성은 미해결이다.
고정 밀도 에너지 교차와 common tangent를 구분하라.
TF 실패·비정형 profile·미정의 composition mask를 지우지 않는다.

Matching의 Euler endpoint, pressure-spline 적분, 다항식 내부 인과성,
감소/상수 P/n, 냉각·receipt·지원구간 계약의 수정과 반례를 보존하라.
Header 질량 기준의 재인코딩은 실제 P, epsilon, mu를 보존하며
epsilon -> epsilon+n*C라는 물리적 에너지 이동과 다르다.
Maxwell 교차를 만들기 위해 상수를 맞추지 않는다.
construct_maxwell_transition의 drop_composition=True는 명시적 thermo-only 복사본을 만들며
원본·물리 열·phase를 유지한다. 원래 abundance를 복원했다는 주장은 하지 않는다.

다음 구현 우선순위는 retained source의 여러 phase 경계를 보존하는 native assembly와
그에 대응하는 CompOSE/CSS stellar intake다. 현재 전체 APR의 Gibbs 경계를 지우는 입력은
거부하며, 원자료의 연속한 단일 branch 선택은 실제 검증됐다.
이 미구현을 APR matching의 물리적 불가능으로 해석하지 않는다.
기존 lower/inner gap·overlap과 source-reference 차이를 임의 보간으로 숨기지 않는다.
각 남은 항목의 최소 계산/구현을 수행하고 실패 원인을 모델·수치·구현·입력범위로 분리하라.
~~~

## 실행 명령

새 output directory를 사용한다. 프로젝트 root에서:

~~~bash
bash fortran_eos/build.sh build/compose_eos_fortran
export PYTHONPATH="$PWD/python_frontend:$PWD/scripts"
export COMPOSE_EOS_FORTRAN_ROOT="$PWD"
export COMPOSE_EOS_FORTRAN_EXE="$PWD/build/compose_eos_fortran/compose_eos_fortran"
export COMPOSE_EOS_NATIVE_BUILD_DIR="$PWD/build/compose_eos_fortran"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python3 -m compose_eos.cli generators
python3 -m compose_eos.cli generate SPG_M5 build/spg_m5_inferred --points 129 --variant public_composition_inferred_20261116
python3 -m compose_eos.cli generate GM_GM1 build/gm1_native --points 129 --n-max 0.5

python3 -m pytest -q python_frontend/tests pipeline/tests
PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg python3 -m pytest -q eos_generator/tests
python3 fortran_eos/tests/run_matching_native_adversarial.py --build-dir build/matching_controls --log-dir build/matching_controls/logs
python3 -m pytest -q fortran_eos/tests/matching_theory_soft_barotrope.py fortran_eos/tests/matching_theory_maxwell_topology.py
~~~

0.5 fm^-3는 이번 matching 통합 캠페인의 명시적 core 상단이다.
원 모델의 고밀도 비인과성을 수리하거나 전체 stellar EOS 범위를 인증한 값이 아니다.
SPG 645행 및 독립 변분 검산은 spg/REPORT.md,
TF 실패·복구와 유한 chord는 tf/ 보고서,
matching 도함수·원전 근거는 data/.../theory/THEORY_AUDIT.md를 따른다.

실제 150조합 재현:

~~~bash
python3 scripts/research_matching_integration_20261117.py --build build/compose_eos_fortran --spg-root data/research_matching_20261117/spg/native_129_repaired --out build/real_matching_pairs
~~~

이 명령의 완료는 모든 조합의 물리적 허용을 뜻하지 않는다.
각 rejection과 source support, 원래의 입력 identity를 함께 보존한다.
정규화 수정 후 Maxwell 후속 실행은 scripts/research_maxwell_integration_20261117.py에 있다.
APR 단일 branch 재현과 전체 phase 보존 미구현 범위는
pipeline_matching/APR_SUPPORT_SCOPE.md와 해당 실행 스크립트에 있다.
