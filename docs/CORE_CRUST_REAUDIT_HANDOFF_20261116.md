# Coding-agent handoff — NS_Cowling, 2026-11-16

```text
NS_Cowling의 적대적 재감사 후속 작업을 이어가라.
Repo: https://github.com/cosmosapjw-quantum/NS_Cowling
Branch: remediation/2026-11-15-m5a-maxwell-eos-contracts
PR: #15, DRAFT / OPEN / UNMERGED 유지.
v6 별 계산 코드 기준: 42950053088fa80752826f4aae56d0a0b013c0e4
추가된 equation-only interface 기준: 594552d7a597dd150be14b9a860a8f5eaa571d7f
(본 문서와 receipt가 추가된 후속 커밋도 읽고, 시작할 때 실제 원격 HEAD를 확인하라.)

먼저 읽을 SSOT:
- docs/CORE_CRUST_FOLLOWUP_INTERFACE_20261116.md
- data/receipts/core_crust_followup_20261116/verification/RECEIPT.json
- docs/CORE_CRUST_ADVERSARIAL_REAUDIT_20261116.md
- data/receipts/core_crust_audit_20261116/AUDIT_RECEIPT.json
- data/receipts/core_crust_audit_20261116/CENSUS_v6.json
- data/receipts/core_crust_audit_20261116/CSS_SURVIVOR_STARS_v6.json
- download/reconstruction_v6_proxy/regeneration_index.json
- data/ns_cowling_cloud_agent_dag.v1.json
- data/eos390_native_fortran_unification_dag.v1.json

보존할 수리·물리 계약:
1. CompOSE reference mass 재표현은 물리적 epsilon,P,mu를 바꾸지 않는다.
   제거된 additive crust energy shift를 복원하지 말라.
2. 모든 source row의 Euler 계약은 matching window 탐색 전에 검사한다.
3. CLD beta 평형은 bulk식이 아니라 전체 cell energy의 조성 변분이다.
   표면/Coulomb 미분과 neutron/proton rest-mass 차이를 포함한다.
4. v6는 exploratory_proxy_CLD, source_calibrated=false이다.
   v6 solver의 Phi와 sigma는 Scurto2024 Eq24–25와 다르다. b의 의미도 다르다.
   scurto_interface.py에는 원전 식과 sigma 미분이 별도 구현되어 있지만
   solver/emitter에 연결되지 않았고 source-calibrated 모델이 아니다.
   함수/모델별 fitted parameters 확보 없이 원전 모델 재현을 선언하지 말라.
5. CSS 구간의 epsilon(P),n(P),mu(P),cs2는 해석식으로 평가한다.
   CSS와 core를 하나의 보간으로 다시 섞지 말라. 접합 양쪽 미분을 보존하라.
6. TOV 인과성 검사는 실제 RHS epsilon(P), 중심 초기화, P(n) 역산과 일치해야 한다.
   다른 보간의 인과성을 대신 근거로 쓰지 말라.
7. 부호 bracket만으로 Cowling 고유근을 채택하지 말고 최종 정규화 잔차를 검사한다.
   유한 주파수 격자의 근 완전성·f/p identity·radial stability를 주장하지 말라.
8. 거절/미생성 표식은 stale EOS 파일이 함께 있어도 우선한다.
   실패 입력/출력 hash와 source scope를 보존하고 실행환경 실패는 exit0으로 숨기지 말라.

현재 데이터:
- v6 proxy crust: SPG_M1(16행),DD2(19),DDME2(14),TW(14).
- M2/M4 placeholder, M3 내부 optimizer-bound, M5 구간 instability로 미생성/거절.
- census80설정/40평가/40source-marked, CSS8, Maxwell0.
- Maxwell20건 mass-convention gate, APR4건 C10, RG8건 beta계약 실패.
  이것을 일반적 Maxwell coexistence 부재로 해석하지 말라.
- 별 계산은 outer crust 없는 유한압력 절단모형이다. DeltaP=0 및 vacuum lapse
  정규화의 범위를 명시하고 완전한 중성자별 반경/주파수로 해석하지 말라.
- 과거 raw/v2–v5 데이터와 실패 receipt는 역사적 근거로 보존한다.

다음 작업 우선순위:
A. equation-only helper를 재구현하지 말라. 다음 최소 작업은 source bulk-model
   입력/identity, 정확한 2408개 질량행, fit 구현/가중치/계수를 식별하는 것이다.
   AME2016 파일 hash만으로 선택행/fit을 식별했다고 보지 말고 SPG surrogate와
   원전 GDFM의 동등성을 가정하지 말라. 입력을 식별하기 전 calibrated CLD는 HOLD.
   식별 후에 전체 energy와 모든 변분 미분을 구현·검증하라. 기존 proxy를 보존하라.
B. 밀도별 전체 변분 잔차, 여러 초기값/전역 최소 후보, 격자 정련을 점검하라.
   현재 finite root scan과 simplex 수렴은 전역 최소 인증이 아니다.
C. source-backed outer crust와 splice를 구성하고 P→0 경계를 갖는 TOV/Cowling
   계산을 수행하라. 계면 slow/rapid조건, mode identity, radial stability를 분리하라.
D. DD2의 질량상수 판본 차이, RG Q5 관례, composition masks와 Maxwellv2,
   EOS390 smooth kernels/전체 campaign의 남은 의무를 기존 DAG로 이어가라.

실행:
bash fortran_eos/build.sh build/compose_eos_fortran
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 -m pytest -q python_frontend/tests pipeline/tests
PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg python3 -m pytest -q eos_generator/tests
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran NS_COWLING_MATCH_OUT_ROOT=build/next_census python3 scripts/match_core_crust.py
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 scripts/run_css_survivor_stars.py build/next_census build/next_stars.json

출력 경로는 새 경로를 사용한다. 실제 최종 검사 수치와 CI identity는 audit receipt와
원격 Actions에서 확인하라. skip은 PASS로 바꾸지 말고, tolerance를 완화하지 말라.
로컬 /home/cosmosapjw/Dropbox/NS_Cowling은 미커밋 상태이므로 덮어쓰지 말라.
격리 checkout/worktree에서 변경을 검증하고 fast-forward push하라.
merge/force-push/원자료 삭제는 하지 말라. 발견한 결함은 RED 반례→수정→GREEN과
독립 물리 검산으로 닫고, source SHA·실행 명령·실패 분류·남은 범위를 인계하라.
```
