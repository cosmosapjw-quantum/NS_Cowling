# 다음 EOS 연구·코딩 작업 인계

```text
PROJECT=cosmosapjw-quantum/NS_Cowling
BRANCH=remediation/2026-11-15-m5a-maxwell-eos-contracts
PR=15
PR_STATE=DRAFT_OPEN_UNMERGED
THIS_LOOP_BASE=c8ff2d6c97a93f1ccf621051861d3091bfe94d81
READ_LIVE_HEAD_FIRST=true
GOAL=원래5core+13crust의 독립적 조성 포함 재생성
STATUS=PARTIALLY_COMPLETE

먼저 원격 실제HEAD/진행 중 변경을 읽고 격리 worktree에서 이어라.
최신 근거:
- docs/research_completion_20261116/REPORT_KO.md
- docs/research_completion_20261116/{ACCEPTANCE_AND_SCOPE,REPRODUCE,PHYSICS_REVIEW}.md
- data/research_completion_20261116/{MODEL_COVERAGE,RESEARCH_LOOP}.json
- data/research_completion_20261116/integration/REGRESSION_RECEIPT.json
- data/research_completion_20261116/gppva_tf/final/SUMMARY.json
- data/research_completion_20261116/tf_density_review/

완료된 계산을 다시 전수조사하지 말고 변경분만 읽어라.
Core5개807행, SPG후보513×5행, GPPVA8개226선택profile,
outer547행/18핵종 층, BSk22원전fit과 별도 경계 고정 재구성이 실행됐다.
최종frontend/pipeline773passed,9skipped; generator142passed.
각 실행이 생성 당시 코드SHA를 갖는다. 최신remote CI는PR본문과Actions에서 별도로 확인하라.

판정:
1. CompOSE100% 일치나 모델의 고밀도c_s²/c²>1만으로 전체를 block하지 않는다.
2. 원 논문functional·species·constraints·calibration·수치지원구간과 원본비교를 분리한다.
3. SPG표는 명시된 candidate surfacefit이며 원저자fit 회복으로 승격하지 않는다.
4. TF는 실제 공간 n/p/e profile과 총분율을 가진다. 미정의 A/Z·free/bound를0으로 만들지 않는다.
5. BSk22analyticfit/first-law수정/경계고정재구성/fullETFSI를 같은모형으로 부르지 않는다.
6. 전체구간global minimum 증명이나nativeport를 모든 생성의 불필요한 선행gate로 만들지 않는다.

다음 실질 작업:
A. SPG원래 finite-nucleus fitting prescription과model별계수를 확인하라.
   Carreau9핵종의48Ni추정값, symmetric n_sat고정면적,
   Scurto2408질량 선택·weight·density prescription을 구분하라.
   이미 실행한20개 대안fit과 저Yp surface불확실성을 활용하고,
   원본EOS압력을 맞추기 위한 숨은refit을 하지 않는다.
B. GPPVA를 필요한 densitygrid로 실제 재생성하고 geometry변화 근처를 정련하라.
   N192/N384radialmesh오차와 densityinterpolation오차를 구분하라.
   fixedN192midpoint/quarter의 직접계산 결과를 이어받고 실패점을 보간으로 메우지 않는다.
C. 공통lowercrust와 각innerfunctional의 질량/에너지/화학퍼텐셜 기준을 정합하라.
   원래common-mass SPG, 종별rest-energy복원SPG, BSkenergy-fit,
   boundary-anchoredBSk는 서로 다른 명시적variant다.
   원전analytic조성의mu_n잔차를 exactbeta-equilibrium으로 숨기지 않는다.
   필요한 경우 fullETFSI구현 또는 원저자입력을 확보하되,
   이미검증된근사표현의 사용가능성과 원전충실도를 별도표시한다.
D. APR/BBB/Skyrme 원본차이는 독립firstlaw/저자코드 근거와 함께 규명한다.
   BBB는 출판SNM/PNM에너지 재구성으로서 사용가능하지만 BHFkernel재해석과는 다르다.

수정은 실제실행·반례보존·독립검토·native조성readback으로 확인하라.
관련테스트와manifest/Gitblob바이트를 검사하고 fast-forward push하라.
PR15는draft/open/unmerged를유지한다. Core-crustmatching과새stellar/mode계산은현재범위밖이다.
```

이 파일은 자기 참조 commit SHA를 포함하지 않는다. 게시된 실제 HEAD는 git 및 PR에서 읽고,
실행 근거의 코드 SHA와 구별한다. 현재 report의 밀도·조성·원전 범위를 유지하면서 실제 계산을 진전시킨다.
