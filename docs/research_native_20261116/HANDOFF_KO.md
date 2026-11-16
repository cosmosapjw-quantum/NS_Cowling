# 로컬 인계

아래 지시를 새 로컬 작업 세션에 전달한다. 정확한 게시 commit SHA는 이 파일을 포함하는 commit을 사용하고, branch tip이 이후 움직였으면 별도로 비교한다.

```text
NS_Cowling의 native EOS 생성 작업을 이어라.
Repository: cosmosapjw-quantum/NS_Cowling
Branch: remediation/2026-11-15-m5a-maxwell-eos-contracts
이번 연구 루프의 baseline: 4de49430560930be84d9a4147d71c3c9841da2db
먼저 git status와 exact HEAD를 확인하고, 사용자의 미커밋 작업을 덮어쓰지 말라.
docs/research_native_20261116/{REPORT_KO.md,REPRODUCE.md,INTEGRATION_REVIEW.md,HANDOFF_KO.md}
data/research_native_20261116/{EXECUTION.json,MODEL_STATUS.json,GOAL_DAG.json,RESEARCH_LOOP.json}
및 public_cohort/SUMMARY.json, source_grid_cohort/SUMMARY.json, regression/REGRESSION.json을 읽어라.

Fortran이 모든 production 물리 계산(functional, 평형, field, gradient, 최소화, phase 선택, EOS 및 native readback)을 수행한다. Python은 frontend/입력 전송/실행/파일 검증/표시에 한정한다. 보존된 Python physics는 명시적 historical oracle/test용이며 production fallback이 아니다.

core5와 crust13의 native 생성 경로를 우선 재사용하라. 전체 재감사나 기존 비싼 계산의 무조건 재실행을 하지 말고 현재 checksum·producer·변경분을 확인하라.
CompOSE와100%일치를 PASS gate로 두지 말라. 원전 충실도·독립 물리 검증·수치 지원구간·CompOSE 비교를 분리하라.
SPG 후보 sigma0,b를 원저자 fit으로 승격하지 말라. 원전의 fixed-density crossing과 별도 common tangent를 구별하라.
BSk source-fit, outer-anchored 처방, ETFSI를 구별하라. lower/inner common tangent를 원전 attachment의 필수조건으로 재도입하지 말라.
미해결/phase변경/서로 다른 segment 사이를 자동 보간하지 말라. 모든 실패·재시도·비정형 profile을 보존하라.
다음 실제 과제는 MODEL_STATUS와 GOAL_DAG의 남은 source/numerical 항목 중 하나를 선택해 실제 계산으로 진전시키는 것이다. 코드 최적화나 형식적인 재심사만 반복하지 말라.
PR15는 draft/open 상태를 유지하며 이번 인계만으로 merge하지 말라.
```

원격 작업 checkout은 /tmp/ns_eos_native_4de4943_20261116이었다. 기존 /home/cosmosapjw/Dropbox/NS_Cowling의 미커밋 변경은 건드리지 않았다. 새 로컬 checkout에서는 위 원격 경로를 사용하지 말고 자신의 저장소 절대경로를 사용한다.
