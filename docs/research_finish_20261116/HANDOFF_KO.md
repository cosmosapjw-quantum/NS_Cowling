# 다음 coding agent 작업 지시

cosmosapjw-quantum/NS_Cowling의 PR #15와 `remediation/2026-11-15-m5a-maxwell-eos-contracts` 브랜치를 이어라. 이 문서를 포함한 최신 커밋의 exact SHA와 원격 변경분부터 확인하라. 이번 루프의 기준은 `573a31485bc9af30915e01fe74f259e4fa5ac54f`였다. 원본 dirty checkout을 건드리지 말고 별도 worktree에서 작업하라.

먼저 `docs/research_finish_20261116/{REPORT_KO.md,INDEPENDENT_REVIEW.md,PUBLIC_GENERATION_CLI.md}`와 `data/research_finish_20261116/{EXECUTION.json,MODEL_COVERAGE.json,GOAL_DAG.json}`을 읽어라. 이미 닫힌 실패의 감사나 전체 census를 무조건 반복하지 말고, 코드·입력 identity가 바뀐 의존 범위만 다시 검사하라.

목표는 원래 core 5종과 crust 13종의 조성을 포함한 독립 재생성이다. CompOSE 100% 일치를 필수 gate로 두지 말고 원 논문 functional·계수·열역학·변분 조건과 밀도별 물리 support로 심사하라. 원 모델의 acausality를 모델 전체의 구현 거절로 바꾸지 마라. 새로운 core–crust 결합이나 stellar 계산은 범위 밖이다.

실행 가능한 18개 경로와 원전 전체 완료를 구분하라. core 5종의 기존 807행은 현 코드와 출력 SHA를 확인해 재사용했고, BBB는 published-energy reconstruction이다. SPG 5모델×5표면 처방의 실제 재최적화·first-law·native composition 검증은 끝났지만 원저자 surface fit은 미확정이다. 후보의 정확한 2,408핵종 목록은 공개 CUTER와 일치한다. CUTER의 네 계수·curvature·0.04 MeV/핵자 처방을 SPG의 두 계수·2% mass 처방과 동일시하지 마라.

후속 물리 작업의 우선순위는 다음과 같다.

1. `geometry/exported`와 branch receipt에서 실제 지원 branch, 비정지 시도, ΔE의 radial-mesh 민감도를 읽고 민감한 전이 구간만 추가 밀도·mesh와 optimizer 정지 허용오차로 정련하라. KKT·radius 잔차와 작은 ΔE를 직접 비교해 energy-error bound인 것처럼 해석하지 마라. 유한 후보의 순서를 global minimum이나 continuum transition으로 부르지 마라. clustered branch 소실만으로 에너지 교차의 root를 만들지 마라. 미해결 family, 더 낮은 feasible 비정지 trial, artificial radius bound를 보존하라.
2. 원저자 SPG의 sigma0,b 또는 정확한 fit protocol을 확보하고 이미 명시한 후보와 구분하여 calibration을 실행하라. 입력 부재는 source-equivalence 주장의 한계이며 모든 후보 계산의 금지가 아니다. source EOS 출력을 역피팅한 값을 독립 원전 입력처럼 쓰지 마라.
3. common lower crust와 13개 inner functional의 정지에너지·chemical-potential convention을 구분해 나머지 crust 구간을 완결하라. 네 TF 부분 crust(.03 fm^-3까지)와 다섯 restored-SPG 별도 후보(.01 fm^-3까지)는 original 13개 full crust가 아니다. BSk analytic 조성/thermodynamic μ 잔차 약 .0255 MeV는 수치 root 허용오차보다 크며 미해결이다. 임의의 상대 energy offset이나 .002 fm^-3 raw append를 쓰지 마라.
4. branch별 조성·phase jump·readback·density convergence를 검증하라. 비균일 TF의 total n/p profile을 free nucleon/cluster fractions로 바꾸지 말고 unavailable native masks와 실제 profile companion을 유지하라. near-uniform 종별 변조를 새로운 물리 phase로 승격하지 마라.

전자장 near-vacuum/dense 실패의 최종 수리와 원래 실패 기록을 재사용하라. 최종 `spatial_tf_equilibrium.py` SHA는 `8be3c6756f1b04a38cb59fdeeacb5c0aed5c98d8337c8ed14c64e7091501bdca`, `spatial_tf.py` band solver SHA는 `2d1f1edc72b170b01194f2705e62905383a2d4292bd84cc58067bc2cc1067c02`이다. dense→band 전환은 같은 Newton 연산자의 순열이며, 원래 source identity와 1,241개 unique profile의 재평가 증거를 보존한다. 계산은 각 초기값 종료 시 checkpoint하라.

RestoredSPGCandidate의 immutable cache/deep copy/critical-cone 수리를 유지하라. public CLI는 fresh output, 요청 행, model/mesh, native/profile hash와 정상 정지·보존 제약을 확인하고 실패를 숨기지 않는다. 실제 core/APR/SPG/TF CLI 실행 증거가 있으므로 dispatcher의 존재만으로 구현을 주장할 필요가 없다.

실제 계산·수정·독립 검토 후 관련 회귀와 checksum을 검증하고 fast-forward push하라. PR #15는 draft/open/unmerged로 유지하고 force push나 merge는 하지 마라. exact SHA, 실행 명령, 원자료·실패, source/numerical/implementation 상태와 남은 최소 행동을 기록하라. 전체 목표는 아직 `PARTIALLY_COMPLETE`다.

## 로컬 인수 절차

기존 dirty checkout을 pull/reset하지 말고 현재 저장소에서 원격 ref를 fetch한 뒤 별도 worktree를 만든다. 게시 메시지 또는 PR #15에 기록된 exact SHA와 아래 `eos_closeout_sha`를 대조한다. 더 최신 HEAD라면 해당 변경분을 먼저 읽는다.

```sh
git fetch origin remediation/2026-11-15-m5a-maxwell-eos-contracts
eos_closeout_sha=$(git rev-parse origin/remediation/2026-11-15-m5a-maxwell-eos-contracts)
git worktree add --detach ../NS_Cowling_finish_review_20261116 "$eos_closeout_sha"
cd ../NS_Cowling_finish_review_20261116
git rev-parse HEAD
sha256sum -c SHA256SUMS_v4.txt
```

회귀 재실행 명령과 환경은 `REPRODUCE.md`를 따른다. 최종 재개 검사 결과는 `data/research_finish_20261116/integration/resume_verification/REGRESSION.json`, geometry 전체 합집합 검사는 `independent_review/GEOMETRY_UNION_RESUME_REVIEW.json`에 있다. 두 anchor 중복을 제외해 156개 계산으로 센다. 이번 인수 작업은 기존 완료 계산 재실행이 아니라, 위 우선순위에 따른 다음 최소 연구 노드 선택에서 시작한다. CI의 실제 성공 여부는 게시 SHA에 대응하는 GitHub Actions run으로 확인한다.
