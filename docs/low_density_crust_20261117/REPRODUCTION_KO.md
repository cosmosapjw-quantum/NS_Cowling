# H 재현과 인계

## 기준 identity

- G base/published commit: 905b845ad0f14d2f1e7f67f9eb802691144850cd.
- H 작업 branch: research/low-density-crust-20261117.
- H의 정확한 local commit/tree 및 배포파일 SHA는 별도 DELIVERY_RECEIPT.json에 기록한다. Commit 자신을 참조하는 hash를 tracked file 안에 넣지 않는다.
- 이 단계는 H local commit과 검토 자료까지 준비한다. H remote push는 사용자의 다음 게시 승인 후 수행한다.

전체 변경 파일과 원시 계산 기록은 H Git commit 및 G를 prerequisite으로 하는 incremental bundle에 있다. Review archive는 보고서와 주요 영수증을 바로 읽을 수 있도록 제공하며, 원시 native STATE/input/stdout 전체는 bundle의 commit에 보존된다.

## 실행 환경

실행된 worktree는 /home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117이며, 고정 G worktree는 sibling crust_refinement_20261117이다. Python은 /home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python, Python 3.12.13, NumPy 2.5.4, SciPy 1.18.1, pytest 9.2.0이었다. OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS는 1이었다.

H는 native TF/model source를 바꾸지 않았다. 재사용한 실행 파일은 다음과 같다.

| 실행 파일 | SHA256 |
|---|---|
| G build/tf_candidate_G_final/tf_candidate_driver | 21f00ecd9cd492e37329e5690d342842c19025d4b649f97f99e84e592522c189 |
| G build/G_verified_general/compose_eos_fortran | 533d55ba2b5264e3e5248395baefc9463b4732b818e3e8af14beb774779f5356 |
| G build/G_verified_general/compose_eos_generate | 65ffd3f3e3a6a8e536104678ae4c56310d71e36654e978821f2ebe997857aa10 |

일반 native build의 실제 G 명령은 bash fortran_eos/build.sh build/G_verified_general이다. TF driver와 생성기 구성은 G source/실행 기록 및 H PLAN의 executable identity를 따른다. 다른 machine에서 rebuild한 바이너리가 같은 SHA일 것이라고 가정하지 않는다. Byte identity가 다르면 별도 실행으로 기록하고 동일한 source/수치 검증을 적용해야 한다.

## 단계별 입력과 실행 기록

1. Native 저밀도 continuation: scripts/extend_low_density_crust_20261117.py 및 scripts/refine_supported_crust_20261117.py. 각 PLAN, PROCESS, STATE receipt가 정확한 argv·branch anchor·mesh·최적화 조건·실행 파일 hash를 기록한다. Failed source를 filter해서 이어붙이지 말고 원래 순서의 실패점에서 support를 끊는다.
2. 기존 Hermite full-EOS admission: scripts/prepare_low_density_stellar_20261117.py INPUTS.json OUTPUT_DIRECTORY. G controls는 stellar/g_baseline, H final inputs는 stellar/h_final에 있다. 초기 h_initial 기록은 최초 공개된 outputs의 역사적 진단이며 최종 source 선택은 h_final을 따른다.
3. Causal opt-in admission/TOV: scripts/prepare_causal_low_density_stellar_20261117.py와 scripts/run_causal_low_density_backgrounds_20261117.py. stellar/causal_extension의 ADMISSION_PROCESS.json, TOV_PROCESS.json 및 사전 PLAN에 실제 실행 인자가 있다. 이미 존재하는 결과를 덮어쓰지 않도록 output directory를 새로 지정한다.
4. DDME-X control matrix: stellar/ddmex_comparison/PLAN.json. Core physics는 exact G sibling source를 사용한 실행이다. 네 source/closure를 같은 solver로 비교했고 production source hash를 기록한다. 보존된 원 실행 runner와 이후 수정을 구분하는 receipt가 있다. Source snapshot은 원래 scripts 경로에 복원해 재현하기 위한 bytes이며, __file__ 기반 ROOT를 사용하므로 snapshot 위치에서 직접 실행하지 않는다. Low-frequency residual failure는 재현 대상에 포함된다.
5. 전체 regression: verification/full_suite_wrapper.py와 FULL_SUITE_RECEIPT.json. 실제 명령은 Python -m pytest -q python_frontend/tests pipeline/tests였으며 고정 G native binary 경로를 환경변수로 전달했다. 1453 passed, 9 skipped, 4 subtests passed, actual exit0이다.

현재 device의 경로를 가진 계획은 재현 가능한 입력의 위치를 고정하기 위한 것이다. 다른 checkout으로 옮길 때에는 data/model/STATE/executable hash를 유지한 채 경로 mapping을 별도로 기록한다. Source input이나 original receipt를 새 경로 값으로 덮어쓰지 않는다.

## 의미론적 주의점

- Strict Maxwell 성공은 선택한 두 source branch와 실제 지원 구간에 한정된다. 모든 geometry의 에너지 최소나 unified EOS 타당성의 증명은 아니다.
- 하부 causal exponential family는 명시적인 phenomenological closure다. 실패 native knot를 수렴한 TF 상태로 바꾸지 않는다.
- 정적 TOV 해는 해당 closure에서의 계산이다. Radial stability, 최대질량, 새 causal 네 모델의 Cowling spectrum은 계산 완료 주장에 포함하지 않는다.
- DDME-X Cowling의 유한 scan, 실제 수용된 root와 거절된 후보를 구별한다. 전체 pytest 통과로 개별 물리 계산의 residual gate를 대체하지 않는다.
