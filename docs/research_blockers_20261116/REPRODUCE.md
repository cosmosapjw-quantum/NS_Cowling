# 재현과 실행 계약

기준 branch: remediation/2026-11-15-m5a-maxwell-eos-contracts. 게시된 exact HEAD에서 실행한다. Python3.12/gfortran13.3, LAPACK/BLAS를 사용했고 native 작업은 BLAS/OMP/MKL 각1 thread, 동시 작업 최대12 이하였다.

~~~bash
python3 -m pip install -e '.[test,plot]'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
bash fortran_eos/build.sh
python3 scripts/validate_blocker_integration.py regression /tmp/ns_blocker_regression_fresh
python3 scripts/validate_blocker_integration.py public /tmp/ns_blocker_public_fresh
~~~

검증 스크립트는 fresh 출력 경로를 사용한다. 반복 시 기존 자료를 덮어쓰지 말고 새 경로를 지정한다.

개별 작업:

- core: bash fortran_eos/tests/core_domain_run.sh /tmp/ns_core_domain_fresh
- BBB: python3 fortran_eos/tests/bbb_source_runner.py --output /tmp/ns_bbb_fresh
- lower: python3 fortran_eos/tests/crust_attachment_checks.py
- SPG inverse·TF: 각 lane 보고서의 source path, input archive, 명령과 resource 제한을 따른다.
- 공개 full paper/source cache는 git에 재배포하지 않았다. URL·hash·재취득 절차·필요한 추출 데이터는 source receipt에 기록한다.

과거 과학 결과의 SHA는 당시 producer identity다. 비관련 진단 파일 변경으로 전체 source집합이 달라도 변경되지 않은 실제 물리 dependency를 확인해 결과를 재사용했다. 실패 시도는 삭제하거나 성공으로 덮어쓰지 않았다.

공개 생성기 추가 출력은 sound_speed.csv와 attachment_diagnostics.json이다. 전자는 c_eq²/c²와 적용 status를, 후자는 실제 source endpoint의 gap/overlap·ΔP·Δμ를 제공한다. 진단의 존재나 hash 성공은 물리적 acceptance와 독립이다. 실패한 미분값은 빈 필드다.

일부 연구 runner는 논문 재현용 고정 data 경로에 결과를 쓴다. 출력 인자를 지원하지 않는 lower/SPG 재현은 별도 worktree에서 수행해 게시된 원자료를 보존한다.
