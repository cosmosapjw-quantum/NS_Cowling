# 다음 coding agent 작업

NS_Cowling의 standalone core/crust composition-aware regeneration을 이어가라. 현재 원격 branch remediation/2026-11-15-m5a-maxwell-eos-contracts와 PR15의 실제HEAD를 먼저 확인하고 격리worktree에서 작업한다. 사용자 목표는 원래5core와13crust 전체의 독립 재생성이다. 전체상태 PARTIALLY_COMPLETE를 유지한다.

먼저 REPORT_KO.md, THERMODYNAMIC_CONTRACT.md, core_literature.md, crust_literature.md 및 data/research_generator_20261116/MODEL_COVERAGE.json, EXECUTION.json을 읽는다. 원본5core+13crust 추출/readback은 실제완료했다. source extraction과 독립미시계산의 완료조건을 혼동하지 않는다.

고밀도에서 c_s^2/c^2>1인 원전모델은 밀도구간별로 판정한다. 원본계산의 유효구간, tabulation domain, extrapolation, numerical support를 기록한다. 원본행과 실패구간을 보존하고 한구간의실패를 전체모델배제로 확대하지 않는다. 물리적오류·수치적미분오류·질량/조성규약 불일치를 원본대조로 구분한다.

현재실행경로
- scripts/build_apr_uniform.py — APR low/high Hamiltonian, npemu, 국소중성 Maxwell 및 전역중성 Gibbs 공존, 밀도지원검사.
- scripts/build_reference_cores.py — GM1원전RMF/npe와 BBB2원전표/parabolic 조성.
- scripts/build_skyrme_uniform.py — KDE0v1/SLy4 전체비대칭Skyrme 및 npemu.
- scripts/generate_ddrmf_core.py — native DD2/DDME2/TW, composition CompOSE출력.
- scripts/generate_gdfm_core.py — 저자Zenodo SPGM1–M5후보3mesonGDFM, NMP재현 및 composition출력.
- scripts/fit_candidate_spg_surface.py — AME2016 후보2408개에서 명시적 vacuum droplet surface fit; 원전과 동등한 calibration은 미확립.
- scripts/build_e390_03a_core_only.py와 build_e390_03b_crust_only.py —5core/13crust 원자료추출.

각스크립트의 현재 인수는 소스를 읽어 확인한다. 생성디렉터리는 새경로를 사용한다. retained source tables와 historical v6proxy를 보존한다.

미완료작업 우선순위
1. candidate_surface_fit.md와 data/research_generator_20261116/candidate_surface를 읽고 실행한5개 fit을 출발점으로 삼는다. 질량RMS15.2–22.4MeV와 Yp>=0.30769의 학습범위는 neutron-rich crust 표면 외삽을 보증하지 않는다. 원전의9magicnuclei calibration과2408massscore의 역할 및 정확한 finite-nucleus 에너지식을 대조한다. source model_par에는 sigma0,b가 없으며 후보fit은 명시적 별도 모델이다.
2. SPGCLD전체energy와각변분미분을 source GDFM+surfacefit으로연결하고 multistart/density정련/조성conservation을검사한다. 기존4mesonproxy M3/M5실패를 source모델실패로이월하지않는다.
3. GPPVA8종은 source-backed nonlinearRMF계수와 spatialmeson/CoulombTF방정식으로 셀크기·geometry를최적화한다. Scurto표면식을GPPVATF의source표면에너지로치환하지않는다.
4. BBB2PCHIP압력은원전Table3대비최대65.18%차이다. pressure_followup.json의 조성경로 미분 가설을 저자 구현/fit계수와 대조한다. 이 대안 압력은15개 원본값을0.6413% 이내로 재현하지만 기본 열역학적 압력으로 채택하지 않았다. APR provider는 Gibbs임을 확인했으며, 여전히 남는 bulk/provider 수치차이와 Skyrme 저밀도 압력 차이를 검증한다.
5. 기존nativeproduction목표에맞춰 새Python물리구현을 Fortran으로이관하고 독립기준값을대조한다.

핵심검증명령
bash fortran_eos/build.sh build/compose_eos_fortran
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 -m pytest -q python_frontend/tests pipeline/tests
PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg python3 -m pytest -q eos_generator/tests

core-crust결합은현재사용자범위밖이다. 실제계산값·source identity·허용오차·compositionmask·미실행항목을정확히보존하라. 테스트로닫힌수정만 commit/fast-forwardpush하고 PR은draft로유지한다. 원본과의비교가실패하면 편차를기록하고물리식/규약/수치표현을수정한다.
