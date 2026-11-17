# Reproduce the 2026-11-16 completion loop

이 문서는 실제 receipt의 인자와 현재 parser를 대조한 재실행 안내다. 기존 실행 성공 여부는 각 receipt가 판정하며, 아래 명령을 적었다는 사실은 새 실행 증거가 아니다. **18개 실행 경로는 5개 core와 13개 whole crust의 완성을 뜻하지 않는다.** 원전 보정, 유효 구간 및 표현의 경계는 [REPORT_KO.md](REPORT_KO.md), [PUBLIC_GENERATION_CLI.md](PUBLIC_GENERATION_CLI.md), [SOURCE_COMPLETION_ADDENDUM_KO.md](SOURCE_COMPLETION_ADDENDUM_KO.md)를 따른다.

## Checkout, native build and regression

아래 명령은 이 문서를 포함한 최종 commit의 repository root에서 실행한다. 먼저 `git rev-parse HEAD`로 재현 대상 commit을 기록한다. 사용한 환경은 Python 3.12.3, NumPy 2.4.2, SciPy 1.17.0, gfortran 13.3.0이다. 새 출력 경로는 아직 존재하지 않아야 한다. 재실행할 때는 기존 증거를 지우지 말고 새 이름을 사용한다.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=python_frontend:scripts
export COMPOSE_EOS_FORTRAN_ROOT=.
export COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran
bash fortran_eos/build.sh build/compose_eos_fortran
python3 -m pytest -q python_frontend/tests pipeline/tests
PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg \
  python3 -m pytest -q eos_generator/tests
```

## Public entrypoint: native integration smoke

소스 checkout에서 `python3 -m compose_eos.cli`는 설치된 `compose-eos` entrypoint와 같은 dispatcher를 호출한다. 다음은 보존된 실제 smoke의 모델·밀도·점 수를 새 출력으로 재현한다.

```sh
python3 -m compose_eos.cli generators
python3 -m compose_eos.cli generate GM_GM1 build/reproduce_cli_gm3 \
  --points 3 --n-min 0.1 --n-max 0.2
python3 -m compose_eos.cli generate APR_APR build/reproduce_cli_apr3 \
  --points 3
python3 -m compose_eos.cli generate SPG_M1 build/reproduce_cli_spg17 \
  --points 17 --variant AME2408_candidate --workers 1
python3 -m compose_eos.cli generate GPPVA_FSU2H build/reproduce_cli_tf3 \
  --densities 0.095 0.0975 0.1 --resolutions 16 24 --workers 1
```

각 `DISPATCH.json`과 자식 receipt를 함께 읽는다. APR는 요청한 3점에 coexistence 끝점 2개를 추가한다. TF smoke의 N16/24 고밀도 3점은 호출·구성·native write/readback의 통합 검사이며, inner crust 전체나 geometry 수렴 검사가 아니다. 자식 exit 0만으로 성공을 인정하지 않는다. Dispatcher는 요청 행, native 출력 해시, profile identity와 유한값 및 수치 상태를 확인한다. 실패 시 `DISPATCH.json`, 자식 receipt와 `execution.log`를 보존한다.

보존된 실행: [GM1](../../data/research_finish_20261116/cli_smoke/core3/DISPATCH.json), [APR](../../data/research_finish_20261116/cli_smoke/apr3/DISPATCH.json), [SPG](../../data/research_finish_20261116/cli_smoke/spg17/DISPATCH.json), [TF](../../data/research_finish_20261116/generation/actual_tf_cli/DISPATCH.json). 별도 [TF evidence adversaries](../../data/research_finish_20261116/independent_review/cli_review/TF_EVIDENCE_GREEN.json)는 실제 복사본의 변조 거절 검사이며 새 TF 최소화가 아니다.

## Core and SPG scientific campaigns

모든 core 경로의 기본 물리 구간에서 161점을 생성하려면 다음을 사용한다. 공통 `--n-min`/`--n-max`를 강제하지 않으므로 모델별 구간 계약이 유지된다. 각 core의 원전 구현 수준과 고밀도 외삽 여부는 별도로 해석해야 한다.

```sh
python3 scripts/generate_core_eos.py build/reproduce_core161 --points 161
```

SPG는 5개 모델 × 5개 표면 처방을 모두 재최적화한다. 두 독립 밀도 격자와, 처방별 공통 밀도에서 실제로 다시 푸는 비교를 구별한다.

```sh
python3 scripts/reoptimize_spg_surfaces.py build/reproduce_spg129 \
  --points 129 --workers 4
python3 scripts/reoptimize_spg_surfaces.py build/reproduce_spg257 \
  --points 257 --workers 4
python3 scripts/compare_spg_reoptimized_surfaces.py \
  build/reproduce_spg129 build/reproduce_spg257 \
  build/reproduce_spg_common129 --points 129 --workers 4
```

이는 보간된 원본 EOS에 맞추는 fit이 아니다. 원저자 M1–M5 표면 보정은 미확정이며 `AME2408_candidate`는 명시적 후보 처방이다. 격자 상한은 각 후보 functional의 **clustered-side common-tangent endpoint**이고, 논문의 동일 밀도 energy crossing은 별도 진단으로 기록된다. 근거와 기존 결과는 [SPG receipt](../../data/research_finish_20261116/spg/RECEIPT.json), [grid129](../../data/research_finish_20261116/spg/grid129), [grid257](../../data/research_finish_20261116/spg/grid257), [common_grid129](../../data/research_finish_20261116/spg/common_grid129)에 있다.

## GPPVA geometry: fresh minimization or validated reuse

전체 새 최소화 재실행은 역사적으로 지정된 27개 밀도 window에서 152개 mesh-density node를 다시 푼다. 각 node의 6개 초기 geometry family와 모든 시도·실패를 보존한다.

```sh
python3 scripts/refine_gppva_geometry.py build/reproduce_geometry_all \
  --meshes 192 384 --subdivisions 2 --workers 4
```

이번 루프의 실제 후반 실행은 dense-to-banded field solver 전환 전에 완료된 **112개 node를 재사용**하고 **40개 node를 새로 계산**했다. 이를 재현하려면 먼저 동일 profile에서 최종 field operator로 재평가한 equivalence receipt를 만든 뒤, 남은 최소화를 실행한다. 기존 재평가는 1,254 attempt / 1,241 unique profile을 검사했으며, 고정 profile의 호환성 확인이지 다시 수행한 전체 최소화는 아니다.

```sh
python3 scripts/validate_retained_tf_profiles.py \
  data/research_finish_20261116/geometry/parallel_final \
  build/reproduce_geometry_compatibility --workers 4
python3 scripts/complete_gppva_geometry.py \
  data/research_finish_20261116/geometry/parallel_final \
  build/reproduce_geometry_remaining \
  --additional data/research_finish_20261116/geometry/fsu2h_recovery \
  --equivalence-receipt build/reproduce_geometry_compatibility/RECEIPT.json \
  --workers 12
```

다음은 위 재사용 경로의 결과를 export한다. 전체 새 최소화 경로를 택했다면 첫 인자를 `build/reproduce_geometry_all`로 바꾼다.

```sh
python3 scripts/export_gppva_geometry.py \
  build/reproduce_geometry_remaining build/reproduce_geometry_export
python3 scripts/summarize_gppva_geometry.py \
  build/reproduce_geometry_remaining build/reproduce_geometry_branches
python3 scripts/audit_gppva_morphology.py \
  build/reproduce_geometry_remaining build/reproduce_geometry_morphology.json.gz
python3 data/research_finish_20261116/independent_review/review_geometry_exports.py \
  --base build/reproduce_geometry_export \
  --output GEOMETRY_EXPORT_REPRO_REVIEW.json
```

마지막 checker의 `--output`은 파일명만 받으며 자기 script 옆에 쓴다. 새 이름을 사용한다. 이는 ledger 선택·profile·hash·native readback 검사이며 field 최소화를 실행하지 않는다. 두 mesh의 차이는 유한 탐색 진단이고 continuum 오차 보증 또는 전역 최소 증명이 아니다. Geometry 변화나 HOLD 공백을 가로지르는 연속 보간도 승인하지 않는다. 상세 명령·가정은 [GEOMETRY.md](GEOMETRY.md), 고정 profile 재평가는 [retained revalidation receipt](../../data/research_finish_20261116/geometry/retained_revalidation/RECEIPT.json), 독립 검사의 구분은 [INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md)를 따른다.

## NL3 upper-density follow-up

The final primary sample still selects a tube at 0.085 fm^-3. Four additional nodes at 0.0875/0.090 and N192/N384 were therefore executed with the same frozen solver. The two 0.085 anchor rows are copied by identity, not recalculated. Reproduce this separate bounded follow-up with fresh outputs:

```sh
python3 data/research_finish_20261116/geometry/extend_upper.py \
  data/research_finish_20261116/geometry/accelerated_final \
  build/reproduce_nl3_upper --model GPPVA_NL3wrL55 \
  --densities 0.0875 0.090 --anchor 0.085 --meshes 192 384 --workers 4
python3 scripts/export_gppva_geometry.py \
  build/reproduce_nl3_upper build/reproduce_nl3_upper_export
```

Primary plus extension contains 156 unique mesh-density nodes; its 18 native tables contain 158 rows because two anchors repeat. Both meshes select uniform profiles at the two new densities. The observed [0.085,0.0875] tube-to-uniform bracket is a finite-profile ordering interval, not a certified energy-crossing root. [Inputs](../../data/research_finish_20261116/geometry/nl3_upper_extension/INPUT.json) and [export receipt](../../data/research_finish_20261116/geometry/nl3_upper_exported/RECEIPT.json) retain the exact source and input identities.

## Partial crust construction and recheck

아래 archival runner들은 output 인자가 없고 정해진 `data/research_finish_20261116/crust_assembly` 경로에 쓴다. **최종 commit을 별도 detached worktree로 연 뒤 그 복사본에서만 실행한다.** 기존 checkout의 보존된 receipt를 작업 출력으로 사용하지 않는다.

```sh
eos_replay_parent="$(mktemp -d)"
git worktree add --detach "$eos_replay_parent/NS_Cowling" HEAD
cd "$eos_replay_parent/NS_Cowling"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=python_frontend:scripts
export COMPOSE_EOS_FORTRAN_ROOT=.
export COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran
bash fortran_eos/build.sh build/compose_eos_fortran
python3 scripts/solve_crust_tf_interfaces.py
python3 scripts/build_restored_spg_crust.py
python3 scripts/assemble_tf_partial_crust.py
python3 scripts/verify_crust_interface_profiles.py
python3 scripts/check_crust_interface_ordering.py
python3 -m pytest -q \
  python_frontend/tests/test_crust_assembly.py \
  python_frontend/tests/test_bsk22_extended_representation.py \
  python_frontend/tests/test_restored_spg_candidate.py
python3 scripts/seal_crust_assembly.py
python3 data/research_finish_20261116/crust_band_recheck/verify_band_profiles.py
```

이 절차는 commit에 포함된 이전 outer-crust 결과, coarse TF profile 및 gap-scan 입력을 명시적으로 재사용한다. 모든 원자료부터 outer crust와 inner EOS를 다시 만드는 단일 파이프라인은 아니다. TF 4개 부분 crust와 rest-energy 처방을 바꾼 SPG 5개 **새 후보 functional**의 부분 crust를 구성·검사하며, 원래 13개 whole crust 또는 core–crust 결합의 완성을 주장하지 않는다. 고정 profile band recheck도 새로운 nuclear 최소화가 아니다. [부분 crust 보고서](crust_assembly/REPORT_KO.md)와 [band recheck receipt](../../data/research_finish_20261116/crust_band_recheck/RECEIPT.json)에 입력 identity, 끝점 및 재개 조건이 있다.
