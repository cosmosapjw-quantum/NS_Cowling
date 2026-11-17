# Native phase 보존 오류 수리 — 2026-11-17

기준 commit은 `4a44333377fc6252ad6bcc0c1024614277199f31`, 원격 CI는 #630 성공이다.
원래 dirty checkout과 기존 연구/실패 결과는 유지하고 `/tmp/ns_cowling_phase_20261117`에서
작업했다. 게시 commit과 그 HEAD의 CI 결과는 최종 인계 메시지를 따른다.

## Validation summary

native CLI가 전체 APR의 Iphase `6→8→7`을 `0`으로 바꾸며 성공을 반환하는 오류를
실제 보존 입력으로 재현했다. 최초의 포괄적 거부 수리는 기존 GPPVA crust의 내부 label도
거부하여 17 failed / 18 passed를 만들었다. 두 raw 결과와 시험 patch를
`data/research_matching_20261117/phase_retention_repair/attempts/overbroad_guard/`에 보존했다.

최종 수리는 retained crust/core 행의 Iphase를 Fortran에서 정확히 복사한다. 두 입력의
어느 label과도 충돌하지 않는 bridge label을 별도로 선택하며, source 숫자에 새로운
물리적 의미를 부여하지 않는다. `eos.compo`에는 phase-only 행을 출력한다.
Npairs=Nquad=0이며 abundance와 closure presence mask를 만들지 않는다.
endpoint 미분은 해당 방향의 연속한 source phase 안에서만 계산한다.
세 점을 확보하지 못하면 기존 최소 stencil 조건으로 거부한다.

native receipt는 전체 output_phase_ids를 기록하며, public API는 retained source와
다시 대조한다. pipeline admission과 CSS intake는 receipt/serialized label 불일치를
거부한다. phase receipt가 없는 과거 all-zero 출력은 읽을 수 있으나 새 admission/stellar solver 진입에서는 거부한다.

CSS intake는 직렬화된 phase 배열을 그대로 반환한다. 현 two-provider solver가 처리하지
못하는 retained source 경계가 있으면 보간 생성 전에
`UNSUPPORTED_RETAINED_SOURCE_PHASE_BOUNDARY`와 실제 경계 index pair를 반환한다.
경계가 없는 경우에만 solver의 private copy에 crust/high 역할 label을 붙인다.

**public API의 전체 APR SUPPORT 거부는 그대로다.** 여러 phase 사이의 허용 보간이나
interface 법칙은 이 수리로 승인되지 않는다. 실제 전체 APR native assembly 호출은
전송 검증을 위한 진단이며 과학적 matching/stellar admission이 아니다.

## Changed files

- `fortran_eos/src/compose_eos_matching.f90`, app: label 보존, phase-only output, stencil 범위, receipt.
- `python_frontend/compose_eos/api.py`: retained source label 검증과 provenance.
- `python_frontend/compose_cowling/css_adapter.py`: intake 보존과 명시적 지원범위 거부.
- `scripts/match_report_admission.py`: native phase receipt와 실제 출력 결합 검증.
- 관련 frontend/pipeline 반례 회귀와 `scripts/research_phase_retention_20261117.py`.

## Commands run

모든 물리 실행은 OMP/OPENBLAS/MKL thread=1, Fortran backend를 사용했다.
raw 명령·결과는 `data/research_matching_20261117/phase_retention_repair/`에 있다.

| Command | Result | Notes |
|---|---|---|
| `bash fortran_eos/build.sh build/compose_eos_fortran` | exit 0 | 실제 최종 executable로 실행 |
| `python -m pytest -q python_frontend/tests pipeline/tests` | 1162 passed, 9 skipped | 최종 리뷰 수리 후 전체 회귀 |
| `python -m pytest -q pipeline/tests/test_matching_pipeline_audit_20261117.py python_frontend/tests/test_css_cowling_adapter.py` | 37 passed | phase receipt 소실/변조, bridge 충돌 포함 |
| `python fortran_eos/tests/run_matching_native_adversarial.py ...` | 3 controls exit 0 | 기존 thermodynamic·interval·extrema 반례 |
| `python -m pytest -q fortran_eos/tests/matching_theory_soft_barotrope.py fortran_eos/tests/matching_theory_maxwell_topology.py` | 6 passed | 기존 이론 반례 |
| `python -m compose_eos.cli generate APR_APR build/phase_retention_actual/APR_APR --points 129 --n-max 0.5` | exit 0, 131 rows | 두 phase endpoint 추가; native readback |
| `python scripts/research_phase_retention_20261117.py --apr build/phase_retention_actual/APR_APR/generated --out build/phase_retention_actual/assembly_003` | PASS_SCOPED_TRANSPORT | native 진단 assembly 다섯 방법; public guard 유지 |

## Numerical impact

SPG_M1은 기존 조건부 public_composition_inferred_20261116 산출물을 hash로 묶어 재사용했다.
APR을 실제 재생성하고 CSS·cubic_F/P·quintic_F/P에 전달했다. CSS 출력 288행, 네 smooth
출력 각289행은 모두 phase sidecar 포함 native CompOSE write/readback을 통과했다.
retained n_B와 Q1은 수치적으로 정확히 같고, retained APR F의 최대 오차는 0 MeV,
crust F의 normalization 재인코딩 오차는 1.137e-13 MeV다. 물리적 에너지 offset은 없다.

별도 Euler 대수 잔차는 최대 4.079e-16(relative), CSS의 독립 analytic ε(P), n(P)
검산 오차는 최대 6.662e-16(relative)다. 기존 2e-12 Euler와 5e-9 CSS gate를 사용했다.
다섯 readback의 native consistency가 통과했다. source phase를 보존하면 해당 경계를
가로지르는 finite-difference stencil은 검사에서 제외되므로 이것을 경계 보간의 승인이나
phase 평형 검증으로 해석하지 않는다. 연속 구간의 이론/미분 반례 검증은 별도로 유지했다.

전체 APR CSS intake에는 source 경계 (221,222), (232,233)가 그대로 전달되고, solver는
명시적 지원범위 오류로 멈춘다. 기존 단일-source-phase CSS의 TOV/Cowling 회귀는 유지했다.

## Reproducibility notes

`APR_APR/`는 새 생성 receipt/input/output, `assembly_003/`는 최종 실제 method별 input internal,
native stdout/stderr, CompOSE output와 RECEIPT.json을 보존한다. 입력 hash와 executable hash는
receipt에 있다. 재현 명령은 새 out directory를 요구하며 기존 시도를 덮어쓰지 않는다.
`SHA256SUMS.json`은 봉인된 산출물과 raw 실패의 byte inventory다.

## Failures / skipped checks

17개 회귀 실패를 만든 포괄 거부는 최종 구현에서 대체했다. 독립 리뷰에서 legacy fallback이 실제 최초 phase 소실 APR을 다시 승인하는 문제가 확인되어 제거했고, 원본 baseline 반례로 거부를 검증했다. 잘못된 비정수 Iphase가 serialization에서 잘리는 경로도 거부한다.
첫 검증 driver 실행의 Python
module shadowing 실패도 `actual_assembly_001.log`에 보존하고 frontend import 경로를 고쳤다.
9 skipped는 전체 회귀의 기존 입력/환경 조건을 따른다. 실패를 통과로 재분류하지 않았다.

게시 후 전체 root `SHA256SUMS_v4.txt` 검사에서 수정 파일 14개의 이전 hash가 남은
패키징 불일치를 확인했다. raw 실패는 `attempts/package_manifest_3ffab7f.log`에 보존했다.
기존 manifest 범위를 유지하며 의도적으로 변경한 파일의 항목과 이 수리의 새 파일만 갱신했다.
물리 source/input, 실행 산출물, 허용오차를 변경하지 않는 포장 목록 수정이다.

## Remaining risks

source fidelity, 독립 물리 검증, 수치 지원범위, CompOSE 숫자 비교는 별도 축이다.
이 작업은 source calibration, 원저자 EOS 동일성, 두 APR 내부 경계의 보간/연속성/전이법칙,
전체 APR stellar evolution을 승인하지 않는다. 18 native 생성 경로의 실행 가능성도 원모형
전체 EOS 재현 완료를 뜻하지 않는다. SPG 조건부 calibration과 기존 고밀도 acausality는
이전 보고서의 제한을 유지한다.

## Next validation step

각 source 경계의 승인된 구간·interface 법칙과 그 native backend/stellar 적용을 마련한 뒤
전체 APR SUPPORT gate를 검토한다. 현재 우선순위의 완료 범위는 label 보존·CompOSE 전달·
stellar intake에서의 정확한 HOLD까지다. 독립 diff 검토는 PASS_SCOPED_REVIEW이며 sibling DECISION_LOG.md에 근거와 해결한 P1/P2를 기록했다.
