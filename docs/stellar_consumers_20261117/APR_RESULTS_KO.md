# D APR 연속 Gibbs 경계와 retained-core 소비자

기준 commit은 C의 `539e6bd8de73f8654d121c0bc30023fb868dc05b`다. 실제 계산은
`data/research_matching_20261117/integration/real_pairs_v2/core/APR_APR/generated`
131행, n_B=0.085–0.5 fm^-3의 native APR 표를 그대로 사용했다. 단계 D의 이 부분은
연속 bulk two-charge Gibbs 경계의 thermodynamic consumer와 남아 있던
SPG_M1–5 × cubic_F/cubic_P/quintic_F/quintic_P/css 25개 조합을 구현·검증했다.
APR 자체의 완전한 별 관측량을 계산한 결과는 아니다.

## 구현과 물리적 범위

새 `python_frontend/compose_eos/apr_phase_consumer.py`의
`APRGibbsConsumer`는 phase 6, 8, 7을 세 구간으로 유지한다. 실제 source endpoint
(마지막 6 행과 첫 7 행)를 양쪽 interpolant에 공유하고, 각 구간에서 epsilon(n_B)를
cubic Hermite master function으로 둔다. knot derivative는 source mu_B다.

```text
mu_B = d epsilon / d n_B
P = n_B mu_B - epsilon
dP/dn_B = n_B d²epsilon/dn_B²
c_s²/c² = n_B (d²epsilon/dn_B²) / mu_B
```

epsilon/P 단위는 MeV fm^-3, mu_B는 MeV, n_B는 fm^-3이다. epsilon은 rest energy를
포함한다. c_s²/c²가 dimensionless이며 c=1의 SI 속도로 해석하지 않는다.
T=0, beta equilibrium, npemu와 기존 APR functional/parameter 및 bulk Gibbs 가정을
그대로 사용한다. 공통 lepton 배경, 전체 전하 중성, surface/Coulomb finite-size 항 생략을
유지한다. source의 `global_phase_minimum_status=NOT_PROVEN`도 그대로다.
기존 source 자체의 sampled tangent-plane 검증은 재사용하며 새 전역 최소 증명을 주장하지 않는다.

Endpoint는 opaque phase tag만으로 승인하지 않는다. 기존
`gibbs_state_at_volume_fraction(0/1)`의 직접 해와 n/P/epsilon/mu를 비교한다.
Source와 직접 endpoint 해 사이의 최대 상대 차이는 2.257e-11이다.
순수 phase 해와 Gibbs endpoint 해의 압력 차이는 5.8e-14 MeV fm^-3 이하,
mu 차이는 3.5e-13 MeV 이하다.

두 endpoint에서 n, epsilon, P, mu가 연속이며, epsilon의 두 번째 미분을 양쪽에서
별도로 반환한다. 전체 mixed interval을 constant-pressure plateau나 density jump로
대체하지 않았다. equilibrium 단일 유체 closure라면 displacement와 Lagrangian
pressure를 연속으로 연결한다. 이 APR 경계에 sharp-Maxwell transfer를 적용하지 않는다.
Frozen-species susceptibility나 composition-dependent buoyancy는 이 consumer가 결정하지 않는다.

Cubic epsilon의 epsilon''는 각 knot interval에서 선형이다. 양 끝에서 epsilon''>0이고
mu>0인 것을 확인하면 구간 전체에서 positive compressibility가 성립한다.
mu−n epsilon''의 미분은 −n epsilon'''이며 구간 내 부호가 일정하므로,
양 끝의 mu−n epsilon''>=0 확인으로 해당 cubic 구간 전체 c_s²/c²<=1을 검사한다.
이는 유한 표로 정의한 interpolant에 대한 검사이며 APR microscopic functional의
모든 상태에 대한 정리가 아니다. 일반 knot에서의 epsilon'' 점프도 유한격자 오차로 남는다.

## 수치 결과

| 경계 | n_B fm^-3 | P MeV fm^-3 | mu_B MeV | 낮은 밀도쪽 c_s²/c² | 높은 밀도쪽 c_s²/c² |
|---|---:|---:|---:|---:|---:|
| 6 → 8 | 0.20439495489345874 | 5.83908413832782 | 989.3170885322884 | 0.0856281067197 | 0.00193408455316 |
| 8 → 7 | 0.23738424893920132 | 5.90285482868961 | 989.6062923414233 | 0.00195227595325 | 0.0726339796526 |

경계 공유 상태의 차이는 부동소수점 반올림 수준이다. P=n mu−epsilon의 cancellation에
의한 endpoint readback 오차를 pressure inverse에서 명시적인 64 epsilon_machine
max(|n mu|,1) 범위로만 snap하여 사용자가 고른 미분 방향을 보존한다.

43개 공통 density sample에서 기존 functional solver와 비교했다.
원표를 수정하지 않고 격자를 절반으로 줄인 66-knot 표와 원래 131-knot 표를 비교했다.

| quantity | 66 knots 최대 상대 오차 | 131 knots 최대 상대 오차 |
|---|---:|---:|
| P | 5.1496943e-5 | 1.3870217e-7 |
| epsilon | 2.5895032e-9 | 1.8444011e-11 |
| mu_B | 5.1677014e-7 | 4.1857428e-9 |

분모는 max(|direct value|,1)이다. 두 격자의 실측 interpolation sensitivity이며
보편적인 asymptotic order를 주장하지 않는다. Endpoint의 독립 one-sided 3-point
finite difference는 h=1e-4, 5e-5, 2.5e-5 fm^-3로 계산했다. 가장 미세한 차분과
consumer dP/dn의 상대 차이는 6→8 낮은쪽 5.55e-7, 높은쪽 1.02e-5,
8→7 낮은쪽 4.93e-7, 높은쪽 1.57e-5다. 차분은 Python source의 정확한 endpoint에서
시작하여 native/Python root 위치의 약 1e-12 fm^-3 차이를 h로 증폭하지 않는다.
각 knot interval의 양 끝에서 측정한 c_s²/c²는 6에서 0.0146200–0.0856281, 8에서
0.00193408–0.00195986, 7에서 0.0726340–0.400713 범위다.
이 숫자는 interval 내부의 c_s² 극값 위치를 별도로 탐색한 결과가 아니다.

## 25개 조합의 실제 생성과 재독해

기존 API SUPPORT guard는 수정하지 않았다. `match_retained_apr`는 정확한
contiguous phase-6 source slice만 native matcher에 전달한다. 기존 native matching의
선택 endpoint 이후 low rows와 source 8→7 suffix를 물리 에너지 기준을 보존하여
연결한다. native의 원래 보고서는 low-branch 결과로 명시하고, 전체 행 수와 phase layout,
source/output row 범위는 별도 assembly 보고서로 기록한다.

25/25 조합은 실제 native kernel 실행 후 전체 CompOSE 표를 쓰고 native reader로
다시 읽었다. retained core의 n/P/epsilon/mu/F 최대 절대 차이는 모두 **0.0**이다.
Source phase와 retained crust/core composition 및 per-row availability mask는
정확히 보존했다. bridge의 abundance와 C3/C4 적용 가능성은 unavailable로 남는다.
원본 source file bytes는 모두 동일하다.

| 범위 | 실제 결과 |
|---|---|
| 20 smooth method × SPG 조합 | CONSTRUCTED_VERIFIED |
| 5 CSS method × SPG 조합 | CONSTRUCTED_VERIFIED |
| 5 CSS native low-branch 독립 admission | ADMITTED |
| 전체 assembled APR EOS의 full-star admission | 주장하지 않음 |
| APR의 full-star TOV/Cowling 관측량 | 이 lane에서 계산하지 않음 |

CSS crust jump는 native의 별도 first-order law로 남고 APR Gibbs endpoint는 continuous다.
기존 CSS two-provider stellar adapter는 retained APR internal boundaries를 거절하므로,
그 adapter를 통과했다고 보고하지 않는다. 전체 assembled 표를 하나의 smooth spline에
넣는 것도 승인하지 않았다. SPG input은 reconstructed inner crust이므로 완전한 표면
연결이 없고 APR upper support는 0.5 fm^-3다. 이 한계는 해당 source의 matching
구성 가능성과 구분한다.

390-cell census에서는 기존 DEFERRED_APR_GIBBS_SUPPORT 25개가 새로
CONSTRUCTED_VERIFIED로 바뀌는 delta다. 5개의 low-branch CSS admission을
전체 assembled-cell 또는 stellar admission으로 더하지 않는다. 기존 C 기록은 수정하지 않았다.

## 검증과 재현

11개 targeted regression이 통과했다. 연속 endpoint, distinct one-sided derivative,
first law, pressure inverse 및 roundoff, 실제 density support, 잘못된 endpoint/phase/Euler,
causality/compressibility rejection, retained energy/masks, source 불변성을 검증한다.
Assembly bookkeeping의 unit test stand-in과 실제 25-cell native 실행 증거를 구분했다.

Native binary SHA256:
`3893fbe061c2564ed4138211f7ea379a973851eb6e785fad22edcec0875152b1`.
실행 interpreter는 `/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python`이다.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONPATH=python_frontend:/tmp/ns_cowling_audit_20261117_pydeps \
COMPOSE_EOS_FORTRAN_EXE=$PWD/build/compose_eos_fortran/compose_eos_fortran \
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/validate_apr_consumers_20261117.py

PYTHONPATH=python_frontend:/tmp/ns_cowling_audit_20261117_pydeps \
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python -m pytest -q \
python_frontend/tests/test_apr_phase_consumer.py
```

두 명령은 exit 0이었다. 최종 test log는 `data/stellar_consumers_20261117/apr/regression.log`,
native campaign log는 같은 경로의 `validation.log`다.
`RUN_CONTEXT.json`은 실제 source/binary/consumer/validator SHA256을 기록하며,
`BOUNDARY_READBACK.json`은 endpoint와 derivative/interpolation sensitivity,
`CELLS.json`과 `pairs/*/CELL.json`은 개별 조합의 결과,
`pairs/*/compose/`는 실제 전체 출력,
`pairs/*/ASSEMBLY.json`은 source row와 전체 phase provenance다.
`SUMMARY.json`은 25/25 constructed/readback와 source 불변 결과를 보존한다.

재생성된 APR source의 composition은 기존 functional/bulk-Gibbs solve가 만든 것이다.
2025-03-14 provider 설명 [CompOSE APR eos.pdf](https://compose.obspm.fr/download/1D/NS/Classical/APR/eos.pdf)의
Gibbs prescription과 연결되지만, 공개 provider table의 미제공 composition을 관측·수록된
원자료처럼 취급하지 않는다. 다른 composition constraints 아래의 SRO Maxwell prescription과도
동일시하지 않는다.

최종 성공 실행은 25개의 native `match-core` 호출이며 nonzero exit를 예외 처리하는
`compose_eos.native.run_native`가 전부 정상 반환했다. 초기 construction pass도 25개를
실행했으므로 이 lane의 두 실행을 합친 matching kernel 수는 50이다. 최종 25개 결과와
기존 C의 delta를 혼동하지 않는다. 최종 실행과 11-test regression은 같은 remote process
PID 1937989에서 순차 실행되어 exit 0으로 끝났고, 그 원시 tool 결과를
`EXECUTION_RECEIPT.json`에 연결했다. 전체 390 상태를 재분류하는 대신
`MATCHING_390_DELTA.json`에 25개 cell ID와 변경 상태를 기록했다.
