# Native TF 압력 감소 구간 정련 — 2026-11-16

기준 HEAD `3f59e3a555fdeac07c3436ed2a17ee243945fa20`. 상태 **PARTIALLY_COMPLETE**. Production 물리·평형·최소화·생성·readback은 모두 Fortran이다. Python은 요청, 파일 identity, 하위 실행 및 결과/진단 표시만 수행했다. 원래 N16 source-grid 결과의 producer 의존 42개 hash와 root manifest 7,231개를 확인한 후 재사용했다. 원 자료와 이력은 변경하지 않았다.

## 실행 범위

원래 N16에서 압력 감소 인접 구간은 FSU2R 18, FSU2H 18, TM1e 17개였다. 각 모델에서 **가장 큰 감소 한 구간**의 정확한 원 grid 양 끝과 그 실제 중간점을 N16, N32, N64, N128에서 별도로 native 생성했다. 총 36 mesh-density 계산이며, 각 density에서 여섯 seed family의 유한 탐색이다. 이는 나머지 50개 원래 N16 감소 구간이나 각 구간 사이 연속 지지의 검증이 아니다.

| 모델 | 원래 N16 최대 감소 bracket [fm⁻³] | 원 N16 양 끝 ΔP | N128 중간점 포함 두 ΔP |
|---|---|---:|---|
| FSU2R | [0.048, 0.048999999] | −0.01929416 | +0.00155537, +0.00158761 |
| FSU2H | [0.050999999, 0.052000001] | −0.02033434 | +0.00180459, +0.00182496 |
| TM1e | [0.052000001, 0.052999999] | −0.04332305 | +0.00217271, +0.00214209 |

압력 단위는 MeV fm⁻³다. N16 양 끝의 P·ε·geometry는 과거 native output을 정확히 재현했다. 추가한 N16 중간점을 포함하면 압력 감소는 각각 1, 1, 2개 인접 구간에 남는다. 세 모델의 N32/N64/N128에서는 각각 세 점 모두 rod 후보가 선택되었고 두 인접 압력 차이가 양수다. 해당 3점 native C6/C10/C11도 통과했으나, 이는 표본 판정이며 이웃 mesh나 전 구간 인과성 증명이 아니다. N64와 N128의 같은 세 점에서 최대 압력 상대 차이는 각각 0.08589%, 0.06614%, 0.07185%다. 독립 최적화된 finite seed 결과이므로 이 차이를 순수 공간 이산화 오차나 continuum bound로 부르지 않는다.

N16에는 선택된 비정형 profile이 존재하며 낮은 압력을 단순히 삭제하지 않았다. 모두 `phase=9`로 부호화되므로 C10/C11의 같은-phase stencil은 서로 다른 비정형 branch를 혼합할 수 있다. 원 profile과 family/geometry 기록을 함께 보존했다. Pressure는 native `P=n mu_B−epsilon` 구성값이므로 Euler identity만으로 독립 압력 미분을 검증했다고 주장하지 않는다.

## 추가 branch 탐색과 실패

세 모델의 bracket 양 끝 각각에 대해 N64, 세 초기 family `(2,droplet),(2,hole),(3,hole)`, 초기 반경 30/45/70 fm의 **54개 Fortran driver 실행**을 추가했다. 54개 모두 process 성공과 native 정지 flag를 반환했다. 각 끝점에서 최소 에너지의 초기 family는 `(2,droplet)`이며, 기존 N64 선택보다 현행 `1e-6 MeV/핵자` 선택 guard를 넘어 낮은 후보는 발견되지 않았다. 같은 spherical-hole 초기 family 안에서도 반경에 따라 다른 정지 후보가 생겼다. 초기 seed를 최종 geometry로 자동 해석하지 않으며 유한 탐색을 전역 최소 인증으로 사용하지 않는다.

본생성기의 224개 시도에는 **7개 vector-field 평가 실패와 1개 `REQUIRES_REFINEMENT`**가 포함되어 있다. 216개 정지 결과와 함께 모든 실패·재시도·비정형 profile을 archive에 보존했다. 출력 또는 native readback 성공을 전체 EOS 물리 승인으로 해석하지 않는다. Nonuniform TF의 free/cluster 분해는 미정의 mask를 유지한다.

## 검증·근거·다음 행동

독립 검토는 producer source hash 162건, 출력 hash 90건, 선택 profile 36개와 54개 branch 실행을 확인했고 **PASS_SCOPED**를 부여했다. Native readback은 18회이나 이 중 6개는 마지막 mesh의 상위 별칭이며 고유 mesh 테이블은 12개다. 현지 검사는 native TF/dispatch 관련 **88 passed**, generator **142 passed**(기존 Matplotlib 경고 1건), native selection **4 checks PASS**였다. 미실행 전체 frontend 회귀는 기존 CI #626 근거와 구별한다.

`data/research_native_refinement_20261116/RECEIPT.json`에 실제 명령·source/executable hash·상태가 있고 `native_mesh_runs.tar.gz`에는 여섯 원 생성 디렉터리, `native_branch_probe.tar.gz`에는 54개 원 출력·오류·시도 receipt가 있다. 각 archive는 추출 후 native receipt와 재실행 스크립트로 확인할 수 있다. 독립 검토와 원 N16 실패 구간·native 검사 결과는 별도 JSON에 있다.

다음은 나머지 원 N16 감소 구간의 우선순위별 실제 mesh/branch 탐색과 이웃 밀도의 독립 압력 미분·연속 support 검증이다. 별도의 원 SPG fit, BSk analytic/ETFSI 및 source attachment 정합은 여전히 미완료다. 본 노드는 원저자 보정, complete crust, 모든 밀도에서의 인과성·압력 단조성 또는 전역 최소를 승인하지 않는다.
