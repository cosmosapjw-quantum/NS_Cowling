# NL3 민감 geometry와 banded RHS 후속

기준 `43b7d8637fefd9bdb8240d9bf8092b921384d794`, 프로젝트 **PARTIALLY_COMPLETE**. 기준 원격 CI #621과 PR draft/open/unmerged를 확인했고 root manifest 5,202건을 검증했다. 이전 156개 unique mesh-density 계산, 1,807개 시도, SPG 후보와 원자료를 source/input identity에 따라 재사용했다. 새 core–crust 결합이나 stellar 계산은 하지 않았다.

## NL3wrL55: 실제 밀도·mesh·정지 조건 탐색

기존 0.085 fm⁻³의 N384 tube·uniform profile을 같은 초기 입력으로 사용해 n=0.085와 0.08625, N384/N576, SLSQP `ftol=1e-13/1e-14`를 각각 실행했다. 두 가지 초기 family만 다룬 **16개 개별 optimizer 실행**이다. 각 종료 시 raw profile이나 예외를 checkpoint했고 입력 및 solver SHA가 실행 중 일정함을 확인했다. 물리 KKT·radius 조건, energy selection tolerance는 변경하지 않았다. 더 작은 ftol은 optimizer의 종료 지시이며 물리 허용오차를 느슨하게 한 것이 아니다.

n=0.085의 tube−uniform 에너지 차이는 N384 **−2.4075385e−5**, N576 **−2.4073727e−5 MeV/핵자**로 차이는 1.6589e−9다. 두 ftol에서 에너지·반복 수·최종 geometry가 같았다. 이는 두 실제 실행의 관측 민감도이지 최적화 오차나 연속극한의 인증 경계가 아니다. n=0.08625에서는 이 두 초기 가지가 uniform으로 수렴하여 별도의 clustered 에너지 가지가 남지 않았다. 가지 소실을 crossing root로 바꾸지 않는다.

같은 새 밀도 0.08625에서 기존 source profile을 초기값으로 여섯 family를 N192/N384 각각 실제 실행했다. 두 mesh 모두 11개 시도가 수렴했고 누락 family는 없다. N192 선택의 `near_uniform_profile`과 N384의 `uniform`은 에너지들이 약 1.25e−12 MeV/핵자 범위에 있는 유한 수치 결과다. near-uniform 종별 변조를 새로운 물리상으로 승격하지 않는다. 기존의 0.085, 0.0875 anchor와 이번 0.08625를 이용해 mesh별 3행을 native write/readback했다. 새 0.08625의 N192/N384 두 mesh-density는 기존 anchor와 구분한다. 추가 N576의 0.085와 0.08625는 두 family 프로브이지 여섯 family 완전 실행이 아니다.

N192 세 행의 native readback은 serialization 및 mask 검사를 통과했지만 물리 검사 overall은 **FAIL(C10,C11)**이다. 위상 구분으로 동일상 미분점이 0개이기 때문이다. 이는 `NO_SAME_PHASE_DERIVATIVE_SUPPORT`이며 초광속 또는 불안정의 검출로 해석하지 않는다. N384 overall은 PASS이나 동일상 미분점이 한 개뿐이다. phase-change bracket을 가로지르는 보간은 허용하지 않는다. finite family 순서로 continuum transition이나 global minimum을 주장하지 않는다.

## SciPy banded RHS 호환성

기존 `_band_solve`는 float64 band 행렬과 `float128` RHS를 SciPy에 전달했다. dtype을 검사하는 RED 테스트를 보존하고, 호출 직전에 RHS를 명시적으로 float64로 변환했다. `screened_electrons`의 `longdouble` potential 상태와 업데이트는 바꾸지 않았다. 실제 원래 실패 입력 8개에서 수정 전후 전자 밀도·potential·chemical potential의 저장된 float64 결과는 모두 같았고, 수정 후 최대 PDE 잔차는 2.5441e−13이다. 이는 확보한 입력과 Linux 환경에서의 한정 결과다. 별도 테스트는 원래 near-vacuum 입력에서 Newton 측 RHS가 `longdouble`, SciPy band solve 입력이 `float64`임을 실제 호출로 확인한다.

계산은 cast 전 source `2d1f1edc…c1067c02`에서 완료한 뒤 수정했다. 이후 source는 `88704abb…229690df`, 전자 equilibrium 코드는 `8be3c675…1501bdca`로 변함없다. 기존 raw 기록은 원래 source와 묶고 변경 후 재평가와 혼동하지 않는다. 독립 검토는 source 변경 전 재평가를 차단한 preflight를 `SOURCE_VERSION_CHANGED_BEFORE_REEVALUATION`으로 기록했다. 이는 계산 실패나 과학적 integrity failure가 아니다.

## 검증과 잔여 범위

최종 local 회귀는 frontend/pipeline **973 passed, 9 skipped**, generator **142 passed**(기존 Matplotlib 경고 1건)이다. 최종 `SHA256SUMS_v4.txt` 5,202건을 통과했다. baseline 실패 입력, RED·GREEN 로그, 전체 실제 시도, 두 mesh native 기록과 source snapshot은 `data/research_next_20261116/`에 있다. 명령·SHA는 `RECEIPT.json`, 파일별 hash는 `SHA256SUMS.txt`를 따른다.

독립 검토의 재평가는 38개 반환 profile(두 가지 16 + 여섯 family midpoint 22)의 밀도·field·KKT와 물리량을 기존 실행 source와 변경 후 평가 source를 구분해 검사한다. 이는 cast의 일반적인 임의 입력에 대한 동등성 또는 전 영역 오차 상한이 아니다.

원 SPG surface fit은 미확정이며 common-tangent와 원전 fixed-density crossing은 별개다. BSk analytic 조성/thermodynamic μ 잔차 약 0.0255 MeV, 네 TF 부분 crust 및 다섯 changed-rest SPG 후보의 범위 제한도 유지한다. 추가 density/mesh/optimizer 수렴과 원전 fit 및 lower/inner 조성 정합은 남아 있다.
