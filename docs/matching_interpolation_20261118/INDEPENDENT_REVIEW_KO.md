# 독립 수치 감사와 범위 판정

- **검증 상태:** H source와 원래 phase doublet를 고정한 보간법 교차검증의 scoped numerical result를 수용한다.
- **증거:** `INDEPENDENT_AUDIT.json`의 **177개 재검증, 실패 0개**. 원래 H의 P/μ doublet, n jump 방향 및 latent heat 항등식 `Δε=μ_tΔn`을 독립적으로 검사했다. 최대 상대적 Maxwell pressure encoding 차이 없이, 절대 pressure 차이는 `1.735×10^-18 MeV fm^-3` 이하, 잠열 항등식 절대잔차는 `5.33×10^-15 MeV fm^-3` 이하.
- **회귀:** focused 108 passed; 전체 1490 passed, 9 skipped, 4 subtests passed, exit 0. 기존 native 바이너리 고정. 소스 185개의 실제 bytes는 실행 전후 변경되지 않았다.
- **구현 재현:** power/rational 보간기 8/8, exponential 8/8, slope-matched Bernstein 8/8. Cubic Hermite 2/8, naive P(n) 원래 energy endpoint 보존 0/8. 단계별 실패는 `RESULTS.json` 등에 별도 보존.
- **Sharp interface 인접 interval:** Hermite/Exponential/Rational 8/8, 순수 power 0/8. Density jump 자체에 interpolation을 적용하지 않았다.
- **보존된 제약:** 고정된 원래 Maxwell 교차점은 재탐색하지 않았다. 이 실험으로 `Maxwell cross-over discovery robust for all interpolation families`라고 주장할 수 없다. 모든 미시적 TF density gap을 채우지도 않았다.
- **해석:** mm 규모 반지름 변화는 선택한 finite shape ensemble의 sensitivity다. Frequency, Love/tidal deformability, radial stability, global TF geometry minimum을 검증했다는 뜻은 아니다.

**감사 중 발견해 수정한 사항:** 초기 보고서에 Maxwell root re-search 미수행을 명시하지 않은 서술상 결함 1개를 추가로 명확히 했다. 수치 계산값은 변경하지 않았다. Scalar array 반환과 NumPy boolean JSON 직렬화 문제의 최초 실패 로그도 `data/matching_interpolation_20261118/`에 보존했다.

**향후 최소 후속:** 원 crust/core branch의 두 source 표를 직접 읽고 same finite bracket에서 `P_-(μ)-P_+(μ)`의 root를 spline family별로 새로 찾는 독립 실험. 이번 검증의 성공률과는 분리해서 집계할 것. 사용자 지정 I 단계는 아직 시작하지 않았다.
