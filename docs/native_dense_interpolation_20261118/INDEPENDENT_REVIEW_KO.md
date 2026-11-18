# Native dense EOS 독립 검토와 게시 범위

**판정:** ACCEPT_SCOPED_NUMERICAL_RESULT (전역 EOS의 물리적 승인은 아님).

- 10개 고정 native 소스, 25,612개 직접 생성 상태, 5개 TF droplet 모델의 새 독립 정지해 30개.
- 독립 출처/물리/홀드아웃 확인 1096개, 미해결 차단 finding 0개.
- RMF·core에서 native 소스 CompOSE 저장 전후 동일성과 물리잔차 확인. APR 6→8→7 개별 phase와 BBB 입력 상한 유지.
- 전체 회귀검사 1,501 passed, 9 skipped, 4 subtests passed, 실행 종료코드 0. 시험 전후 코드 해시 동일.
- 원 EOS 25,612행의 상대 Euler identity 최대잔차 5.7e-14.
- 도표는 nonlinear 8개 방법만 포함하며 linear 보간법은 전부 제외. 별도의 원시 로그에는 음속 미정의 diagnostic을 보존.
- 실패한 TF 24개 신규 밀도점은 수용된 profile로 승격하지 않았으며 추가 geometry를 검증하지 않았다.
- 이전 단계의 sharp Maxwell 교차근은 고정했으며 새로운 interpolation family로 교차근을 다시 찾지 않았다.
- 동일 physical model을 정밀 계산한 결과에 대한 interpolation holdout error이지 실제 핵물질 EOS의 진리값이나 모델 간 물리적 우월성 판정이 아니다.

## 다음 실행 최소 작업

같은 source의 두 원 phase table에서 각 family로 독립적으로 Maxwell P₋(μ)=P₊(μ) 교차근을 재탐색하고, pressure preference, n jump, latent heat와 지원 구간을 source identity별로 비교한다.
