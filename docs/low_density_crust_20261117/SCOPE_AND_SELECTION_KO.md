# H 실행 범위와 계산 선택

G commit 905b845ad0f14d2f1e7f67f9eb802691144850cd의 원격 push가 확인된 뒤 H를 시작했다. H의 원격 게시 승인은 아직 없으며, 이 단계는 실제 native 저밀도 계산·whole-EOS 구성·조건부 stellar 비교까지 수행한다.

1. G에서 고른 다섯 crust model을 N128에서 저밀도로 연장하고 NL3, DDME-X, NL3wrL55, TM1e의 전이창을 재계산한다. DDME2의 전이창은 exact G N128 상태를 재사용한다. 에너지 영점, 모델 벡터, source acceptance 기준은 고정한다.
2. 원래 순서의 실패 밀도에서 native support를 끊는다. 한 차례의 반대 방향 accepted-neighbor 초기화로 구체적인 support gap을 점검하고, 같은 실패의 무제한 restart는 하지 않는다. 이는 optimizer 실패를 물리적 비존재로 해석하지 않기 위한 유한 계산이다.
3. 고정 E outer/surface에 대해 기존 direct Hermite 연결의 최초 및 다음 admissible knot를 관측량 계산 전에 선택한다. DDME-X 비교는 N64/.007, N128/.007, N128/.003, N128/.004의 네 case로 mesh/source 효과와 lower-join 효과를 구분한다.
4. 같은 endpoint가 일반 causal 조건을 만족하나 cubic 표현만 실패하는 경우, 유도된 exponential n/μ 연결을 명시적인 opt-in 현상론적 표현으로 구현한다. κ는 압력 적분과 endpoint만으로 정한다. microscopic TF source나 failed native gap의 복구로 세지 않는다.
5. 고정 1.4 M_sun background와 finite-frequency Cowling 비교를 계산한다. Slow/rapid conversion은 이상화한 경계조건이다. 실제 설정, 원시 residual, 거절 결과와 유한 repair를 기록한다. 전체 spectrum, radial stability, 모든 geometry의 global minimum은 이 단계의 결론이 아니다.

새 독립 H review는 한 번만 실시한다. 구체적인 오류의 수정과 그 재확인은 같은 review 안에서 수행한다. 재현 가능한 결과와 남은 실패를 확정한 뒤 H local commit을 만들고 게시 승인을 요청한다.
