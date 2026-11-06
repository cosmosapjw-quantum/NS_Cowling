#  scan·root·정규화 계약

기준 물리식과 차원은 cowling.md, cutoff와 lapse는 tov.md를 따른다. 주파수의 근 찾기와 모드 식별을 구별한다.

- 공개 입력 frequency_hz를 omega_rad_s=2pi*frequency_hz로 변환한다. Bisection의 xtol 단위는 rad/s이다.
- scan은 성공한 인접 trial만 bracket으로 사용한다. 독립 표본의 작은 잔차만으로 root를 확정하지 않는다.
- coordinatebracket과 정규화 잔차는 동시에 확인한다. 미해결 noise/stagnation은 rootfailure로 남긴다.
- x=r/Rs, w=W/Rs^3, v=V/Rs^2; integral(w²+6v²)dx=1. 부호는 중심의 A0방향을 유지한다. 이 norm은 단순 진단이며 물리 관성 정의가 아니다.
- 기본 smooth node구간 [0.001,0.999], 상대 zero1e-7; W,V 각각 검토. 0/0→f,1/1→p1의 규칙은 이번 smoothbarotrope 가족에 한정하며 변화·해상도 검사를 동반한다.
- basis는 `scaled_WV_L2_dx_v1`, 식별은 `smooth_l2_WV_nodes_v1`, 배경은 `구간별 선형 q 보간`이다.
-  adaptive는 모든 q knot에서 종료한 뒤 다음cell로 이동한다. RK4의 substeps는 cell당 고정이며 equation/BC오류에 대한 독립정답은 아니다.

실행 근거: fixtures/modes/smooth_l2/results.json, refinements.csv, process_evidence.zip. 기준 제작계획의 root/identity구분을 따른다.
