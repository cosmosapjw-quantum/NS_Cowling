# 공존 후보와 SPG 미분 후속 인계

기준 `a0b158f3dc4af9825192bf6a185a6331619868a5`. 전체 5 core와 13 원본 crust의 독립 재생성은 **PARTIALLY_COMPLETE**다. production 물리는 Fortran이고 Python은 입출력 및 검사에 사용한다. 본 작업은 기존 원자료, 실패 및 후보 처방을 보존한다.

- TF의 두 공존 후보 각 중간 밀도에서 N128 여섯 초기 구조를 계산했다. 이 유한 표본에서 기존 공존선 아래에 검출된 후보는 없다. 이는 전역 최소 또는 연속 전이점 증명이 아니다.
- N256 TM1e는 두 branch root와 반간격 잔차 및 세 내부 공존선 표본을 계산했다. 결과와 한계는 `TF.md` 및 실행 receipt를 따른다.
- N256 FSU2R은 주 root와 반간격 잔차가 목표 이내였으나 첫 공존선 표본의 tube solve가 `REQUIRES_REFINEMENT`로 종료되어 **HOLD**다. 상태 거부의 정확한 원인은 미확인이며 과거 실패를 덮어쓰지 않는다. 다음 최소 행동은 해당 상태의 optimizer 종료 코드, 제약 활성 및 프로필을 기록하는 별도 native 진단과 초기값/branch 정련이다. 정지 허용오차 완화는 금지한다.
- SPG 760개 원 조성 상태에서 5개 step의 실제 전체 에너지 조성 미분을 계산하고 기존 해석 미분과 비교했다. `SPG.md`는 해석 미분의 제한된 검산이다. `(1-f)`의 잔차 감소는 관측 패턴이며 원저자 처방이나 보정 fit의 입증이 아니다. 원저자 코드·입력 및 다른 변분 방향은 여전히 별도 의무다.

원전 충실도, 독립 물리 검증, 수치 지원 밀도·mesh 범위, CompOSE 비교는 별도로 판정한다. 조성 미정의 부분을 0으로 채우거나 실패점을 보간하지 않는다. 전역 최소, 전체 원본 crust, 물리적 표면, radial stability, mode identity와 spectrum completeness는 해결하지 않았다. 본 작업에서 core–crust 접합이나 별 계산을 새로 추가하지 않았다.

재현과 source dependency 차이는 `REPRODUCE.md`, lane별 명령·실패·입력과 출력 hash는 `TF.md`, `SPG.md`, `data/research_coexistence_20261116/`에 있다. 기존 전체 체크섬과 원 자료는 그대로다. 정확한 출판 commit과 CI는 외부 git/CI 기록에서 확인한다.
