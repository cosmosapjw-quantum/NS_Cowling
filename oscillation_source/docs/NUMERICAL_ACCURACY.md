# 같은 문제의 작은 차이를 계산하기

전체 유도는 강의 제3·5·6·10장, 사용 안내는 `documents/Precision_Study.pdf`, 두 코드의 대조 결과는 `documents/Code_Comparison.pdf`에 있다.

같은 EOS nodal 변형과 목표 질량에서 두 모드의 변화가 수치적으로 분해된다. 예제의 Δf는 약 -0.0104485 Hz, Δp1은 약 -0.469742 Hz다. `examples/precision/reference_summary.json`은 읽기 자료이며 계산 실패의 대체값이 아니다.

`python tools/stellar_methods.py --output ../study`를 기본 실행으로 사용한다. 이 명령은 main, method, boundary, mass, null의 21개 단일별 계산과 절대 주파수 수렴, 고유함수 식별, 비순환적인 scalar 적분 검사를 연결한다. `variational_checks.json`의 마지막 차이가 앞서 정한 수치 여유보다 크거나 감소하지 않으면 성공으로 판정하지 않는다. 낮은 수준의 `precision_study.py`는 개별 실험의 도구이며 전체 검증을 대신하지 않는다.

중심 다음 항은 `Model(r0_m=100.0, center_order=2)`로 선택한다. 이 설정을 기본 1 m 설정과 한 refinement 수열에 섞지 않는다. 별도 출력 경로에서 시작점 민감도를 확인한다.

수치적 재현과 실제 핵물질의 정확도는 다르다. 수치 여유는 엄밀한 상한이나 신뢰구간이 아니며, 원표·보간 모형·섭동 closure의 물리적 불확실성을 포함하지 않는다. Scalar 검사는 연속 cold barotrope에 한정된다. 물질적 밀도 jump에는 별도의 경계항이 필요하다.
