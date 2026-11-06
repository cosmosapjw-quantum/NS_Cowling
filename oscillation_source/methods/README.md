# 수치 방법과 EOS 표현

이 라이브러리는 적분·근 찾기와 한 열역학 함수의 값·미분·역함수를 제공한다.
TOV와 Cowling의 미분방정식, 중심·표면 조건은 `calculation/src/physics/`에 따로 있다.
`use ns_numerical_methods`로 적분과 근 찾기를 호출한다. EOS에는
`init_table_eos`, `cell_state`, `table_pressure_coordinate`를 사용한다.
소스를 읽을 수 있지만 내부 구현을 다시 작성하는 것은 필수 과제가 아니다.

중심 급수는 `use ns_center_series`의 `tov_center_second`와 `cowling_center_second`로 제공한다. 작은 반지름 조건과 같은 EOS 구간 안에 머무는 조건을 확인하며, 범위를 벗어나면 실패를 반환한다. 자세한 유도는 강의 제3·5장에 있다.
