# T05 고정된 smooth 폴리트로프 배경

`background.csv`는 `tools/calibrate_tov.py`의 실제 Fortran TOV 실행에서 얻은 고정 입력이다.
`metadata.json`과 `provenance/`의 실행·입력·출력 해시를 `tools.background_io.read_background`로 확인한다.
행의 물리량은 SI이고 q만 무차원이다. 원점의 행은 없으며 첫 행은 r0=1 m에서의 중심 급수 시작값이다.
중심압력 Pc와 첫 행의 P_start는 같지 않다. 원점 부근 접근은 이후 배경 adapter가 정본 중심 급수를 사용해야 한다.

본 데이터의 검증 범위는 TOV 평형뿐이다. Cowling 고유주파수, 고유함수 및 그 보간 허용오차는 아직 교정하지 않았다.
표면은 P_cut=1e22 Pa의 유효 절단면이며 정확한 진공 접합이 아니다. 마지막 nu는 그 절단면에서 Schwarzschild 시간값으로
정규화한 모델의 시간 기준이다. 표면을 외삽하거나 density/pressure를 0으로 덮어쓰지 않는다.

압력·밀도 등은 T04의 동일 해석적 EOS에서 계산했다. 791개 행의 저장 해상도 자체를 후속 진동의 정확도 보증으로 사용하지 않는다.
명시된 EOS 식과 중심·표면·단위·수치 설정은 provenance/input.nml 및 `docs/T05_TOV.md`에서 확인한다.
CSV의 재배열이나 열 변경 뒤에는 기존 metadata와 같은 파일이라고 부르지 않는다.
