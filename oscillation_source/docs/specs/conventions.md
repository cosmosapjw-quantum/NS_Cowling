# 변수와 계산 조건

우리는 차갑고 비회전인 구대칭 완전유체의 평형과 선형 비방사 진동을 계산한다. Metric signature는 (-,+,+,+)이며 Cowling 근사에서 metric perturbation은 0이다. 시간 의존성은 exp(i omega t), omega=2 pi f, f는 Hz이다.

별의 커널은 SI를 사용한다. r은 m, m은 kg, P와 epsilon은 J m^-3, n은 m^-3, mu는 J이다. epsilon은 정지에너지를 포함한다. chi=dP/d epsilon은 무차원이며 압축응답은 평형 barotrope와 동일하다. 회전, 탄성, 자기장, 초유체, 조성 동결과 반응률은 기본 식에 포함하지 않는다.

평형식은 m'=4 pi r^2 epsilon/c^2, P'=-(epsilon+P) nu', nu'=G(m+4 pi r^3 P/c^2)/(c^2 r^2(1-2 G m/(r c^2)))이다. 양의 중심 시작 반경과 양의 압력 끝점을 명시한다. 끝점 lapse는 nu(R)=log(1-2GM/(Rc^2))/2로 정한다. 이는 유한 압력에서 자른 문제의 시간 기준이지 정확한 진공 접합을 뜻하지 않는다.

W의 단위는 m^3, V의 단위는 m^2이다. Q=(omega/c)^2 r^2 exp(lambda-2nu) V+nu' W, W'=Q/chi-ell(ell+1) exp(lambda)V, V'=2nu'V-exp(lambda)W/r^2를 사용한다. 정칙 중심은 W=A r^(ell+1), V=-A r^ell/ell이고 표면에서는 Delta P=-(epsilon+P)exp(-lambda)Q/r^2=0을 요구한다.

모드의 근과 이름은 별도로 판정한다. 물질적 계면에서는 W와 Delta P가 연속이며 양쪽 상태를 분리한다. 주파수 구간의 유한한 탐색으로 스펙트럼 전체를 찾았다고 결론내리지 않는다.

Fortran 입출력은 NAMELIST와 name/type/value TSV를 사용한다. 명시적으로 성공한 결과만 물리값으로 읽는다. learner와 solutions는 같은 함수 이름을 사용하지만 동시에 링크하지 않는다. 어느 쪽을 실행하든 상대쪽으로 자동 전환하지 않는다.
