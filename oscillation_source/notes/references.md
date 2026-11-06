# 참고문헌과 읽기 범위

[C1] R. R. Silbar and S. Reddy, “Neutron Stars for Undergraduates,” *American Journal of Physics* **72**, 892-905 (2004); erratum **73**, 286 (2005). DOI: 10.1119/1.1703544. arXiv:nucl-th/0309041. 구조와 물질 모형을 계산으로 연결하는 읽을거리다. 동반 핵물질 식은 본문에 표시한 가정과 매개변수로 직접 계산하며, 모든 수치를 이 문헌의 결과로 간주하지 않는다.

[C2] S. Yoshida and U. Lee, “Nonradial Oscillations of Neutron Stars with a Solid Crust: Analysis in the Relativistic Cowling Approximation,” *Astronomy & Astrophysics* **395**, 201-208 (2002). DOI: 10.1051/0004-6361:20021270. arXiv:astro-ph/0210591. 유체와 고체 지각을 포함하는 문제의 차이와 모드의 공간적 성질을 비교하는 데 사용한다. 탄성 지각의 spectrum을 이 강의 유체 예제의 정답으로 쓰지 않는다.

[C3] G. Montefusco, M. Antonelli and F. Gulminelli, “Frozen and $\beta$-equilibrated $f$ and $p$ modes of cold neutron stars: nuclear metamodel predictions,” *Astronomy & Astrophysics* **694**, A150 (2025). DOI: 10.1051/0004-6361/202452727. arXiv:2410.08008. 평형 압축응답과 고정 조성 응답의 물리적 구분에 사용한다. 이 강의의 두 방정식에 고정 조성 효과가 자동 포함된다는 근거가 아니다.

[N1] E. Hairer and G. Wanner, DOPRI5, Université de Genève, <https://www.unige.ch/~hairer/prog/nonstiff/dopri5.f>. Embedded Runge-Kutta 계수의 확인에 사용한다. 동반 controller의 모든 선택을 해당 프로그램과 동일시하지 않는다.

[U1] BIPM, SI defining constants, <https://www.bipm.org/en/measurement-units/si-defining-constants>.

[U2] NIST, 2022 CODATA adjustment, <https://physics.nist.gov/cuu/Constants/Table/allascii.txt>. 실제 사용한 상수는 `fixtures/constants/codata2022.json`에 있다.

[U3] IAU 2015 Resolution B3, arXiv:1510.07674. 명목 태양 mass parameter의 정의에 사용한다.

[F1] Fortran Package Manager, Manifest reference, <https://fpm.fortran-lang.org/spec/manifest.html>. `fpm.toml`의 target과 local dependency의 의미에 사용한다.

본문의 간단한 수치 예제와 두 코드의 대조값은 동반 계산 자료에서 나온다. 위 문헌의 수치 결과를 그대로 전사한 것이 아니다. 코드의 비교는 같은 물리 가정·단위·표면을 먼저 맞춘 경우에만 유효하다.


[A1] F. D. Swesty, “Thermodynamically Consistent Interpolation for Equation of State Tables,” *Journal of Computational Physics* **127**, 118–127 (1996), doi:10.1006/jcph.1996.0162. 하나의 thermodynamic potential에서 값과 미분을 얻는 일관성 원칙의 참고다. 논문의 bi-quintic Helmholtz scheme을 여기서 사용했다는 뜻은 아니다.

[A2] R. H. D. Townsend and S. A. Teitler, “GYRE: An open-source stellar oscillation code based on a new Magnus Multiple Shooting Scheme,” *MNRAS* **435**, 3406–3418 (2013), arXiv:1308.2965, doi:10.1093/mnras/stt1533. Stellar shooting의 한계를 다루는 더 고급의 대안이다. 동반 주 계산은 DOPRI/RK4 양방향 방법이며 GYRE를 호출하지 않는다. 별도로 대조한 연구 구현은 GL4 Magnus 전달을 사용하지만 GYRE 그 자체는 아니다.

[A3] E. Hairer and G. Wanner, *Fortran Codes*, nonstiff integration, DOPRI5, Université de Genève. <https://www.unige.ch/~hairer/software.html>. 이 패키지는 자체 DOPRI 단계·controller와 RK4를 사용하며, 문헌 코드를 전체적으로 복제하거나 그 모든 기능을 지원하지 않는다.


[A4] R. H. D. Townsend, R. V. Kuenzi and J. Christensen-Dalsgaard, “Tools for Characterizing the Numerical Error of Stellar Oscillation Codes,” arXiv:2505.20212v2 (2025). 코드 간 일치와 개별 코드의 오차 특성화가 다른 질문이라는 점을 읽는다. 이 강의의 경험적 여유를 그 논문의 변분적 오차 측정식이나 엄밀한 상한으로 동일시하지 않는다.
