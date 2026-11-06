# 공통 표기·단위·상수 수식표

이 수식표는 SI 변수와 동반 함수의 표기를 정리한다.

## 물리량의 의미

| 기호 | 의미 | SI 단위 | 대표 코드 이름 |
|---|---|---|---|
| \(r\) | areal radius | m | `radius_m` 또는 `r_m` |
| \(m\) | 중력질량 | kg | `mass_kg` |
| \(P\) | 압력 | Pa = J m\(^{-3}\) | `pressure_pa` |
| \(\epsilon\) | 정지질량을 포함한 총 에너지밀도 | J m\(^{-3}\) | `energy_j_m3` |
| \(n_B\) | baryon 수밀도 | m\(^{-3}\) | `n_b_m3` |
| \(\mu_B\) | baryon당 화학퍼텐셜 | J | `mu_j` |
| \(\rho_0=m_b n_B\) | 지정한 기준 baryon 질량의 정지질량 밀도 | kg m\(^{-3}\) | `rho0_kg_m3` |
| \(\rho_E=\epsilon/c^2\) | 총 에너지에 대응하는 질량밀도 | kg m\(^{-3}\) | `rho_e_kg_m3` |
| \(\omega\) | 각주파수 | s\(^{-1}\) | `omega_rad_s` |
| \(f=\omega/(2\pi)\) | 시간당 cycle 수 | Hz | `frequency_hz` |
| \(k=\omega/c\) | Cowling 식의 역길이 주파수 변수 | m\(^{-1}\) | `k_per_m` |
| \(\chi=dP/d\epsilon\) | 지정한 barotropic 압축응답 | 1 | `chi` |

\(k\)는 단순히 \(\omega/c\)로 정의한 변수이며, 별 안의 국소 음향 radial wavenumber와 자동으로 같지 않다. \(\rho_0\)와 \(\rho_E\)도 서로 바꾸어 쓸 수 없다. \(m_b\)는 함수 인수로 명시하고, 물리상수 파일에서 모델의 공통 nucleon 질량 선택까지 조용히 결정하지 않는다.

시공간 signature는 \((- + + +)\), 시간 인자는 \(e^{i\omega t}\)이다. 단위구면 조화함수에 대해 \(\xi^r=e^{-\lambda}WY_{\ell m}/r^2\), \(\xi^a=-VD^aY_{\ell m}/r^2\)이면 \([W]=L^3,[V]=L^2\)이다. 각좌표는 무차원이고 \(\xi^r\)는 길이 성분이므로 이 정의에서 차원을 확인한다.

## 정확한 정의와 측정값

`constants_id = SI2019_CODATA2022_IAU2015`을 사용한다.

- \(c=299792458\;\mathrm{m\,s^{-1}}\), \(h=6.62607015\times10^{-34}\;\mathrm{J\,s}\), elementary charge \(e=1.602176634\times10^{-19}\;\mathrm C\)는 현 SI 정의상 exact이다.[U1]
- \(G=6.67430\times10^{-11}\;\mathrm{m^3\,kg^{-1}\,s^{-2}}\)의 CODATA 2022 표준불확실성은 \(1.5\times10^{-15}\) SI이다. 전자·중성자·양성자 질량도 같은 CODATA 판본으로 고정한다.[U2]
- \(\hbar=h/(2\pi)\)는 기본 상수에서 계산한다. 유한 이진 실수로 저장할 때의 반올림은 SI 정의의 불확실성과 다른 개념이다.
- IAU 명목 태양 mass parameter는 \((GM)_\odot^{\rm N}=1.3271244\times10^{20}\;\mathrm{m^3\,s^{-2}}\)이다. 이는 실제 태양의 최신 측정 질량이 아니라 고정된 환산량이다.[U3]

테스트의 \(10^{-15}\) 수준 비교 허용치는 이 상수와 산술의 구현을 확인하는 기준이다. 측정 불확실성을 그 수준까지 줄였다는 뜻은 아니다.

## 독립 환산 기준

\[
1\;\mathrm{MeV}=10^6(1.602176634\times10^{-19})\;\mathrm J,
\quad 1\;\mathrm{fm}=10^{-15}\;\mathrm m,
\]
\[
1\;\mathrm{MeV\,fm^{-3}}=1.602176634\times10^{32}\;\mathrm{J\,m^{-3}},
\quad 1\;\mathrm{fm^{-3}}=10^{45}\;\mathrm{m^{-3}}.
\]

압력과 총 에너지밀도는 단위의 차원이 같지만 의미는 다르다. 이 환산 함수는 두 의미를 추론하지 않으며 어느 열을 입력했는지 호출자가 기록한다.

질량 출력의 정의를 \(q_\odot=G_{\rm ref}M/(GM)_\odot^{\rm N}\)로 고정하면
\[
M=q_\odot\frac{(GM)_\odot^{\rm N}}{G_{\rm ref}},\qquad G_{\rm ref}=G_{\rm CODATA2022}.
\]

`nominal_solar_mass_to_kg`와 역변환은 이 정의를 함께 사용한다. 별의 target mass와 출력에서 다른 \(G\)를 적용하지 않는다.

## 변환과 물리적 허용 판정의 경계

단순 배율 변환은 음의 변화량에도 정의된다. 예를 들어 음의 \(\delta\epsilon\)을 SI로 바꾸는 일은 음의 평형 에너지밀도를 가진 EOS를 허용한다는 뜻이 아니다. 따라서 일반 변환은 유한 signed scalar를 허용하며 EOS의 양성·causality·유효구간은 후속 EOS 계층에서 검사한다. \(\rho_0\)의 생성 함수는 정의에 따라 \(n_B\ge0,m_b>0\)를 별도로 요구한다.

모든 변환 결과는 `status_t`와 함께 반환한다. 실패에서는 status가 우선이고 출력 숫자는 사용하지 않는다. NaN/Inf, overflow 및 비영 입력이 정확히 0으로 사라지는 전체 underflow를 거절한다. 표현 가능한 subnormal 값은 금지하지 않지만, 상대 정확도 보증은 아니다.

Fortran의 `real(real64)`는 차원 검사 자료형이 아니다. \(2\pi\), 전체 부호, rest/internal-energy 의미는 차원만으로 검출할 수 없다. 왕복 변환, 독립 기준값, 의미가 고정된 입력과 deliberate-bug fixture를 함께 확인한다.

## 출처

[U1] BIPM, *SI defining constants*: <https://www.bipm.org/en/measurement-units/si-defining-constants>

[U2] NIST, *Fundamental Physical Constants, 2022 CODATA adjustment*, 전체 표: <https://physics.nist.gov/cuu/Constants/Table/allascii.txt>

[U3] IAU 2015 Resolution B3, arXiv:1510.07674, resolution의 nominal solar conversion 표: <https://arxiv.org/abs/1510.07674>

출처에서 전사한 값·단위·표준불확실성은 `fixtures/constants/codata2022.json`에 기계 판독형으로 보존한다. 상수표와 실제 측정 불확실성의 의미는 구별한다.
