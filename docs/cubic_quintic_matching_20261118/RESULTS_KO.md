# Cubic vs Quintic Hermite: dense native EOS, Maxwell/CSS, Cowling

대상: NS_Cowling, 기준 commit 2b32672da2b9b602942a362ff7208510708be459, 연구 branch research/cubic-quintic-matching-modes-20261118.

## 정의와 중요한 구분

T=0, beta equilibrium에서 n_B=dP/dmu_B, epsilon=mu_B*n_B−P, c_s^2/c^2=n_B/(mu_B*dn_B/dmu_B)를 사용했다.
Cubic P(mu)는 훈련 노드의 (P,n)을, quintic P(mu)는 (P,n,chi=dn/dmu)를 보존한다. Quintic의 chi는 정확한 미시모델 미분값이 아니라 오직 훈련에 사용된 n(mu)에 대한 PCHIP 미분으로 구했다. 따라서 높은 차수가 작은 물리오차를 보장하지 않는다.
모든 native EOS 오차는 동일 physical functional의 직접 계산한 withheld node에 대한 RMS 상대오차율(%). f/p1 및 R은 실제 참값을 알 수 없으므로 같은 1.4 M_sun, closure, 표면/상전환 가정에서 quintic−cubic의 상대 변화율(%)이며 observational error가 아니다.
Cubic_F/quintic_F는 F를 보간하는 native matcher, cubic_P/quintic_P는 P를 보간하는 native matcher이며 동일 parameterization 안에서만 차수를 비교한다. CSS는 sharp phase transition의 별도 prescription, Maxwell은 P,mu 교차에 근거한 first-order construction이다.

## 1. 원 EOS에서 직접 계산한 core/crust 물리량 오차율

| 원 EOS·source stride | 방법 | n_B RMS [%] | P RMS [%] | epsilon RMS [%] | cs2 상대 RMS [%] | sampled physical admission [%] |
|---|---|---:|---:|---:|---:|---:|
| RMF 5개 / 16 | cubic | 5.2272e-05 | 1.44152e-06 | 5.25047e-05 | 0.0127804 | 100 |
| RMF 5개 / 16 | quintic | 7.69128e-05 | 1.85268e-06 | 7.71015e-05 | 0.0152723 | 100 |
| RMF 5개 / 8 | cubic | 7.19738e-06 | 9.39237e-08 | 7.23678e-06 | 0.00396498 | 100 |
| RMF 5개 / 8 | quintic | 1.13857e-05 | 1.35007e-07 | 1.14401e-05 | 0.00610601 | 100 |
| Core 5개 / 4 | cubic | 0.000441804 | 3.10523e-06 | 0.00049419 | 0.139816 | 100 |
| Core 5개 / 4 | quintic | 0.000497679 | 3.55189e-06 | 0.00055253 | 0.275804 | 100 |
| Core 5개 / 8 | cubic | 0.00151194 | 1.66079e-05 | 0.00163855 | 0.4744 | 100 |
| Core 5개 / 8 | quintic | 0.00152938 | 1.83415e-05 | 0.00166411 | 0.521794 | 100 |
| TF crust 5개 / H_fixed_knots | cubic | 0.00744232 | 0.000157273 | 0.00745366 | 0.209806 | 100 |
| TF crust 5개 / H_fixed_knots | quintic | 0.00659227 | 0.000164876 | 0.00660178 | 0.691836 | 100 |

각 sparse stride는 같은 native source에 대한 별도 검증이므로, 2개 stride를 서로 독립된 물리 모형으로 합산하지 않았다. APR의 서로 다른 Gibbs phase는 일절 넘나들지 않고, native source 자체의 acausal high-density 부분은 interpolation failure에서 분리했다.

### 모델별 동일 소스 n_B RMS [%]

| Model | cubic [%] | quintic [%] | quintic / cubic |
|---|---:|---:|---:|
| core/APR_APR | 7.50692e-05 | 4.03971e-05 | 0.538132 |
| core/BBB_BHF2 | 0.00337637 | 0.00341346 | 1.01098 |
| core/GM_GM1 | 8.31959e-05 | 4.41587e-05 | 0.53078 |
| core/RG_KDE0v1 | 9.53899e-05 | 0.000151607 | 1.58933 |
| core/RG_SLY4 | 7.51709e-05 | 0.000118829 | 1.58078 |
| rmf/DDMEX_ORIGINAL_TABLE_RHO0 | 2.04568e-05 | 3.88988e-05 | 1.90151 |
| rmf/DDMEX_PROVIDER_ROUNDED_N0 | 2.16195e-05 | 4.02817e-05 | 1.86321 |
| rmf/NL3_PROVIDER_PRECISION | 7.6616e-05 | 6.67564e-05 | 0.871312 |
| rmf/PK1_CANONICAL | 5.96894e-05 | 9.66717e-05 | 1.61958 |
| rmf/TM1_CANONICAL | 5.78194e-05 | 0.000112429 | 1.94449 |
| tf/DDME2 | 0.0188579 | 0.0172882 | 0.916764 |
| tf/DDMEX | 0.000312619 | 0.000535638 | 1.71339 |
| tf/NL3 | 3.03981e-05 | 6.97626e-05 | 2.29497 |
| tf/NL3wrL55 | 0.000201472 | 0.000244843 | 1.21527 |
| tf/TM1e | 0.0069026 | 0.00460008 | 0.666427 |

TF N128는 실제 stationary droplet geometry·선택한 native branch에서 새로 수용된 밀도점 30개에 한정된다. Failed optimizer mode 8은 실제 데이터가 아니며, 전 geometry global minimum은 미검증이다.

## 2. Native crust-core matching의 source-fixed 형상 비교

SPG_M1/M3 × PK1/SLy4에서 cubic_F, quintic_F, cubic_P, quintic_P, CSS 각각을 실제 Fortran native kernel로 재실행했다. 총 20개 CompOSE 생성·재읽기 완료. 연결 구간의 소리속도 극값은 다음과 같다.

| crust/core | method | bridge n_left, n_right [fm^-3] | min(cs2/c2) | max(cs2/c2) | CSS P_t [MeV fm^-3] / latent heat [MeV fm^-3] |
|---|---|---|---:|---:|---|
| SPG_M1__PK1_CANONICAL | cubic_F | 0.0765392, 0.102085 | 0.00977759 | 0.225068 | none |
| SPG_M1__PK1_CANONICAL | quintic_F | 0.0765392, 0.102085 | 0.00967182 | 0.350156 | none |
| SPG_M1__PK1_CANONICAL | cubic_P | 0.0765392, 0.102085 | 0.00803415 | 0.202823 | none |
| SPG_M1__PK1_CANONICAL | quintic_P | 0.0765392, 0.102085 | 0.00786168 | 0.298404 | none |
| SPG_M1__PK1_CANONICAL | css | 0.0765392, 0.102085 | 0.110342 | 0.110342 | 0.512779 / 14.6851 |
| SPG_M1__RG_SLY4 | cubic_F | 0.0765392, 0.11462 | 0.00827647 | 0.0247583 | none |
| SPG_M1__RG_SLY4 | quintic_F | 0.0765392, 0.11462 | 0.00813315 | 0.0342838 | none |
| SPG_M1__RG_SLY4 | cubic_P | 0.0765392, 0.11462 | 0.00749917 | 0.0298384 | none |
| SPG_M1__RG_SLY4 | quintic_P | 0.0765392, 0.11462 | 0.00556601 | 0.0314897 | none |
| SPG_M1__RG_SLY4 | css | 0.0765392, 0.121381 | 0.012214 | 0.012214 | 0.512779 / 0.259633 |
| SPG_M3__PK1_CANONICAL | cubic_F | 0.0811502, 0.102085 | 0.0107344 | 0.27003 | none |
| SPG_M3__PK1_CANONICAL | quintic_F | 0.0811502, 0.102085 | 0.010574 | 0.417529 | none |
| SPG_M3__PK1_CANONICAL | cubic_P | 0.0811502, 0.102085 | 0.00896621 | 0.247906 | none |
| SPG_M3__PK1_CANONICAL | quintic_P | 0.0811502, 0.102085 | 0.00878466 | 0.38495 | none |
| SPG_M3__PK1_CANONICAL | css | 0.0811502, 0.102085 | 0.133927 | 0.133927 | 0.565705 / 12.3958 |
| SPG_M3__RG_SLY4 | cubic_F | 0.0811502, 0.115942 | 0.00857806 | 0.0234625 | none |
| SPG_M3__RG_SLY4 | quintic_F | 0.0811502, 0.115942 | 0.00834581 | 0.0275807 | none |
| SPG_M3__RG_SLY4 | cubic_P | 0.0811502, 0.115942 | 0.00775309 | 0.0249263 | none |
| SPG_M3__RG_SLY4 | quintic_P | 0.0811502, 0.115942 | 0.0059079 | 0.0321867 | none |
| SPG_M3__RG_SLY4 | css | 0.0811502, 0.121381 | 0.0121998 | 0.0121998 | 0.565705 / 0.145534 |

Native match 생성/readback의 성공은 full stellar interface admission과 별개다. 실제 quintic_F 출력 중 PK1의 두 pair는 stellar 소비자의 upper Hermite 물리조건 검사에서 거절됐다. 원 source나 native matching 성공 상태를 임의로 수정하지 않았다.

## 3. Original-source Maxwell 교차근에 대한 보간 오차율

Reference는 동일한 native crust/core의 전체 dense source를 이용한 native Maxwell root. Cubic/quintic은 원 source를 각각 희소화해 P_low(mu)=P_high(mu) crossing을 새로 찾았다.

| crust/core | sparse level | method | mu_t error [%] | P_t error [%] | Δn error [%] | latent heat error [%] |
|---|---|---|---:|---:|---:|---:|
| SPG_M1__PK1_CANONICAL | sparse_4x | cubic | 5.17086e-07 | 0.000548894 | -0.00126076 | -0.00126024 |
| SPG_M1__PK1_CANONICAL | sparse_4x | quintic | -9.75184e-07 | -0.00103557 | 0.0022714 | 0.00227043 |
| SPG_M1__PK1_CANONICAL | sparse_8x | cubic | 4.97806e-06 | 0.00528456 | 0.0267884 | 0.0267934 |
| SPG_M1__PK1_CANONICAL | sparse_8x | quintic | -1.02872e-05 | -0.0109223 | -0.053523 | -0.0535333 |
| SPG_M3__PK1_CANONICAL | sparse_4x | cubic | 5.97232e-07 | 0.000617365 | -0.000901323 | -0.000900725 |
| SPG_M3__PK1_CANONICAL | sparse_4x | quintic | -1.03143e-06 | -0.00106728 | 0.00139263 | 0.0013916 |
| SPG_M3__PK1_CANONICAL | sparse_8x | cubic | 5.69645e-06 | 0.0058885 | 0.0314316 | 0.0314373 |
| SPG_M3__PK1_CANONICAL | sparse_8x | quintic | -1.00685e-05 | -0.0104193 | -0.0529682 | -0.0529783 |

SPG_M1/M3 → PK1 두 pair에는 full source native Maxwell crossing 하나씩이 있으며, SPG_M1/M3 → SLy4 두 pair에는 검사된 공통 support에서 교차가 없었다. 이는 특정 source window의 수치 사실이지 다른 geometry/phase를 포함한 일반적인 no-go가 아니다.
CSS 전이압력과 잠열을 Maxwell 값과 비교할 수는 있지만 두 방법의 다른 물리 가정 때문에 그 수치를 interpolation relative error로 정의하지 않는다.

## 4. TOV/relativistic Cowling 비교 (M=1.4 M_sun)

H의 실제 Maxwell source 4개에서 lower/outer/surface closure는 고정하고 upper phase 내부 P(mu)만 cubic→quintic으로 바꿨다. Native matcher 자체를 재실행한 SPG→PK1 두 pair의 stellar EOS도 별도로 계산했다. 지표는 모두 같은 EOS source와 같은 mass 조건의 quintic−cubic 상대 변화율이다.

| paired same-input EOS | interface conversion | ΔR/R [%] | Δn_c/n_c [%] | Δf/f [%] | Δp1/p1 [%] |
|---|---|---:|---:|---:|---:|
| H_fixed_maxwell / H_DDME2_N128__primary | slow | -3.65605e-06 | -1.18025e-05 | -1.23513e-05 | 1.41503e-05 |
| H_fixed_maxwell / H_DDME2_N128__primary | rapid | -3.65605e-06 | -1.18025e-05 | 1.02644e-05 | 6.73733e-06 |
| H_fixed_maxwell / H_DDMEX_ORIGINAL_TABLE_RHO0_N128__primary | slow | -6.26611e-07 | -4.34154e-06 | 6.05543e-06 | 2.68966e-06 |
| H_fixed_maxwell / H_DDMEX_ORIGINAL_TABLE_RHO0_N128__primary | rapid | -6.26611e-07 | -4.34154e-06 | 1.69335e-06 | 1.68568e-06 |
| H_fixed_maxwell / H_GPPVA_NL3wrL55_N128__primary | slow | -3.05387e-06 | -3.24415e-06 | -2.25555e-06 | 1.30237e-05 |
| H_fixed_maxwell / H_GPPVA_NL3wrL55_N128__primary | rapid | -3.05387e-06 | -3.24415e-06 | 7.8736e-06 | 1.74851e-05 |
| H_fixed_maxwell / H_NL3_PROVIDER_PRECISION_N128__primary | slow | -7.92284e-07 | -9.11121e-07 | 4.55464e-06 | 2.96474e-06 |
| H_fixed_maxwell / H_NL3_PROVIDER_PRECISION_N128__primary | rapid | -7.92284e-07 | -9.11121e-07 | -3.62904e-07 | -4.363e-06 |
| native_matching / SPG_M1__PK1_CANONICAL | slow | 0.000151795 | -0.000668331 | -0.000293245 | -0.00895962 |
| native_matching / SPG_M3__PK1_CANONICAL | slow | 0.000166787 | -0.000737164 | -0.000697348 | -0.00856559 |

미승인/보류된 mode는 결코 주파수 0 또는 mode 부재라는 뜻이 아니다. fixed frequency window에서 다중 root, matching residual, interface junction 또는 (W,V)=(1,1) node 분류를 통과하지 못한 경우다.

### CSS/Maxwell vs smooth Hermite: 서로 다른 물리 constructions의 별 관측량 (차이이지 오차 아님)

| EOS source | method | conversion | R [km] | f [Hz] | p1 [Hz] |
|---|---|---|---:|---:|---:|
| SPG_M1__PK1_CANONICAL | cubic_F | slow | 14.4156 | 1851.9 | 5267.2 |
| SPG_M1__PK1_CANONICAL | cubic_P | slow | 14.4156 | 1851.86 | 5266.95 |
| SPG_M1__PK1_CANONICAL | quintic_P | slow | 14.4156 | 1851.86 | 5266.47 |
| SPG_M1__PK1_CANONICAL | css | slow | 14.4158 | 1852.1 | 5335.36 |
| SPG_M1__PK1_CANONICAL | css | rapid | 14.4158 | 1851.6 | 5263.73 |
| SPG_M1__PK1_CANONICAL | native_maxwell | slow | 14.4204 | 1842.31 | 5147.46 |
| SPG_M1__PK1_CANONICAL | native_maxwell | rapid | 14.4204 | 1842.28 | 5143.13 |
| SPG_M3__PK1_CANONICAL | cubic_F | slow | 14.4155 | 1851.93 | 5266.45 |
| SPG_M3__PK1_CANONICAL | cubic_P | slow | 14.4155 | 1851.91 | 5266.31 |
| SPG_M3__PK1_CANONICAL | quintic_P | slow | 14.4156 | 1851.9 | 5265.85 |
| SPG_M3__PK1_CANONICAL | css | slow | 14.4157 | 1852.19 | 미승인 |
| SPG_M3__PK1_CANONICAL | css | rapid | 14.4157 | 1851.74 | 5264.07 |
| SPG_M3__PK1_CANONICAL | native_maxwell | slow | 14.4204 | 1842.3 | 5147.35 |
| SPG_M3__PK1_CANONICAL | native_maxwell | rapid | 14.4204 | 1842.27 | 5142.83 |

Cowling은 ell=2, equilibrium Gamma1=Gamma_eq, perfect-fluid bulk, δg=0, 무한대 관측자 기준이고, slow/rapid conversion을 분리했다. Frozen-species, 탄성 crust, 방사 감쇠 또는 full GR mode를 확인한 결과가 아니다.
같은 조건에서 계산한 주파수 차이가 numerical refinement shift와 비슷하면 물리적으로 분리 가능한 차이라고 주장하지 않는다.

## 원시 재현 근거 및 미해결 항목

Native source 25,612 densities + TF N128 30 independent withheld source, 20 native matching outputs, 4 pairs independent Maxwell re-search, 1.4M_sun H/native matching TOV plus Cowling modes. 모든 최초 실패 로그와 root 미승인 사례는 새 파일에 남겼다.
미확정: 모든 pasta geometry 전역 최소, 미시 susceptibility 직접 재계산, 보간 함수와 source parameter의 모델 불확실성, low-frequency mode 완전 스펙트럼, 모든 20 native matching에 대한 stellar admission, full GR/radial stability. 원래 H/APR/BBB support 밖 외삽은 금지했다.

결론: Quintic은 cubic보다 항상 오차를 줄이지 않았으며, Maxwell 교차근의 정밀도에도 자동 이득이 없었다. 더 높은 차수를 채택하려면 P,n,chi의 진짜 source fidelity와 c_s²/도함수 수렴을 먼저 보여야 한다.
