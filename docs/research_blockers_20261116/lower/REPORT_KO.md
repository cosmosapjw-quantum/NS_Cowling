# BSk22 하부 crust: 공개 표현, 열역학적 적분, 원전 연결의 해소 범위

기준은 `a49aaf9f82ddfda484e0f965aefdfe73aa7918f2`다. 이번 범위는 R2–R9와 C2–C8의 하부 crust blocker에 한정한다. 기존 native lower kernel과 source-fit 출력은 변경하지 않았다. 새 `compose_eos_crust_attachment.f90`가 원전 경계 진단과 명시적 압력 적분 대안을 제공한다. 부모 작업의 공개 CLI 통합은 별도 검증 경계다.

## 판정

**완전한 ETFSI 재계산은 공개 BSk22 분석적 EOS를 구현하기 위한 필수 조건이 아니다.** 저자 배포 페이지는 개별 fit 출력을 모두 그대로 쓰는 대신 열역학적 항등식으로 밀도 또는 에너지밀도를 계산하여 일관성을 회복하라고 명시한다. 따라서 압력 적분은 원전에서 정당화된 표현 선택이다. 다만 새 microscopic ETFSI 최소화, 원저자 실제 joining algorithm의 재현, 조성의 완전한 변분적 일관성은 다른 주장이다. 이번 구현은 이 세 주장을 하지 않는다.

원전은 처음부터 HFB outer crust와 ETFSI inner crust의 근사 차이를 인정한다. 그러므로 full ETFSI 실행이 outer seam을 반드시 없애리라는 전제도 성립하지 않는다. 원전 density attachment를 내보내는 작업과 독립적인 Gibbs/Maxwell 평형을 구성하는 작업은 분리한다. 이 구분은 물리적 불일치를 은폐하지 않으며 모든 jump와 cross-seam 미지원 상태를 기록한다.

## 직접 확인한 공개 원전

1. [Pearson et al. 2018, 2019 정오표 반영판](https://arxiv.org/html/1903.04981v1), Appendix C.1–C.4, §3.3, Table13. HTML SHA256 `e1b9a6d5f5b60726dc3178e3397c90881b7fcb68e1d4f3ef59cd9f216d963240`.
2. [저자 분석 fit 배포 안내](https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html), [bskfit18.f](http://www.ioffe.ru/astro/NSG/BSk/bskfit18.f), 2023-02-13 수정본. 코드 SHA256 `4bd3b716f04e40c69165fa83b5cd8ecb03c30aa1f2b45342fb070cef09a3a99c`. 기존 SOURCE.json에 기록된 코드와 byte identity가 확인된다.
3. [Potekhin et al. 2013](https://www.aanda.org/articles/aa/full_html/2013/12/aa21697-13/aa21697-13.html), §3.1 및 결론. 압력 표현을 유지하면서 열역학적 적분으로 n 또는 rho를 교정하는 선행 구성이다. 이 논문의 BSk19–21 수치 오차를 BSk22 오차로 전용하지 않았다.
4. [BSk22 제공자 PDF](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK22/eos.pdf)는 원 ETFSI/HFB 계산과 별도 Ioffe analytic-fit 배포를 구별한다. 최종 CompOSE 표는 원 계산 결과이지 재생성 solver 코드가 아니다.
5. [Shchechilin et al. 2023](https://arxiv.org/html/2309.01591v1)는 깊은 crust에서 Fortran 수렴 문제를 개선하기 위해 Python/SciPy로 재작성했다고 기술한다. 이번 Exa/SciSpace/Parallel 공개 검색과 저자 배포 경로에서는 전체 BSk22 ETFSI 소스 저장소를 식별하지 못했다. 이는 한정된 미확인 상태이며 비공개라고 증명한 것도, 유료 장벽으로 접근 실패한 것도 아니다. 저자에게 연락하지 않았다.

원본 획득 내역은 `data/research_blockers_20261116/lower/SOURCE_INTAKE.json`에 있다. Ioffe HTTPS 인증서 체인 실패는 공개 HTTP 주소로 해소했고, 다른 서비스의 검색 결과와 기존 pinned hash로 identity를 교차 확인했다. 원본 전체 논문/코드는 기존 프로젝트 정책에 따라 ignored `build/blocker_lower_sources`에 보존하고 새 추적 파일에는 URL·hash·해석만 기록했다. 원문의 추가 업로드를 사용자에게 요구할 이유가 없다.

## 원래 근사와 fit 오차

| 대상 | 원전에서 확인한 범위 | 해석 |
|---|---|---|
| e_eq(n) fit, C1/HTML75 | 대략 1e-9–3 fm^-3; 전형적 상대 오차 약0.5% 이하, phase boundary 최대3% 미만 | neutron rest energy를 제외한 e_eq의 경험적 fit 오차; 전체 epsilon의 3%라는 뜻이 아님 |
| P(rho) fit, C4/HTML78 | log10(rho/[g cm^-3])≈6–16; 전형적 약1%, 경계 최대 약4% 이하 | 유한 원 데이터 비교 결과; 엄밀한 연속구간 상계나 P(n) 동일밀도 상계를 자동 제공하지 않음 |
| 저자 INCRST의 Ye | 코드 주석 BSk22 max.error 0.0047 | 무차원 조성 fit 오차; %.0047로 바꾸지 않음 |
| 저자 INCRST의 neutron mu | 코드 주석 수치0.045 /0.072MeV | 바로 위 주석은 rms/max.abs.for CMUn>0.2/max.abs의 세 표기를 두 수치에 연결함. 두 값의 통계적 의미를 확정한 엄밀 상계로 사용하지 않음 |
| HFB 대 ETFSI at drip | BSk22 e_eq=-1.85 대-1.78MeV; P=4.99e-4 대5.10e-4MeV fm^-3; Z=38 대40 | Table13, §3.3: 내부에너지 약4%, 압력 약2% 차이. 동일 nuclei 비교에서는 약5%. 수치 반복 허용오차가 아님 |

마지막 차이의 원전 설명은 ETF kinetic/spin-current 근사, 섭동적 proton shell 처리와 neutron shell 생략, 제한된 density profile, pairing 근사, HFB collective/Wigner 항 생략이다. 이를 gfortran 반올림 또는 적분 tolerance 실패로 분류하면 안 된다.

기존 17점 native source-fit seam의 deltaP=-3.349914e-6MeV fm^-3와 deltaMu=-0.05081565MeV, density overlap=-2.90064e-8fm^-3는 그대로 남는다. 표의 ETFSI/HFB 차이와 현재 smooth analytic fit seam은 서로 다른 비교이며 숫자가 같아야 하는 대상도 아니다. SPG/GPPVA 제공자의 `.002` attachment 및 TW PDF(BSk22) 대 live page(BPS) 충돌도 기존 SOURCE_ATTACHMENT_CONTRACT를 그대로 따른다.

## 직접 유도: 같은 에너지 기준의 압력 적분

T=0, 국소 전하중성, 한 개의 보존 baryon density n을 갖는 barotrope를 가정한다. epsilon은 rest energy를 포함하고, e=epsilon/n, P는 MeV fm^-3, mu는 MeV이다. 첫 법칙과 Euler 관계는

\[
 d\epsilon=\mu\,dn,\qquad P=n\mu-\epsilon,
 \qquad \frac{de}{dn}=\frac{P}{n^2}.
\]

공개 에너지 fit을 e_s(n), 압력 fit의 합성을 P_s(n)=P_fit(rho_s(n)), rho_s=epsilon_s/c^2로 정의한다. 압력을 적분할 때 rho의 argument는 원래 rho_s로 고정한다. 적분된 epsilon을 다시 P_fit에 대입하면 다른 암시적 EOS가 되므로 그렇게 하지 않는다.

고정 anchor n0에서 **원래 source energy를 그대로** e0=e_s(n0)로 선택하면

\[
 I(n)=\int_{n_0}^{n}\frac{P_s(u)}{u^2}\,du,
 \quad e_I(n)=e_s(n_0)+I(n),
 \quad \epsilon_I=ne_I,\quad \mu_I=e_I+P_s/n.
\]

경계에 맞추려고 추가한 Cn 항은 없다. source neutron mass939.5654MeV와 원 fit의 에너지 영점을 유지한다. 차이를 직접 쓰면

\[
 e_I(n)-e_s(n)=\int_{n_0}^{n}
 \frac{P_s(u)-u^2 e_s'(u)}{u^2}\,du.
\]

이는 source fit 사이 미분 불일치의 누적 보정이다. n=n0에서 정확히0이며 상수 offset이 아니다. 첫 법칙과 dmu_I/dn=P_s'(n)/n은 직접 미분으로 성립한다. 따라서

\[
 \Gamma_I=\frac{nP_s'}{P_s},\qquad
 \frac{c_s^2}{c^2}=\frac{dP_s}{d\epsilon_I}
 =\frac{\Gamma_I P_s}{\epsilon_I+P_s}.
\]

수치적으로 x=ln(n/n0)를 써서 I=integral P_s(n0 exp x)/(n0 exp x) dx를 적분한다. native 기존 adaptive Simpson을 재사용했고, API 기본 상한은0.002fm^-3, 별도 extended=true에서만0.01fm^-3다. 이 식은 평형 barotrope의 미분이다. 동결 조성 perturbation 또는 독립 species chemical derivatives가 필요하면 해당 비평형 functional이 추가로 필요하다.

기존 outer-anchored 대안은 다른 anchor 계약이다. 물리적으로 계산한 outer 상태 (P_L,mu_L)와 P_s(n0)=P_L를 먼저 정한 뒤 e0=mu_L-P_L/n0를 사용한다. 이 경우 n0 위치의 일차 변동은 압력 일치 때문에 상쇄되어

\[
 \delta e_I(n)=\delta\mu_L-\delta P_L/n_0
 +\int_{n_0}^{n}\delta P_s(u)/u^2\,du
\]

가 된다(최종 n 고정, anchor root의 일차 변동 포함). 원 source anchor와 outer 물리 anchor를 같은 representation 이름으로 섞지 않는다. 이번 신규 구현은 전자를 제공하고 기존 후자의 출력은 변경하지 않았다.

## 오차 전파와 실제 구현

동일 n에서 |P_s-P_true|<=r P_true, 0<=r<1의 **구간 전체 envelope를 별도 가정**하고, 고정 n0의 에너지/입자 anchor 오차를 sigma0라 하자. 그러면 a=r/(1-r)에 대해

\[
 |\delta e_I|\le\sigma_0+a I,\qquad
 |\delta\mu_I|\le\sigma_0+a I+a P_s(n)/n,
 \qquad |\delta\epsilon_I|\le n(\sigma_0+a I).
\]

P_s>0, n>=n0를 사용했다. 단위는 각각 MeV, MeV, MeV fm^-3다. 적분값의 수치 추정오차 q는 별도 기록하며, model-bound 추정 자체에 미치는 a q까지 포함하면 합산 수치기여 추정은 (1+a)q 및 n(1+a)q다. adaptive Simpson의 q를 엄밀 interval-arithmetic 상계라고 부르지 않는다. n0 자체가 불확실하면 source anchor의 sigma0에 placement 효과도 포함해야 한다.

논문의 약4% P(rho) residual을 위 r에 바로 대입하면 동일 n의 rho_s(n) 오차 전파가 빠진다. 예를 들어 |delta ln rho|<=eta와 |d ln P/d ln rho|<=L이 추가로 확보되어야 좌표합성 오차 약L eta를 통제할 수 있다. 논문 error plot만으로 이 강한 전제를 확보했다고 선언하지 않는다. 압력 함수값의 envelope만으로 P' 또는 sound-speed 오차도 제한할 수 없다.

`bsk22_reference_integral(n0,n,anchor_error,pressure_envelope,state,budget,status,epsabs,extended)`는 이 조건부 계산을 실행한다. `state`에는 원 `source_epsilon`, `source_pressure`, source thermo residual과 composition mu residual이 유지된다. budget은 anchor, 적분값, 모델 상계, quadrature estimate, 원 source energy와의 차이를 별도 필드로 반환한다. 원자료를 fit해서 error envelope를 축소하지 않았다.

`.002`에서 예시 가정 sigma0=.01MeV,r=.04를 넣으면 e 오차 상계0.119044MeV, mu 상계0.164332MeV, epsilon 상계2.38088e-4MeV fm^-3가 나온다. 이는 입력 가정의 결과이며 BSk22 원전으로 인증된 불확실도 수치가 아니다.

## source attachment 구현

`source_attachment_diagnostics(left,right,boundary,diagnostic,status)`의 endpoint vector는 `[n,epsilon,P,mu]`다. 양의 유한값과 nominal boundary의 양측 endpoint 사이 포함 여부만 확인한다. 원전 density overlap도 허용해 delta_n<0로 기록하며, density gap은 delta_n>0다. deltaP, deltaMu, deltaEpsilon, 각 endpoint Euler residual을 반환한다. `source_prescription_valid`는 이 자료 구조의 유효성일 뿐 물리 평형 판정이 아니다.

`cross_interval_supported=false`, `equilibrium_established=false`는 항상 유지한다. P나 mu가 우연히 가까워도 미계산 구간의 interpolation, phase coexistence, 전역 안정성을 생성하지 않는다. `source_attachment_supports`는 해당 두 endpoint만 지원한다고 답하며 전 구간 EOS support query가 아니다. 원 lower/inner source 표와 upper functional을 자동 이어 붙이지 않는다.

## 실제 검산과 원전 chemical-fit 차이

`fortran_eos/tests/crust_attachment_checks.py`가 gfortran13.3, F2018, -O2, -fcheck=all로 세 모듈과 전용 driver를 컴파일했다. BLAS/OMP/MKL thread는 각각1이다. 41개 log-spaced n0=2.69403e-4–.002fm^-3에서 독립 Python 압력 구현 및 SciPy quadrature와 비교했다.

| 검산 | 실제 최댓값 |
|---|---:|
| 독립 적분의 e 차이 |1.51144e-11MeV|
| finite-difference d epsilon/dn -mu |1.20524e-8MeV|
| finite-difference n²de/dn -P |1.56659e-11MeV fm^-3|
| sound-speed 미분 비교 |1.57576e-12 (c_s²/c²)|
| source energy와 적분 energy의 차이 |0.0641632MeV/입자|
| source-energy-anchored mu와 author-code composition mu 차이 |0.0572895MeV|

Euler·baryon·charge closure, 양의 compressibility, 0<c_s²/c²<1을 검사했고 NaN/음의 error/잘못된 support 입력을 거부했다. 1e-6→1e-12MeV quadrature tolerance refinement에서 마지막 두 e 차이는 약1.4e-12MeV였다. 이는 모델 error가 수치 tolerance보다 훨씬 크다는 구분을 지지한다. source-fit17점의0.086913MeV 및 outer-anchored 확장 결과의약0.0255MeV와 이번0.0572895MeV를 같은 통계로 섞지 않는다.

추가로 정오표 반영 논문의 C18/HTML92는

\[
 \mu_n^*=\frac{\hbar^2(3\pi^2 n_f)^{2/3}/(2m_n)}
 {1+9.41\sqrt n+7.62n+259n^2},
\]

반면 배포된 INCRST와 기존 native 구현은

\[
 \mu_n^*=\frac{0.211\,(m_nc^2)n_f^{2/3}}
 {1+9.41\sqrt{n_f}+7.62n_f-259n_f^2}
\]

를 사용한다(계수는 명시된 fm,MeV 관례). total n 대 free n_f와 마지막 항 부호가 모두 다르다. 이 원전 간 차이에 대해 공개 코드가 우선인 현재 representation을 유지하며 논문식으로 몰래 바꾸지 않았다. 해당 식 오류의 궁극적 저자 의도를 결정했다고 주장하지 않는다.

`crust_attachment_fit_variants.py`가 같은41점에서 두 식을 직접 계산했다. max|mu_code-mu_paper|=0.0915368MeV, max|mu_paper-mu_Euler_source|=0.101153MeV, max|mu_code-mu_Euler_source|=0.0874686MeV였다. 이41점의 마지막 수치는 과거17점0.086913과 표본이 다르다. `CHEMICAL_FIT_VARIANTS.json`은 n,n_f,Ynf, 두 chemical fit, source/integrated Euler mu를 전부 보존한다. 두 fit의 차이는 구현 regression이 아니라 실제 code-paper provenance 차이다.

## 완료 범위와 재현

Derived: 고정 source anchor와 물리 outer anchor에 대한 적분식·일차 오차전파. Literature-supported: analytic consistency route와 HFB–ETFSI 근사 차이. Implementation-verified/numerically checked: 새 native helper의 독립 적분·미분·support·closure 검산. Unresolved: 원저자 lower/inner joining 세부, microscopic source-equivalence 및 composition stationarity. Bounded-not-found: 공개 전체 BSk22 ETFSI solver 저장소. Public source-prescription export의 진행은 이 마지막 항목을 기다릴 필요가 없다.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 fortran_eos/tests/crust_attachment_checks.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 fortran_eos/tests/crust_attachment_fit_variants.py
```

결과는 `data/research_blockers_20261116/lower/{TEST_RECEIPT.json,reference_integral_rows.csv,CHEMICAL_FIT_VARIANTS.json,test.log}`에 있다. 구현 파일의 frozen SHA256은 `c4dd6f0762fad107997f31c0ae3b1ea39c1cb10ed6f2b80813c13ba38b8955c2`다. 원래 lower kernel SHA는 receipt에 따로 기록되어 있으며 변경되지 않았다.
