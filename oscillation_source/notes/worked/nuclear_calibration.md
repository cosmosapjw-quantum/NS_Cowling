# 핵물질 포화조건과 한 점의 npe 평형

이제 함수형을 실제로 정하고 포화조건과 한 점의 조성평형을 계산하자. 아래의 baryon당 에너지는 정지에너지를 제외한 MeV이며, SI 변환에는 제공 상수를 쓴다. 입력 `binding_mev`가 뜻하는 것은 양의 결합에너지 크기가 아니라 $E_b=-16\,\mathrm{MeV}$라는 **음의 에너지**다. 폴리트로프의 $C$와 혼동하지 않도록 regulator는 $C_R$로 표기한다.

## 1. 어떤 자유도를 먼저 선택하는가

$u=n/n_0$, $n_0=0.16\,\mathrm{fm}^{-3}$라 두고

\[
e_0(u)=T_0u^{2/3}+\frac A2u+\frac{B}{\sigma+1}F(u),
\qquad F(u)=\frac{u^\sigma}{1+C_Ru^{\sigma-1}}
\]

를 선택한다. 이는 함수형의 가정이다. $A,B,T_0$의 단위는 MeV이고 $u,C_R,\sigma$는 무차원이다. 공통 nucleon 질량 $m_N=(m_n+m_p)/2$를 쓰는 대칭 핵물질의 두 Fermi sea는

\[
T_0=\frac35\frac{\hbar^2}{2m_N}
\left(\frac{3\pi^2 n_0}{2}\right)^{2/3}
\]

를 준다. 이 식에서 $n_0$는 $\mathrm{m}^{-3}$로 변환하고 J를 MeV로 변환한다. 제공 CODATA 2022 질량과 SI 정의 상수로 $T_0\simeq22.10752741\,\mathrm{MeV}$이다.

포화점에서는 $e_0(1)=E_b$, $P(n_0)=0$, $K_0=240\,\mathrm{MeV}$를 요구한다. 고정 조성에서 $P=n^2de_0/dn=n_0u^2e_0'(u)$이므로, 단위 환산을 제외한 영압력 조건은 $e_0'(1)=0$이다. $K_0=9n_0^2d^2e_0/dn^2|_{n_0}$는 $9e_0''(1)=K_0$로 바뀐다. 함수형과 $C_R$를 먼저 선택하면 남은 세 제약으로 $A,B,\sigma$를 구할 수 있다. 제약 세 개가 함수형과 regulator까지 선택해 주는 것은 아니다.

## 2. 계수 문제를 한 개의 방정식으로 줄인다

$F_1=F(1)$과 도함수를 직접 구하면

\[
F_1=\frac1{1+C_R},\quad
F'_1=\frac{\sigma+C_R}{(1+C_R)^2},\quad
F''_1=\frac{(\sigma-1)[\sigma+C_R(2-\sigma)]}{(1+C_R)^3}.
\]

따라서 세 조건은

\[
T_0+\frac A2+\frac{BF_1}{\sigma+1}=E_b,
\qquad \frac23T_0+\frac A2+\frac{BF'_1}{\sigma+1}=0,
\]
\[
-2T_0+\frac{9BF''_1}{\sigma+1}=K_0
\]

이다. 세 식을 한꺼번에 풀기 전에 공통으로 들어가는 계수를 찾아보자. 앞의 두 식을 빼면 그 계수 하나가 사라진다. $D=T_0/3-E_b$라 두면

\[
B=\frac{(\sigma+1)D}{F'_1-F_1},\qquad
A=2\left(E_b-T_0-\frac{BF_1}{\sigma+1}\right).
\]

단, $F'_1-F_1=(\sigma-1)/(1+C_R)^2$이므로 $\sigma=1$에서는 이 분해가 특이해진다. 이를 무조건 작은 분모로 나누는 구현은 허용하지 않는다. 마지막 조건에 대입하면

\[
K_0+2T_0=9D\frac{\sigma+C_R(2-\sigma)}{1+C_R}.
\]

따라서 $C_R\ne1$, $D\ne0$인 이 풀이의 범위에서

\[
\boxed{\displaystyle
\sigma=\frac{(1+C_R)(K_0+2T_0)/(9D)-2C_R}{1-C_R}}
\]

를 얻는다. 코드에서는 선언한 $\sigma$ 구간의 bracketed solve로 같은 문제를 푼다. 방금 구한 폐형식은 그 결과를 다른 연산으로 확인할 기준이다. 그런데 $C_R=1$에서는 곡률 조건이 $\sigma$를 더 이상 정하지 못한다. 단순히 분모를 작은 값으로 바꿔 실행할 문제가 아니다. 어떤 제약이 독립성을 잃었는지 돌아보아야 한다.

$C_R=0$이면 $F=u^\sigma$이므로

\[
\sigma=\frac{K_0+2T_0}{3T_0-9E_b},\quad
B=\frac{\sigma+1}{\sigma-1}\left(\frac{T_0}{3}-E_b\right),\quad
A=E_b-\frac53T_0-B
\]

로 줄어든다. 마지막 $A$ 식과 앞의 일반식은 위 포화조건 아래에서 대수적으로 같다.

## 3. 계산 가능한 worked example

동일한 $E_b=-16\,\mathrm{MeV}$, $K_0=240\,\mathrm{MeV}$에서 독립적인 60자리 scalar 계산은 다음을 준다. 표시 자릿수는 물리적 불확실도가 아니라 계산 예제의 재현을 위한 것이다.

| 입력 $C_R$ | $\sigma$ | $A$ [MeV] | $B$ [MeV] |
|:--|--:|--:|--:|
| 0 | 1.3513292382 | -209.24796675 | 156.40208773 |
| 0.2 | 1.5269938573 | -182.64138436 | 161.36320861 |

대입 후 반드시 세 residual $e_0(1)-E_b$, $n_0e_0'(1)$, $9e_0''(1)-K_0$를 **서로 다른 단위로** 기록한다. 코드 출력은 각각 `binding_residual_mev`, `pressure_residual_mev_fm3`, `k0_residual_mev`이다. 계수 하나가 표와 비슷하다는 이유만으로 포화조건을 검증하지 않는다. 유한차분으로 $K_0$를 확인할 때에는 여러 $h$를 사용하여 절단오차와 cancellation을 구분한다.

두 행이 같은 포화점을 만족해도 고밀도 물질까지 같아진 것은 아니다. 선택한 밀도 구간의 $P$, $\epsilon$, $dP/d\epsilon$를 비교해야 한다. “Regulated”라는 이름은 그 비교를 대신하지 않는다. 인과성도 실제 도함수로 확인한다.

## 4. 대칭에너지의 convention을 정확히 읽는다

채택한 코드는

\[
e(n,\alpha)=e_0(u)+S(u)\alpha^2,\qquad
\alpha=\frac{n_n-n_p}{n},
\]
\[
S(u)=S_{\mathrm{kin},0}u^{2/3}+
(S_0-S_{\mathrm{kin},0})u,\qquad
S_{\mathrm{kin},0}=(2^{2/3}-1)T_0
\]

를 사용한다. 여기서 $S_{\mathrm{kin},0}$는 두 Fermi sea의 **순수 neutron endpoint와 대칭물질 사이의 차이**에 대응한다. 이것을 $\alpha=0$에서의 정확한 Taylor 계수라고 부르면 안 된다. 실제 자유기체식

\[
T_{\mathrm{exact}}(u,\alpha)=T_0u^{2/3}
\frac{(1+\alpha)^{5/3}+(1-\alpha)^{5/3}}2
\]

의 Taylor 전개는 $T_0u^{2/3}[1+(5/9)\alpha^2+O(\alpha^4)]$다. 따라서 $5T_0/9$와 $(2^{2/3}-1)T_0$를 같은 계수로 부를 수 없다. 하나는 중심 근방의 전개이고 다른 하나는 endpoint 차이를 반영한 포물선형 근사다. `exact_kinetic_excess_mev`는 이 구별을 살피는 진단값이다. 이미 $S\alpha^2$에 kinetic 기여가 있으므로 exact kinetic을 다시 더하면 같은 기여를 중복해 넣는다.

## 5. 한 밀도에서 beta equilibrium을 구한다

공통 nucleon rest mass와 npe 조성의 가정 아래 $x=n_p/n$, $n_e=nx$, $\alpha=1-2x$이다. 고정 $n$에서 parabolic baryon 에너지를 조성으로 미분하면

\[
\mu_n-\mu_p=4S(n)(1-2x).
\]

전자 화학퍼텐셜에는 rest energy가 포함된다.

\[
\mu_e(n_e)=\sqrt{m_e^2c^4+
\hbar^2c^2(3\pi^2n_e)^{2/3}}.
\]

따라서 MeV 단위로 계산할 scalar residual은

\[
F(x)=4S(n)(1-2x)-\frac{\mu_e(nx)}{1\,\mathrm{MeV}},
\qquad 0\le x\le\frac12.
\]

여기서 분모의 $1\,\mathrm{MeV}$는 J 단위 값에서 MeV 수치로 환산하는 에너지다. $S>0$, $n>0$이면 내부 구간에서 첫 항은 감소하고 전자 항은 증가하므로 $F$는 감소한다. 양 끝이 반대 부호인 경우 유일한 내부 root가 있다. 저밀도에서 endpoint 부호 조건 자체가 실패할 수 있으므로 모든 밀도에서 내부 root가 있다고 가정하지 않는다.

$n=n_0$, $S_0=30\,\mathrm{MeV}$를 대입한 독립 고정밀 계산은

\[
x_p\simeq0.03755722652,\quad
\alpha\simeq0.92488554696,\quad
\mu_e\simeq110.98626564\,\mathrm{MeV}
\]

를 준다. 단계별로 전하중성, alpha 변환, beta residual을 별도 확인한다. 이 예제는 **muon을 제외한 모델의 답**이다. 현실적 중성자별 조성을 예측했다는 뜻이 아니다.

총 에너지밀도에서는 baryon rest energy를 복구하고 electron total energy를 더한다.

\[
\epsilon=n\,[m_Nc^2+e(n,\alpha)]+
\epsilon_e(nx),
\]

모든 항을 J와 $\mathrm{m}^{-3}$로 통일한다. 다음으로 밀도를 조금 바꾸어 제1법칙을 확인한다고 하자. $n\pm\delta n$에서 조성도 다시 풀어야 할까? 평형 곡선을 따라 미분한다면 그렇다. 조성을 그대로 고정하면 frozen-composition 응답을 구하는 것으로, 여기서 선택한 equilibrium response와 다른 질문이 된다.

## 6. 학생이 남길 증거와 해석의 경계

풀이를 확인한 뒤에는 계수표를 덮고, regulator가 바뀌었을 때 무엇을 다시 계산해야 하는지 설명해 보자. 제출물에는 함수형·에너지 영점·세 포화 residual·beta residual·밀도구간을 연결한다. 표의 숫자를 재현하는 것과 모델의 선택을 이해하는 것은 다르다. 이 계산만으로 many-body 핵이론이나 crust 전이점까지 유도한 것은 아니다.

관련 함수는 `ns_nuclear_toy.f90`와 `ns_beta_equilibrium.f90`이다. 본 풀이의 포화조건과 한 점의 조성 검사는 별 전체의 수렴 시험과 별개이다.
