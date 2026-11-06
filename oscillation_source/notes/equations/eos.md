# 차가운 단일 폴리트로프의 수식표 수식과 계산 계약

이 수식표는 사용하는 지정된 barotrope를 정의한다. 미시적인 β평형 계산이나 실제 중성자별의 관측값을 제시하는 자료가 아니다. SI가 기본이며, 모든 에너지밀도에는 명시한 선형 항을 포함한다. 

## 정의와 제1법칙

$n>0$는 baryon 수밀도 $[n]=\mathrm{m}^{-3}$, $P$는 압력 $\mathrm{J\,m^{-3}}$, $\epsilon$은 총 에너지밀도, $\mu$는 baryon당 에너지 J이다. 차가운 단일 보존량의 barotrope에서

\[
P=n\frac{d\epsilon}{dn}-\epsilon,
\qquad \frac{d}{dn}\left(\frac{\epsilon}{n}\right)=\frac{P}{n^2}.
\]

$P=Kn^\Gamma$, $K>0$, $\Gamma>1$를 대입하고 적분하면

\[
\epsilon=Cn+\frac{P}{\Gamma-1},\qquad
\mu=\frac{d\epsilon}{dn}=C+\frac{\Gamma}{\Gamma-1}\frac{P}{n},
\qquad \epsilon+P=n\mu.
\]

$[C]=\mathrm{J}$이고 $[K]=\mathrm{J\,m}^{3\Gamma-3}$이다. $C$에는 선택한 baryon 정지에너지와 선형 에너지 기준을 포함한다. 이 과제는 $C>0$인 cold branch만 지원한다. $C$를 변경하는 것은 중력원이 될 $\epsilon(P)$와 $\mu$를 바꾸므로 단순한 출력 원점 변경이 아니다. $C$의 예제값을 측정된 nucleon 질량이나 원자질량 상수라고 부르지 않는다.

## 기본 저장 표현

수치적으로는 $K$와 $n^\Gamma$를 따로 형성하지 않고 양의 기준점 $(n_a,P_a)$를 사용한다.

\[
P(n)=P_a(n/n_a)^\Gamma,\qquad K=P_a/n_a^\Gamma.
\]

이는 동일한 함수의 표현이다. $n^\Gamma$나 $K$가 중간에 overflow/underflow할 때 기준점 표현이 유리하다. `init_polytrope_k`는 유한한 SI $K$를 받는 보조 생성자이며 로그에서 기준 압력을 계산한다. 정의되지 않은 $K=0$이나 계산 결과 자체의 비유한 값을 보정하지 않는다.

## 압축률과 분석적 역함수

\[
\chi=\frac{dP}{d\epsilon}=\frac{\Gamma P}{\epsilon+P},\qquad c_s^2=c^2\chi,
\qquad n(P)=n_a(P/P_a)^{1/\Gamma}.
\]

$q_a=\Gamma P_a/[(\Gamma-1)n_a]$를 에너지 J로 정의하면

\[
n(\mu)=n_a\left[\frac{\mu-C}{q_a}\right]^{1/(\Gamma-1)},\qquad
\frac{dn}{d\mu}=\frac{n^2}{\Gamma P}=\frac{n}{\mu\chi}.
\]

역함수는 root solver가 아니라 위 식으로 평가한다. $n^2$의 중간 overflow를 피하도록 마지막 도함수를 계산한다. 단위는 $\mathrm{m^{-3}\,J^{-1}}$이다. 반환 state는 복원된 $n$에서 동일한 forward 함수로 평가하며 독립적인 $P,\epsilon,\chi$ 보간은 사용하지 않는다.

## 유효구간

생성자에 $0<n_{\min}<n_{\max}$, $n_a\in[n_{\min},n_{\max}]$를 명시한다. 모든 query는 닫힌 구간 안이어야 한다. $n=0$, $P=0$ 또는 $\mu=C$를 이 과제에서 진공 상태로 연장하지 않는다. 이후 TOV의 표면 처방은 별도 작업이다.

$C>0$일 때 $q=\Gamma P/[(\Gamma-1)Cn]>0$라 놓으면

\[
\chi=(\Gamma-1)\frac{q}{1+q}
\]

이며 $q$와 $\chi$는 $n$에 따라 증가한다. 따라서 양의 압축률과 **전체 명시 구간**의 인과성은 상한의 $\chi\le1$로 확인할 수 있다. $1<\Gamma\le2$의 수학적 고밀도 극한은 $\chi\to\Gamma-1\le1$이다. $\Gamma>2$는 제한된 밀도 구간만 허용될 수 있다. 코드는 초광속 구간을 자르거나 $\chi$를 1로 clipping하지 않는다. 부적합 구간은 초기화 자체를 거절한다.

이는 해석적인 단일 폴리트로프에만 적용되는 단조성 논증이다. 일반 표 EOS의 인과성을 끝점 두 개로 검사해도 된다는 뜻이 아니다. 부동소수점의 경계 판정은 interval arithmetic에 의한 엄밀 인증이 아니다.

## 역함수의 조건수와 표현 한계

\[
\left|\frac{d\ln n}{d\ln\mu}\right|
=\frac{\mu}{(\Gamma-1)(\mu-C)}.
\]

$\mu-C$가 작으면 상대오차가 증폭된다. 실제 압축응답과 별개인 수치적 조건화 문제이다. 구간의 두 $\mu$가 같은 real64 값으로 붕괴하거나 하한이 $C$와 구별되지 않으면 `stagnation`으로 초기화를 거절한다. 모든 query가 동일 정밀도의 $n$을 복원한다고 보장하지 않는다.

먼저 $P$ 또는 $\mu$ 입력이 저장된 endpoint 범위 안에 있는지 **엄격하게** 검사한다. 한 ULP라도 입력이 밖이면 `out_of_domain`이다. 안쪽의 정당한 입력을 역변환할 때만 반올림으로 넘어간 $n$을 32 machine epsilon 상대폭 안에서 endpoint로 되돌릴 수 있다. 그보다 큰 오차는 `stagnation`이다. 이것은 자료 밖 외삽이나 물리적 domain 확장이 아니다.

## 상태와 호출 계약

`eos_state_t`는 `pressure_pa, energy_j_m3, n_b_m3, mu_j, chi, side, valid`를 갖는다. 생성되지 않은 모델, 실패한 query와 도함수 범위 초과에서 `valid=false`이고 값은 NaN으로 지운다. 외부 TSV에는 실패한 상태의 물리값을 내보내지 않는다. 0이라는 정상 물리값과 실패 상태를 구별한다.

`SIDE_REGULAR=0, SIDE_INNER=1, SIDE_OUTER=-1`은 일반 조건이다. 이 smooth 단일 branch에서 허용된 세 label은 같은 값이며 압력 query가 요청한 label을 보존한다. 다른 정수는 거절한다. Jump의 양쪽 EOS를 결정하는 기능은 /T10에 속한다.

- `init_polytrope(model, n_anchor, p_anchor, gamma, c_j, n_min, n_max, status)`
- `init_polytrope_k(model, k_si, gamma, c_j, n_anchor, n_min, n_max, status)`
- `eos_at_n(model, n_b_m3, state, status)`
- `eos_at_pressure(model, p_pa, side, state, status)`
- `eos_at_mu(model, mu_j, state, dn_dmu, status)`
- `polytrope_domain(model, lower, upper, status)`
- `polytrope_ready(model)`

Model은 explicit read-only 인수이다. 상위 TOV callback은 caller의 생존 범위에서 이 model을 bind할 예정이며 전역 mutable EOS를 공유하지 않는다. 두 model을 번갈아 평가해도 상태가 섞이지 않아야 한다.
