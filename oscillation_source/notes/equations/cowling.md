# Cowling 단일-shot 수식표: 제5장 / 

## 모델 범위와 변수

이번 단계는 차갑고 비회전인 구대칭 완전유체 배경 위의 relativistic Cowling 근사에서, `ell=2` polar fluid perturbation의 **한 trial angular frequency**를 적분한다. Metric perturbation, elastic shear, superfluid relative motion, rotation, magnetic field는 포함하지 않는다. 배경의 압축응답은  polytropic barotrope의

\[
\chi=\frac{dP}{d\epsilon},\qquad 0<\chi\le1
\]

를 그대로 사용한다. 이 선택은 frozen-composition response를 자동으로 포함하지 않는다.

시간 의존성은 `exp(i omega t)`, `omega>0`이며

\[
k=\frac{\omega}{c}\quad [k]=L^{-1}.
\]

단위구면 harmonic의 정규화와 TOV lapse convention은 최종 spec을 따른다. 변위는

\[
\xi^r=e^{-\lambda}\frac{W(r)}{r^2}Y_{\ell m},\qquad
\xi^a=-\frac{V(r)}{r^2}D^aY_{\ell m}.
\]

각좌표가 무차원이므로

\[
[W]=L^3,\qquad [V]=L^2.
\]

## radial Cowling RHS

`ell=2`, `L_ell=ell(ell+1)=6`에서

\[
Q=k^2r^2e^{\lambda-2\nu}V+\nu'W,
\]

\[
\frac{dW}{dr}=\frac{Q}{\chi}-6e^\lambda V,
\qquad
\frac{dV}{dr}=2\nu'V-\frac{e^\lambda}{r^2}W.
\]

따라서 `Q`의 차원은 `L^2`, `dW/dr`은 `L^2`, `dV/dr`은 `L`이다. `omega`에 Hz를 넣어 자동 변환하지 않는다. Hz를 사용할 때에는 제공 단위 함수로 먼저 `omega=2*pi*f`를 만든다.

`cowling_rhs_point`는 한 배경점에서 위 radial 식만 평가한다. `chi<=0`, `chi>1`, `r<=0`, `omega<=0`, 비유한 입력과 표현 불가능한 지수항은 명시적으로 거절한다. smooth 경로에서는 물리적 jump를 smoothing하거나 crossing하지 않는다.

## 중심 정칙성

중심의 regular branch는 leading order에서

\[
W=A_0r^{\ell+1}+\cdots,\qquad
V=-\frac{A_0}{\ell}r^\ell+\cdots.
\]

`ell=2`이면 `W=A0*r^3`, `V=-A0*r^2/2`이다. fixed background는 정확한 `r=0` 행이 아니라 `r0=1 m`에서 시작하므로, 이 leading regular state를 첫 점에 적용한다. 이번 단계에서는 higher-order center series를 도입하지 않는다. `r0` 수렴은 mode calibration에서 다시 확인한다.

선형성 때문에 `A0`의 크기는 고유조건을 바꾸지 않는다. 구현은 adaptive error control이 임의 출력 진폭에 끌려가지 않도록 내부적으로 `abs(A0)=1`의 방향만 적분하고, 성공한 출력에 `abs(A0)`를 복원한다. `A0<0`은 해의 전체 방향을 뒤집는다.

##  고정 배경과 수치 좌표

배경은

`fixtures/fixed_background/smooth_poly_t05/background.csv`

이며 `q=ln(P_start/P)` 좌표로 저장되어 있다. `r,P,epsilon,nu,lambda,dnu/dr,chi`는 모두 SI/무차원 물리량이다. 배경 접근자는 현재 `구간별 선형 q 보간`, 즉 **각 저장 열을 q에서 piecewise-linear하게 보간**한다. 이 선택은 계산 representation의 일부이며 정확한 EOS 미분을 생성하지 않는다.

Cowling 물리식은 radial derivative 식을 유지한다. 전체 shot의 적분 좌표만 와 맞추어 q로 두고

\[
\frac{dr}{dq}=\frac{P}{(\epsilon+P)\nu'}>0,
\qquad
\frac{d(W,V)}{dq}=\frac{d(W,V)}{dr}\frac{dr}{dq}
\]

를 사용한다. 이는 low-pressure cutoff 부근의 매우 작은 `chi`에서 radial step이 불필요하게 경직되는 것을 완화하고 고정 배경의 q representation을 직접 사용하기 위한 **수치 어댑터**이다. 새로운 물리식이나 surface event를 뜻하지 않는다.

는 배경 interpolation을 바꾸어 ODE convergence를 억지로 단조롭게 만들지 않는다. ODE tolerance를 충분히 줄이면 piecewise-linear background representation의 error/noise floor에 도달할 수 있다. 이 현상은 실패를 smoothing으로 숨길 이유가 아니라 background-resolution convergence를 별도로 측정할 이유이다.

## 표면 조건과 정규화 잔차

Lagrangian pressure perturbation은 채택 변수에서

\[
\Delta P=-\frac{(\epsilon+P)e^{-\lambda}}{r^2}Q.
\]

실제 자유표면에서는 `Delta P=0`이나 / fixture는 `P_cut>0`의 유효 절단 문제이다. 따라서 마지막 저장점의 조건은 **finite-cutoff approximation에서의 residual diagnostic**이며 정확한 진공 표면 calibration으로 부르지 않는다.

유효 표면 `R_s`에서

\[
z=\left(\frac{W}{R_s^3},\frac{V}{R_s^2}\right),\qquad
b=\left(R_s\nu'_s,(kR_s)^2e^{\lambda_s-2\nu_s}\right)
\]

를 정의한다. 그러면

\[
\frac{Q(R_s)}{R_s^2}=b\cdot z
\]

이고 기본 잔차는

\[
\mathcal R(\omega)=\frac{b\cdot z}{\|b\|_2\|z\|_2}.
\]

두 norm 중 하나가 0이거나 비유한 경우에는 root처럼 취급하지 않고 실패한다. `W` 또는 `V` 한 성분이 0이라는 이유로 pole이 생기지 않는다. 양의 전체 amplitude 배율은 residual을 보존하고 음의 전체 배율은 residual의 부호만 뒤집는다. 한 frequency scan에서 임의로 해의 부호를 다시 선택하면 안 된다.

## 단일 trial의 범위

단일 trial에서 확인할 것은

1. 고정 배경을 명시된 representation으로 읽고,
2. 한 trial `omega`에서 regular center state를 시작하며,
3. Cowling IVP를 cutoff까지 적분하고,
4. signed normalized surface residual을 반환하는 것

이다.

단일 trial만으로는 spectrum scan, bracket, root refinement와 모드 식별을 완료할 수 없다. `mode_calibrated=false`, `identity_status=not_attempted`는 의도한 상태이다. 검증한 scan/bracket 구조를 실제 Cowling residual에 붙이는 일은 범위다.
