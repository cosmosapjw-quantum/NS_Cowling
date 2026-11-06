# SI TOV 수식표와 양의 압력 표면

기준은 `notes/chapters/ch03_tov.md`이다. Signature는 (−,+,+,+),
\(ds^2=-e^{2\nu}c^2dt^2+e^{2\lambda}dr^2+r^2d\Omega^2\)이다. \(m\)은 kg 단위의 중력질량이며
\(\epsilon\)은 정지에너지를 포함한 J m⁻³이다. \(P\)는 Pa, \(r\)는 m이다.

\[
m'=4\pi r^2\epsilon/c^2,\quad
\nu'=\frac{G(m+4\pi r^3P/c^2)}{c^2 r^2(1-2Gm/(rc^2))},\quad
P'=-(\epsilon+P)\nu'.
\]

`dm_dr`의 차원은 kg m⁻¹, `dp_dr`는 Pa m⁻¹, `dnu_dr`는 m⁻¹이다. \(r=0\), \(P\le0\),
\(m<0\) 또는 \(1-2Gm/(rc^2)\le0\)에서 보통의 RHS를 호출하지 않는다. EOS callback은 질의한 압력과 같은 상태를 반환해야 한다.

## 중심 시작값

\(\epsilon_c=\epsilon(P_c)\)를 고정한 급수의 선도항으로
\[
m(r_0)=\frac{4\pi\epsilon_c}{3c^2}r_0^3+O(r_0^5),\qquad
P(r_0)=P_c-\alpha r_0^2+O(r_0^4),
\]
\[
\alpha=\frac{2\pi G}{3c^4}(\epsilon_c+P_c)(\epsilon_c+3P_c).
\]
\(\nu(r_0)\)의 additive constant는 자유로 두고, 마지막에 시간 기준을 정한다. 이 시작값은 원점의 정확한 해가 아니므로
\(r_0\)를 줄이는 시험을 독립적으로 수행한다. 시작 압력 보정이 중앙 압력의 1% 미만이고 실수형에서 구별되는 구간만
현재 시작값의 admitted 범위로 사용한다. 이 1%는 유효성 guard이지 전체 정확도 보증이 아니다.

## 명시적인 수치 표현 변경

EOS는 \(P>0\)의 유한 폐구간만 제공한다. 표면을 찾으려고 그 아래를 외삽하거나 음수를 clipping하지 않는다.
이번 smooth 적분은 \(P_s=P(r_0)\), \(q=\ln(P_s/P)\)를 독립변수로 사용한다. 이는 수식표 radial TOV를 바꾸는 것이 아니라 연쇄법칙에 따른 adapter이다.

\[
\frac{dr}{dq}=-\frac{P}{P'},\quad
\frac{dm}{dq}=m'\frac{dr}{dq},\quad
\frac{d\nu}{dq}=\frac{P}{\epsilon+P}.
\]

초기 상태는 \((r_0,m(r_0),\nu(r_0))\)이며 \(q=0\)에서 시작한다.
\(q_s=\ln P_s-\ln P_{\rm cut}\)가 유효 표면의 정확한 좌표이다. 끝점 압력은 이 좌표의 정의로 결정되며,
마지막 양압 grid cell의 반경을 표면으로 선택하지 않는다. 결과로 계산한 \(r(q_s)\)가 표면 반경이다.
Stage의 물리량과 출력은 여전히 SI이고, \(dq\)만 무차원이다. 동일 adaptive routine을 사용하며, controller를 다시 구현하지 않는다.

`h_initial`, `h_min`, `h_max`는 이 호출에서 **dq**를 뜻한다. 오류 성분은 `(r[m],m[kg],nu[1])`이다.
중심의 \(dr/dq\sim1/r\)를 해상하기 위해 초기 dq를 `0.1*(Pc-Ps)/Pc` 이하로 제한한다. 이후 step 증가는  controller가 담당한다.
이 경로는 압력이 단조 감소하는 smooth 문제용이다. 일반 radial event-root solver, 진공 \(P=0\), 또는 sharp phase restart의 완성이 아니다.
제시한 `eos_events`는 smooth 압력 stop만 받는다. 사용자가 불연속 callback을 숨겨 넣는 것을 자동 검출한다고 주장하지 않는다.

## 유효 cutoff와 lapse

이번 fixture는 `finite_cutoff_schwarzschild_time`을 사용한다.
\[
\nu(R_s)=\tfrac12\ln[1-2GM_s/(R_sc^2)].
\]
양의 압력에서 자른 별이 진공에 정확히 접합됐다는 의미가 아니라, 생략 외층을 가진 유효 계산의 시간 기준이다.
실제 진공 모델은 외층 연장 또는 실제 표면 극한을 따로 정해야 한다. \(P_{\rm cut}\) 변화는 ODE tolerance 변화와 달리 계산 문제의 경계를 바꾼다.
\(\nu\mapsto\nu+C\), \(\omega\mapsto e^C\omega\)가 같은 \(\omega^2e^{-2\nu}\)를 주므로 이후 진동의 시간 기준을 보존해야 한다.

## 같은 질량

입력으로 주어진 Pc 폐구간을 로그 간격으로 실제 적분해 질량이 증가하는지 검사한다.
이것은 유한 표본의 branch guard이며 방사형 안정성의 증명도 전구간 단조성 증명도 아니다.
\(x=\ln(P_c/P_{c,\min})\)에서 \(F=[M(P_c)-M_*]/M_*\)를  bisection으로 푼다.
마지막 재적분에서 \(|M-M_*|\le M_{\rm atol}+M_{\rm rtol}M_*\)를 별도로 확인한다.
로그 Pc root의 좌표 오차만 작다고 같은 질량에 도달했다고 보고하지 않는다.

## 독립 기준식

균일 총 에너지밀도 \(\epsilon_0\)의 Schwarzschild 내부해를 형식적 정적 기준문제로 사용한다.[1]
\[
s=\frac{\epsilon_0+P_c}{\epsilon_0+3P_c},\qquad
u_s=\frac{s(\epsilon_0+3P_{\rm cut})}{\epsilon_0+P_{\rm cut}},\qquad
R_s^2=\frac{3c^4}{8\pi G\epsilon_0}(1-u_s^2).
\]
여기서 \(u_s\)는 제곱근 형태의 대수적 보조량이며 metric 함수 \(\nu\)와는 다른 양이다.
\(M_s=4\pi\epsilon_0R_s^3/(3c^2)\). 이 상수밀도 모델은 인과적인 EOS나 Cowling 입력이 아니며,
시험용 record의 `chi` placeholder는 압축응답으로 해석하지 않는다. 실제 배포 배경은 causal polytrope만 사용한다.

폴리트로프에서는 열역학과 hydrostatic equilibrium을 합쳐 \(d\ln\mu=-d\nu\), 즉 \(\mu e^\nu=\mathrm{constant}\)를 얻는다.
이 항등식으로 lapse와 EOS의 관계를 별도로 검사하되, 질량 방정식의 정확성까지 이 하나로 증명하지 않는다.

[1] K. Schwarzschild (1916), *On the gravitational field of a sphere of incompressible fluid according to Einstein's theory*,
S. Antoci 번역, arXiv:physics/9912033. 기준식은 본 수식표의 SI·중심압력·cutoff 표기로 재작성했다.
