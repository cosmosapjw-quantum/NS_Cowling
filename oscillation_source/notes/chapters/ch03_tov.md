# 상태방정식에서 평형 별로

EOS가 준비되었다. 그런데 같은 물질로 만든 별이라도 중심 상태가 다르면 질량과 반경은 달라진다. 우리는 하나의 중심 상태를 선택하고, 차갑고 비회전인 구대칭 완전유체를 정해 둔 양의 압력 cutoff까지 적분한다. 이미 TOV를 작성했다면 익숙한 식이다. 이번에는 **그 코드의 중심·표면·시간 기준이 진동 계산에 무엇을 전달하는지**에 주목하자.

## 가정

- metric signature는 \((-,+,+,+)\)이고
  \[
  ds^2=-e^{2\nu(r)}c^2dt^2+e^{2\lambda(r)}dr^2+r^2d\Omega^2
  \]
  를 쓴다.
- \(m(r)\)은 kg 단위의 중력질량, \(\epsilon\)은 정지에너지를 포함한 J m\(^{-3}\), \(P\)는 Pa이다.
- EOS는 한 압력에서 \(\epsilon,n,\mu,\chi=dP/d\epsilon\)를 같은 열역학 convention으로 반환한다.
- 첫 예제용 별은 매끄러운 단일 barotrope이며, sharp material interface는 이 장의 적분구간 안에 숨겨 넣지 않는다.
- EOS가 \(P=0\)까지 정의되지 않으면 surface를 외삽하지 않고 양의 \(P_{\rm cut}\)에서 문제를 끝낸다.

## 평형식의 유도와 SI 형태

정적 구대칭 metric에서 질량함수는
\[
e^{-2\lambda}=1-\frac{2Gm}{rc^2}
\]
로 정의한다. Einstein 방정식의 시간-시간 성분은 구 내부의 총 에너지가 질량함수를 어떻게 증가시키는지 주어
\[
\boxed{\frac{dm}{dr}=\frac{4\pi r^2\epsilon}{c^2}}
\]
가 된다. 따라서 \(dm/dr\)의 차원은 kg m\(^{-1}\)이다. 여기서 \(\epsilon/c^2\)가 들어가야 하며, 정지질량밀도와 총에너지 대응 질량밀도를 혼동하면 안 된다.

정적 유체의 운동량 보존 \(\nabla_\mu T^{\mu r}=0\)은
\[
\frac{dP}{dr}=-(\epsilon+P)\frac{d\nu}{dr}
\]
를 준다. 남은 Einstein 방정식으로 lapse gradient를 없애면
\[
\boxed{\frac{d\nu}{dr}=\frac{G\left(m+4\pi r^3P/c^2\right)}{c^2r^2\left(1-2Gm/(rc^2)\right)}}
\]
이고, 따라서
\[
\boxed{\frac{dP}{dr}=-(\epsilon+P)
\frac{G\left(m+4\pi r^3P/c^2\right)}{c^2r^2\left(1-2Gm/(rc^2)\right)}}.
\]
이 세 식이 이 패키지의 SI TOV system이다. \(d\nu/dr\)의 차원은 m\(^{-1}\), \(dP/dr\)의 차원은 Pa m\(^{-1}\)이다.

## 중심 정칙해

보통 RHS에는 \(r^{-2}\)가 있으므로 \(r=0\)에 그대로 대입하지 않는다. 정칙 중심에서는 \(m=O(r^3)\), \(P-P_c=O(r^2)\)이다. \(\epsilon_c=\epsilon(P_c)\)라 두고 leading order를 대입하면
\[
\boxed{m(r)=\frac{4\pi\epsilon_c}{3c^2}r^3+O(r^5)},
\]
\[
\boxed{P(r)=P_c-\alpha r^2+O(r^4)},\qquad
\alpha=\frac{2\pi G}{3c^4}(\epsilon_c+P_c)(\epsilon_c+3P_c).
\]
따라서 실제 적분은 작은 양의 \(r_0\)에서 이 값을 초기조건으로 사용한다. \(r_0\)를 줄이는 것은 center-series 근사의 민감도를 보는 시험이지, 표면 cutoff를 낮추는 것과 같은 수치 refinement가 아니다.

## 중심 급수로 첫 적분 구간을 줄인다

중심 오차를 줄이려면 시작점을 무조건 원점에 가깝게 잡아야 할까? 꼭 그렇지는 않다. 정칙해의 다음 계수를 알면 $1/r$ 항을 포함한 식을 아주 작은 반경부터 반복해 적분하는 부담을 줄일 수 있다. 중심을 포함한 하나의 매끄러운 EOS 구간에서
\[
m=m_3r^3+m_5r^5+O(r^7),\qquad
P=P_c+P_2r^2+P_4r^4+O(r^6)
\]
로 둔다. $\epsilon_c=\epsilon(P_c)$, $\chi_c=(dP/d\epsilon)_c$이며 $\epsilon_2=P_2/\chi_c$다. TOV 우변의 같은 차수끼리 비교하면
\[
\begin{aligned}
m_3&=\frac{4\pi\epsilon_c}{3c^2},&
g_1&=\frac{4\pi G}{c^4}\left(\frac{\epsilon_c}{3}+P_c\right),\\
P_2&=-\frac{(\epsilon_c+P_c)g_1}{2},&
m_5&=\frac{4\pi\epsilon_2}{5c^2},\\
g_3&=\frac{G}{c^2}\left[m_5+\frac{4\pi P_2}{c^2}
+\left(m_3+\frac{4\pi P_c}{c^2}\right)\frac{2Gm_3}{c^2}\right],\\
P_4&=-\frac{(\epsilon_c+P_c)g_3+(1+\chi_c^{-1})P_2g_1}{4}
\end{aligned}
\]
를 얻는다. 여기서 $\nu'=g_1r+g_3r^3+O(r^5)$다. 새로운 물질항을 추가한 것이 아니다. 같은 TOV 식의 국소해를 한 항 더 알고 시작하는 것이다.

계수의 평가는 제공 함수 `tov_center_second`에 맡긴다. 대신 시작점이 같은 EOS 구간 안에 있고 상대 보정이 작은지 확인해야 한다. 도함수가 불연속인 knot를 하나의 Taylor 급수로 건너갈 수는 없다. 이 조건은 급수를 사용하는 국소 기준이지 생략항의 엄밀한 상한은 아니다. 제5장의 중심 변위에도 같은 판단을 적용한다.


## 압력 좌표와 유효 표면

매끄러운 별에서 압력이 단조 감소하므로
\[
q=\ln\frac{P_s}{P},\qquad P_s=P(r_0)
\]
를 독립변수로 쓸 수 있다. 연쇄법칙으로
\[
\frac{dr}{dq}=-\frac{P}{dP/dr},\qquad
\frac{dm}{dq}=\frac{dm}{dr}\frac{dr}{dq},\qquad
\frac{d\nu}{dq}=\frac{P}{\epsilon+P}.
\]
이는 TOV의 물리식을 바꾸는 것이 아니라 이미 정의한 radial equations를 다른 좌표로 적분하는 adapter다. 선언한 \(P_{\rm cut}>0\)이면
\[
q_s=\ln P_s-\ln P_{\rm cut}
\]
가 끝 좌표이고, 그 위치의 areal radius를 \(R_s\)라 부른다. 마지막으로 저장된 양압 점을 임의로 surface로 선택하거나, domain 아래로 EOS를 extrapolate하지 않는다.

## 시간 정규화와 열역학적 불변량

\(\nu\)를 적분하고 나면 상수가 남는다. 질량과 반경만 보았다면 지나치기 쉽지만, 진동 주파수를 비교하려면 시간 기준을 정해야 한다. 이 유한 cutoff 문제에서는
\[
\nu(R_s)=\frac12\ln\left(1-\frac{2GM_s}{R_sc^2}\right)
\]
로 고정한다. 이는 같은 좌표시간에서 결과를 비교하기 위한 정규화다. 양의 압력에서 자른 유체가 진공 Schwarzschild 해와 물리적으로 완전 접합되었다는 뜻은 아니다.

차갑고 단일 보존량으로 기술되는 barotrope에서는 \(dP=n\,d\mu\)와 \(\epsilon+P=\mu n\)이므로
\[
\frac{d\ln\mu}{dr}=\frac{1}{\epsilon+P}\frac{dP}{dr}=-\frac{d\nu}{dr},
\qquad
\boxed{\mu e^\nu=\mathrm{constant}}.
\]
이 관계는 EOS와 lapse를 함께 점검하는 독립 잔차다. 이것 하나가 \(dm/dr\)이나 모든 profile 열의 정확성까지 증명하는 것은 아니다.

## 같은 질량의 두 별

**두 선택을 비교하자.** EOS를 바꾼 뒤 중심압력을 그대로 두는 실험과 최종 질량을 다시 맞추는 실험은 무엇이 다른가? 아직 주파수를 계산하지 말고, 두 실험에서 바깥 근 찾기가 필요한 쪽과 그 잔차를 먼저 설명한다.

EOS가 달라지면 같은 \(P_c\)가 같은 질량을 주지 않는다. 따라서 비교 연구에서는 각 EOS에 대해 지정한 안정한 branch 구간에서 \(P_c\)를 scan하고, 그 안에서
\[
F(P_c)=\frac{M(P_c)-M_*}{M_*}
\]
의 root를 찾는다. 마지막에는 반드시 별을 다시 적분해 \(|M-M_*|\)가 선언한 질량 허용오차를 만족하는지 확인한다. 좁은 \(P_c\) bracket만으로 same-mass가 보장되는 것은 아니다.

## 유도 확인

1. \(dm/dr\)에 \(\epsilon/c^2\)가 들어가는 이유를 단위와 질량함수 정의에서 설명할 수 있는가?
2. \(P'=-(\epsilon+P)\nu'\)와 lapse 식을 결합해 TOV 압력식을 다시 얻을 수 있는가?
3. \(m=O(r^3)\), \(P-P_c=O(r^2)\)를 가정해 중심 leading coefficients를 직접 얻을 수 있는가?
4. \(q=\ln(P_s/P)\)의 연쇄법칙에서 \(dr/dq>0\)가 되는 이유를 설명할 수 있는가?

## 코드 맵

- 직접 구현 경로: `learner/src/structure/ns_tov.f90`의 SI point RHS, center series, log-pressure adapter, 적분 연결.
- 직접 구현 경로: `learner/src/structure/ns_target_mass.f90`의 same-mass residual과 bracketed solve.
- 제공 기반: EOS record, 상태/오류 record, adaptive integrator, profile I/O와 plotting helper.
- compact reference: `notes/equations/tov.md`는 이 장의 수식·단위·failure contract를 빠르게 대조하는 카드다.

## 독립 검산

- 균일 총에너지밀도 Schwarzschild interior solution의 \(M,R\)를 형식적 정적 기준으로 사용한다. 이것은 인과적 EOS나 진동 입력이 아니다.
- 약한 중력 \(\Gamma=2\) polytrope에서는 Lane-Emden scale을 비교한다.
- \(r_0\), ODE tolerance, \(P_{\rm cut}\)을 **서로 다른 실험**으로 변화시킨다.
- \(\mu e^\nu\)의 공간적 일정성을 확인하되 이것을 다른 식의 검증으로 대체하지 않는다.

## 해석 질문

- 왜 \(r_0\)를 줄이는 것과 \(P_{\rm cut}\)을 낮추는 것은 같은 convergence test가 아닌가?
- \(\nu\)의 additive constant가 \(M,R\)에는 보이지 않지만 진동 주파수에는 중요한 이유는 무엇인가?
- same-mass root가 여러 branch를 가로지르면 왜 물리 비교가 모호해지는가?
- \(\mu e^\nu\)가 잘 보존되는데도 질량 profile이 틀릴 수 있는 반례를 어떻게 만들 수 있는가?

## 범위 한계

이 장은 cold, nonrotating, spherical perfect-fluid equilibrium과 finite positive-pressure cutoff를 다룬다. 정확한 진공 표면, 방사형 안정성의 전역 증명, rotation, magnetic field, elastic crust, sharp-interface restart는 여기서 구현하지 않는다. 다음 진동 장에서는 이 고정 배경을 입력으로 사용하되, 배경이 존재한다는 사실과 고유주파수가 검증됐다는 사실을 구분한다.
