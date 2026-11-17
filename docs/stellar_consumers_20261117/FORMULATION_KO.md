# D: EOS·TOV·Cowling의 물리 계약

## 가정과 계산 대상

기준은 C commit 539e6bd8de73f8654d121c0bc30023fb868dc05b의 B primary,
data/core_support_20261117/primary_v1/REPRESENTATION.json이다.
surface, outer polytrope, lower bridge, upper P(mu)의 실제 연속식을 그대로 평가한다.
BSk22 lower segment·SPG·SLy4와 명시한 surface/bridge로 구성된 reference-assisted hybrid이며,
native unified crust와 구분한다.

별은 정적·구대칭·비회전·비자성 완전유체다.
배경과 섭동 모두 beta 평형 응답을 선택한다: Gamma_1=Gamma_eq.
전단 탄성·초유체·frozen-species 부력은 이번 계산의 자유도가 아니다.
APR 소비자 검증과 matching은 별도 경로이며 SLy4 whole-star 결과를 APR 별로 해석하지 않는다.

## TOV, enthalpy 및 단위

계량과 물리적 식은
\[
ds^2=-e^{2\Phi}c^2dt^2+(1-2GM/c^2r)^{-1}dr^2+r^2d\Omega^2,
\]
\[
M'=4\pi r^2\epsilon/c^2,\quad
P'=-(\epsilon+P)\Phi',\quad
\Phi'=\frac{G}{c^2}\frac{M+4\pi r^3P/c^2}{r(r-2GM/c^2)}
\]
이다. epsilon은 rest energy를 포함한다.
원전: [Oppenheimer–Volkoff 1939](https://doi.org/10.1103/PhysRev.55.374),
[Tolman 1939](https://doi.org/10.1103/PhysRev.55.364).

내부 길이 단위 km에서 \(m=GM/c^2,\ \bar\epsilon=G\epsilon/c^4,\ \bar P=GP/c^4\)를 쓴다.
뒤의 두 양은 \(\mathrm{km}^{-2}\) 단위다.
\(G=6.67430\,10^{-11}\,\mathrm{m^3kg^{-1}s^{-2}}\), \(c=299792458\,\mathrm{m/s}\),
\(1\,\mathrm{MeV\,fm^{-3}}=1.602176634\,10^{32}\,\mathrm{Pa}\)로 변환한다.
태양질량은 native와 같은 nominal \((GM)_\odot^N=1.3271244\,10^{20}\,\mathrm{m^3s^{-2}}\)에 묶인다.
이는 실제 태양질량의 무오차 관측값이 아니라 변환 convention이다.
[IAU 2015 B3](https://arxiv.org/abs/1510.07674)

Cold first law \(d\epsilon=\mu_Bdn_B\), Euler 관계 \(P=n_B\mu_B-\epsilon\)에서
\(dP=n_Bd\mu_B\)이고,
\[
h=\int_0^P\frac{dP}{\epsilon+P}=\ln(\mu_B/\mu_{\rm vac}),\qquad
\frac{dh}{dr}=-\Phi'
\]
를 얻는다. 따라서
\[
\frac{dr}{dh}=-\frac{r(r-2m)}{m+4\pi r^3\bar P},\qquad
\frac{dm}{dh}=4\pi r^2\bar\epsilon\,\frac{dr}{dh}.
\]
Upper Hermite \(P(\mu_B)\)에서는 \(n_B=P'\), \(\epsilon=\mu_BP'-P\)를 같은 다항식에서 얻는다.
별 적분을 위해 원 EOS의 knot를 PAV 재투영하거나 energy offset을 추가하지 않는다.

Regular center는
\[
m(r)=\frac{4\pi}{3}\bar\epsilon_c r^3+O(r^5),\qquad
h_c-h(r)=2\pi(\bar\epsilon_c/3+\bar P_c)r^2+O(r^4)
\]
이다. 이 식으로 유한 \(r_0\)에서 출발하고 \(r_0\) 변화로 생략항을 검사한다.
배경은 analytic vacuum \(h=0\)까지 적분한다.
유한 표면 절단 비교에서는 해당 절단 모델의 \(m_s,R_s,h_s\)로
\[
\Phi(h)=\tfrac12\ln(1-2m_s/R_s)-h+h_s
\]
를 정규화한다. 원래 vacuum 모델의 시간 정규화와 혼용하지 않는다.

검산용 \(u=r^2,\ v=m/r\) 식
\[
u_h=-\frac{2u(1-2v)}{v+4\pi u\bar P},\qquad
v_h=-\frac{(1-2v)(4\pi u\bar\epsilon-v)}{v+4\pi u\bar P}
\]
은 원래 식의 변수변환이다. CAS로 두 잔차가 정확히 0임을 확인했다.

B의 실제 지원은 \(0\le n_B\le1.0\,\mathrm{fm}^{-3}\)다.
지원 밖 입력은 거절하고 extrapolation하지 않는다.
표본 수열의 양의 \(dM/dn_c\)는 최대질량이나 모든 종류의 안정성을 증명하지 않는다.
최대질량에는 turning point를 포함하는 추가 native EOS 지원이 필요하다.

## Cowling 방정식과 정칙성

Cowling은 \(\delta g_{\mu\nu}=0\)을 선택하므로 방사 감쇠시간을 계산하지 않는다.
[McDermott, Van Horn, Scholl 1983](https://doi.org/10.1086/161006)은 이 근사와
정칙 중심·표면의 Lagrangian pressure 조건을 준다.
평형과 frozen response 구분 및 식 비교는
[Montefusco, Antonelli, Gulminelli 2025](https://arxiv.org/html/2410.08008v2)를 이용했다.

\[
\xi^r=e^{-\Lambda}W Y_{\ell m}e^{i\omega t}/r^2,\quad
\xi^\theta=-V\partial_\theta Y_{\ell m}e^{i\omega t}/r^2,\quad
\xi^\phi=-V\partial_\phi Y_{\ell m}e^{i\omega t}/(r^2\sin^2\theta).
\]
\(e^{2\Lambda}=(1-2m/r)^{-1}\), \(k=\omega/c\), \(f=\omega/(2\pi)\)다.
\(\omega\)는 무한대의 시간에 대한 angular frequency다.
\(\chi_{\rm eq}=dP/d\epsilon=c_{\rm eq}^2/c^2\),
\(\Gamma_{\rm eq}=(\epsilon+P)\chi_{\rm eq}/P\)이면 선택한 평형 식은
\[
W'=\chi_{\rm eq}^{-1}[k^2r^2e^{\Lambda-2\Phi}V+\Phi'W]
-\ell(\ell+1)e^\Lambda V,\qquad
V'=2\Phi'V-e^\Lambda W/r^2.
\]
\(r\)과 \(k^{-1}\)의 길이 단위는 일치한다.

\(W=Ar^{\ell+1}\), \(V=Br^\ell\)를 직접 대입하면 \(A=-\ell B\)다.
따라서 정칙 중심은
\[
W=Ar^{\ell+1}+O(r^{\ell+3}),\qquad
V=-Ar^\ell/\ell+O(r^{\ell+2}).
\]
위 2025 논문의 중심 근방 prose에 있는 동일 계수 표기를 복사하지 않고,
ODE 자체를 만족하는 부호와 계수를 사용했다.

Scaled 변수 \(w=W/r^{\ell+1}\), \(v=V/r^\ell\)는
\[
w'=[\Phi'/\chi_{\rm eq}-(\ell+1)/r]w+
[k^2re^{\Lambda-2\Phi}/\chi_{\rm eq}-\ell(\ell+1)e^\Lambda/r]v,
\]
\[
v'=-e^\Lambda w/r+(2\Phi'-\ell/r)v
\]
를 만족한다. 중심의 \(1/r\) 행렬 고유값은 \(0,-(2\ell+1)\)이며
정칙 방향은 \((1,-1/\ell)\)다. 원래 식과의 동일성은 CAS 및 수치 식 비교로 확인한다.

표면 \(\Delta P=0\)은
\[
Q=k^2r^2e^{\Lambda-2\Phi}V+\Phi'W=0
\]
이다. 유한 \(h_{\rm cut}>0\)에서 이 조건으로 outer solution을 안쪽으로 적분한다.
중심에서 출발한 solution과 내부 지점에서 normalized determinant를 0으로 맞춘다.
표면에서 부과한 조건의 작은 잔차만으로 고유해를 인증하지 않는다.
root bracket, matching 잔차, 결합 불일치, node 수, 중심·표면·오차·matching 지점 수렴을 함께 기록한다.

## 평형 closure의 물리적 의미

일반적인 perturbative 음속 \(\chi_{\rm ad}\)에 대해
\[
\mathcal A=\epsilon'/(\epsilon+P)-P'/(\Gamma_1P)
=-\Phi'(\chi_{\rm eq}^{-1}-\chi_{\rm ad}^{-1}),
\]
\[
N^2=-c^2e^{2(\Phi-\Lambda)}\Phi'\mathcal A
=c^2e^{2(\Phi-\Lambda)}(\Phi')^2
(\chi_{\rm eq}^{-1}-\chi_{\rm ad}^{-1}).
\]
선택한 \(\chi_{\rm ad}=\chi_{\rm eq}\)에서는 \(N^2=0\)이다.
이는 실제 frozen-composition 별의 부력이 0이라는 결론이 아니다.
그 예측에는 fixed-species off-equilibrium derivatives와 일관된 전체 섭동계가 필요하다.
평형 경로의 조성값만으로 그런 도함수를 복원하지 않는다.
문헌의 다른 모델에 대한 fractional frequency difference를 이 hybrid EOS의 오차로 대입하지 않는다.

## APR continuous Gibbs consumer

[APR provider sheet 2025-03-14](https://compose.obspm.fr/download/1D/NS/Classical/APR/eos.pdf)는
고밀도 pion-condensed phase에 bulk Gibbs 처방을 명시한다.
온도·조성 제약이 다른 APR 구현의 Maxwell 처방과 구분한다.

실제 generated table의 retained \(6\rightarrow8\rightarrow7\) endpoint에 묶인
phase별 cubic Hermite \(\epsilon(n)\)에서
\[
\mu_B=\epsilon'(n),\quad P=n\epsilon'(n)-\epsilon(n),\quad
\chi=n\epsilon''(n)/\epsilon'(n)
\]
를 얻는다. 경계의 \(n,\epsilon,\mu_B,P\)는 연속이고
\(\epsilon'',\chi\)는 양쪽 값을 유지한다.
Cubic interval에서는 \(\epsilon''\)가 선형이고
\((\mu_B-n\epsilon'')'=-n\epsilon'''\)의 부호가 고정되므로 endpoint 검사로
positive compressibility와 causal inequality를 확인할 수 있다.
이는 binary64 polynomial 검사이며 outward-rounded interval proof는 아니다.

이 경계에는 density jump나 delta-function surface layer가 없다.
선택한 평형 완전유체에서는 \(W,V\), radial displacement 및 Lagrangian pressure를
연속으로 전달하고 계수의 one-sided derivative만 달리 적용한다.
Sharp Maxwell의 slow/rapid phase-conversion transfer를 삽입하지 않는다.

25개 APR matching은 low homogeneous branch에서 native matcher를 호출한 다음
source-bound mixed/high-density suffix의 물리량과 composition mask를 보존한다.
native low-branch admission과 assembled output 검증을 구분한다.
CSS 접합의 별도 sharp law를 내부 APR Gibbs 경계에 적용하지 않는다.

## 증거 연결

실제 numerical acceptance는 main RESULTS_KO.md와 lane별 결과를 따른다.
SOURCE_LEDGER.json은 원전의 지지 범위를, SYMBOLIC_CHECK.json은 실제 CAS 입력·출력 및
수정한 symbolic 치환 오류를 보존한다. CAS 통과는 실행 오차나 source 정확성을 대체하지 않는다.
Native table 경로와 exact-contract reference의 차이는 별도 보간 비교 receipt에 남긴다.
