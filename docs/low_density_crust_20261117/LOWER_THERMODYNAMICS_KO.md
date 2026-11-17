# H 하부 연결: endpoint 판정과 구성적 대안

기준 G commit은 `905b845ad0f14d2f1e7f67f9eb802691144850cd`다. 이 문서는 H의 native 저밀도 연장에 사용할 열역학 판정과 G에서 거절된 DDME2 lower bridge의 원인을 정리한다. H의 새 TF 상태나 stellar 결과의 완료 보고서는 아니다.

**결론:** G의 DDME2 여섯 endpoint 조합은 고정 E outer endpoint와 **안정하고 인과적인 연속 barotrope를 구성할 수 있다**. 거절 원인은 양끝값을 모두 고정한 단일 cubic \(P(\mu)\)의 추가 제약이다. 실제 같은 endpoint를 보존하는 결정적 exponential-\(n/\mu\) 연결을 구성하여 검산했다. 이 예는 미시적 TF 재생성이 아닌 현상론적 표현이며, H 생산 contract에는 채택하지 않았다. 기존 outer의 더 이른 217개 저장 knot를 검사하는 것만으로는 이 여섯 DDME2 cubic을 수용할 수 없었다.

## 1. 정의와 근거의 범위

찬 평형 barotrope에서 \(n\)은 baryon 수밀도, \(P\)는 압력, \(\epsilon\)은 정지에너지를 포함한 에너지밀도다. 이 문서의 chemical potential은 저장 Q3를 별도로 조절한 값이 아니라

\[
\mu=\frac{\epsilon+P}{n},\qquad
d\epsilon=\mu\,dn,\qquad dP=n\,d\mu
\]

를 만족하는 baryon Euler chemical potential이다. 단위는 \(n\,[\mathrm{fm}^{-3}]\), \(P,\epsilon\,[\mathrm{MeV\,fm}^{-3}]\), \(\mu\,[\mathrm{MeV}]\)다. \(c_s^2/c^2=dP/d\epsilon\)로 쓰며 \(c=1\)을 물리 단위로 가정하지 않는다. 이 열역학 유도에는 \(\hbar\)가 직접 등장하지 않으며, native 계산의 기존 \(\hbar c\), 질량 및 결합상수 입력은 바꾸지 않았다. 시공간 부호는 프로젝트의 \((-+++)\)를 따르지만 아래 판정은 metric에 독립적이다.

- **derived:** endpoint 필요조건, cubic의 필요충분조건, exponential 연결의 존재 구성.
- **numerically checked / implementation-verified:** G의 실제 열 개 join trial, E outer 217개 knot를 사용한 2,170개 유한 조합, 실제 `HermiteInterval.validate()`와 해석식의 일치, 구성예의 독립 quadrature 및 first-law 차분.
- **literature-supported, reused:** TF의 npe 구성·중성미자 없는 beta 평형·cell neutrality와 source convention은 [G source physics](../crust_refinement_20261117/SOURCE_PHYSICS_KO.md)의 이미 읽은 Xia et al. arXiv:2208.12893v1 및 Niu et al. arXiv:2506.11492v2 근거를 재사용한다. 새로운 원전 검색이나 새 parameter calibration을 하지 않았다.
- **unresolved / outside scope:** lower bridge 내부의 핵종·neutron-drip 조성, 완전한 microscopic free-energy 최소, 모든 geometry의 global minimum, G에서 미해결인 DDME-X 원전 J/L. 이 결과는 그 문제들의 해결 증거가 아니다.

Outer 표현의 SSOT는 [E fine contract](../../data/stellar_extension_20261117/eos/fine/REPRESENTATION.json)다. 이 표현은 공개 outer의 \(n,P\)와 첫 energy anchor를 사용한 pressure-primary 적분이며, 모든 raw \(\epsilon,\mu\) knot를 그대로 보존한 원저자 microscopic EOS와 같지 않다. 이 구분과 isotope 공백의 유효 표현은 [A source audit](../whole_crust_20261117/SOURCE_AUDIT_KO.md), [A physics derivation](../whole_crust_20261117/PHYSICS_REVIEW_KO.md)에 이미 기록되어 있다. 여기서는 E 표현 자체를 고정해 재사용한다.

## 2. 일반 smooth causal connector의 정확한 endpoint 조건

양끝 상태를 \(L,R\)라 하고

\[
a=\mu_L>0,\quad b=\mu_R>a,\quad d=b-a,\quad
\Delta P=P_R-P_L>0,\quad m=\frac{\Delta P}{d}
\]

로 둔다. \(n_L,n_R,\epsilon_L,\epsilon_R>0\)이며 각 endpoint는 위 Euler 관계를 만족한다. 안정한 연속 연결은 \(n'(\mu)>0\)를, 인과적인 연결은

\[
0<\frac{c_s^2}{c^2}
 =\frac{n}{\mu n'}\le1
\]

을 만족해야 한다. 안정성만으로도 \(n_L<m<n_R\)가 필요하다. 따라서 양끝 \(n,P,\mu\)의 증가와 평균밀도 조건은 첫 검사지만 충분한 인과성 검사는 아니다.

\(y(\mu)=n(\mu)/\mu\)를 도입하면

\[
\mu n'-n=\mu^2y',\qquad
\frac{c_s^2}{c^2}=\frac{y}{y+\mu y'}.
\]

양의 \(n,\mu\) 영역에서 인과성은 \(y'\ge0\)와 동치다. 따라서

\[
y_L=\frac{n_L}{a}\le y_R=\frac{n_R}{b},
\qquad
y_L\mu\le n(\mu)\le y_R\mu
\]

이며 적분하면

\[
\boxed{\quad
\frac{n_L}{2a}(b^2-a^2)
\le \Delta P\le
\frac{n_R}{2b}(b^2-a^2).
\quad} \tag{1}
\]

\(y_L<y_R\)이고 \(n(\mu)\)가 연속이면 (1)의 두 부등식은 **strict**해야 한다. 왼쪽 또는 오른쪽 endpoint 근방에서 각각 적분차가 양수이기 때문이다. \(y_L=y_R\)일 때는 유일하게 \(n(\mu)=y_L\mu\), \(\Delta P=y_L(b^2-a^2)/2\)인 luminal 연결이 가능하다. \(\mu_L=\mu_R\)인 sharp coexistence doublet는 이 smooth-interval 정리의 대상이 아니다.

다음 절의 구성은 \(y_L<y_R\)와 (1)의 strict 부등식이 이 클래스에서 **충분조건이기도 함**을 보인다. \(\epsilon(\mu)=\mu n-P\)로 두면

\[
\frac{d\epsilon}{d\mu}=\mu n'>0.
\]

따라서 \(\epsilon_L>0\)에서 전 구간의 양의 에너지밀도가 따라오며, \(\epsilon_R\)도 고정 endpoint Euler 값과 일치한다. 이것은 source energy에 \(Cn\)을 추가하는 조절이 아니다.

차원상 \(y\)는 \(\mathrm{fm}^{-3}\mathrm{MeV}^{-1}\), (1)은 \(\mathrm{MeV\,fm}^{-3}\), \(\mu n'-n\)은 \(\mathrm{fm}^{-3}\)다. \(c_s^2/c^2\)는 무차원이다.

## 3. 단일 cubic Hermite가 요구하는 더 좁은 창

기존 `HermiteInterval`은 \(x=(\mu-a)/d\in[0,1]\)에 대해

\[
P(x)=P_L+d n_Lx
+d(3m-2n_L-n_R)x^2
+d(-2m+n_L+n_R)x^3
\]

를 사용한다. 따라서

\[
n(x)=n_L+(6m-4n_L-2n_R)x
 +(3n_L+3n_R-6m)x^2, \tag{2}
\]

\[
n'_L=\frac{6m-4n_L-2n_R}{d},\qquad
n'_R=\frac{-6m+2n_L+4n_R}{d}. \tag{3}
\]

Causal margin \(C=\mu n'-n\)의 미분은 \(dC/d\mu=\mu n''\)다. (2)의 \(n''\)는 상수이고 \(\mu>0\)이므로 \(C\)의 최솟값은 양끝 중 하나에 있다. \(C_L,C_R\ge0\)이면 (3)의 양끝 susceptibility도 양수이며, \(n'\)가 선형이므로 전 구간에서 양수다. 따라서 \(n,\epsilon>0\)도 유지된다.

결국 양의 ordered endpoint에 대해 단일 cubic의 안정성과 인과성은 정확히

\[
\boxed{\quad
\frac{4n_L+2n_R+n_Ld/a}{6}
\le m\le
\frac{2n_L+4n_R-n_Rd/b}{6}.
\quad} \tag{4}
\]

에 동치다. 순수 안정성의 lower bound \((2n_L+n_R)/3\)보다 (4)의 lower bound가 조금 더 높다. 한 endpoint가 equality인 경우 \(c_s=c\)가 허용되며 susceptibility가 0이 되는 것은 아니다.

이는 endpoint와 해당 cubic에 국한된 판정이다. (4)가 실패하고 (1)이 통과한다면 다른 causal barotrope가 가능하다. \(\mu,n\)가 양수라는 가정을 제거한 일반 polynomial 정리가 아니며, exact arithmetic 식 (4)와 구현의 명시적 floating-point roundoff 처리를 구별한다.

## 4. G DDME2 실패의 실제 숫자

E의 고정 왼쪽 endpoint는

\[
(n_L,P_L,\epsilon_L,\mu_L)
=(0.000269021,\;
0.0004981910045986,\;
0.2520862865914808,\;
938.9024559275276).
\]

G [full-contract admission](../../data/crust_refinement_20261117/stellar/FULL_CONTRACT_ADMISSION.json)의 모든 DDME2 trial을 그대로 읽었다.

| mesh | \(n_R\) | 실제 \(m\) | cubic causal lower bound | cubic 결과 |
|---|---:|---:|---:|---|
| N64 | 0.008 | 0.00267420708087 | 0.00284619722065 | 거절 |
| N64 | 0.009 | 0.00297482543632 | 0.00317954052715 | 거절 |
| N128 | 0.006 | 0.00205096597611 | 0.00217950740636 | 거절 |
| N128 | 0.007 | 0.00236106424590 | 0.00251285275855 | 거절 |
| N128 | 0.008 | 0.00266879229738 | 0.00284619707056 | 거절 |
| N128 | 0.009 | 0.00297585055478 | 0.00317954059241 | 거절 |

\(n_R=0.006\), N128의 \(\mu_R=942.2544534730021\)이고, (1)의 평균밀도 허용구간은

\[
0.0002695012190717406
<m<
0.005989327731379396.
\]

실제 \(m\)은 이 구간 안에 있지만 (4)의 lower bound보다 작다. 이 때문에

\[
n'_L=-0.00022979973371364018
\quad[\mathrm{fm}^{-3}\mathrm{MeV}^{-1}],
\qquad
C_L=-0.21602855535522864
\quad[\mathrm{fm}^{-3}]
\]

가 된다. 이는 roundoff 규모가 아니며 기존 구현이 거절해야 하는 실제 cubic 불안정성이다. 여섯 DDME2 trial 모두 같은 분류다. 반면 G의 네 TM1 control은 (4)를 통과하고 기존 구현 결과와 일치한다. 모든 열 개 trial은 일반 조건 (1)을 통과한다.

## 5. 전 endpoint를 보존하는 결정적 구성예

다음 family는 일반 조건 (1)의 strict 범위 전체를 실제로 구성한다.

\[
\Phi_k(x)=
\begin{cases}
\dfrac{e^{kx}-1}{e^k-1}, &k\ne0,\\
x, &k=0,
\end{cases}
\]

\[
y(\mu)=y_L+(y_R-y_L)\Phi_k(x),\qquad
n(\mu)=\mu y(\mu). \tag{5}
\]

모든 유한 실수 \(k\)에 대해 \(\Phi_k(0)=0,\Phi_k(1)=1\)이며

\[
\frac{d\Phi_k}{dx}
=\frac{k e^{kx}}{e^k-1}>0
\]

이다. 따라서 \(n,P,\epsilon\)가 양수인 채 \(0<c_s^2/c^2<1\)를 전 구간에서 보장한다. 압력은

\[
P(\mu)=P_L+\int_a^\mu s\,y(s)\,ds,\qquad
\epsilon(\mu)=\mu n(\mu)-P(\mu)
\]

로 정한다. 유일한 \(k\)는 \(\int_a^b s\,y(s)\,ds=\Delta P\)로 고정한다. Stellar 질량·반지름·모드 주파수는 입력되지 않는다.

\(0<x<1\)에서 \(\Phi_k(x)\)는 \(k\)의 엄격한 감소함수다. 실제로
\(\Phi_k(x)=\int_0^x e^{kt}dt/\int_0^1e^{kt}dt\)이고, \(k\)가 증가하면 양의 weight가 큰 \(t\)로 이동한다. 그 미분은 indicator \(\mathbf 1_{t\le x}\)와 \(t\)의 음의 covariance로 쓸 수 있다. 따라서 압력 적분도 연속이며 엄격하게 감소한다. \(k\to+\infty\)와 \(k\to-\infty\)에서 각각 (1)의 lower/upper bound로 접근한다. 이로써 strict 범위 전체에서 유일한 유한 root가 존재한다.

계산용 primitive는

\[
J_k(x)=\int_0^x(a+dt)\Phi_k(t)\,dt
=\frac{
a[(e^{kx}-1)/k-x]
+d[((kx-1)e^{kx}+1)/k^2-x^2/2]
}{e^k-1}
\]

이고, \(J_0(x)=ax^2/2+dx^3/3\)다. 그러면

\[
P(x)=P_L+d\{y_L(ax+dx^2/2)+(y_R-y_L)J_k(x)\}.
\]

실제 진단은 pressure root를 독립 adaptive quadrature로 구하고, 이 primitive와 다시 비교했다. 큰 양의 \(k\)에서 overflow를 피하도록 \(\Phi_k=e^{k(x-1)}(1-e^{-kx})/(1-e^{-k})\)를 사용했다. 이 유한 진단에서 쓰인 primitive의 \(k\)는 0에서 충분히 떨어져 있으며, 일반 production 구현으로 채택하려면 \(k\approx0\)의 series 및 큰 \(|k|\)의 안정적 primitive 처리를 별도로 넣어야 한다. 현재 진단 script는 범위를 assert하여 그 바깥을 조용히 사용하지 않는다.

N128 DDME2 \(n_R=0.006\)의 경우

\[
k=2.4902682205225353
\]

이며, 1,001개 표본에서 \(c_s^2/c^2=0.0007468317140\)–\(0.0013699720674\)였다. 이 숫자는 표본 extrema이고, continuum causality의 근거는 표본의 촘촘함이 아니라 위 \(y'>0\) 유도다. 열 개 G 구성예 전체에서 analytic pressure와 독립 quadrature의 최대 절대차는 \(5.5511\times10^{-17}\,\mathrm{MeV\,fm}^{-3}\), first-law 중앙차분의 최대 상대오차는 \(2.1426\times10^{-10}\)였다.

이 family는 자유로운 microscopic EOS 복원이 아니다. \(k\)는 주어진 family 안에서 endpoint가 결정하지만, family 선택 자체는 phenomenological model choice다. Lower endpoint와 upper native endpoint의 실제 sound-speed 미분을 동시에 맞추는 조건은 부과하지 않았고, 접합점의 양측 미분은 일반적으로 다르다. 다른 가족과의 비교는 표현 민감도 비교가 되며, 관측값에 맞춰 \(k\)를 다시 조절하면 현재의 endpoint-preserving 계약을 벗어난다. 생성한 중간점의 composition은 미정이며 native TF state로 세지 않는다.

## 6. 더 이른 outer knot의 지위와 유한 검사

E outer의 prefix를 더 일찍 끝내고 그 knot의 이미 적분된 \(n,P,\epsilon,\mu\)를 그대로 쓰는 것은 수학적으로 정의 가능한 **다른 lower closure**다. 첫 energy anchor나 retained prefix를 바꾸지 않으면 상대 energy shift가 생기지 않는다. 그러나 기존 outer의 유효구간 일부를 새 phenomenological bridge로 대체하므로, 그 선택 자체를 microscopic 재생성이라고 부를 수 없다.

물리적으로 선택하려면 retained outer와 실제 native TF support의 적용영역, composition 및 energy convention을 먼저 정하고, stellar 관측량과 무관한 규칙을 사용해야 한다. 단지 cubic 통과 여부가 유리해진다는 이유만으로 nuclear shell·pairing·neutron-drip의 일관성이 복구되지는 않는다. TF의 낮은 밀도 해가 수치적으로 수렴하는 것과 그 근사가 outer 핵종 physics를 재현하는 것도 구별한다. 이 H support는 source가 마지막 knot를 어떤 정확한 neutron-drip 조건으로 정했는지 새로 결정하지 않았다.

진단에서는 E의 **217개 저장 knot 전부**를 각 G trial에 적용했고, 기존 cubic validator와 식 (4)를 비교했다. DDME2 여섯 trial의 총 1,302개 조합에서 cubic 수용은 **0개**였다. 여섯 경우 모두 (4)에 가장 가까운 knot도 원래 마지막 **row 216**이었다. 따라서 이 유한 G 집합에서 earlier-outer cut은 direct cubic 실패를 해결하는 방법이 아니다. 새 H native lower endpoint나 다른 source EOS에 대한 보편적 불가능성으로 확장하지 않는다. 전체 열 개 trial의 2,170개 조합은 [scan CSV](../../data/low_density_crust_20261117/thermo_construction/EARLIER_OUTER_SCAN.csv)에 보존했다.

## 7. H 실제 저밀도 연장에 적용할 유한 절차

1. Frozen F parameter identity와 G에 연결된 실제 native continuation을 유지한다. 실패점 사이를 채운 보간 상태를 accepted TF point로 세지 않는다.
2. 실제 새 knot를 낮은 밀도부터 검사한다. 고정 E endpoint에 대해 양의 \(n,P,\epsilon\), 증가하는 \(n,P,\mu\), Euler 관계를 확인한 뒤 (1)과 (4)를 별도로 계산한다.
3. (4)를 통과하고 그 위 native support 및 sharp core–crust contract가 모두 인증되면, 관측량 계산 전에 정한 leftmost-admissible 규칙으로 lower bridge를 짧게 한다. 낮은 \(n_R\)가 언제나 admission을 보장하는 것은 아니므로 실제 \(P_R,\epsilon_R\)가 필요하다.
4. (1)을 통과하지만 (4)가 실패한 경우는 **단일 cubic 표현의 실패**로 남긴다. 현 문서의 (5)는 필요한 경우 검토할 명시적 대안이며, 새 whole-EOS contract와 소비자 검증 없이 자동으로 production에 넣지 않는다.
5. (1) 자체가 실패하면 **그 고정 endpoint pair와 매끄러운 causal single-barotrope 클래스**에 한정된 incompatibility다. 파라미터 변경, energy offset 또는 전체 sharp-Maxwell 후보의 no-go로 바꾸지 않는다.

이번 support에서는 공유 physics 구현, E outer/surface, F parameter, source EOS 또는 production contract를 수정하지 않았다. 새 native TF 결과와 TOV/Cowling 비교는 H production 보고서를 따른다.

## 8. 실행 영수증과 CAS 상태

- [분석 script](../../data/low_density_crust_20261117/thermo_construction/analyze_endpoints.py)의 실제 실행은 **exit 0**, `2026-11-17T14:29:17.486308+00:00`–`2026-11-17T14:29:18.666921+00:00`다. [RUN_01.json](../../data/low_density_crust_20261117/thermo_construction/RUN_01.json), stdout, stderr를 보존했다.
- [ENDPOINT_CONSTRUCTION.json](../../data/low_density_crust_20261117/thermo_construction/ENDPOINT_CONSTRUCTION.json)에 E contract, G admission, frozen parameter JSON, 기존 구현과 source-note의 SHA256 및 네 G matched table의 파일별 재검증을 기록했다. Script 전후 입력 hash는 같았다.
- Wolfram context 호출은 MCP `-32603 Internal error`로 실패했다. 첫 evaluator 입력에는 괄호 문법 오류가 있어 `$Failed`를 반환했고, 한 번 수정한 실제 두 번째 호출에서는 (2), (3), \(C'_x=(a+dx)n''_{xx}/d\), 양끝 \(C\)에 대한 네 symbolic zero check가 모두 `True`였다. Primitive와 \(k\to0\) limit도 반환했다.
- 두 번째 CAS 호출에는 `Symbol::undefined2`, `Power::infy`, `Infinity::indet` 메시지가 동반되었다. 따라서 warning-free CAS run이라고 보고하지 않는다. 정확한 입력과 두 출력, context failure를 [WOLFRAM_RECEIPT.json](../../data/low_density_crust_20261117/thermo_construction/WOLFRAM_RECEIPT.json)에 모두 남겼다. 직접 유도와 독립 수치 적분이 본 결론의 근거이며 CAS는 보조 검산이다.
