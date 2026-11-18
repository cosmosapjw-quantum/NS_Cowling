# H EOS matching 보간법별 교차검증 (2026-11-18)

**근거:** 공개 확인된 H commit `634a586fdd0c214bf7fa9718b02ebafa3e927008`의 8개 strict N128 whole-EOS contract. 원 TF source 값, 고정 outer/surface와 sharp Maxwell doublet는 수정하지 않았다. I 연구 단계는 시작하지 않았다. 이 결과는 별도 로컬 branch에서 실제 실행된 유한 보간법 실험이다.

## 두 matching 문제의 분리

- **Sharp coexistence:** \(P_-=P_+,\mu_-=\mu_+,n_-\ne n_+\). 실제 Maxwell density jump는 어떤 family에서도 보간·제거하지 않는다.
- **Smooth lower bridge:** 공개 E outer와 H의 실제 stationary TF source 시작점 사이의 미시적으로 미해결인 구간. 대체 barotrope는 명시적인 phenomenological closure이고 실패한 native density knot를 새 stationary state로 승격시키지 않는다.

## cold barotrope와 방법의 동등 비교

\(n=dP/d\mu,\;\epsilon=\mu n-P,\;c_s^2/c^2=n/(\mu\,dn/d\mu)\). \(y=n/\mu\), \(a=\mu_L,b=\mu_R,d=b-a\)라 두면 \(y_R>y_L\)이고 서로 다른 endpoint에서 일반 smooth causal connector의 필요충분조건은 \(y_L(b^2-a^2)/2<\Delta P<y_R(b^2-a^2)/2\)이다. 충분성은 증가하는 \(y(\mu)\)의 적분값을 연속 조정하여 구성한다. cubic은 이보다 좁은 shape 조건을 요구한다.

| 보간법 | 8개 중 endpoint/물리 수용 | 구분 |
|---|---:|---|
| Cubic Hermite \(P(\mu)\) | 2/8 | 6개는 cubic 내부 안정성 또는 인과성 불충족 |
| Exponential monotone \(n/\mu\) | 8/8 | endpoint 압력 적분이 단일 \(\kappa\)를 고정 |
| Power-law \(\Phi=x^\alpha\) | 8/8 | bridge 왼쪽에서 \(c_s=c\); 이웃과 큰 미분 불일치 |
| Rational \(\Phi=x/[q+(1-q)x]\) | 8/8 | endpoint 보존·인과성, 근 적분 독립 검산 |
| Nonnegative Bernstein slope-matched | 8/8 | 이웃 양끝 음속과 pressure area를 동시에 고정 |
| 단순 log-log \(P(n)\) + first-law reconstruction | 0/8 | 원래 \(n,P\) endpoint는 보존하지만 원래 \(\epsilon_R\) 불일치 |

Log-log polytrope의 0/8은 polytrope 일반의 물리적 불가능성이 아닌 *동일한 세 상태량 endpoint를 모두 고정한 계약*의 실패다. 그로 인한 화학퍼텐셜의 endpoint 차이는 아래 원자료에 각 조합별로 남겼다.

## 실제 1.4 태양질량 TOV 및 음속 비교

표의 반지름 변화는 **같은 모델·같은 첫 연결점에서 exponential을 기준**으로 한 \(\Delta R\)이다. 단위는 mm이며, \(c_s^2/c^2\)의 최댓값은 해당 phenomenological lower bridge에서만 평가했다.

| Crust → core | \(n_R\) [fm⁻³] | \(R_\mathrm{exp}\) [km] | \(\Delta R_\mathrm{power}\) [mm] | \(\Delta R_\mathrm{rational}\) [mm] | \(\Delta R_\mathrm{Bern}\) [mm] |
|---|---:|---:|---:|---:|---:|
| DDME2 → PK1_CANONICAL | 0.006 | 14.420416664 | -2.7477 | +1.5135 | -2.4757 |
| DDMEX_ORIGINAL_TABLE_RHO0 → GM1 | 0.003 | 13.875173981 | -0.1818 | +0.0739 | -0.1808 |
| NL3_PROVIDER_PRECISION → GM1 | 0.005 | 13.875197801 | -0.7066 | +0.6325 | -0.1086 |
| GPPVA_NL3wrL55 → PK1_CANONICAL | 0.001 | 14.420421582 | -0.0120 | +0.0683 | +0.0181 |

8개 전체 조합의 수용된 method 차이 중 최대 \(|\Delta R|\)는 **3.7060 mm**다. 이 수치는 방법 간 유한 차이이지 전통적인 관측 오차 또는 보간 불확실성의 엄밀한 상한이 아니다. 같은 중심밀도로 수행한 tighter TOV control의 최대 \(|\Delta R|\)는 0.0216 mm여서 일부 매우 작은 방법 차이는 이 실험의 수치 감도와 같은 규모다. 최대 target mass residual은 4.43e-08 M_sun이다.

| 모델(primary) | exponential 최대 cₛ²/c² | rational 최대 | power 최대 | Bernstein 최대 / 차수 |
|---|---:|---:|---:|---:|
| DDME2 | 0.00136997 | 0.00143403 | 1 | 0.00263275 / D=4 |
| DDMEX_ORIGINAL_TABLE_RHO0 | 0.00121231 | 0.00114447 | 1 | 0.00263275 / D=4 |
| NL3_PROVIDER_PRECISION | 0.00110192 | 0.00097369 | 1 | 0.00263275 / D=4 |
| GPPVA_NL3wrL55 | 0.114647 | 0.00748124 | 1 | 0.0316131 / D=16 |

Exponential·power·rational 간 원시 pressure area의 독립 quadrature 최대 절대잔차는 4.22e-16 MeV fm⁻³이었다. 별의 반지름이 거의 같더라도 \(c_s^2\) 프로파일은 크게 달랐다. Power-law 방법은 이 8개 lower bridge에서 causal gate는 통과하지만 왼쪽 endpoint에서 \(c_s=c\)로 귀결되어 압축률 경계 불연속이 크다.

## 양끝 음속을 보존하는 Bernstein family

\(y=y_L+(y_R-y_L)\Phi(x),\;g=\Phi'=\sum_{i=0}^{D}b_i B_i^{D}(x)\)로 두고, 모든 \(b_i\ge0\)를 제한했다. 따라서 \(g\ge0\)이고 \(c_s^2/c^2\le1\)가 전 구간에서 보장된다. \(\int_0^1 gdx=1\), endpoint slope 두 조건, pressure area 한 조건은 Bernstein 계수에 **선형 제약**을 준다. 검증에는 비음수 linear programming feasible certificate를 사용했으며, 가능한 경우 이산 기울기 변화를 줄이는 추가 constrained minimization을 적용했다.

양끝 \(c_s^2/c^2\)는 E outer의 마지막 polytrope interval과 N128 source upper의 첫 cubic Hermite interval에서 정해진 one-sided 값이다. 이 값들은 직접 계산된 microscopic sound-speed derivative가 아니다. 따라서 이를 맞춘 결과는 **고정된 표현에 대한 2차 압력 도함수의 접합**이지, 원전 microscopic response의 검증은 아니다.

8개 모두 해를 구했고, 최소 feasible Bernstein 차수의 분포는 {4: 6, 16: 1, 6: 1}였다. NL3wrL55의 첫 연결에서 D=16을 필요로 하는 것은 매우 단순한 저차 smooth connector로는 endpoint slope와 area를 동시에 만족하기 어렵다는 유한-degree 사실이지 일반적인 no-go가 아니다.

**실험 범위의 제한:** 기존 H의 Maxwell 교차근과 sharp doublet를 고정해 주변 보간을 검증했다. 원래 crust/core source에서 서로 다른 보간식을 적용하여 Maxwell crossing root search를 처음부터 다시 탐색하지는 않았다. 따라서 아래 전이 보존 8/8은 *이미 발견된 전이의 소비자 보존성*을 뜻하며, 보간법별 새로운 Maxwell root 발견률이나 그 이동량을 측정하지 않는다. 해당 질문은 별도 원자료 교차근 재탐색 실험으로 남긴다.

## Sharp Maxwell 주변 interval만 변경한 별도 검증

전이 \(j\)의 **밀도 불연속 interval \(j\)**은 그대로 두고, 바로 인접한 smooth interval \(j-1\)과 \(j+1\)의 형식만 바꿨다. 8개 모두 Hermite/Exponential/Rational은 수용됐다. Power-law는 8개 모두 한쪽 이상에서 \(\alpha<1\)이어서 endpoint susceptibility가 발산하고 \(c_s^2\to0\)인 비정칙 표현으로 거절됐다. 이는 Maxwell의 물리적 부재 판정이 아니다. 모든 통과 case에서 transition \(P_t\), \(\mu_t\), \(n_-,n_+\)와 source 저장값은 동일했다.

## 재현·검증 및 잔여 범위

- H 원본 branch는 원격 push 완료 후 변경하지 않았다. 실험은 로컬 `research/matching-interpolation-crosscheck-20261118` worktree에서 실행했다.
- 새 source의 bridge 8개, one-sided transition adjacent 8개, slope-matched Bernstein 8개에 대해 실제 TOV를 계산했다. 성공·거절의 기계 기록과 각 밀도구간의 수치 결과는 `data/matching_interpolation_20261118/`에 보존했다.
- 집중 회귀는 **108 passed, exit 0**, 전체 frontend/pipeline은 **1490 passed, 9 skipped, 4 subtests passed, exit 0**이었다. 전체 검사에서 source 파일 185개의 실행 전후 SHA가 같았고, 기존 native 실행 파일을 고정해 재사용했다. SymPy 1.14.0으로 rational 및 power 압력 원시함수 미분 항등식을 각각 확인했다.
- Wolfram 연결 시도가 서버 오류로 실패하여 이번 symbolic 결과를 Wolfram 검증으로 표기하지 않는다.
- Native optimizer가 실패한 6개 밀도점, TM1e N128 transition bracket, 모든 geometry global minimum, 새 causal closure의 low-frequency Cowling mode는 여전히 해결되지 않았다. 별의 radial stability와 complete mode spectrum도 주장하지 않는다.

## 현재 결과에 근거한 방법 선택

1. **기준선:** H에서 이미 검증된 exponential \(n/\mu\).
2. **독립적 형상 교차검증:** rational \(n/\mu\). 양끝 1차 법칙을 동일하게 만족하되 내부 압축률 함수가 다르다.
3. **미분 연결의 조건부 후보:** slope-matched Bernstein. 실제 TF source의 미분 응답을 더 직접 확인할 수 있을 때 선택의 물리적 신뢰도가 올라간다.
4. **Stress test:** power-law. causal PASS와 endpoint derivative artifact의 차이를 보여준다.
5. **거절의 올바른 해석:** cubic failure를 general smooth causal no-go로 해석하지 않으며, log-log polytrope의 source energy 불일치를 허용 EOS로 처리하지 않는다.

## 주요 문헌

- [Haensel & Potekhin (2004), Analytical representations of unified equations of state of neutron-star matter](https://www.ioffe.ru/astro/Stars/Paper/hp04.pdf), first-law consistency and interpolation.
- [Fortin et al. (2016), Neutron star radii and crusts: Uncertainties and unified equations of state](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.94.035804).
- [Ferreira & Providência (2020), Neutron Star Properties: Quantifying the Effect of the Crust–Core Matching Procedure](https://www.mdpi.com/2218-1997/6/11/220).
- [Canullán-Pascual et al. (2025), Consistent Crust-Core Interpolation and Its Effect on Non-radial Neutron Star Oscillations](https://scispace.com/papers/consistent-crust-core-interpolation-and-its-effect-on-non-2fysc659doqb).

**재개 지점:** I는 시작하지 않는다. Matching 범위에서 이어갈 가장 작고 물리적으로 중요한 작업은 실제 원 미시모델로부터 인접 구간의 \(dn/d\mu\)를 더 신뢰성 있게 추정하고, 해당 불확실성을 포함해 독립 family 간 EOS/별 mode 차이가 수치 오차를 넘는지를 판정하는 것이다.
