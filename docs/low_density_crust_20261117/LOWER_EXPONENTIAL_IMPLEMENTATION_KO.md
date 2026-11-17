# H opt-in causal lower bridge 구현과 검증

이 addendum은 [첫 열역학 진단](LOWER_THERMODYNAMICS_KO.md)을 freeze한 **뒤** parent의 H 실행 범위 안에서 추가한 생산 코드 변경을 기록한다. 첫 문서의 “아직 생산 contract에 채택하지 않은 구성예”는 그 진단 시점의 상태다. 여기서는 해당 family를 명시적 opt-in으로 구현하고 소비자와 함께 검증했다. 실제 어떤 H source knot와 stellar 모델에 채택했는지는 H admission 및 stellar 영수증을 따른다.

## 구현 범위와 API

변경한 생산 파일은 `python_frontend/compose_eos/whole_crust.py` 하나다. 기존 `HermiteInterval` 클래스 본문은 그대로이며, 기존 lower contract의 기본 선택도 Hermite다. `CausalExponentialInterval`과 lower-contract dispatch를 추가하고 pressure inverse의 type annotation만 넓혔다. `stellar_background.py`, `stellar_sharp.py`, Cowling 소비자는 수정하지 않았다.

```python
from compose_eos.whole_crust import CausalExponentialInterval, WholeCrustBarotrope

bridge = CausalExponentialInterval(
    n_left, p_left, epsilon_left,
    n_right, p_right, epsilon_right,
)
assert bridge.validate()["status"] == "PASS"
record["lower_bridge"] = bridge.contract()
barotrope = WholeCrustBarotrope.from_contract(record)
```

Constructor의 optional `kappa=`는 저장된 contract를 다시 읽을 때 사용한다. 새 construction은 이를 생략해 양끝 \(n,P,\epsilon\)만으로 유일한 pressure-integral root를 구한다. 저장된 \(\kappa\)를 읽을 때도 endpoint root를 다시 계산하여 일치 여부를 검사한다. 자유롭게 주파수·질량·반지름에 맞추는 입력이 아니다.

Contract의 식별자는

- `kind = "causal_exponential_n_over_mu"`
- `representation = "causal_exponential_density_over_chemical_potential"`

이다. 저장 \(\kappa\), source endpoint의 정확한 \(n,P,\epsilon\), \(\mu\) support, 명시적 미정 composition, `source_energy_shift="NONE"`와 analytic validation을 포함한다. `kind`가 없으면 새 family를 자동 추측하지 않는다. 모르는 kind/representation, 누락·변조된 \(\kappa\), endpoint와 불일치하는 저장 \(\mu\) support는 거절한다.

새 클래스는 기존의 `evaluate_mu`, `pressure(t)`, `number(t)`, `susceptibility(t)` interface를 제공한다. `table(count, template=...)`는 미시적 조성을 만들지 않는 derived thermo table을 출력한다. 양끝 \(n,P,\epsilon\)는 입력값을 정확히 반환하고, 그 사이에서는 하나의 명시적 first-law 표현을 평가한다.

## 해석적 수용조건과 수치 표현

첫 문서의 notation \(a=\mu_L\), \(d=\mu_R-\mu_L\), \(x=(\mu-a)/d\), \(y=n/\mu\)를 그대로 사용한다. Positive ordered endpoint, \(y_R\ge y_L\), 그리고 서로 다른 \(y\) endpoint에 대한 strict causal pressure-integral bounds를 먼저 검사한다.

\(y_L=y_R\)의 degenerate domain은 \(n=y_L\mu\), \(c_s=c\)인 유일한 연결이다. 이 경우 canonical \(\kappa=0\)만 허용하고, pressure integral의 일치는 명시한 floating-point roundoff bound 안에서 확인한다. 일반 \(y_L<y_R\)의 strict integral bounds는 tolerance로 넓히지 않는다. 에너지에 \(Cn\)을 추가하거나 실패 endpoint를 움직이는 처리는 없다.

Overflow와 cancellation을 피하기 위해 \(\Phi_\kappa\)의 primitive를 두 적분으로 나눈다.

\[
F_0(x)=\int_0^x\Phi_\kappa(t)\,dt,\qquad
F_1(x)=\int_0^x t\Phi_\kappa(t)\,dt,
\]

\[
P(x)=P_L+d\left[
y_L\left(ax+\frac{d x^2}{2}\right)
+(y_R-y_L)(aF_0+dF_1)\right].
\]

\(|z|=|\kappa x|\le1/2\)에서는 \(r(\kappa)=\kappa/\operatorname{expm1}(\kappa)\), \(r(0)=1\)를 사용한 series

\[
F_0=r(\kappa)x^2
\sum_{j=0}^\infty\frac{z^j}{(j+1)!(j+2)},\qquad
F_1=r(\kappa)x^3
\sum_{j=0}^\infty\frac{z^j}{(j+1)!(j+3)}
\]

를 계산한다. 구현은 \(j=18\)까지 사용한다. 이 유한 구간에서 잘린 항의 크기는 double-precision roundoff보다 작다. \(\kappa=0\)은 정확히 \(F_0=x^2/2,F_1=x^3/3\)로 처리한다.

양의 큰 \(\kappa\)에는 \(e^{\kappa(x-1)}\)와 \(e^{-\kappa}\)만 사용하고 \(e^{+\kappa}\)를 만들지 않는다. 음의 큰 \(\kappa=-q\)에는 \(1-e^{-qx}\), \(1-(1+qx)e^{-qx}\)의 형태를 사용한다. 작은 \(qx\)는 위 series가 맡으므로 두 큰 항의 차감으로 작은 integral을 잃지 않는다. 실제 테스트에는 \(\kappa=\pm800\), \(\pm10^{-9}\), 0을 포함했다.

압력 root는 bounded expanding bracket과 Brent solver로 구한다. Endpoint pressure integral을 다시 평가한 잔차와 roundoff bound를 contract에 남긴다. 유한 float64 입력으로 극한 causal bound 근처의 root를 분리할 수 없으면 `ValueError`로 보고한다. 그런 수치적 미해결을 물리적 no-go로 해석하지 않는다.

## 연속구간 검증의 근거

\(y'>0\)이므로 \(n'=y+\mu y'>0\), \(\mu n'-n=\mu^2y'>0\)이다. 이 항등식이 일반 family의 continuum 안정성·인과성 근거다. 보고하는 extrema는 단순 grid 표본의 최솟값이 아니다.

\[
n''=y'(2+\kappa\mu/d),\qquad
\frac{d}{d\mu}(\mu^2y')=\mu y'(2+\kappa\mu/d)
\]

이므로 susceptibility와 causal margin은 양끝 및 support 안의 \(\mu=-2d/\kappa\)만 확인하면 된다. Sound speed의 stationary condition은

\[
F(\mu)=y+\mu\kappa A/d=0,\qquad
A=y_L-\frac{y_R-y_L}{e^\kappa-1}.
\]

여기서 \(F'=\kappa y/d\)이므로 내부 root는 최대 하나다. 구현은 안정적인 \(\kappa A=\kappa y_L-(y_R-y_L)r(\kappa)\)를 사용해 양끝과 부호가 바뀌는 유일한 내부 root를 평가한다. \(\kappa=0\)의 극한도 유한하다. 양의 source energy와 \(\epsilon'=\mu n'>0\)에서 energy positivity가 따른다.

이 검증은 선언한 barotrope의 검증이다. Native optimizer 실패가 있었던 구간을 이 family가 덮는다면 그 구간은 **실패를 보존한 상태에서 추가한 phenomenological closure**이며, 성공한 TF interpolation이나 microscopic 재생성으로 표기해서는 안 된다.

## 실제 테스트와 freeze

`test_causal_lower_bridge.py`의 새 26개 parametrized case와 기존 whole-crust/background/sharp tests를 함께 실행했다. 실제 결과는 **71 passed, exit 0**, runtime 2.78 s다. [TEST_01.json](../../data/low_density_crust_20261117/thermo_construction/TEST_01.json), stdout, stderr에 명령과 exact hash가 있다.

새 검증에는 다음이 포함된다.

- 독립 quadrature로 만든 endpoint에서 \(\kappa\) 복원, 안정적 primitive 및 exact source endpoint 반환.
- First law와 susceptibility의 독립 중앙차분, analytic extrema와 probe 일치.
- Strict causal bounds와 \(y_L=y_R\) degenerate domain, nonphysical endpoint, nonfinite/boolean parameter 거절.
- 저장된 \(\kappa\), kind/representation 및 \(\mu\) support 변조 거절과 JSON roundtrip.
- Pressure inverse와 chemical-potential/enthalpy 평가의 일치, 변경되지 않은 outer 및 upper source 상태.
- G DDME2 고정 endpoint의 cubic failure와 새 causal connector 수용을 같은 입력으로 비교한 회귀.
- Sharp-interface enthalpy TOV와 독립 pressure-coordinate TOV의 질량·반지름 비교, 인터페이스 밀도 공백 유지.

이번 code freeze:

| 파일 | SHA256 |
|---|---|
| `python_frontend/compose_eos/whole_crust.py` | `7bb94e914019045a636f2008d5ed4ce2b6e346c8f79e19f0100416c5e924b513` |
| `python_frontend/tests/test_causal_lower_bridge.py` | `ebcba5313a21ea178dd75ec594f17fd827daa3de786f13fba481e0c221914a93` |

이는 opt-in 생산 경로의 구현 검증이며, H의 새 stellar 결과나 모든 microscopic 후보의 타당성을 미리 승인하는 영수증이 아니다. H admission은 실제 connected native source support, endpoint 순서, 이 connector, upper smooth intervals와 sharp coexistence contract를 각각 확인해야 한다.
