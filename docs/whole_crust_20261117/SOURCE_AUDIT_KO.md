# Whole-crust A: 입력의 독립적인 수치 진단

## 범위와 재현

기준 자료는 세 계열 Maxwell 산출물이 저장된 `8c917baf3f5975a6109c3b6e59177079c285c4e0`의
기존 입력 및 `docs/matched_eos_reference_20261117/PLAN_KO.md`의 A 계획이다.
이번 진단은 EOS를 생성하거나 수정하지 않는다. 기존 CoMPOSE 행을 명시적인 density
index로 읽고 물리적인 (n,P,\epsilon,\mu_B)를 복원했다. 원본 파일 SHA256과 선택한
prefix의 semantic hash는 `data/whole_crust_20261117/source_audit/SOURCE_AUDIT.json`에 있다.

```bash
python scripts/diagnose_whole_crust_sources_20261117.py
```

실제로 21개 입력을 읽고 2,480개의 유효한 **새 cubic bridge 후보**를 기록했다.
후보 수는 미시 모형 또는 Maxwell 평형 해의 수가 아니다. 후보 표는
`BRIDGE_CANDIDATES.csv`에 있으며 모든 값의 단위는 (n\,[\mathrm{fm}^{-3}]),
(P,\epsilon\,[\mathrm{MeV\,fm}^{-3}]), (\mu\,[\mathrm{MeV}])다.

## 물리량 복원과 표현 진단

복원식은 (P=nQ_1\), (\epsilon=n m_n^{\rm ref}c^2(1+Q_7)\),
(\mu_B=m_n^{\rm ref}c^2(1+Q_3)\)이다. 이번 보간 표현 검사는

\[
\mu_E=(\epsilon+P)/n
\]

를 사용하는 별도 진단이다. 원본 Q3를 고친 것으로 기록하지 않는다.
찬 한-charge barotrope에 대해 endpoint의 (P,n,\mu_E\)를 고정한 cubic (P(\mu_E)\)를
만들어, 각 interval 전체에서 다음 다항식의 극값을 계산했다.

\[
n=P'>0,\qquad n'=P''>0,\qquad \mu_E n'-n\ge0.
\]

마지막 조건은 (c_s^2/c^2=n/(\mu_E n')\le1\)이다. 표본점만의 검사가 아니다.
원본 (P,\epsilon,n\)의 유지와 정확한 연속 보간의 안정성은 별도 조건이다.
비교할 phase가 바뀌는 구간에는 별도 경계 또는 명시적 bridge가 필요하다.

## 공개 BSk22 lower block

공개 SPG_M1–M5의 첫 217개 행은 열역학적으로 동일한 lower block이다.
공급자 Iphase는 모두 0이지만, composition에는 **18개 동위원소 구간과 17개 경계**가 있다.

| 값 | 첫 행 | 마지막 행 |
|---|---:|---:|
| (n) | (4.67962\times10^{-10}) | (2.69021\times10^{-4}) |
| (P) | (9.000000079422\times10^{-12}) | (4.981910045986\times10^{-4}) |
| 저장된 (\mu_B) | 929.80107634726 | 938.89623004093 |
| (\mu_E) | 929.80109230585 | 938.89754032057 |

원본의 (\mu_E-\mu_B\) 절댓값 최댓값은 **0.00132422950 MeV**이며, 마지막 몇 행에
한정된 현상이 아니다. 동위원소별 구간 안에서도 증가하는 체계적인 잔차가 있다.
Q3만 Euler 값으로 투영해도 열역학 표현 문제가 모두 해결되지는 않는다.

실제로 이 217개 행의 인접 knot를 전부 유지하는 cubic (P(\mu_E)\) Hermite는
**216개 interval 중 40개에서 (dn/d\mu_E>0\) 조건을 위반**한다. 여기에는 여러 동위원소
구간 내부의 interval도 포함된다. 이는 해당 원자료의 물리적 EOS 전체에 대한 불가능성
명제가 아니라, 지정한 knot와 미분값을 동시에 유지하는 해당 cubic 표현의 실패다.
그 실패 interval을 통과하는 전체 표를 같은 보간법으로 인증할 수 없다.

공개 마지막 행의 밀도는 모든 native SPG의 첫 밀도보다 높지만, chemical potential은
native 첫 행보다 낮다. Native 첫 행의 (\mu_E\)는 모두 약 938.900000000 MeV다.
또한 공통 (\mu\) 지지구간이 없다. 따라서 마지막 public 행과 첫 native 행의 직접 append나
단일 Maxwell 접합은 현재 제공된 이 endpoint만으로 성립하지 않는다.

더 낮은 public endpoint를 선택하면 명시적인 새 smooth bridge는 구성할 수 있다.
단순 cubic의 최소 폭 후보는 다음과 같다. 이것만으로 앞의 전체 public lower block을
사용할 수 있다는 뜻은 아니며, 40개 표현 실패와 17개 isotope 경계 문제는 함께 해결해야 한다.

| Native SPG | Public 왼쪽 행 | (n_L) | Native 오른쪽 행 | (n_R) |
|---|---:|---:|---:|---:|
| M1 | 171 | 0.000159336 | 0 | 0.0002331736909364 |
| M2 | 153 | 0.000123830 | 0 | 0.0002416219329338 |
| M3 | 175 | 0.000178627 | 0 | 0.0002355410299484 |
| M4 | 179 | 0.000198954 | 0 | 0.0002475602850152 |
| M5 | 147 | 0.000112224 | 0 | 0.0002361127426875 |

행 번호는 모두 0부터 센다. Endpoint 선택에 앞서 필요한 convexity 조건

\[
n_L<\frac{P_R-P_L}{\mu_R-\mu_L}<n_R
\]

을 검사했고, 후보 표에는 전체 polynomial-extrema 결과도 포함했다.

## 기존 독립 outer-crust의 재사용 가능성

`data/research_finish_20261116/crust_assembly/restored_spg/SPG_M1/compose`의 첫 547개 행은
당시 `ASSEMBLY.json`에서 명시한 **independent_outer** 구간이다. 뒤의 restored-SPG
variant와 이 prefix를 혼동하지 않는다. 해당 prefix는 다음 범위를 갖는다.

| 값 | 첫 행 | neutron-drip 쪽 마지막 행 |
|---|---:|---:|
| (n) | (1.00000000000027\times10^{-9}) | (2.69432006401090\times10^{-4}) |
| (P) | (2.993814854565\times10^{-11}) | (4.992059167690\times10^{-4}) |
| (\mu_E) | 930.489177292586 | 939.565420520000 |

**17개의 명시적인 공통 (P,\mu\) 동위원소 전이**가 있고, 나머지 **529개 매끄러운 interval
전부**가 위 cubic 표현의 안정성·causality 검사를 통과한다. 동위원소 전이 사이의 density gap을
보간하지 않아야 한다. 동일한 의미의 prefix가 과거 TF 부분-crust 자료에도 보존되어 있다.

이 outer의 물리적 chemical potential은 public lower와 다르다. 이 차이를 header 재인코딩으로
제거할 수 없으며, 공개본과 동일한 EOS라고 부르지 않는다. 에너지에 (C n\)을 추가하지 않고,
현재 native SPG와의 새 bridge를 정의할 때 실제로 가능한 첫 endpoint는 다음과 같다.

| Native SPG | Native 오른쪽 행 | (n_R) | (\mu_R) |
|---|---:|---:|---:|
| M1 | 47 | 0.001957033667502 | 940.6485098294 |
| M2 | 45 | 0.001930680458687 | 940.6645290892 |
| M3 | 47 | 0.002012347936895 | 940.6565360904 |
| M4 | 46 | 0.002130507270560 | 940.7059026377 |
| M5 | 47 | 0.001938971209431 | 940.6075919176 |

이 경로는 이미 존재하는 검증된 outer와 명시적 lower-inner bridge를 사용하는 대안이다.
Bridge 내부의 미시적 조성·평형은 별도 모델링 가정으로 남는다. 이 표는 어떤 경로를 최종
baseline으로 채택했는지를 선언하는 생산 receipt를 대신하지 않는다.

## 기존 SPG–SLy4 core 연결의 표현 차이

기존 `campaign/pairs/SPG_Mi__RG_SLY4__quintic_P/compose`는 (n=0.5\)까지의 표다.
M1/M3/M5의 표는 모든 인접 knot를 사용하는 cubic (P(\mu_E)\) 표현 검사를 통과한다.
M2/M4에서는 **interval 128**, 즉 native crust 마지막 행에서 기존 quintic bridge 첫 행으로
가는 구간이 해당 cubic으로 표현될 때 (dn/d\mu_E\) 조건을 위반한다.

이것은 기존 C2 quintic 자체가 실패했다는 판정이 아니다. 서로 다른 보간 표현을 같은 것으로
취급하면 생기는 문제다. M2/M4를 사용하려면 원래 quintic coefficient를 소비자가 사용하거나,
그 quintic에서 추가 knot를 평가한 뒤 새 표현을 독립적으로 확인해야 한다.

CSS 표 5개는 interval 128의 기존 공통 (P,\mu\) density jump를 제외한 모든 interval이
매끄러운 cubic 표현 검사를 통과한다. 이때 sharp interface의 identity와 출력 provenance를
계속 유지해야 한다.

## 이번 진단의 결론

M1/M3/M5의 기존 smooth upper 표와, 정확한 isotope 전이가 보존된 independent outer는
연속적인 whole-crust 후보를 만드는 데 즉시 사용할 수 있는 입력이다. 둘 사이에는 표에
명시한 새 bridge가 필요하며, (n<10^{-9}\)에는 별도의 surface prescription이 필요하다.
공개 lower를 그대로 선택하는 경로에는 40개의 표현 실패와 17개 isotope 경계에 대한
추가 처리가 필요하다. 최종 production에서 선택한 경로·보간·surface 처방의 검증 상태는
별도 생산 기록으로 판단한다.

이 진단 이후 A의 생산 계획은 원래 공개 lower 출처를 유지하기 위해 별도의
**pressure-primary 표현**을 채택했다. 이 표현은 원본 (n,P\) knot와 첫 energy/baryon을
고정하고 제1법칙으로 나머지 energy를 적분하며, 원본의 모든 (\epsilon\) knot를 동시에
고정하는 앞의 실패한 cubic과 다르다. 동위원소 switch의 빈 구간은 명시적인
regularization으로 표기한다. 위 independent-outer 경로는 비교 가능한 대안으로 남으며,
실제 baseline과 구성 민감도 결과는 `PLAN_KO.md`와 생산 receipt를 따른다.
