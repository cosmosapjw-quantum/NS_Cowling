# NS_Cowling 개선 B — 실제 core 지원 검증과 primary EOS 확장

2026-11-17 KST. 기준: A branch `research/whole-crust-20261117`, commit `4069b258c906790d315e1f4760008bea9e41b598`.
작업 branch: `research/core-support-20261117`. B 결과의 게시에는 별도 사용자 승인이 필요하다.

## 결과

A의 primary whole-crust EOS를 실제 생성된 SLy4 core로 확장했다. 기존 704행에 69행을 추가한 **773행**이며, 연속 표현의 실제 지원은 **표면부터 1.0 fm⁻³**다. A의 모든 기존 행 필드와 아래쪽 연속식을 보존했다. 표의 유한 양의 최저 밀도 아래는 A의 동일한 analytic surface prescription이 P=0까지 이어진다.

GM1, SLy4, KDE0v1의 기존 385행씩을 재사용하고, 내부점과 끝점에서 총 1,158개 native 진단을 추가했다. 세 모델의 실제 생성 범위 0.02–1.0 fm⁻³에서 정의한 연속 P(μ) 표현은 양의 압축률과 인과성을 만족했다. 이 결과는 해당 표현의 구간 검사와 native functional의 유한 수치 검증이다. 모든 실수 밀도에서 microscopic functional의 ground state를 증명한 것은 아니다.

B의 source-variant 결정은 공통 fitted SLy4 kinetic coefficient를 유지하는 것이다. 새로운 particle-dependent variant와 에너지 이동은 필요하지 않았다. GM1의 공개 provider 차이는 기존 unresolved 상태와 수치를 보존했다.

## B1–B5 완료 범위

| 항목 | 실제 산출물 | 판정 |
|---|---|---|
| B1 source identity | A와 extended SLy4의 parameter/species/mass/energy convention ledger | 동일 physics dictionary, 서로 다른 실제 sampling |
| B2 실제 core domain | 세 모델 각 384개 연속 구간, native 내부·끝점, 독립 Skyrme 구현 비교 | 실제 0.02–1.0 fm⁻³의 명시된 범위 검증 |
| B3 primary 상단 연결 | 704+69행, 연속 소비자, native CompOSE write/read | 완료 |
| B4 variant 정책 | 공통 K 유지, 별도 variant 불필요, GM1 residual 보존 | 완료 |
| B5 후속 지원 계약 | 실제 n/P 상한, 상한 밖 입력 거절, C→D 인계 | 완료 |

A–D 중 A와 B의 실행 결과가 준비되었다. C와 D는 이번에 실행하지 않았다. B의 실행 완료와 GitHub 게시는 별도 상태다.

## 물리 정의와 검산

온도 T=0, neutrino-free cold beta equilibrium, 균일하고 spin-saturated한 nuclear functional을 사용한다. SLy4와 KDE0v1은 n,p,e,μ, GM1은 원래 정의한 n,p,e다. ε에는 rest energy를 포함한다. 단위는 n_B: fm⁻³, ε와 P: MeV fm⁻³, μ: MeV다. 음속은 c_s²/c²=dP/dε라는 무차원 비율로 쓴다.

중성인 beta-equilibrium trajectory에서 μ=μ_n이며, 각 연속 조각을 P(μ)로 표현한다.

```text
n = dP/dμ
ε = μ n − P
dε/dμ = μ dn/dμ
dP/dn = n / (dn/dμ)
c_s²/c² = n / (μ dn/dμ)
```

두 native endpoint의 n,P,ε를 cubic Hermite P(μ)에 보존한다. μ=(ε+P)/n에서 시작하므로 endpoint를 맞추기 위한 ε→ε+C n 보정은 없다. 구간 전체에서 n>0, dn/dμ>0와 μ dn/dμ−n≥0를 검사했다. 다항식의 stationary candidate를 평가하고 음속 유리함수의 극점 후보도 검사했다. Wolfram Language에서 first-law/Euler/음속 항등식의 잔차 세 개가 정확히 0이었다. 이 기호 검산은 물리 source fidelity를 별도로 증명하지 않는다.

구간 certificate는 binary64로 정의된 보간식에 대한 수치 검사다. Outward-rounded interval arithmetic의 엄밀한 enclosure로 부르지 않는다. 압축률의 양성은 비균일 phase, spin, pairing, 유한파장 불안정성까지 모두 배제하는 조건이 아니다.

## 세 core의 실제 지원

| 모델 | 입자 | 실제 생성 행 | 실제 n_B 범위 [fm⁻³] | 연속 보간식 c_s²/c² 최소–최대 |
|---|---|---:|---|---|
| GM1 | npe | 385 | 0.02–1.0 | 0.0006333208332–0.7564927757 |
| SLy4 | npeμ | 385 | 0.02–1.0 | 0.002977075640–0.8318776264 |
| KDE0v1 | npeμ | 385 | 0.02–1.0 | 0.003426604429–0.7490182825 |

Metadata에 남은 SLy4 1.4 및 KDE0v1 1.8 fm⁻³는 이번 실제 지원으로 인정하지 않았다. 1.0 fm⁻³는 실제 생성한 자료의 끝이다. 이번에 관측한 인과성 파괴점이나 별의 최대질량 경계가 아니다. 이후 D에서 필요한 central density가 이 범위를 넘으면 그 필요를 근거로 다음 확장을 수행한다.

기존 native knot의 최대 Euler 잔차는 4.55×10⁻¹³ MeV fm⁻³, 전하중성 잔차는 2.54×10⁻¹⁵ fm⁻³ 이하다. 추가 native derivative 진단의 최대 상대 first-law 잔차는 GM1 2.033×10⁻¹¹, SLy4 1.197×10⁻¹¹, KDE0v1 1.786×10⁻¹¹이다. 끝점은 실제 지원 안에 머무는 one-sided 5-point stencil을 사용했다.

SLy4와 KDE0v1은 각각 771개 native state 기록을 독립 Python 구현과 비교했다. 최대 |ΔP|는 각각 1.307×10⁻¹², 9.095×10⁻¹³ MeV fm⁻³이고 최대 |ΔY_p|는 약 1.2×10⁻¹⁵ 이하다. 두 구현이 같은 functional/계수에 따르는지 검증한 결과이며 공개 provider와의 독립적 calibration 일치 주장과 구별한다.

Muon onset은 SLy4 약 0.1184892381, KDE0v1 약 0.1109181666 fm⁻³다. 연속적인 lepton 조성 변화로 취급하며 nuclear first-order density jump를 부여하지 않았다. GM1의 npe 정의는 유지했다. 저밀도 homogeneous branch를 A의 crust 대신 사용하지 않았다.

## A primary의 확장 결과

| 항목 | 값 |
|---|---:|
| 보존한 A 행 | 704 |
| 추가한 native 행 | 69 (source indices 316–384) |
| 전체 행 | 773 |
| 기존 upper anchor n_B | 0.5 fm⁻³ |
| 첫 추가 n_B | 0.5001965919793008 fm⁻³ |
| 최종 n_B | 1.0 fm⁻³ |
| 최종 P | 506.42133443967515 MeV fm⁻³ |
| 최종 ε | 1212.9757155668863 MeV fm⁻³ |
| 최종 μ | 1719.3970500065614 MeV |
| 최종 연속 소비자 c_s²/c² | 0.8318776264232889 |
| 최종 native analytic c_s²/c² | 0.8318619126114971 |

보간식과 native derivative의 마지막 차이는 약 1.57×10⁻⁵이다. 양자를 같은 수치라고 쓰지 않는다. 별 관측량의 허용 오차는 D의 수렴 검증에서 정해야 한다.

A의 마지막 n_B=0.5 상태와 원래 native core row128은 n,P,ε,μ가 정확히 같다. A의 704행에 대한 27개 decoded ROW_FIELDS는 확장 표를 native로 읽은 후에도 array-exact로 같았다. 이것은 전체 serialized EOS 파일의 byte identity를 뜻하지 않는다. 파일에 새 행을 추가했으므로 파일 bytes는 바뀐다.

Surface, outer, lower bridge 계약 및 기존 upper n/P/ε 배열과 모든 기존 smooth interval 계수를 보존했다. 추가 tail은 header 재인코딩 후에도 물리 ε를 유지하며 native source와의 최대 차이는 2.27×10⁻¹³ MeV fm⁻³다.

기존 anchor에서 반환되는 음속은 one-sided interpolation derivative에 따라 왼쪽 0.3497481753549825, 오른쪽 0.34973230083369256이다. 이 약 1.59×10⁻⁵ 차이는 ordinary interpolation knot의 2차 미분 차이이며 phase transition이 아니다. n,P,ε,μ는 연속이고 기존 내부 구간의 음속은 바뀌지 않았다.

A의 마지막 행은 composition masks가 false인 thermo-only 행이다. 기존 zero placeholder와 mask를 보존했고, 새 native 69행에만 source가 실제 제공한 fractions/masks가 있다. 이 availability 변화는 물리적인 조성 불연속으로 해석하지 않는다. 연속 소비자는 equilibrium thermodynamics를 제공하며 fully frozen-species closure를 제공하지 않는다. 기존 native 코드가 frozen이라고 부르는 진단 또한 proton fraction을 고정하되 e–μ 재분배는 평형화한 응답이므로 후자의 closure를 대체하지 않는다.

## 검증 결과

- 새 extension 동작 검사: 구현 전 7개 실패를 보존했고 구현 후 7개 통과.
- 독립 analytic domain certificate 검사: 6개 통과.
- 연관 regression: **155 passed, 0 failed**, 실제 exit code 0.
- 전체 frontend/pipeline CI 범위: **1,370 passed, 9 skipped, 0 failed**, exit code 0, 95.42초. 9개 skip은 현재 환경에 준비되지 않은 기존 raw CompOSE fixture 관련이며, 목록은 FRONTEND_PIPELINE_REPAIRED.log에 있다.
- 최초 전체 실행의 6개 실패는 mpmath 부재로 인한 환경 오류였다(1,364 passed/6 failed/9 skipped). 기존 mpmath 1.3.0 경로를 scoped PYTHONPATH로 연결해 해소했으며 원래 실패와 수정 실행을 모두 보존했다.
- Native CompOSE write/read 및 전체 773행 consumer-knot 일치: 최대 상대오차 n 3.006×10⁻¹⁵, ε 3.098×10⁻¹⁵, μ 1.092×10⁻¹⁵.
- A의 32개 source/input 파일 hash를 실행 전 기준 commit과 대조하고 실행 후 불변을 확인했다.
- 독립 검토는 실제 표·contract·polynomial·source hash·상한 밖 입력 거절을 다시 확인했고 material defect 없이 종료했다.
- 실제 상한 다음 floating-point pressure는 ValueError로 거절한다.

## source-variant 결정과 남은 과학적 문제

SLy4는 K=20.73553 MeV fm², m_n c²=939.5654 MeV, m_p c²=938.272 MeV, 기존 lepton/mass convention을 유지했다. A와 B의 source dictionary가 같다. 역사적 Fortran core 파일의 차이는 lower computational bound .08→.02 및 주석이며, 그 파일의 functional/solver/kinetic 수식은 같다. 역사적 executable 자체의 byte identity까지 확정하지 않았다.

기존 split-K control은 공개 trajectory 차이를 상당히 줄였으나 provider가 실제 사용한 convention의 attribution은 unresolved다. 이번 지원 확장에 variant가 불필요하므로 기본 functional을 바꾸지 않았다. GM1은 rounded-parameter reconstruction이고 exact public nodes에서의 최대 pressure residual 0.4863857320%를 유지했다. 새 fitting으로 provider identity를 만들지 않았다. 자세한 원전·경로·hash는 SOURCE_LEDGER_KO.md와 SOURCE_LEDGER.json에 있다.

TOV 질량–반지름 곡선, maximum mass, radial stability, Cowling frequencies 및 composition-dependent buoyancy는 아직 계산하지 않았다. 이번 결과는 그 계산에 전달할 EOS와 실제 지원을 완성한 것이다.

## 재현 및 사용

필요 환경: repository의 기존 Python test dependencies(NumPy/SciPy/pytest/mpmath), GNU Fortran. Python orchestration은 native 생성 상태를 재사용하며 CompOSE IO도 native backend를 호출한다. 새 physics generator를 Python으로 대체하지 않았다.

```bash
python scripts/audit_core_support_20261117.py --output data/core_support_20261117/domain_reproduction_NEW
python scripts/extend_core_support_20261117.py --output data/core_support_20261117/reproduction_NEW
PYTHONPATH=python_frontend python -m pytest -q python_frontend/tests/test_core_extension.py python_frontend/tests/test_core_support_domain_audit.py
```

확장 실행기는 기존 output directory를 거절한다. 재현할 때 새로운 directory를 사용한다. Domain audit도 위 명령처럼 새로운 output 경로를 지정해 기존 실행 근거를 보존한다. Test dependencies는 pyproject.toml의 test extra를 따른다.

연속 EOS 사용:

```python
import json
from compose_eos.whole_crust import WholeCrustBarotrope
with open("data/core_support_20261117/primary_v1/REPRESENTATION.json") as f:
    eos = WholeCrustBarotrope.from_contract(json.load(f))
state = eos.state_at_pressure(100.0)  # MeV fm^-3
```

## 다음 실행: C, 이후 D

C의 정식 목적은 **TF와 나머지 core/method의 필요한 확대**다. 기존 DD2 reference seed의 실제 n_B/Ye/Ymu와 유한 TF 해를 출발점으로 continuation을 수행한다. 제공되지 않은 meson fields, cell radius, cluster 의미를 만들어 채우지 않는다. 허용 geometry별 에너지·잔차·mesh/cell-size 수렴과 phase 선택을 비교하고 실제 지원이 있는 나머지 조합만 추가한다. Requested/eligible/attempted/admitted/diagnostic/unsupported census는 구분한다.

C 이후 D는 APR retained 6→8→7 Gibbs boundary를 소비자까지 보존하고, 실제 whole-crust/core 범위의 TOV와 명시한 perturbation closure의 Cowling 관측량을 계산하는 단계다. B 자체는 D의 최대중심밀도 필요를 선결정하지 않는다.

구체적 다음 입력과 종료조건은 `data/core_support_20261117/NEXT_STEP_C.json`에 있다. B의 게시 승인 후 현재 branch의 이 결과를 게시하며, C는 다음 작업으로 남겨 둔다.
