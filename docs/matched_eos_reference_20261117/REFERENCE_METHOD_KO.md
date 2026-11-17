# 생성 EOS와 공개 표의 물리량 비교

이 비교는 현재 생성된 유한 밀도 구간의 EOS가 공개 표와 얼마나, 어떤 방식으로 다른지 측정한다.
공개 표와의 수치 일치를 생성 EOS의 공통 합격조건으로 정하지 않는다. 대신 표현상의 차이,
명시된 모델 선택, 계산으로 설명된 차이, 아직 설명되지 않은 잔차를 구분한다.
실제 matching 성립 및 native admission은 별도 campaign 결과로 판정한다.

## 입력과 비교 범위

`scripts/compare_generated_reference.py`는 생산 reader와 독립적인 한정된 raw CompOSE
reader를 사용한다. 단일 T=0, Yq slice만 허용하고 파일의 밀도축 시작 인덱스를 실제로 적용한다.
각 thermo/composition 행의 인덱스, 필수 Q1–Q7, 중복·결측·비유한 값, composition의 존재 여부를 검사한다.
일부 cold provider가 T축을 0으로 선언하고 모든 행을 1로 표기한 경우에는 유일한 scalar
label임을 확인하고 선언값과 실제 행 label을 둘 다 기록한다. 서로 다른 행 label을 섞은 경우는 거절한다.

| 모델 | native phase | 공개 phase | 요청 밀도 범위 / fm⁻³ | 비교 조성 |
|---|---:|---:|---|---|
| GM_GM1 | 6 | 3 | 0.08–0.5 | homogeneous free proton fraction |
| RG_SLY4 | 6 | 3 | 0.08–0.5 | proton, electron, muon fractions |
| SPG M1–M5 | 1 | 7 | 0.002 이상, 실제 공통 지원까지 | electron, muon fractions |

서로 다른 phase ID 자체가 서로 다른 물리상을 뜻하지는 않는다. 위 mapping은 제공자 phase legend와
생성 request에 근거한다. 동일 phase ID가 비연속 구간에서 재등장해도 중간의 다른 phase를 가로질러
보간하지 않는다. 정확한 비교 밀도는 요청 범위, 두 표의 실제 phase 지원의 교집합으로 제한한다.

SPG 공개본 하부의 BSk22 segment는 생성 spherical CLD와 다른 모델이며 이 비교에 포함하지 않는다.
공개 nuclear group tuple의 비표준 의미를 생성 표의 표준 cluster number fraction과 동일시하지 않는다.
이번 비교는 그룹 A, Z와 그룹 abundance를 사용하지 않는다. 제공되지 않은 particle fraction은
결측으로 유지하며, 결측 구간을 가로지르는 조성 보간을 하지 않는다.

## 표현을 복원한 물리량

각 파일 header의 neutron reference rest energy를 \(m_n^{\rm ref}c^2\)라 하면,

\[
P=n_B Q_1,\qquad
\epsilon=n_B m_n^{\rm ref}c^2(1+Q_7),\qquad
\mu_B=m_n^{\rm ref}c^2(1+Q_3).
\]

\(n_B\)는 fm⁻³, \(P,\epsilon\)은 MeV fm⁻³, \(\mu_B\)는 MeV다.
\(\epsilon\)에는 rest energy가 포함된다. `mu_b`와
\(h_B=(\epsilon+P)/n_B\)를 별도 관측량으로 비교한다. 공개 Q5가 0이 아닌 경우를 포함하여
Q3를 조용히 Euler enthalpy로 대체하지 않는다. \(m_n^{\rm ref}\)가 다르다는 사실만으로
실제 nuclear kernel의 질량 convention이 같거나 다르다고 단정할 수 없다.

`energy_per_baryon`은 \(\epsilon/n_B\)다. 출력을 맞추기 위한 \(\epsilon\to\epsilon+C n_B\)
보정이나 평균 차이 제거는 하지 않는다. 서로 다른 header로 같은 물리량을 재인코딩한 analytic
control은 차이가 0으로 복원되어야 한다.

공개 표의 실제 밀도 node에서 생성 표를 PCHIP로 보간한다. 같은 node의 linear interpolation도
계산하여 표현 선택의 민감도를 기록한다. 둘의 차이는 엄밀한 오차 상한이 아니다. 상대차는
해당 참조 관측량 최대 절댓값의 10⁻¹²보다 작은 분모에서 정의하지 않고 절대차를 유지한다.

## 압력 차이의 에너지 미분 진단

추가 진단은 같은 연속 phase에서

\[
P_E(n_B)=n_B^2\frac{d}{dn_B}\left(\frac{\epsilon(n_B)}{n_B}\right)
\]

를 CubicSpline과 PCHIP 두 표현으로 계산한다. 두 방법 모두 각 phase의 처음·마지막 두 node를
미분 비교에서 제외한다. 상전이 경계나 이웃 phase를 통과하는 stencil을 사용하지 않는다.
각 표의 실제 제공 \(P\)와 \(P_E\)의 차이, 그리고

\[
(P_{\rm gen}-P_{\rm ref})-(P_{E,\rm ref}-P_{\rm ref})
=P_{\rm gen}-P_{E,\rm ref}
\]

를 따로 남긴다. 공개 ε trajectory와 공개 P가 다르게 구성되어 있다면 그 차이가 생성/공개 P
잔차를 얼마나 설명하는지 직접 측정할 수 있다. 이는 보존 표에 대한 수치적 분해이며,
공개 저자의 내부 구현 절차나 오류 원인을 확인했다는 뜻은 아니다. 유한 표의 미분 오차도 남는다.
검증용 analytic control에서는 ε는 같은 채 제공 P만 20% 바꾼 표를 만들어 위 분해가 정확히
그 차이를 회수하는지 확인한다.

## 이미 계산된 모델 차이 근거

GM1은 원 논문의 rounded coupling ratios와 공통 939 MeV nuclear mass, npe composition을
사용한 reconstruction이다. 공개 thermo header 939.56533/938.272 MeV만으로 저자의 nuclear
kernel mass를 확정할 수 없다. 기존 포화물성 계산과 원 밀도 grid 비교는
[`reference_cores.md`](../research_generator_20261116/reference_cores.md)에 있다.
새 잔차에 단순히 “질량 차이”라는 설명을 붙이지 않고, 별도 controlled calculation이 재현한
영향만 수치적으로 설명된 것으로 분류한다.

SLy4는 \(\alpha=1/6\), 발표된 rounded coefficients, fitted kinetic coefficient
20.73553 MeV fm²와 별도 neutron/proton rest masses, npemu equilibrium을 사용한다.
과거 mass·kinetic·higher-precision constants·lepton \(\hbar c\)·source composition sensitivity
계산은 저밀도 공개 pressure 잔차의 대부분을 설명하지 못했다.
그 결과와 source energy derivative 검산은
[`skyrme_implementation.md`](../research_generator_20261116/skyrme_implementation.md)에 보존되어 있다.
이 근거를 재사용하며 동일한 sensitivity 계산을 새 승인 gate로 반복하지 않는다.

SPG cohort의 variant는 `public_composition_inferred_20261116`이다.
사용한 \(\sigma_0,b\)는 공개 조성에서 추론한 candidate inputs이며 원저자 핵질량 fit의 복구가 아니다.
`UNRESOLVED_SOURCE_BETA_SURFACE_DERIVATIVE_PATTERN`과 `original_mass_calibration=UNRESOLVED`는
현재 생성 receipt의 실제 상태다. 별도 AME2408 candidate calibration과 이 cohort를 혼동하지 않는다.
모델 정의·calibration 선택과 공개 표의 미분압력 차이를 각각 기록한다.

## 실행과 출력

최신 다운로드한 reference가 준비되면 다음처럼 실행한다.

```sh
python scripts/compare_generated_reference.py \
  --reference-root data/matched_eos_reference_20261117/references \
  --retained-reference-root eos_generator/eos_data \
  --output data/matched_eos_reference_20261117/comparison \
  --plots
```

`--reference-map`은 model별 raw CompOSE directory를 지정하는 JSON mapping이다.
출력 경로는 새 디렉터리나 비어 있는 디렉터리여야 한다.

`COMPARISON.json`에는 각 입력의 axis/row labels, reference masses, 파일 SHA-256, actual phase
지원, 모든 관측량의 최대·중간·RMS 차이와 위치, 보간·미분 민감도, Python/NumPy/SciPy 버전을
기록한다. 모델별 CSV는 비교한 모든 node와 결측값을 보존한다. 그림은 출처나 원인 판정을
대신하지 않는다. `--retained-reference-root`를 지정하면 과거 git reference와 최신 reference의
byte identity 및 물리량 차이도 별도 subsection/CSV로 기록한다.

## 원 제공자와 재사용 근거

- [CompOSE GM1](https://compose.obspm.fr/eos/54), [cold SLy4](https://compose.obspm.fr/eos/134).
- [SPG M1 crust](https://compose.obspm.fr/eos/311), [Scurto et al. 원 논문](https://arxiv.org/html/2402.15548v1), [저자 bulk dataset](https://zenodo.org/records/11084025).
- [Glendenning–Moszkowski 원 manuscript](https://escholarship.org/uc/item/6sw424k1), [SLy4 parameter table](https://arxiv.org/abs/astro-ph/0006135), [Gulminelli–Raduta unified framework](https://arxiv.org/abs/1504.04493).
- [`CALIBRATION_SOURCE_CONTRACT_KO.md`](../research_finish_20261116/CALIBRATION_SOURCE_CONTRACT_KO.md): 원저자 calibration, 공개 bulk 자료, 별도 candidate 및 lower segment의 구분.

최종 원인 판정은 현재 reference intake와 실제 `COMPARISON.json`을 함께 사용한다.
설명되지 않은 차이는 크기·밀도 범위와 함께 `UNRESOLVED_SOURCE_DIFFERENCE`로 유지한다.

## 2026-11-17 최신 공개 자료로 실행한 결과

최초 분석은 Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0으로 실행했다.
게시한 동일 명령의 재현 산출물은 Python 3.12.13, NumPy 2.5.4, SciPy 1.18.1에서 생성했고,
그 실제 환경과 실행 상태는 `comparison/COMPARISON.json`과 `validation/EXECUTION.json`에 있다.
이 절의 수치는 새 다운로드의 immutable bytes와 현재 sealed 생성 표의 비교다.
`comparison/COMPARISON.json`의 파일 hashes, 실제 비교 구간과 node별 CSV가 재현 근거다.

| 모델 | 비교 node | 최대 \(|\Delta P/P_{ref}|\) | 최대 \(|\Delta\epsilon/\epsilon_{ref}|\) | 최대 \(|\Delta\mu_B|\) / MeV |
|---|---:|---:|---:|---:|
| GM_GM1 | 132 | 0.4863865% | 0.0276846% | 0.782798 |
| RG_SLY4 | 20 | 0.1827865% | 0.0042571% | 0.0743081 |
| SPG_M1 | 148 | 15.6648311% | 0.0000173% | 0.193478 |
| SPG_M2 | 160 | 18.6600169% | 0.0000312% | 0.243184 |
| SPG_M3 | 149 | 15.3726822% | 0.0000172% | 0.188483 |
| SPG_M4 | 158 | 16.5083026% | 0.0000233% | 0.207226 |
| SPG_M5 | 146 | 16.9184899% | 0.0000199% | 0.207797 |

각 관측량의 최댓값이 나타나는 밀도는 서로 다를 수 있다. 압력 최대 상대차는 GM1의
0.0829066667 fm⁻³, SLy4의 0.08317638 fm⁻³, SPG 각 공개 inner-crust의 첫 node에서 발생한다.
최댓값만으로 전체 밀도 범위의 차이를 대표하지 않는다.

### 과거 reference와 최신 reference의 차이

과거 git reference와 최신 공개 archive의 byte identity가 달라도 물리량이 달라졌다고 단정할 수 없다.
GM1의 retained/current 압력 최대 상대차는 1.03×10⁻¹⁰, ε 차이는 5.64×10⁻¹²다.
SPG M1–M5의 thermo 물리량은 비교 node에서 0 또는 약 10⁻¹⁶ 상대차다.
이 경우는 인코딩·출력 정밀도의 차이와 물리량 동등성을 구분한다.

반면 **SLy4 retained/current 압력 차이는 6.71857%**다. 과거 git 표를 비교 기준으로 쓰면
현재 생성물의 최대 압력 차이가 6.12431%였지만, 최신 공개 표를 사용하면 0.1827865%다.
따라서 과거 source-fidelity 문서의 큰 저밀도 잔차를 현재 공개본의 문제로 재사용해서는 안 된다.
이것이 provider revision인지 저장소 내 과거 변환인지의 원인 귀속은 별도 provenance 조사에 맡긴다.
이번 비교가 직접 확인한 것은 두 reference artifact의 물리량 차이와 그 영향이다.

### SPG pressure 차이의 제한된 수치적 설명

SPG의 ε 차이는 매우 작지만 압력 차이는 크다. 동일 표의 에너지 미분 진단으로 이를 분해했다.
아래 값은 각 phase 양끝 두 node를 제외한 내부 범위에서의 최대 상대차다.

| 모델 | 내부 공개 node | 공개 \(|P_E-P|/|P|\) 최대 | \(|P_{gen}-P_{E,ref}|/|P_{ref}|\) 최대 |
|---|---:|---:|---:|
| SPG_M1 | 144 | 15.17948% | 0.0162320% |
| SPG_M2 | 156 | 18.11425% | 0.0394070% |
| SPG_M3 | 145 | 14.89433% | 0.0214368% |
| SPG_M4 | 154 | 16.00810% | 0.0318997% |
| SPG_M5 | 142 | 16.40530% | 0.0238185% |

이 내부 구간에서 공개 표의 \(P_E-P\)가 생성/공개 압력 잔차 대부분을 수치적으로 설명한다.
CubicSpline/PCHIP 간 해당 미분압력 민감도는 공개 압력으로 나누면 최대 0.00315–0.00524%다.
생성 ε의 미분압력과 생성 P의 상대차는 비교 가능한 node에서 최대 3.03×10⁻⁵다.
이 결과는 공개 압력을 그대로 목표로 forcing/retuning할 근거를 제공하지 않는다.
현재 생성물의 열역학 미분을 유지하면서 차이를 데이터 수준에서 설명할 수 있다.

이 수치적 분해가 원저자 pressure code의 특정 오류나 derivative 생략을 증명하지는 않는다.
기존 receipt의 beta/surface derivative 원인 미확정 상태와 surface calibration의 출처 상태를 유지한다.
또한 공개 조성에서 parameter를 추론한 variant이므로 SPG ε의 높은 일치를 독립 blind prediction의
성공률로 해석하지 않는다.

### 차이 ledger

| 항목 | 현재 근거 상태 | 판정 |
|---|---|---|
| header reference-mass 차이 | implementation-verified analytic control | 물리량 복원으로 제거; 물리 energy offset을 fit하지 않음 |
| GM1 압력 약 0.49% 잔차 | numerically checked / attribution unresolved | 선언된 kernel/parameter 선택의 영향은 별도 controlled 계산 근거와 함께만 귀속 |
| 과거 SLy4 큰 잔차 | implementation-verified artifact comparison | 현재 공개 표에서 대부분 해소; retained/current 차이 원인은 별도 조사 |
| 현재 SLy4 압력 약 0.18% 잔차 | numerically checked controlled explanation | 입자별 kinetic 계수 control에서 0.01127%로 감소; 정확한 provider 선택 귀속과 남는 잔차는 미확정 |
| SPG 큰 pressure 차이 | numerically checked decomposition | 내부 범위에서 공개 ε 미분압력과 제공 P의 차이가 대부분 설명 |
| SPG 원저자 pressure construction | unresolved | 위 데이터 분해를 내부 kernel 원인 확정으로 승격하지 않음 |
| SPG calibration | established receipt identity / original-source unresolved | 공개 조성 추론 candidate로 명명; 원저자 mass-fit 복구 주장 없음 |
| composition group semantics | established scope exclusion | 확실한 particle fractions만 비교, 누락 값은 결측 |

analytic regression 9개가 실제 통과했다. tests는 reference-mass 재인코딩, 임의 axis origins,
scalar T alias 및 혼합 label 거절, phase 경계, 결측 조성 보존, Q3/enthalpy 구분,
0 분모, 미분의 내부 지원, 알려진 pressure mismatch 분해를 검사한다.
최초 pytest 부재는 격리된 `/tmp/ns_cowling_comparison_pydeps` 설치로 해결했으며 생산 환경을 변경하지 않았다.

후속 원밀도점 계산은 [`DIFFERENCE_DIAGNOSIS_KO.md`](DIFFERENCE_DIAGNOSIS_KO.md)에 있다.
현재 SLy4의 공통 kinetic coefficient를 입자별 질량으로 계산하는 control이 작은 잔차의 대부분을
재현했다. 이는 생산 functional을 변경한 결과가 아니다. GM1의 header-mass 및 printed-precision
control은 전체 차이를 설명하지 못했으므로 그 원인 미확정 상태는 유지한다.

원격 저장소의 `git log --follow`에서 보존 SLy4 thermo 파일은
`5d56f1bca333abbf9cc1f7500e69fde9f7791ac2`(2026-11-06, baseline import)에 도입됐으며,
이후 그 경로의 수정 이력은 없었다. 이 사실만으로 import 이전의 변환과 제공자 revision을
구분할 수는 없다. 새 분석의 reference authority는 이번 다운로드와 그 SHA-256이다.
