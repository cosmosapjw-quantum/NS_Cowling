# 조밀한 직접 생성 EOS와 보간법의 정량 교차검증

**대상:** NS_Cowling, 2026-11-18. 기준 Git commit: a0fa5fb93851a9647030f3589c24a62ddea527b1. 결과 branch: research/native-dense-interpolation-benchmark-20261118.

## 핵심 결론

이전 H 단계의 보간식 간 비교에서 한 걸음 더 나아가, 동일한 물리 모형의 방정식을 Fortran native 계산기로 밀도 격자에서 재계산했어. 학습에 사용하지 않은 실제 계산점으로 보간 오차와 미분 관계를 검증한 결과, **균질 RMF에서는 Rational·Exponential, core에서는 Hermite, TF crust에서는 Bernstein·Hermite**가 수용률과 정확도 사이의 좋은 균형을 보였어. 이는 특정 범위의 유한 계산값이며 보편적인 최적 보간법의 증명은 아니야.

원래 EOS의 phase와 물리 지원 범위를 고정했으며, APR의 Gibbs 6→8→7 구간은 각 phase 내에서만 보간했어. 도표 5개와 아래 순위표에서는 요청대로 linear를 제외했어.

## 직접 계산한 reference의 범위

| 범주 | 모델 | 새 native 상태 |
|---|---|---:|
| 균질 RMF | NL3, PK1, TM1, DDME-X 두 매개변수 입력 | 5 × 4,097 = 20,485 |
| Core | APR, BBB-BHF2, GM1, KDE0v1, SLy4 | APR 1,027 + 4 × 1,025 = 5,127 |
| 합계 | 서로 구별되는 원 입력 10개 | **25,612** |

RMF 입력 밀도는 0.015–1.0 fm⁻³, core 입력은 해당 source의 지원 범위야. BBB-BHF2는 원전 microscopic support 상한 0.8 fm⁻³을 넘기지 않았고, APR의 phase 7은 고밀도에서 자체적인 인과성 제한이 나타났어.
TF 5개 crust 모델에서는 기존 N128 정지해를 학습점으로 고정하고 새 밀도점 57개를 요청했어. 최초 26개와 양쪽 수용 profile 재초기화 후 추가 7개가 수용됐지만, 이 중 3개는 기존 훈련점에서의 반복 검산이었기 때문에 **독립 홀드아웃은 30개**야. 나머지 24개 실패는 원래 수용 조건을 완화하지 않은 채 보존했어.

## 지표·방법별 물리적 검증

사용한 기본 관계는 nB=dP/dμB, ε=μB·nB−P, cₛ²/c²=nB/(μB·dnB/dμB)야. 압력만 맞추는 보간법도 미분하여 재구성한 nB가 직접 계산한 native 값을 놓칠 수 있기 때문에, 밀도·압력·에너지/바리온·음속 오차를 각각 비교했어.
표의 RMS는 보간에 사용하지 않은 native 점에서 평가했어. 상대오차는 ppm(10⁻⁶), 에너지 오차는 MeV/baryon, 음속제곱 오차는 cₛ²/c²의 절대값이야. Coverage와 native-causal 물리 수용률이 **둘 다 99.5% 이상**인 방법만 전체 순위에 들어가며, 수용된 일부 점의 RMS가 작더라도 수용률이 낮으면 순위에서 제외해.

### 균질 RMF 5개 입력

| 방법 | n RMS [ppm] | P RMS [ppm] | ε/n RMS [MeV] | cₛ² RMS | Coverage | 물리 수용 | 순위 가능 |
|---|---:|---:|---:|---:|---:|---:|---|
| Cubic Hermite P(μ) | 0.3792 | 0.01039 | 0.0003595 | 6.275e-06 | 100.00% | 100.00% | 예 |
| Exponential n/μ | 0.3206 | 0.008287 | 0.0003041 | 5.242e-06 | 100.00% | 100.00% | 예 |
| Rational n/μ | 0.3089 | 0.007817 | 0.0002928 | 4.549e-06 | 100.00% | 100.00% | 예 |
| Bernstein (slope) | 0.7062 | 0.01925 | 0.0006678 | 1.028e-05 | 100.00% | 100.00% | 예 |
| Power-law | 1.689 | 0.0383 | 0.001623 | 3.273e-05 | 11.19% | 11.19% | 아니오 |
| PCHIP | 60.72 | 1.186 | 0.05792 | 0.001758 | 100.00% | 100.00% | 예 |
| Akima | 32.32 | 0.3978 | 0.03339 | 0.003421 | 100.00% | 100.00% | 예 |
| Natural cubic | 396 | 12.92 | 0.4026 | 0.008271 | 100.00% | 99.84% | 예 |

**현재 수용 기준에서의 1위:** Rational n/μ (n RMS 0.3089 ppm).

### Core EOS 5개 입력

| 방법 | n RMS [ppm] | P RMS [ppm] | ε/n RMS [MeV] | cₛ² RMS | Coverage | 물리 수용 | 순위 가능 |
|---|---:|---:|---:|---:|---:|---:|---|
| Cubic Hermite P(μ) | 11.49 | 0.1237 | 0.01243 | 0.0008702 | 100.00% | 100.00% | 예 |
| Exponential n/μ | 10.68 | 0.1079 | 0.01166 | 0.0008458 | 98.83% | 99.96% | 아니오 |
| Rational n/μ | 10.41 | 0.102 | 0.01143 | 0.0008442 | 98.83% | 99.96% | 아니오 |
| Bernstein (slope) | 11.81 | 0.1389 | 0.01289 | 0.001147 | 98.80% | 99.93% | 아니오 |
| Power-law | 5.202 | 0.0598 | 0.005061 | 0.0001309 | 16.68% | 16.87% | 아니오 |
| PCHIP | 33.73 | 0.4977 | 0.03911 | 0.002156 | 100.00% | 99.99% | 예 |
| Akima | 61.76 | 0.9263 | 0.06969 | 0.003974 | 100.00% | 99.99% | 예 |
| Natural cubic | 272.2 | 5.598 | 0.3679 | 0.01243 | 100.00% | 99.71% | 예 |

**현재 수용 기준에서의 1위:** Cubic Hermite P(μ) (n RMS 11.49 ppm).

### Spatial TF N128 crust 5개

| 방법 | n RMS [ppm] | P RMS [ppm] | ε/n RMS [MeV] | cₛ² RMS | Coverage | 물리 수용 | 순위 가능 |
|---|---:|---:|---:|---:|---:|---:|---|
| Cubic Hermite P(μ) | 74.42 | 1.573 | 0.07017 | 6.174e-06 | 100.00% | 100.00% | 예 |
| Exponential n/μ | 77.04 | 1.639 | 0.07264 | 6.026e-06 | 100.00% | 100.00% | 예 |
| Rational n/μ | 78.35 | 1.672 | 0.07387 | 5.952e-06 | 100.00% | 100.00% | 예 |
| Bernstein (slope) | 70.67 | 1.742 | 0.06662 | 2.318e-05 | 100.00% | 100.00% | 예 |
| Power-law | 41.27 | 3.191 | 0.03891 | 1.949e-05 | 63.33% | 63.33% | 아니오 |
| PCHIP | 1297 | 61.74 | 1.223 | 0.0004562 | 100.00% | 100.00% | 예 |
| Akima | 1505 | 70.11 | 1.419 | 0.0003161 | 100.00% | 100.00% | 예 |
| Natural cubic | 1795 | 195.4 | 1.692 | 0.0006808 | 100.00% | 100.00% | 예 |

**현재 수용 기준에서의 1위:** Bernstein (slope) (n RMS 70.67 ppm).

## 모델별 보간법 결과

| 직접 생성한 모델 | 기준에서 가장 낮은 n RMS를 보인 수용 방법 | n RMS [ppm] | PCHIP 대조 [ppm] |
|---|---|---:|---:|
| rmf/NL3_PROVIDER_PRECISION | Rational n/μ | 0.3426 | 65.57 |
| rmf/PK1_CANONICAL | Rational n/μ | 0.3981 | 60.72 |
| rmf/TM1_CANONICAL | Rational n/μ | 0.3999 | 59.52 |
| rmf/DDMEX_PROVIDER_ROUNDED_N0 | Rational n/μ | 0.1469 | 58.76 |
| rmf/DDMEX_ORIGINAL_TABLE_RHO0 | Rational n/μ | 0.14 | 58.77 |
| core/APR_APR | Cubic Hermite P(μ) | 0.5629 | 54.9 |
| core/BBB_BHF2 | Rational n/μ | 23.12 | 36.12 |
| core/GM_GM1 | Rational n/μ | 0.2744 | 21.67 |
| core/RG_KDE0v1 | Rational n/μ | 0.6917 | 21.09 |
| core/RG_SLY4 | Rational n/μ | 0.5096 | 21.39 |
| tf/DDME2 | Power-law | 80.53 | 3201 |
| tf/DDMEX | Rational n/μ | 3.124 | 113.7 |
| tf/NL3 | Cubic Hermite P(μ) | 0.304 | 32.73 |
| tf/NL3wrL55 | Bernstein (slope) | 1.673 | 105.2 |
| tf/TM1e | Bernstein (slope) | 42.38 | 1366 |

## TF crust의 native 계산과 제한

| 모델 | 신규 계산 요청 | 최초 수용 | 재초기화 후 추가 수용 | 독립 native 홀드아웃 |
|---|---:|---:|---:|---:|
| DDME2 | 15 | 5 | 2 | 4 |
| DDMEX | 12 | 5 | 2 | 7 |
| NL3 | 8 | 8 | 0 | 8 |
| NL3wrL55 | 12 | 6 | 0 | 6 |
| TM1e | 10 | 2 | 3 | 5 |

이번 crust 결과는 **Wigner–Seitz N128 droplet 하나의 stationary profile**에서 확인한 값이야. Rod/slab/tube/bubble 간 자유에너지 비교, geometry 전역 최소, 물리적 crust–core 전이의 존재 여부 또는 연속체 limit를 확정하지 않았어. DDME2의 수용된 power-law 결과는 신규 점이 4개뿐이며 좌측 endpoint 음속의 인공적 특징 때문에 기본법으로 승격하지 않았어. 실패 native 상태는 보간값으로 채우지 않았어.

## 물리적으로 따로 분리한 실패

- APR phase 7은 약 n_B>=0.870 fm⁻³에서 직접 native 해 자체가 acausal이므로 보간법에 그 원인을 전가할 수 없어.
- BBB-BHF2 생성기는 원 논문의 에너지 표를 재구성한 kernel로부터 계산한 결과이며, 원래 microscopic BHF many-body 문제를 다시 풀었다는 뜻은 아니야.
- PCHIP/Akima는 압력 데이터에는 적합할 수 있어도 dP/dμ에서 복원한 baryon density가 원 실제 계산과 달라질 수 있어.
- Natural cubic spline은 압력의 매끄러움만으로 안정성·인과성이 보장되지 않아. 보고된 원자료 점의 유한차분 cₛ²와 interpolation 자체의 미분 오차도 분리해야 해.
- Sharp Maxwell density jump를 평활화하지 않았고, 서로 다른 interpolation으로 교차근을 새로 찾아 전이 존재성의 변화를 측정한 연구는 이번 범위 밖이야.
- 적은 TF 홀드아웃과 source 자체의 매개변수 오차가 있는 상황에서 작은 RMS 차이를 모델 간 물리적 우열로 해석하지 않아.

## 도표 (linear는 표시하지 않음)

1. [RMF_NATIVE_HOLDOUT_NO_LINEAR](figures/RMF_NATIVE_HOLDOUT_NO_LINEAR.png)
2. [CORE_NATIVE_HOLDOUT_NO_LINEAR](figures/CORE_NATIVE_HOLDOUT_NO_LINEAR.png)
3. [TF_NATIVE_HOLDOUT_NO_LINEAR](figures/TF_NATIVE_HOLDOUT_NO_LINEAR.png)
4. [NATIVE_HOLDOUT_COVERAGE_NO_LINEAR](figures/NATIVE_HOLDOUT_COVERAGE_NO_LINEAR.png)
5. [NATIVE_MODEL_METHOD_HEATMAP_NO_LINEAR](figures/NATIVE_MODEL_METHOD_HEATMAP_NO_LINEAR.png)

## 실제 검증 상태

- 독립 source/readback/phase/holdout 검사: **1096건, 실패 0건**.
- 전체 frontend/pipeline: **1,501 passed, 9 skipped, 4 subtests passed, 실제 exit 0**.
- 소스 파일 208개는 suite 전후 byte identity 유지.
- 직접 native 25,612행의 재계산 Euler 항등식 상대 잔차 최댓값: **5.7e-14**.
- 실패 로그와 생성기·컴파일러·표·실행파일 해시는 기계 판독 가능한 JSON과 native 상태 파일에 보존했어. 검증되지 않은 각종 EOS 물리적 가설은 별도로 unresolved로 남겼어.

## 다음 연구 순서 (아직 실행하지 않음)

동일한 원 crust/core 두 branch를 **각기 다른 보간법으로 근사해 Maxwell P₋(μ)-P₊(μ)=0의 교차근을 독립적으로 재탐색**하고, 존재 여부·전이압력·밀도점프·잠열을 dense native 기준으로 직접 비교해야 해. 이번에 고정한 H Maxwell doublet와 주변 보간을 바꾼 결과를 혼합해서는 안 돼.

## 주요 문헌

- Haensel & Potekhin (2004): https://www.ioffe.ru/astro/Stars/Paper/hp04.pdf
- Fortin et al. (2016): https://doi.org/10.1103/PhysRevC.94.035804
- Ferreira & Providência (2020): https://www.mdpi.com/2218-1997/6/11/220
- Verma & Mallick (2025): https://arxiv.org/html/2501.06453
