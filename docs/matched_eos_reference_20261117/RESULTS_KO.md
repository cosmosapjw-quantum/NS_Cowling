# 실제 core–crust EOS 생성과 공개자료 진단 결과

기준일: 2026-11-17. 시작점은 local 작업이 끝난 PR16의
`32074d48793e265ab21f7b7db0553f6b3d1283dc`다.
작업 branch는 `research/matched-eos-reference-20261117`이며 기존 PR16과 종료된 one-shot
CLAIM은 보존했다. 이번 사용자 요청으로 별도 campaign을 실행했다.

## 1. 이번에 실제로 확보한 것

**SPG M1–M5 × GM1·SLy4의 10개 모델 쌍에 대해 20개 matched EOS 표를 확보했다.**
이번에 새로 matching한 표는 19개이고, 기존 SPG M1–GM1 CSS 한 개는 같은 입력·방법·출력
identity를 확인하고 현재 reader, native consistency, 독립 CSS admission으로 재검증했다.
실행 요청은 30개이며 신규 native kernel 호출은 29회다.

| 방법 | 요청 | 새 kernel 호출 | 실제 결과 | 검증 범위 |
|---|---:|---:|---|---|
| quintic_P | 10 | 10 | 10 CONSTRUCTED_VERIFIED | native construction, 수치·메타데이터 전체 readback, serialized consistency |
| CSS | 10 | 9 | 9 ADMITTED + 1 REUSED_ADMITTED | 위 검증과 독립 CSS interface admission |
| 직접 Maxwell | 10 | 10 | 10 REJECTED | 실제 native 호출 후 현재 지원구간의 실패 원인 보존 |

EOS 파일은 `data/matched_eos_reference_20261117/deliverable/matched_eos_20.zip`에 있다.
CATALOG.json/CSV에는 성공한 표와 거절된 경우를 모두 기록했다. ZIP 안의 EOS 100파일을
원 output inventory와 비교하고 다시 압축 해제해 SHA-256을 확인했다.

```text
ZIP bytes: 574088
ZIP SHA256: bf57c864bd29785442f6f002a3a63e952c354b4eca1add47ff076ddea10f8401
```

이 표들은 **유한 밀도 지원의 inner-crust+core, thermo-only EOS**다. Core 상단은
`n_B=0.5 fm^-3`이고 표면/outer-crust 연결은 아직 포함하지 않았다. 원 source의 composition은
보존되지만 matched barotrope는 명시적으로 abundance를 제거한 출력이며 결측을 물리적 0으로
해석하면 안 된다. 현재 결과는 whole-star EOS나 TOV/Cowling 계산의 완료를 뜻하지 않는다.

SPG의 variant는 `public_composition_inferred_20261116`이다. 원저자 mass calibration과의
동일성은 확인되지 않았다. SPG와 GM1/SLy4의 조합은 서로 다른 functional을 연결한 후보이므로
각 public unified EOS의 재현으로 명명하지 않는다. Smooth/CSS bridge 자체도 명시적 모델 선택이다.

## 2. Matching 결과의 물리적 해석

직접 Maxwell은 동일한 baryon chemical potential에서 두 source의 압력이 같은 점을 찾는다.
현재 입력에서 8개는 공통 μ_B 구간에 압력 교점이 없었고, SPG M2–SLy4 및 M5–SLy4는
μ_B 지원 자체가 겹치지 않았다. 이는 현재 유한 source 지원에 대한 계산 결과다.
다른 밀도 지원이나 모델에서도 Maxwell이 불가능하다는 명제로 확대하지 않는다.

CSS는 crust 쪽 density jump와 constant-sound-speed bridge를 구성한다. 여기서
`chi=c_s^2/c^2=dP/dε`이며 ε는 rest energy를 포함한다. Jump에서

\[
P_+=P_-,\qquad \mu_{B,+}=\mu_{B,-},\qquad
\Delta\epsilon=\mu_B\Delta n_B
\]

를 검사한다. 이 인공 bridge가 원래 두 미시모델의 직접 Maxwell 교점이라는 뜻은 아니다.
두 경우를 서로 다른 method로 보존한다.

| Crust | GM1 CSS chi | GM1 Δε / MeV fm⁻³ | GM1 Seidov ratio | SLy4 CSS chi | SLy4 Δε / MeV fm⁻³ |
|---|---:|---:|---:|---:|---:|
| SPG M1 | 0.113806 | 42.3559 | 1.14343 | 0.0121184 | 0.117757 |
| SPG M2 | 0.142747 | 28.5627 | 0.664680 | 0.0291520 | 0.720017 |
| SPG M3 | 0.121030 | 40.4666 | 1.02907 | 0.0129895 | 1.19307 |
| SPG M4 | 0.172912 | 25.5810 | 0.535495 | 0.0314492 | 0.188671 |
| SPG M5 | 0.101297 | 36.7208 | 1.04415 | 0.0200664 | 0.505427 |

GM1에 붙인 M1·M3·M5의 `seidov_small_core_stable=false`를 그대로 보존했다. Seidov 기준은
작은 high-density core가 발생하는 특정 sharp transition의 국소 안정성 진단이다.
이 flag가 true여도 별 전체의 안정성이나 radial/nonradial mode 안정성을 판정한 것은 아니다.
현재는 outer-crust/표면과 충분한 고밀도 지원을 가진 별 모형도 아직 계산하지 않았다.

10개 smooth 출력의 `diagnostic_only`는 모두 false였다. Smooth의 현재 검증 범위를
독립 CSS/Maxwell admission과 같은 것으로 표기하지 않았다.

## 3. 최신 공개본을 별도로 확보한 이유와 실제 identity

CompOSE의 실제 download href를 따라 GM1, cold SLy4, SPG M1–M5의 작은 archive를 받았다.
DD2 crust와 unified는 아래 고정 Zenodo 버전을 받고 게시 MD5와 실제 bytes를 대조했다.

| Reference | 원 제공자 | 고정 archive/판정 |
|---|---|---|
| GM1 | https://compose.obspm.fr/eos/54 | cold npe, core 비교만 수행 |
| SLy4 | https://compose.obspm.fr/eos/134 | cold npemuN, phase 3의 homogeneous 구간 |
| SPG M1–M5 crust | https://compose.obspm.fr/eos/311 , /316 , /317 , /318 , /319 | lower BSk22를 제외한 공개 CLD phase 7 |
| DD2 crust | https://compose.obspm.fr/eos/210 | https://zenodo.org/records/14508936 , v1, MD5 9a70c2aa57bcf415fde43b46c533f763 |
| DD2 unified | https://compose.obspm.fr/eos/217 | https://zenodo.org/records/14808942 , v2, MD5 d8e258273081fe73ac4ca665756a4de5 |

각 source의 URL, 받은 시각, 원 archive, SHA-256, 실제 table 파일, provider page와 README를
`references/`에 보존했다. DD2 archive의 CC-BY-4.0 및 원 EOS 논문 인용 요구도 보존했다.
DD2의 원 functional은 Typel et al. (2010), TF/pasta 계산은 Grill et al. (2014), BSk22 lower
segment는 Pearson et al. (2018)에 해당한다. Zenodo depositor 명단을 원 EOS 저자로 대체하지 않는다.

저장소의 과거 table과 최신 archive는 byte identity가 모두 같지는 않았다.
GM1과 SPG의 주요 thermo 물리량은 재인코딩·출력 정밀도 수준에서 동등했다.
반면 **SLy4는 실제 물리량이 달랐다.** 현재 campaign의 공통 밀도점에서 과거/현재 참조 표의
압력 최대 상대차는 6.71857%였다. 원격 Git 이력에는 그 보존 파일이 2026-11-06 baseline
`5d56f1bca333abbf9cc1f7500e69fde9f7791ac2`에 들어온 뒤의 변경이 없다.
Import 이전 변환인지 제공자 revision인지의 귀속은 미확정이다. 이번 보고서는 새 archive를 기준으로 한다.

## 4. 공개본과 다른 이유: 측정, 대조 계산, 남은 잔차

단위는 `n_B: fm^-3`, `P,ε: MeV fm^-3`, `μ_B: MeV`다. 각 파일의 reference mass로
`P=n_B Q1`, `ε=n_B m_n^ref c²(1+Q7)`, `μ_B=m_n^ref c²(1+Q3)`를 복원했다.
Q3와 enthalpy `(ε+P)/n_B`는 별도로 비교했다. Phase 경계를 가로질러 보간하지 않았고,
공개 밀도점에서 생성 표의 PCHIP 값을 비교하며 linear control을 함께 남겼다.
보간 표현 차이는 엄밀한 오차 상한이라고 주장하지 않는다.

| 모델 | 공통 node | 최대 압력 상대차 | 최대 총 ε 상대차 | 판정 |
|---|---:|---:|---:|---|
| GM1 | 132 | 0.4863865% | 0.0276846% | 내부 matching 성립; provider 차이의 완전한 귀속 미확정 |
| SLy4 | 20 | 0.1827865% | 0.0042571% | kinetic 선택의 효과를 별도 원밀도점 계산으로 대부분 설명 |
| SPG M1 | 148 | 15.6648311% | 0.0000173% | 공개 ε 미분압력과 제공 P의 차이로 내부 잔차 대부분 설명 |
| SPG M2 | 160 | 18.6600169% | 0.0000312% | 위와 같음 |
| SPG M3 | 149 | 15.3726822% | 0.0000172% | 위와 같음 |
| SPG M4 | 158 | 16.5083026% | 0.0000233% | 위와 같음 |
| SPG M5 | 146 | 16.9184899% | 0.0000199% | 위와 같음 |

최댓값이 나타나는 밀도는 관측량마다 다르다. 전체 밀도별 값과 absolute/RMS/median residual은
`comparison/COMPARISON.json`과 CSV에 있다. 결과를 단일 최대값으로만 판정하지 않는다.

### SLy4: 과거 큰 잔차를 바로잡고 kinetic convention을 계산으로 분리

과거 git table을 기준으로 한 현 생성물의 압력 차이는 6.12431%였으나 최신 공개본에서는
0.1827865%다. 따라서 과거 큰 저밀도 잔차를 현재 provider의 문제로 재사용하지 않는다.

현재 생산 모델은 공통 `K=20.73553 MeV fm²`를 사용한다. Interaction, rest mass, lepton
처리를 유지하고 kinetic term만 `K_q=(ħc)²/(2m_q c²)`로 바꾼 독립 control을 계산했다.
`K_n=20.7212489912`, `K_p=20.7498130573 MeV fm²`다. 같은 공개 원밀도점 20개에서
β평형을 다시 풀면 최대 압력 상대차가 **0.1827748% → 0.0112669%**, 최대 proton-fraction
차이는 **8.94044e-5 → 5.75593e-6**으로 감소한다. E/A 최대차도 0.0430150 → 0.00288629 MeV다.

고정 조성에서 `δε=δK_n τ_n+δK_p τ_p`, `δP=2δε/3`, `δμ_q=δK_q(3π²n_q)^(2/3)`를
유도했고 finite difference로 최대 상대 residual 3.40854e-11을 확인했다.
이 계산은 선언된 kinetic 선택의 효과를 재현했다. 제공자의 실제 입력을 확보했다는 뜻은 아니며
생산 functional을 조용히 변경하지 않았다. 남는 작은 잔차와 정확한 provider 귀속은 분리해서 남긴다.

### SPG: 큰 압력 차이는 데이터 수준에서 설명되지만 저자 코드 원인은 미확정

같은 연속 phase에서 `P_E=n_B² d(ε/n_B)/dn_B`를 CubicSpline과 PCHIP로 계산하고 phase 양끝
각 두 node를 제외했다. 공개 표 자체의 `P_E-P`는 내부 구간에서 최대 14.89–18.11%다.
반면 현재 생성 압력과 공개 에너지 미분압력의 차이는 다음과 같다.

| 모델 | 내부 node | 최대 abs(P_gen-P_E,ref)/abs(P_ref) |
|---|---:|---:|
| M1 | 144 | 0.0162320% |
| M2 | 156 | 0.0394070% |
| M3 | 145 | 0.0214368% |
| M4 | 154 | 0.0318997% |
| M5 | 142 | 0.0238185% |

두 미분 보간법의 민감도는 공개 압력 대비 최대 0.00315–0.00524%다. 따라서 해당 내부 구간의
생성/공개 압력 차이 대부분은 공개 ε trajectory와 제공 P 사이의 차이로 수치적으로 분해된다.
공개 압력에 강제로 맞춰 생성기의 열역학 미분을 바꿀 이유는 없다.
이것이 원저자 pressure kernel의 특정 항 누락을 증명하지는 않는다.
`UNRESOLVED_SOURCE_BETA_SURFACE_DERIVATIVE_PATTERN`과 원 mass calibration 미확정은 유지한다.
Public composition을 이용한 candidate이므로 ε의 높은 일치를 독립 blind prediction으로 부르지 않는다.

### GM1: 아직 완전히 설명되지 않은 차이

현재 GM1은 원전 rounded parameters와 공통 kernel mass 939 MeV의 명시적 reconstruction이다.
Header neutron reference mass를 kernel mass로 바꾸는 control은 최대 압력 차이를
0.486386%에서 0.634622%로 키웠다. 따라서 header 값을 근거로 생산 질량을 교체하지 않았다.
다섯 원전 계수를 printed half-digit 범위 안에서 변화시킨 제한된 control은 0.178919%까지
줄였지만 전체 차이를 설명하지 못했다. Held-out 114개 밀도에서도 E/A 잔차가 남았다.

GM1은 `UNRESOLVED_SOURCE_DIFFERENCE`다. 이 상태를 “정당한 이유가 완전히 확인됐다”로
표시하지 않는다. 내부 열역학·matching 결과는 독립적으로 유효한 근거로 보존하고, 제공자 복제
동일성이나 완전히 해명된 공개본 차이를 요구하는 분석에는 이 미확정 항목을 함께 적용한다.
추가 offset fitting이나 실제 자료가 없는 원인 추정으로 닫지 않는다.

## 5. 공개 DD2 자료를 다음 계산에 연결한 실제 준비

두 Zenodo archive의 공통 283행을 비교했다. Density 정밀도 차이는 최대 8.19e-8 relative,
압력은 1.41e-7 relative, ε는 8.20e-8 relative이며 Ye/Ymu는 동일하다.
각 archive에서 다섯 밀도 좌표의 Ye, Ymu, total Yp=Ye+Ymu, P, ε, μ_B를 추출했다.
`reference_support/SEED_CATALOG.json`이 TF initial-ansatz 좌표와 branch 진단의 입력이다.
완전한 공간 field profile, cell radius 또는 nuclear group 의미를 공개 표에서 임의로 복원하지 않았다.

실제 하부 217행의 끝은 약 `2.69021e-4 fm^-3`, inner-crust 66행의 시작은
`2.0000001e-3 fm^-3`다. Unified core의 첫 tabulated point는 `0.071795955 fm^-3`다.
웹 설명의 threshold와 실제 node 지원을 구분했고, 이 경계·간격을 자동 spline으로 통과시키지 않는다.
Whole-crust 작업은 명시적인 연결법과 density origin을 갖춰야 한다.

## 6. 실행·검증과 남은 작업

- 현재 native executable build 성공; 모든 matching은 이 대화에서 제어한 별도 worktree에서 실행.
- 새 campaign 회귀 9개, reference analytic 회귀 9개 포함 **frontend/pipeline 1260 passed, 9 skipped**.
  Skip은 기존 `inputs/raw_compose`의 M0/E390 자료가 stage되지 않은 항목이다.
- 30 cell의 실행/재사용/실패와 각 output readback, consistency, method admission을 보존.
- 최신 reference 비교, 원밀도점 control, DD2 seed catalog의 실제 실행 exit 0.
- Wolfram Language에서 kinetic/first-law/상수 offset의 7개 대수 residual이 모두 0.
  첫 호출의 protected-symbol 경고와 수정한 호출을 함께 보존했다.
- ZIP의 EOS 100파일이 입력 inventory 및 ZIP readback SHA와 일치.

새 전체 generator 물리 suite를 반복 실행하지 않았다. 생산 Fortran 및 generator functional은
이번 작업에서 바뀌지 않았고 시작 head의 기존 CI 근거를 재사용했다. 이번 branch의 GitHub CI가
별도로 실행됐다는 주장은 하지 않는다. 실제 이번 실행의 로그와 source hash가 현재 근거다.

다음 작업 순서와 종료 조건은 [`PLAN_KO.md`](PLAN_KO.md)에 있다. 우선순위는 usable whole-crust
closure, 필요한 core 지원 확장, source 차이의 명시적 variant 정책, TF continuation, APR boundary
지원, 이후 TOV/Cowling이다. 나머지 390칸을 모두 성공으로 만드는 것을 목표로 삼지 않는다.

## 주요 재현 명령

기존 산출물을 덮어쓰지 않도록 출력에는 새로운 경로를 사용한다.

```sh
bash fortran_eos/build.sh build/compose_eos_fortran
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python scripts/generated_pair_campaign.py --build build/compose_eos_fortran \
  --out data/new_matched_campaign
python scripts/compare_generated_reference.py \
  --reference-root data/matched_eos_reference_20261117/references \
  --retained-reference-root eos_generator/eos_data --output build/new_reference_comparison --plots
python scripts/diagnose_reference_differences_20261117.py \
  --source-root data/matched_eos_reference_20261117/references --output build/new_reference_diagnosis
python -m pytest -q -ra python_frontend/tests pipeline/tests
```

동일한 질문의 원인 추정을 반복하기보다 위에 보존한 근거를 재사용하고, 다음 단계에서 입력·모델·
source 지원이 실제로 바뀌는 지점만 새로 계산한다.
