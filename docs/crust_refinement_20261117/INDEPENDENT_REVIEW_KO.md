# G 독립 검토: 범위를 한정한 승인

검토 기준은 F commit `d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a`에서 시작한 G 작업이다. 생산 담당과 분리된 검토자가 신규 TF registry/driver/adapter, source support 재생성, sharp stellar 소비자와 결과 보고서를 읽고 필요한 독립 계산을 수행했다. **Fresh independent review는 한 번이며, 아래에서 발견한 구체적 문제의 제한된 수리와 재확인을 같은 검토 패스에 포함했다.** 검토자가 생산 구현이나 matching 입력을 수정하거나 commit/push하지 않았다.

## 판정

**ACCEPT_SCOPED. 정의된 G 산출물에는 남은 blocking finding이 없다.**

승인 범위는 계산한 native branch와 유한 support에서의 23개 verified Maxwell 변형 출력, 명시적으로 선택한 lower closure 아래의 TM1–GM1 전체 EOS, 그리고 기록된 perfect-fluid/Cowling/phase-conversion 가정 아래의 검출된 수치 모드다. 23개 출력은 여섯 pair 유형의 mesh/anchor 변형이다. 전 geometry의 microscopic global minimum, continuum convergence, source identity가 미해결인 공개 DDME-X 표의 복구, 완전한 진동 스펙트럼이나 radial stability의 승인을 뜻하지 않는다.

G의 기술 검토 승인은 게시에 대한 사용자 승인을 대신하지 않는다. F push의 exact remote-ref receipt는 보존되어 있고 G는 별도 게시 승인 대상이다.

## 1. 직접 수행한 독립 검사

| 대상 | 독립 검사 또는 확인한 실행 근거 | 결과 |
|---|---|---|
| TF 표현 수정 | 중간/최종 API의 함수 AST 비교 | `export_compose`만 변경. `run_state`, 초기 profile 보간, build 함수 불변 |
| Native state | 776개 실제 STATE 기록에서 accepted 529개 profile의 cell weights로 baryon/charge 평균 재계산 | 최대 오차 1.3635e−15 / 2.1684e−19 fm⁻³. 저장된 압력의 Euler 산술 잔차 0 |
| 실제 출력 | TF 7개 + support 16개 출력의 파일 SHA, ADMITTED 판정, 직렬화 후 native verify 확인 | 23개 모두 기록과 byte identity 일치 |
| Support 집계 | POINT/attempt 및 matching 원 기록을 직접 재집계 | 요청 462, accepted 453, native attempt 674, Maxwell 29, admitted variant 16 |
| F 입력 보존 | 기준 commit 대비 기존 F/source data의 git diff 확인 | 수정 없음 |
| Sharp TOV/Cowling | 별도의 압력 좌표 TOV, unscaled W,V 방정식, 직접 2×2 접합 선형계 사용 | 아래 여섯 finite-jump 비교 통과 |
| 최종 build/suite | 새 G general native build 및 전체 frontend/pipeline receipt/log 확인 | build exit 0; 1427 passed, 9 skipped, 4 subtests passed; pytest 실제 exit 0 |
| 실행 identity | 전체 suite 전후 361개 code hash와 최종 bytes 확인 | 변경 없음; 검토 freeze receipt에 재기록 |

Native state의 위 검사는 직접적인 constraint/산술 검사다. 압력이 코드에서 \(P=n\mu_B-\epsilon\)로 계산되므로 Euler 산술 잔차 0만으로 에너지 미분 first law가 검증되었다고 해석하지 않았다. 생산 담당이 **다시 최적화한 인접 밀도**로 구한 DDME2 N128/TM1 N64 derivative checks를 별도로 읽었으며, 상대 차이는 각각 3.261e−7, 1.962e−9이다.

독립 stellar control은 안정하고 인과적인 합성 EOS의 유한 density jump를 사용했다. 생산 코드의 scaled variables나 transfer matrix를 독립 적분에 재사용하지 않았다. 물리적 압력 좌표에서 unscaled W,V를 적분하고, \(\Delta P\)와 slow의 W 또는 rapid의 V를 보존하는 2×2 선형계를 직접 풀었다. Slow/rapid 각각 100, 300, 900 Hz에서 생산 enthalpy/scaled 방식과 비교한 normalized determinant의 최대 절대 차이는 **1.6698e−8**이었다. 같은 유한 mode surface에서 R와 M의 상대 차이는 **3.0799e−11, 1.2936e−12**였다. 실제 process exit는 0이다.

이 control은 finite-jump 구현의 교차검산이며 실제 TM1 model spectrum을 다른 독립 물리 코드로 재현했다는 주장이 아니다. 근거는 `review/INDEPENDENT_SHARP_CONTROL.json`과 동일 디렉터리의 재현 script다.

## 2. 발견한 구체적 문제와 수리

### R1. 중심 Taylor 영역이 interface를 삼키는 입력

초기 구현에서는 합성 control의 \(n_c=5+10^{-8}\) 또는 \(5+10^{-10}\), \(r_0=10^{-3}\) km에서 \(h_0<h_t<h_c\)가 되어도 별을 반환하고 interface 목록을 비웠다. 고밀도 중심 phase의 일부를 미해상도 Taylor 영역에 묻는 구현 경계 오류였다. 생산 TM1의 중심 밀도는 약 0.3192 fm⁻³이고 해당 조건과 충분히 떨어져 있었다.

수리는 \(h_0\le h_t<h_c\)인 경우 작은 \(r_0\)를 요구하며 거절하는 guard 하나다. 문제 입력 두 개의 거절을 직접 재확인했고, \(r_0=10^{-4}\) km로 줄이면 interface가 \(r=8.7232\times10^{-4}\) km에서 실제로 검출됐다. 최종 module에서 guard 하나를 제거한 AST가 원 생산 snapshot과 같음을 독립 확인했다.

원 생산 source SHA는 `3cbac6c1a5417d4fe4f5870a6b64e41f0c51980cf45110d778883bcd57315cc4`, 최종 source는 `d5d11e595a298541abe660a514d9275dc7152bf46cdd44a5078331c8f8111c82`다. 생산을 최종 source로 다시 수행했다고 주장하지 않았으며, 모든 생산 입력에서 guard가 false인 증거와 원 snapshot을 보존했다. 새 regression은 실제 RED 뒤 stellar 7 tests GREEN이었고 최종 전체 suite도 통과했다. **CLOSED.**

근거: `CENTER_INTERFACE_PROBE.json`, `CENTER_GUARD_CLOSEOUT.json`, 생산 측 `stellar/CENTER_GUARD_REPAIR.json`.

### R2. Brent 수렴과 수치 모드 admission의 구분

초기 primary slow 저주파 root는 177.186042603 Hz에서 Brent가 수렴했지만 normalized determinant와 join mismatch가 약 **6.93e−5**였다. N32 mesh control의 slow 저주파와 rapid 첫 root도 각각 **2.39e−6, 1.60e−6**으로 정한 1e−6 join 기준을 넘었다. `root_converged=true`만으로 실제 접합조건을 통과했다고 처리하면 안 되는 사례였다.

원 출력은 그대로 남기고 이 세 근에만 같은 배경/closure/surface cutoff의 tighter ODE와 root controls를 적용했다.

| 제한된 재계산 | 선택된 frequency (Hz) | join mismatch |
|---|---:|---:|
| N64 primary, slow 저주파 | 177.185059402 | 1.45225e−7 |
| N32 primary, slow 저주파 | 177.204991069 | 8.56541e−8 |
| N32 primary, rapid 첫 root | 1904.184188330 | 3.51344e−9 |

세 결과의 surface/interface 조건도 기준을 통과했다. 별도의 joint numerical refinement 다섯 root 역시 최대 determinant/join **2.54e−7**로 통과했다. 세 original 실패를 보존한 선택 방식으로 **세 모델 control의 15개 root variant와 다섯 joint refinement root**를 확인했다. 판정은 determinant/join<1e−6, surface/traction<1e−10이며 kinematic residual은 모두 0이다.

177 Hz profile의 W 부호 변화 세 개는 중심의 roundoff 잡음에만 있지 않았다. Lower bridge 및 outer crust에 있고, lower bridge의 음속 변화에 따라 local Lamb frequency가 해당 주파수를 가로지른다. 부호 변화 자체를 구현 오류 또는 특정 mode family의 증거로 삼지 않았다. Final report의 finite-scan 및 미분류 root 서술이 적절하다. **CLOSED.**

근거: `LOW_ROOT_INSPECTION.json`, `MODE_FINDINGS_CLOSEOUT.json`, 생산 측 `stellar/MODE_ADMISSION.json`과 세 tight-control artifact/receipt. 독립 closeout의 20 variant는 main 15 + refinement 5이며 서로 다른 20개 물리 mode species를 뜻하지 않는다.

### R3. TF 보고서의 수식 문자열

TF 보고서의 최초 저장에서 문자열 escape가 일부 TeX backslash와 제어문자를 손상시켰다. Literal/raw 방식으로 문서를 바로잡은 뒤 수식과 C0 control/tab 부재를 확인했다. 최종 보고서 SHA는 `e56e93d13cbde14f39d882f2a863ce922c992c34a74cfed8baa5d00147dffb44`이다. 수치 source나 계산 결과는 바뀌지 않았다. **CLOSED.**

검토용 center probe에서도 최초 JSON 직렬화가 NumPy boolean 때문에 exit 1이었다. 그 진단 script의 boolean만 명시적으로 변환했고 이후 실제 결과를 얻었다. 이 오류를 생산 물리 실패로 분류하지 않았으며 `EXECUTION_RECEIPTS.json`에 실패와 수리 실행을 구별했다.

## 3. 물리·수치 해석의 경계

TF registry의 새 TM1/NL3/DDME-X case는 기존 case를 덮어쓰지 않는 additive 변경이다. 고정 parameter vector, rest-mass convention, scalar/isovector mapping과 uniform fixed-profile regression을 확인했다. 공개 NL3/DDME-X 원 파일은 그대로이고 새 native counterpart를 별도 source로 부르는 것이 정확하다. DDME-X supplemental J/L identity는 unresolved로 남겨야 한다.

DDME2의 F-selected basin과 G 첫 basin은 N16에서 달랐다. 정확한 F profile로 실행한 네 anchor control에서 N32/64/128이 생산 anchor와 일치한다는 기록은 초기 branch ambiguity를 제한적으로 해소한다. 전체 밀도 영역이나 모든 geometry의 최소 에너지 증명으로 확장할 수 없다. Support의 두 anchor도 고해상도 공통 accepted 상태에서 같은 branch와 양립하지만 실제 실패점은 support gap으로 남겼다.

마지막 mesh 단계의 transition-pressure 변화는 약 0.017–0.824% 범위이고 latent-energy 변화가 더 큰 경우도 있다. 이것은 측정한 mesh sensitivity다. 보고서가 임의로 continuum-converged EOS나 엄밀한 오차 상한으로 바꾸지 않은 점을 확인했다.

전체 stellar EOS는 inherited E outer/surface 및 explicit Hermite lower closure에 조건부다. TM1의 두 join을 관측량 계산 전에 선택했고 원 source 물리량을 이동시키지 않았다. DDME2의 검사한 direct joins가 거절된 것은 그 closure의 유한 실패이며, 전 범위 no-go가 아니다. Support에서 새로 얻은 네 pair도 surface까지의 전체 EOS로 자동 승격하지 않았다.

Slow와 rapid의 접합조건은 표면층 없는 perfect-fluid interface의 두 이상화다. Bulk \(\Gamma_1=\Gamma_{\rm eq}\), Cowling, electron-only native crust 적용 범위 등을 명시했다. Cold EOS 표로 실제 conversion rate, elastic crust modes, frozen-composition spectrum 또는 radial stability를 결정한 것으로 읽으면 안 된다. 원전의 nonradial 조건과의 대응은 Sotani–Tominaga–Maeda [Sec. II D](https://arxiv.org/abs/gr-qc/0108060), Tonetto–Lugones [Sec. III](https://arxiv.org/abs/2003.01259)를 직접 확인했다.

## 4. 최종 결과와 재현성의 일치

`RESULTS_KO.md`의 여섯 finest transition 값, 776 TF 호출/539 요청/529 accepted 집계, 1.4 M⊙ TM1 star와 다섯 refined frequencies는 원 기록과 일치한다. Stellar primary 값은 \(n_c\simeq0.3192215091\) fm⁻³, \(R\simeq13.875442603\) km이며 interface 반경은 약 12.886507311 km다. Primary mass crossing과 여섯 background samples를 최대질량 또는 radial stability 증거로 쓰지 않았다.

Production matching은 이전 D native matcher, 신규 TF는 G 전용 driver, 전체 suite는 실제 다시 build한 G general executable을 사용했다. 각 binary의 SHA와 source/compiler receipt를 구분했다. Historical interim API PLAN의 SHA를 최종 API로 조용히 바꾸지 않았고, 새로운 output directory를 요구하는 prospective final-API reproduction plan은 준비 검증만 수행했다. 그 준비를 전체 재생산 실행으로 표기하지 않았다.

전체 suite는 새 G native executable로 **1427 passed, 9 skipped, 4 subtests passed in 145.37 s**, 실제 child exit 0이었다. Build source signature와 test 전후 361개 source hashes가 일치했다. Final report의 게시 상태도 F remote push 완료와 G 승인 대기를 구분한다.

검토 결과의 exact artifact/file hashes, 최종 보고서 및 검사한 source identity는 `review/REVIEW_RECEIPT.json`에 둔다. 이 검토 뒤에는 추가한 계산을 다시 일반 감사 대상으로 돌리는 루프를 요구하지 않는다. 다음 작업은 보고서에 지정한 낮은 밀도 source의 실제 연장, 필요한 mesh/geometry 비교와 mode 해석으로 이어갈 수 있다.
