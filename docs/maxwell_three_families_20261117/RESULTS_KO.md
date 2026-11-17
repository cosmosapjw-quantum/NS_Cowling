# Sharp Maxwell: 세 계열 195개 모델 조합의 실제 생성 결과

2026-11-17에 등록된 crust 13종 × core 5종을 `native_generated`, `public_original`,
`public_reconstructed` 세 계열에서 각각 평가했다. 다른 계열을 섞은 접합은 없다.
최종 결과는 **5개 모델 쌍의 유한 지지구간 Maxwell 표 6개**다. 공개 원본과 공개본의
최소 열역학 재구성 계열에서는 현재 제공된 유효 구간 안의 단일 접합이 생성되지 않았다.
이를 해당 물리 모형 전체에서 Maxwell 해가 존재하지 않는다는 정리로 해석하지 않는다.

생산 기준 코드 커밋은 `71e7929e65d46801528b7b331d4ccd89624c6ff2`다.
출발점 `59183fcd62c91aa1506ad1608965d2c52950e2ee`는 2026-11-11의 sharp Maxwell
cohort 분기 `f46fdd5ea10ee83582f6be8b4ec851b234895ea1`보다 뒤의 관련 수정들을 포함한다.
두 기준의 실제 파일을 비교한 뒤 이 출발점에 이번 수정을 적용했다.

## 최종 결과와 범위

| 같은 계열끼리의 조합 | 모델 조합 수 | 실제 native 구간 조합 계산 | 출력 있는 모델 쌍 | 검증된 접합표 |
|---|---:|---:|---:|---:|
| 코드 생성 × 코드 생성 | 65 | 427 | 5 | 6 |
| 공개 원본 × 공개 원본 | 65 | 60 | 0 | 0 |
| 공개본 최소 재구성 × 공개본 최소 재구성 | 65 | 81 | 0 | 0 |
| 합계 | **195** | **568** | **5** | **6** |

모델 쌍과 유효 연결 구간 쌍을 구분한다. source에 불연속 phase 경계나 불량 interval이
있으면 각 유효 연결 구간을 먼저, 상대 EOS와 무관하게 정한다. 최소 3개 source knot를
갖는 모든 연결 구간 조합을 실제 kernel에 제출한다. 같은 모델 쌍에서 여러 유한 범위
접합표가 나올 수 있으므로 이를 여러 독립적인 미시 모형으로 세지 않는다.

195개 모델 조합 중 143개는 적어도 한 구간 쌍이 native kernel에 제출되었다.
나머지 52개는 공개 RG 두 core의 미해결 Q5 규약 때문에 제출할 유효 구간 쌍이 없었다.
자료 취득 실패나 실행환경 실패로 남긴 모델 조합은 없다.

### 실제 생성된 6개 표

아래 밀도 단위는 fm^-3, 압력·에너지밀도는 MeV fm^-3, 화학퍼텐셜은 MeV다.
같은 core를 붙여도 서로 다른 source component는 개별 표로 보존했다.

| Crust | Core | Candidate | n_low | n_high | P_t | mu_t | Delta epsilon | 출력 밀도 범위 |
|---|---|---|---:|---:|---:|---:|---:|---|
| GPPVA_FSU2 | GM_GM1 | 002 | 0.03679499789 | 0.04573703635 | 0.07774785550 | 944.14549887 | 8.44258536 | 0.017999999–1.0 |
| GPPVA_FSU2R | RG_KDE0v1 | 003 | 0.02363931299 | 0.02952937460 | 0.08119127103 | 948.29478432 | 5.58551471 | 0.02–1.0 |
| GPPVA_FSU2R | RG_SLY4 | 003 | 0.02495360326 | 0.03450184806 | 0.08958085975 | 948.64015298 | 9.05784841 | 0.02–1.0 |
| GPPVA_FSU2R | RG_SLY4 | 004 | 0.02720685808 | 0.03487660688 | 0.09085099541 | 948.67676783 | 7.27611250 | 0.026000001–1.0 |
| GPPVA_FSU2H | RG_KDE0v1 | 004 | 0.03061639591 | 0.03653283848 | 0.11073808235 | 949.19033112 | 5.61583008 | 0.028000001–1.0 |
| GPPVA_FSU2H | RG_SLY4 | 004 | 0.03088977904 | 0.04092133942 | 0.11261408408 | 949.25133332 | 9.52247206 | 0.028000001–1.0 |

두 FSU2R–SLy4 표는 불량 interval로 분리된 서로 다른 유효 crust 구간을 사용한다.
그 사이를 보간하거나 두 표를 하나의 전역 EOS로 합치지 않았다. 각 표의 독립적인
coexistence 조건과 유한 지지범위가 검증된 것이며, 하나의 전역 평형 해가 둘이라는
뜻도 아니다.

모든 산출물에서 새 접합은

\[
P_-(\mu_t)=P_+(\mu_t),\qquad n_+>n_-,\qquad
\epsilon_+-\epsilon_-=\mu_t(n_+-n_-)
\]

를 만족한다. 각 매끄러운 source 구간에서는

\[
n_B=\frac{dP}{d\mu_E}>0,\qquad
\frac{dn_B}{d\mu_E}>0,\qquad
0<\frac{c_s^2}{c^2}=\frac{n_B}{\mu_E\,dn_B/d\mu_E}\le1
\]

을 native 다항식의 구간 극값으로 검사했다. 성긴 sample의 PASS를 전 구간의
안정성으로 대신하지 않았다. 6개 접합의 정규화된 잠열 항등식 잔차 최댓값은
`9.489372789568138e-16`이다. 이는 미시 모형의 정확도나 mesh 오차의 크기를 뜻하지 않는다.

**출력은 유한한 inner-crust/core 지지구간의 접합표다.** neutron-star 표면까지의
outer crust, 모든 기존 내부 phase 전이, TOV/Cowling 결과 또는 최대질량을 함께
인증한 것은 아니다. GPPVA 입력은 N16 유한 mesh의 stationary native 후보이며,
고해상도 수렴과 global minimum은 별도로 미해결이다.

## 세 계열의 정의와 원본 보존

`native_generated`는 이 저장소의 생성기가 계산한 결과다. GM1·KDE0v1·SLy4는 이번에
각각 385개 밀도점, 0.02–1.0 fm^-3에서 새로 계산했다. 다른 native 입력은 이전에
실제로 생성된 결과의 receipt와 출력 해시를 검증해 재사용했다. 기존 native 18종의
238개 봉인 파일이 모두 일치했다.

이 분류가 원저자의 microscopic solver를 모두 다시 실행했다는 뜻은 아니다.
BBB는 출판된 SNM/PNM energy table에 기반한 구현이며, SPG는 공개 조성으로 추정한
surface parameter variant를 사용한다. 이 SPG의 원저자 mass calibration과 surface
equilibrium 등가성은 미해결이다. TF는 유한 mesh·geometry·seed 탐색이다. 해당
`source_contract`, variant, receipt와 제한을 각 계열의 source bundle에 포함했다.

`public_original`는 실제 공급자 원본이다. 18종 모두 확보했고 10종은 신규 다운로드,
8종은 같은 날 확보한 archive의 모든 해시를 검증한 재사용이다. 원본 ZIP, 공급자
페이지, PDF/README, CoMPOSE 파일과 provenance를 보존했다. APR_APR은 등록된
CompOSE 68번 자료이며, 다른 unified release인 328번으로 모델을 바꾸지 않았다.

`public_reconstructed`는 공개 원본의 에너지·밀도를 유지하고, 진단된 압력·Euler
chemical-potential 불일치만 수선한 계열이다. 코드 생성 EOS나 공개본에 맞춘 핵
매개변수 재적합 결과를 이 계열에 섞지 않았다. 이미 일관되거나 규약이 미해결인
부분은 바꾸지 않았다는 사실까지 receipt에 기록했다.

공개 GPPVA의 Nadd 열은 각 모델 PDF가 정의한 pasta geometry 코드다. 원본 row와
값은 `SOURCE_AUXILIARY.json`에 보존하고, working phase는 원래 Iphase와 geometry의
조합으로 나눈다. 그 경계를 넘어 spline을 만들지 않는다. KDE의 미정의 추가 열은
이름을 추측하지 않고 보존하며, matching의 thermodynamic-only projection에서만
명시적으로 제외한다. Q1–Q7을 바꾸는 작업과 이 metadata projection을 구분한다.

## 실제 수정한 결함

1. **Retained source row의 묵시적 변경.** 기존 native 함수는 유지되는 원본 행에서도
   Q3를 Euler 값으로 교체하고 Q2/Q5를 0으로, Q6를 Q7로 덮었다. 이제 유지되는 행은
   물리량 그대로 re-encode하며, 새 coexistence 두 행만 계산한다. 서로 다른
   CoMPOSE reference mass를 사용하는 경우에도 물리적 에너지와 화학퍼텐셜은 보존한다.
2. **v2의 tail 누락.** boundary phase만 반환하며 앞·뒤 source 행을 버리던 문제를
   고쳤다. inherited phase 경계와 원 source row/phase 대응표를 따로 기록한다.
   이 보존이 미검증 내부 경계를 새 Maxwell 전이로 인증한다는 뜻은 아니다.
3. **증거·실패 보고.** 실제 계산 전에 admitted로 기록하던 동작, 실패 category 유실,
   무시되던 `return_report=False`, 없는 조성에서 neutrality residual을 0으로 쓰던
   문제를 수정했다. segment가 support seal, beta flag, masks, axes를 유지한다.
4. **계산 하한과 기본 sampling 하한의 혼동.** 해석적 GM1/Skyrme functional의
   기본값 0.08을 계산 금지 경계로 재사용하던 문제를 고쳤다. 계산 지원 하한은
   세 모델만 0.02로 분리하고 default는 0.08로 유지했다. APR/BBB의 제한은 바꾸지
   않았다. 저밀도 9개 상태에서 독립적인 두-step 유한차분 first law를 검산했다.
5. **새 접합점의 1 ULP 거짓 실패.** FSU2R–SLy4 candidate_004에서 두 동일 압력을
   `Q1=P_t/n`으로 저장한 뒤 복원할 때 `-1.3877787807814457e-17 MeV fm^-3`의
   차이가 생겼다. 정확한 quadratic EOS로도 재현했다. 새 coexistence 두 행에만
   `64*machine_epsilon*pressure_scale`의 표현 오차 한계를 적용한다. 원래 압력값은
   바꾸지 않으며, phase 내부의 압력·밀도·에너지 순서 검사는 유지한다.

수정 전 실패 로그와 수정 후 검산 로그를 모두 `validation/`에 보존했다.

## 재구성의 수학과 판정

CoMPOSE의 `mn`을 정지 에너지 단위 `m_n^ref c^2`로 읽고

\[
F(n_B)=\frac{\epsilon}{n_B}=m_n^{\rm ref}c^2(1+Q_7),
\qquad x=\ln(n_B/n_0)
\]

로 둔다. 찬 한 conserved-charge barotrope의 관계는

\[
P=n_B\frac{dF}{dx},\qquad
\mu_E=F+\frac{dF}{dx}
\]

다. 각 문헌상 phase 안에서 source energy knots를 통과하는 cubic을 사용하며,
다른 좌표 `F(n_B)`의 미분과 인쇄 자리수의 민감도를 대조한다. 차이가 그 진단
기준을 넘는 경우에만 Q1/Q3를 수정한다. 에너지 이동, smoothing, 핵 모형 재적합,
조성 재평형 또는 교차점을 만들기 위한 tuning은 없다.

54개 준비 입력과 공개 원본/재구성의 18개 대응을 직접 비교했다. 밀도, Q7 및 물리적
epsilon, 질량 기준, Q2/Q4/Q5, working phase와 axes가 정확히 보존되었다. 압력 변경
행은 모두 각 segment의 불일치 진단과 일치한다. 최초·최종 생산 입력 54개도 internal
바이트와 source 바이트 및 유효 component가 동일하다.

SPG5의 source pressure와 energy derivative의 큰 불일치는 재구성 후 해당 inner-crust
구간의 native 형상 검사를 통과했다. 하지만 thermodynamic consistency의 복구만으로
상대 core와의 공통 압력 교차점이 생겨야 할 이유는 없다. Energy-primary 선택 자체도
명시적 선택이며, 표만으로 어느 열이 원저자의 의도에 가까운지 항상 유일하게 결정할
수는 없다. 세부 수식과 판정 기준은 `THERMO_RECONSTRUCTION_KO.md`에 있다.

## 공개 원본과 재구성 계열에서 출력이 0개인 이유

우선 공개 core 표에 붙어 있는 기존 crust를 새 core로 오인하지 않았다.

| Core | 실제 선택한 공개 core 시작 밀도 [fm^-3] |
|---|---:|
| APR_APR | 0.085 |
| BBB_BHF2 | 0.16 |
| GM_GM1 | 0.048 |
| RG_KDE0v1 | 0.04582369 |
| RG_SLY4 | 0.05248075 |

최신 GM1 공개본은 attached crust까지 모두 Iphase=3으로 쓰므로 phase 번호만으로
core를 고르면 안 된다. 공급자의 lower-crust cutoff 위 첫 실제 grid point를 썼다.
새로 계산한 native FSU2–GM1 전이의 core-side 밀도 `0.045737...`은 공개 GM1의
첫 core 점 `0.048`보다 낮다. 서로 다른 source의 전이가 같다는 뜻은 아니지만,
코드 생성과 공개 표의 실제 지원 범위 차이를 분명히 보여준다.

공개 SLy4와 KDE0v1은 provider가 beta-equilibrium이라고 기술하지만 저장된 Q5가
표준 한-charge beta-equilibrium 계약에 맞지 않는다. SLy4는 대체로 neutron–proton
정지 에너지 차이, KDE는 전자 정지 에너지에 가까운 수치를 보인다. 이는 export
convention의 단서이고, 그 정확한 변환이 입증된 것은 아니다. 단순히 0으로 바꾸지
않았고 재구성에서도 보존했다. native category `BETA_EQUILIBRIUM_FAIL`은 여기서
**입력 규약 미해결**을 뜻하며 실제 물질이 beta equilibrium을 위반한다는 확정 판정은
아니다. 이 때문에 각 공개 계열에서 13 crust × 2 core = 26개 모델 쌍이 차단된다.

나머지 공개 조합의 native 결과는 다음과 같다. 이 표는 모델 쌍이 아니라 구간 쌍 수다.

| 계열 | 공통 mu 구간 없음 | 공통 구간에서 P 교차 없음 | 합계 |
|---|---:|---:|---:|
| 공개 원본 | 40 | 20 | 60 |
| 공개 재구성 | 54 | 27 | 81 |

재구성 후 유효 source component가 늘어 계산 가능한 구간 조합도 늘었지만, 그 안에서
단일 Maxwell 접합은 나오지 않았다. 공개 에너지에 `C*n_B`를 더해 상대 chemical
potential을 맞추거나, 공개 표를 지원 영역 밖으로 외삽해 성공 건수를 늘리지 않았다.

## SPG_M4–BBB_BHF2의 두 교차점

독립 다항식 계산으로 17개 union interval의 모든 실근을 조사했다.

| 값 | 첫 교차 | 둘째 교차 |
|---|---:|---:|
| mu [MeV] | 952.927698608 | 954.989712061 |
| P [MeV fm^-3] | 0.399432670499 | 0.578971049365 |
| n_core−n_crust [fm^-3] | +0.004827506537 | −0.005629295196 |
| 높은 압력의 phase 전환 | crust → core | core → crust |

둘째 교차 이후에도 유효한 원래 crust 자료가 남아 있어 첫 근만 고르고 core를 계속
사용할 수 없다. 두 근은 단순근이며 압력 잔차는 2.22e-16 MeV fm^-3 이하이고,
재진입한 crust의 압력 우위는 약 1.15e-4 MeV fm^-3까지 커져 roundoff가 아니다.
둘째 교차에서 실제 평형 envelope를 따라가는 밀도 점프는 양수다. 따라서 이를
열역학적 불가능이라고 부르지 않는다. 단일 crust→core 출력과 다른 phase-chain 및
continuation 표현이 필요한 경우로 분류하고, 두 근과 그림을 그대로 제공했다.

## 검증 산출물과 사용

최종 생산 receipt는 `data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json`,
최종 source manifest는 `prepared_v2/MANIFEST.json`이다. 최초 실행은 `campaign/`에,
최종 재실행은 `campaign_v2/`에 따로 남겼다. 1 ULP 결함 발견으로 다시 계산한 것이며
실패했던 결과를 삭제하거나 같은 실행에서 성공한 것으로 바꾸지 않았다.

최종 runtime은 Python 3.12.13, NumPy 2.5.4, SciPy 1.18.1, gfortran 13.3.0이다.
실제 matching binary SHA256은
`2efbe0dcc6b037df151eeb7f582c7955ee42feb216895bbc95ad36627f242837`다.
각 결과를 native CoMPOSE writer로 쓴 뒤 다시 읽어 모든 열과 metadata를 비교했고,
직렬화된 값으로 독립 admission 및 native consistency 검사를 다시 수행했다.

최종 코드에서 실제 실행한 전체 frontend/pipeline 및 두 matching-theory suite는
**1,301 passed, 9 skipped**였고 실행시간은 147.62초였다. 로그는
`validation/final_full_tests.log`다. 9개 skip은 기존 M0/E390 경로
`inputs/raw_compose`를 요구하는 별도 legacy integration fixture의 미배치 때문이다.
이번 source intake나 195개 조합을 실행하지 않았다는 뜻은 아니다. 최종 sharp Maxwell
집중 회귀 45개, 재구성 대조 10개, 패키징 계약 4개도 해당 로그와 개별 증거에 남겼다.
이 수치는 실제 해당 실행의 결과이며 GitHub Actions 실행 성공을 주장하는 것은 아니다.

최종 배포본은 `data/maxwell_three_families_20261117/deliverable_v2/`에 있다.

- `sharp_maxwell_native_generated.zip`: 실제 채택된 6개 CoMPOSE EOS, 65개 모델 쌍 기록,
  source 및 generation receipt.
- `sharp_maxwell_public_original.zip`: 원본 18종, 65개 모델 쌍 기록 및 실패 근거.
  채택된 matched EOS는 0개다.
- `sharp_maxwell_public_reconstructed.zip`: 18종의 재구성 또는 명시적 pass-through 표,
  수정량·불변량·65개 조합의 실패 근거. 채택된 matched EOS는 0개다.
- `CATALOG.json`, `CATALOG.csv`, `MAXWELL_CENSUS.png`, `MAXWELL_CENSUS.pdf`.

원본 prefix/tail을 버리지 않고 확인할 수 있도록 별도 `source_retaining_thermo.csv`도
동봉했다. 이는 원 source 앞부분 + 검증된 유한 접합표 + 원 source 뒷부분의 물리량을
행 provenance와 함께 보여 주는 **진단용 view**다. 경계 보간이나 전 영역 검증을 하지
않았으므로 채택된 EOS 수에 포함하지 않았다. 원본의 composition/Nadd와 phase 의미는
별도 source bundle에 보존되어 있다.

3개 ZIP의 1,973개 member를 다시 열어 CRC·바이트 수·SHA256을 모두 검증했다.
실패한 후보의 compose 파일은 `matched_eos/`에 포함하지 않는다. 상세 archive hash는
`PACKAGING_RECEIPT.json`에 있다.

새 경로에서 재현하는 명령은 다음과 같다. 기존 고정 입력·receipt를 보존하므로 출력
디렉터리가 이미 있으면 덮어쓰지 않는다.

```bash
bash fortran_eos/build.sh build/compose_eos_fortran
python scripts/maxwell_three_family_campaign.py prepare \
  --output data/maxwell_reproduction_01/prepared
python scripts/maxwell_three_family_campaign.py execute \
  --manifest data/maxwell_reproduction_01/prepared/MANIFEST.json \
  --output data/maxwell_reproduction_01/campaign
python scripts/package_maxwell_three_families_20261117.py \
  --campaign data/maxwell_reproduction_01/campaign/CAMPAIGN.json \
  --output data/maxwell_reproduction_01/deliverable
```

## 원자료와 추적 경로

- CoMPOSE manual: https://arxiv.org/abs/2203.03209
- APR: https://compose.obspm.fr/eos/68
- BBB: https://compose.obspm.fr/eos/69
- GM1: https://compose.obspm.fr/eos/54
- KDE0v1: https://compose.obspm.fr/eos/92
- SLy4: https://compose.obspm.fr/eos/134
- 나머지 crust 13종의 공급자 URL, 실제 다운로드 URL, PDF, archive 및 모든 hash는
  `intake/INTAKE.json`, `intake/SOURCE_INVENTORY.json`과 계열별 source bundle에 있다.
- 이전 코드 생성–공개본 차이의 정량 분석은
  `docs/matched_eos_reference_20261117/DIFFERENCE_DIAGNOSIS_KO.md`에 보존되어 있다.

에너지-primary 선택, source 모형·calibration·mesh 범위, 공개 RG chemical-potential
규약, 미해결 inherited 경계는 결과와 함께 유지해야 하는 과학적 조건이다. 이를
제거한 채 파일 생성 성공을 완성된 별 EOS 또는 원저자 EOS의 복원으로 해석해서는 안 된다.
