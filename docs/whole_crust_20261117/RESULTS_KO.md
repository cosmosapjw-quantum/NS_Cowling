# A. Whole-crust 연결 결과

2026-11-17. 생산 코드 commit: `5c0e94649cd586d63323a9ec3fe908fd91586d6c`.
생산 원자료 checkpoint: `f7a4037bb089a167098265f4c7f3820c734b9dec`.
시작점: `8c917baf3f5975a6109c3b6e59177079c285c4e0`.

## 완료 범위

원래 A 계획의 최소 완료 조건을 충족했다. 공개 BSk22 outer와 코드 생성 SPG_M1–SLy4
결과를 명시적인 하단 bridge 및 analytic surface로 연결한 **reference-assisted hybrid**를
실제 생성했다. 주 baseline 하나와 구성 민감도 일곱 개, 총 8개 CoMPOSE 표·5,638행이다.
성공한 생산의 새 Fortran piece 실행은 17회다. 원본 source, 사용하지 않은 native prefix,
조성 형식 보정량, 열역학 재구성량, 행별 origin, 실행 입력·출력·hash를 함께 보존했다.

연속 consumer는 정확한 P=0 표면과 상단 n_B=0.5 fm^-3까지를 제공한다. 유한 CoMPOSE 표의
첫 밀도는 4.67962e-18 fm^-3이며 그 아래는 명시한 surface 식을 사용한다. 모든 표의 최대
압력은 68.33246940869708 MeV fm^-3다. CSS의 density jump 내부는 보간하지 않는다.

이는 기존 native-generated/public-original/public-reconstructed 세 계열의 sharp-Maxwell
195칸 census와 별도인 A 산출물이다. 그 census와 원본 자료는 변경하지 않았다. 이번
8개를 독립적인 microscopic EOS 8종이나 새로운 sharp-Maxwell 성공 8개로 집계하지 않는다.

| 상단 구성 | 하단 연결 native row | surface gamma | 행 수 | 역할 |
|---|---:|---:|---:|---|
| quintic_P | 40 | 5/3 | 704 | 주 baseline |
| quintic_P | 36 | 5/3 | 708 | 하단 연결점 민감도 |
| quintic_P | 37 | 5/3 | 707 | 하단 연결점 민감도 |
| quintic_P | 40 | 4/3 | 704 | 표면 처방 민감도 |
| CSS | 40 | 5/3 | 702 | core–crust 구성 민감도 |
| CSS | 36 | 5/3 | 706 | 두 구성 선택의 대조 |
| CSS | 37 | 5/3 | 705 | 두 구성 선택의 대조 |
| CSS | 40 | 4/3 | 702 | 두 구성 선택의 대조 |

주 baseline ID는 `BSk22pp__SPG_M1__RG_SLY4__quintic_P__baseline`이다.
실제 표와 계약은 `data/whole_crust_20261117/production_v2/models/`에 있다.

## 물리량·단위·구성식

n_B는 fm^-3, P와 epsilon은 MeV fm^-3, mu_B는 MeV다. epsilon은 rest energy를 포함한다.
현재 구성은 T=0의 평형 barotrope다. c_s^2/c^2=dP/d epsilon이며 c를 1로 둔 속도와
혼동하지 않는다. 내부 composition을 지정하지 않은 derived piece가 있으므로 frozen
composition derivative나 buoyancy는 이 barotrope만으로 정해지지 않는다.

### 1. 공개 outer의 압력을 유지한 재구성

공개 outer의 n/P/epsilon을 모두 고정한 cubic P(mu_E) 표현은 216개 구간 중 40개에서
양의 압축률을 만족하지 못했다. 이것은 해당 표현의 실패이며, 원래 microscopic EOS가
존재하지 않는다는 결론은 아니다. A에서는 n/P knot와 첫 epsilon/n을 고정하고 다음
관계로 종속 열역학량을 재구성했다.

\[
P(n)=P_i(n/n_i)^{\Gamma_i},\quad
\Gamma_i=\frac{\ln(P_{i+1}/P_i)}{\ln(n_{i+1}/n_i)},\quad
F(n)=\frac{\epsilon}{n},\quad \frac{dF}{dn}=\frac{P}{n^2},\quad
\mu=F+\frac{P}{n}.
\]

Gamma=1에서는 로그 적분 극한을 사용한다. 이 적분은 source의 첫 에너지에 고정되며
상대 core energy를 C*n만큼 이동하거나 matching에 맞춰 적분상수를 조정하지 않는다.
217개 n/P knot는 유지했고, 최대 |Delta F|는 **0.004915606957638374 MeV/baryon**,
최대 |Delta epsilon|/epsilon은 **5.2458563072179814e-6**다. 원본과의 차이는 공개
압력과 제1법칙을 동시에 만족시키는 이 선택에서 생겼으며 행별 correction을 기록했다.
이 pressure-primary 변형은 앞선 energy-primary reconstructed family와 구분된다.

원본의 17개 isotope switch는 공통 P/mu를 갖춘 두 endpoint로 제공되지 않는다.
이 17개 구간은 **명시적 tabulation-gap regularization**으로 기록했다. 그 안의
composition은 미정이며, microscopic sharp nuclear transition을 복구했다는 주장을
하지 않는다. 같은 핵종 내부의 보간과 isotope gap의 처방은 contract에서 구분한다.

### 2. 공개 조성의 형식 보정

최초 생산은 CoMPOSE writer의 조성 검사에서 중단됐다. 원본 첫 217개 phase-0 행은
free neutron/proton/muon이 없고 단일 핵종 A,Z 및 Ye≈Z/A를 제공한다. 그러나 nuclear
group의 마지막 값은 1이다. CoMPOSE의 number-fraction 규약으로 이를 읽으면 내부
baryon fraction이 Yh=A가 되어 정규화가 맞지 않는다. 이 원인과 실패 경로는
`production/CAMPAIGN.json`에 FAILED로 남겼다.

공식 매뉴얼의 Yi=n_i/n_B와 group의 정의, raw field inventory, presence masks,
source hash, phase, 정수 A/Z 및 전하 관계를 확인한 뒤 **해당 outer 217개 행에만**
내부 Yh=1, 표준 출력 Yi=1/A를 적용했다. 이는 주어진 순수 핵종 inventory의 바리온
보존으로 결정된다. 원본 Ye의 최대 인쇄 잔차 약 5e-6은 그대로 보존했다. 이 adapter는
열역학량을 바꾸지 않으며 energy 재구성과 별도 ledger를 갖는다. Inner 조성의 의미는
이 규칙으로 추정하지 않는다. 원본 abundance를 보존한 공개 source snapshot도 포함했다.

근거: [CoMPOSE manual, section 4.2.3, Eq. 4.22](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf).
URL 파일명은 3.00이고 확인한 PDF 표지는 version 3.01이다.

### 3. Outer–inner의 명시적 연결

주 baseline의 하단 bridge는 다음 두 실제 상태를 잇는 cubic P(mu)다.

| endpoint | n_B [fm^-3] | P [MeV fm^-3] | mu_B [MeV] |
|---|---:|---:|---:|
| 재구성한 outer 마지막 행 | 2.69021e-4 | 4.981910045986e-4 | 938.9024559275272 |
| native SPG_M1 row 40 | 1.4255858527497875e-3 | 1.4430127457960573e-3 | 940.3247772576805 |

\[
n=P_{,\mu},\qquad \epsilon=\mu P_{,\mu}-P,\qquad
\frac{c_s^2}{c^2}=\frac{P_{,\mu}}{\mu P_{,\mu\mu}}.
\]

양끝 n/P/epsilon을 동시에 맞추고 구간 전체에서 n>0, P_{,mu mu}>0,
mu P_{,mu mu}-P_{,mu}>=0을 다항식 극값으로 검사했다. 주 bridge의 최대
c_s^2/c^2는 **0.006974341835140935**다. row 36과 37로 연결점을 바꾼 대조본의
최대값은 각각 0.0338087582965788, 0.016125555528160033이다. 세 선택 모두
허용되므로 row 40은 유일한 microscopic 전이점이라는 뜻이 아니다.

사용하지 않은 native row 0–39는 주 baseline에서 삭제·수정된 source가 아니라
채택하지 않은 prefix로 별도 CSV에 보존했다. 기존 native upper knot는 유지한다.
**그 knot 사이의 연속 표현은 별도로 인증한 piecewise cubic P(mu_E)다.** 원래
quintic/CSS 생성기의 off-knot 함수와 모든 점에서 같다고 주장하지 않는다.

### 4. 표면과 cutoff

첫 outer knot 아래는 다음 처방이다.

\[
P=Kn^\gamma,\qquad
\epsilon=a n+\frac{P}{\gamma-1},\qquad
a=\frac{\epsilon_a-P_a/(\gamma-1)}{n_a},\qquad
h=\int_0^P\frac{dP}{\epsilon+P}=\ln\frac{\mu}{a}.
\]

주 gamma=5/3의 a는 929.7530114783499 MeV다. n→0에서 P,epsilon,h→0이고
mu→a다. 1<gamma<=2와 a>0에서 이 처방은 안정적이고 인과적이다. Wolfram Language로
surface 구성식의 네 항등식 잔차가 0임을 확인한 receipt를 보존했다.

| n_cut/n_anchor | 누락 h, gamma=5/3 |
|---:|---:|
| 1e-4 | 1.1141345673754571e-7 |
| 1e-8 (현재 유한 표) | 2.400330294743814e-10 |
| 1e-12 | 5.171354855149837e-13 |

gamma=4/3 대조에서는 현재 cutoff의 누락 h가 약 1.78267e-7이다. 이는 EOS surface
처방과 유한 cutoff의 수렴을 보여 주며, 별의 질량·반지름 수렴 결과는 아니다.
이 surface는 현상론적 저밀도 연장이며 원자물질의 완전한 ground-state 계산이 아니다.

## 구현 수정과 실제 검증

Fortran kernel이 새 outer/surface/bridge 상태를 계산하고, 독립 Python 표현이 같은
상태와 구간 내부의 식·극값을 대조했다. 새 코드와 기존 native CoMPOSE binary의 hash는
`production_v2/CAMPAIGN.json`에 있다. 출력의 readback은 n/P/epsilon만이 아니라
조성 값·presence·행 상태·header까지 대조했다.

실제 회귀에서는 luminal 경계의 cancellation roundoff를 발견해 causal margin의
연산 전 스케일로 오차 한계를 계산하도록 수정했다. 허용 오차를 넘는 acausality는
여전히 거절하며 c_s^2를 clip하지 않는다. 또한 CSS의 거의 같은 두 endpoint 압력이
허용 band 안에 있지만 scalar root solver의 bracket 밖이 되는 실패를 수정했다.
signed pressure band 전체를 명시적 low/high 상태로 처리하며 source 행은 보존한다.
Derived phase-only 행의 PARTIAL metadata도 native reader의 규약과 일치시켰다.

| 실제 확인 | 결과 |
|---|---|
| Python 회귀 + 실제 gfortran native control | 62 passed |
| 생성 표의 native CoMPOSE write/readback | 8/8 PASS |
| 각 구간의 analytic first law·안정성·인과성 | 허용된 모든 smooth 구간 PASS |
| 독립 검토자의 실제 표 전체 재평가 | 5,638행, 최대 n/epsilon 상대차 3.1086244689504383e-15 |
| 독립 검토자의 입력·구현·출력 hash 재검사 | 68개 일치 |
| native Euler closure 상대 잔차 최대 | 4.839952794897757e-16 |
| 연속 upper 표현의 최대 c_s^2/c^2 | 0.3497481753549825 |
| P=0 상태 | 8개 모두 n=P=epsilon=h=0 |
| 실제 CSS density jump | Delta n=1.2336473282301141e-4 fm^-3 |
| 실제 CSS latent identity 잔차 | 1.942890293094024e-16 MeV fm^-3 |

Native C7/C8/C9는 유한차분 진단이다. 최대 상대 잔차 **5.27475%, 1.88479%,
0.00209292%**를 그대로 기록했다. 이를 0이나 엄밀한 first-law 오차로 바꿔 보고하지
않는다. 새 consumer의 piece별 analytic identity 및 경계 검증은 별도 근거다. 단일
전역 spline이나 무차별 수치 미분을 적용하면 지정한 EOS와 다른 표현이 된다.

독립 최종 판정과 그 범위는 `FINAL_REVIEW_KO.md`에 있다. 광범위한 기존 전체 suite를
이번에 다시 실행했다고 주장하지 않으며, 62개 집중 회귀와 실제 생성·readback·독립
artifact 검증을 이번 변경의 근거로 사용한다.

## 사용과 산출물

사용할 때는 `compose/`와 `REPRESENTATION.json`을 함께 취급한다. 기본 consumer는
P 입력으로 n/P/epsilon/mu/c_s^2/c^2/h를 제공하고 지원 범위 밖 입력은 거절한다.
CSS에서는 `transition_side="low"` 또는 `"high"`를 명시할 수 있다.
포터블 묶음의 `evaluate.py`는 NumPy/SciPy만으로 평가하며 실행 예는 `USAGE_KO.md`에 있다.

주 묶음은 `whole_crust_primary.zip`, 전체 대조본은 `whole_crust_all_variants.zip`이다.
각 묶음은 멤버별 SHA256·크기·CRC와 독립 압축 해제 후 CLI 실행을 확인하도록 패키징한다.
실제 패키징 완료 여부는 `deliverable_v2/PACKAGING_RECEIPT.json`이 최종 근거다.

핵심 원자료: `production_v2/CAMPAIGN.json`, `OUTER_RECONSTRUCTION.json`,
`OUTER_COMPOSITION_NORMALIZATION.json`, 모델별 `RESULT.json`, `ROW_PROVENANCE.json`,
`REPRESENTATION.json`. 요약은 `SUMMARY.json`, 재개 지점은 `NEXT_STEP_B.json`이다.

## 다음 순서: B. 필요한 core 지원 확장과 source-variant 정책

다음 실행은 이미 생성된 native core를 A의 upper consumer에 연결하는 작업이다.
GM1, KDE0v1, SLy4는 앞선 실행에서 **각 385행, 실제 0.02–1.0 fm^-3**까지 생성돼 있다.
KDE0v1/SLy4 receipt의 계산 허용범위 metadata가 각각 1.8/1.4를 적더라도 실제 생성
상단은 1.0이다. 이 차이를 구분하고 기존 자료를 재사용한다.

1. A 주 baseline에 필요한 SLy4부터 기존 upper와 확장 source의 parameter, 입자,
   rest-mass/kinetic convention을 대조하고 차이 ledger를 고정한다.
2. 실제 생성된 core 구간에서 beta equilibrium, 제1법칙, 양의 압축률, causality,
   phase 경계를 확인해 사용 가능한 실제 endpoint를 정한다.
3. 기존 0.5 endpoint와 그 위 native 구간을 일관된 표현으로 연결한다. 동일 density
   anchor가 필요하면 그 지점만 native 계산에 추가하며 근거 없는 외삽은 하지 않는다.
   Whole-crust contract의 연속성·readback을 다시 확인한다.
4. SLy4의 fitted common kinetic coefficient는 기본값으로 유지한다. 입자별 kinetic
   prescription이 필요할 때만 별도 이름·입력·beta-equilibrium 해를 갖는 native variant로
   추가한다. GM1의 provider/calibration 잔차는 해소된 만큼만 귀속하고 불명 부분은 남긴다.

B의 종료 조건은 확장된 whole-crust+core와 실제 지원구간·source 차이 ledger, 또는
물리/수치/domain 실패의 정확한 기록이다. 목표 stellar central density가 그 범위에
들어오는지는 D의 실제 stellar 계산에서 확인한다. 이후 원래 순서는 **C. TF와 나머지
core/method의 필요한 확대**, **D. APR boundary 및 TOV/Cowling 소비자 연결**이다.

현재 산출물은 평형 barotrope와 구성 민감도다. TOV M–R, 최대 질량, radial stability,
Cowling mode, frozen-composition buoyancy, microscopic nuclear phase 복원을 이번 A의
결과로 포함하지 않는다.
