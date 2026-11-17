# 독립 EOS 재생성 연구·구현 결과 — 2026-11-16

시작점은 PR15의 실제 HEAD `c8ff2d6c97a93f1ccf621051861d3091bfe94d81`이다.
원래 목표인 core5종과 crust13종에 대응하는 **18개 명시적 생성 경로를 실제 실행**했다.
이 수에는 SPG 후보 surface calibration과 BBB의 출판 에너지 재구성이 포함된다.
따라서 이를 원 저자의 모든 미시적 모형·보정계수를 그대로 회복한18종이라고 부르지 않는다.
전체 목표 상태는 **PARTIALLY_COMPLETE**이며, 아래의 구현·검증 완료 범위는 실제로 사용할 수 있다.

이번 판정은 원 논문의 functional·종·제약·변분 미분, 독립 수치 계산과 밀도별 지원 범위를 기준으로 한다.
CompOSE100% 일치는 무조건적 pass gate가 아니다. 고밀도 비인과성도 그 모델의 저밀도 생성 경로를 폐기하지 않는다.
원본과의 차이, 원전 모형의 근사, 수치 구현 오류를 서로 다른 항목으로 보존했다.

## 실제 생성과 검증

|계열|실행 결과|해석 범위|
|---|---|---|
|Core5종|독립 밀도격자807행, 조성 포함 native write/readback|APR·GM1·KDE0v1·SLy4의 명시된 functional, BBB 출판 에너지·parabolic asymmetry 재구성|
|SPG M1–M5|513점씩2,565행 및257점씩1,285행 정련 비교|원전 GDFM 후보와 명시적 AME2408 surface fit의 구형 CLD|
|GPPVA8종|226개 선택 프로파일; 각23/23개 원천 사용 밀도점 확보|공간 의존 meson/Coulomb장, 국소 n/p/e, cell radius와 geometry 최적화|
|공통 outer crust|547행,18개 핵종 층, neutron-drip root|실제 HFB22/AME2016/Cu 핵질량을 사용하는 독립 Gibbs 최소화|
|BSk22 lower-inner|원전 pressure-fit/energy-derivative 표현 각각129행|출판 analytic ETFSI 표현과 조성; ETF kernel 자체를 다시 최소화한 결과는 아님|
|선택적 BSk22 경계 재구성|129/257행, native readback|원전 압력을 유지하고 외각 경계에서 first law 적분; 새로운 명시적 근사 표현|

Core의 독립 first-law 미분 오차는 최대5.184e-8MeV이다.
SPG는 반경을 virial 관계로 제거하고 gas proton을 독립 변수로 포함했다.
전이 부근 M1/M3에서 순수 neutron gas가 양성자 삽입에 불안정하다는 반례를 실제로 재현했으며,
최종 표에 각각23/20개의 양성자 gas 상태가 나타난다.
257→513점 압력 보간 차이는0.0583–0.0628%이고, 모든 행의 최대 projected gradient는1.329e-8이다.
기존 M3 optimizer 실패점은4.67e-13 residual의 정지해로 수정됐다.

TF의 조성은 공간 neutron/proton/electron profile과 총입자 분율로 제공한다.
Diffuse nucleus와 rod/slab에서 정의하지 않은 discrete A/Z나 bound/free 분해를0으로 위장하지 않는다.
Native table의 availability mask와 companion profile reader가 이 구분을 유지한다.

## GPPVA 물리와 수치 범위

8종의 bulk2,400상태, beta80상태, native DD비교12상태를 별도로 검증했다.
DD rearrangement, rho half-vertex, FSU/NL3 scalar normalization 및 양의 유효질량 가지를 확인했다.
국소 전자는 screened Poisson 문제와 일정한 electrochemical potential로 제거하며,
이 제거 전후 functional gradient와 전하 제약을 독립 검산했다.

|모형|선택한 전체 행|마지막 clustered n_B|첫 uniform n_B|원천 사용 구간의 최대 radial-mesh 압력 변화|
|---|---:|---:|---:|---:|
|DD2|28|0.070|0.075|0.1922%|
|DDME2|27|0.070|0.075|0.1729%|
|TW|28|0.075|0.080|0.1724%|
|FSU2|28|0.050|0.055|0.1780%|
|FSU2R|29|0.090|0.095|0.3385%|
|FSU2H|29|0.095|0.100|0.1704%|
|NL3wrL55|29|0.080|0.085|0.1000%|
|TM1e|28|0.095|0.100|0.2101%|

밀도 단위는fm^-3이다. 모든 모형에서0.002–0.1의23개 요청점이 계산됐으며 그 구간에 보류된 점은 없다.
전이 위치는 위 표의 **표본 bracket**이며 정밀 root가 아니다.
N192와 일부 저밀도 N384 계산을 비교했다. 독립 검토는 최종226개 프로파일 전체를 재계산했고,
에너지 차이6.90e-12MeVfm^-3, 밀도 미분 압력 차이4.89e-8MeVfm^-3,
최대 chemical KKT residual0.007282MeV를 확인했다.
한편0.002미만 TF 확장의 일부 실패점은 유지했다. 원래 공급자의 그 구간은 BSk22이고,
이 TF 외삽을 BSk22 대체물로 해석하지 않는다.

Radial-mesh 정련은 밀도 보간 수렴과 다르다.
고정 N192의 실제 midpoint와 추가 quarter-density 계산은
`data/research_completion_20261116/tf_density_review/`에서 독립적으로 다룬다.
고정 N192의88개 실제 midpoint를 검토하고, 가장 큰 비전이 보간 오차 구간에서
모형별2개 quarter-density를 추가해16개 상태를 세초기값으로 계산했다.
같은 추가 상태에 대한 log-log PCHIP 압력 오차는8모형 모두 감소했다:
DD2 0.687→0.196%, DDME2 0.370→0.0499%, FSU2 3.497→0.3948%,
FSU2H 0.587→0.2464%, FSU2R 0.823→0.2265%, NL3wrL55 0.442→0.0950%,
TM1e 0.792→0.1992%, TW 0.749→0.2434%.
FSU2의 큰 오차는 추가 초기값에서 다른 가지가 발견된 결과가 아니라 거친격자의 곡률 문제였다.
이는 선택한 최악구간에서의 실제 정련 검증이며 전체 연속구간의 동일오차 상한을 뜻하지 않는다.
Finite multi-start 최저 후보를 전역 최소의 증명으로 부르지 않으며,
거친 밀도표의 보간 정확도를 공간 격자 오차로 대신하지 않는다.

## 밀도별 물리적 지원

Core 직접 미분으로 구한 c_s^2/c^2=1의 국소 경계는 APR0.8689531516,
KDE0v1 1.4109658933, SLy4 1.2650056715fm^-3이다.
세 finite-difference step에서 root spread는 각각1.49e-12,2.67e-12,7.96e-12fm^-3이다.
이는 수치 root의 정련이며 실제 물리 모형 불확실성은 아니다.
BBB2와 GM1의 실행 구간에는 비인과적 sampled interval이 없었지만,
BBB2의0.8fm^-3 위 원 논문 값부터 이미 extrapolation이라는 제한은 남는다.

인접 표 행의 secant는 직접 미분 경계와 같지 않다.
예를 들어 첫 acausal secant가 pointwise crossing을 bracket한다고 가정하지 않았다.
전체 고밀도 진단 행은 보존하고 물리적 사용 범위를 별도로 기록했다.
자세한 생성 명령과 모델별 가정은 CORE_GENERATION.md와 MODEL_COVERAGE.json에 있다.

## 원전 입력과 공통 저밀도 crust

Outer crust의 최종 neutron drip은128Sr,
n_B=2.6943200640e-4fm^-3, P=4.9920591677e-4MeVfm^-3,
g=939.5654205200MeV이다.
257/513점 조사에서18개 핵종 층이 같았고 전이 압력 차이는최대5.414e-11relative였다.
Pearson2018 표의 반올림 값과는 압력0.367%,밀도0.456%까지 다르지만 이를 실패 기준으로 사용하지 않았다.

BSk22 저자는 에너지·압력·조성에 대한 별도 analytic fit과 Fortran 코드를 공개했다.
[저자 안내](https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html) 자체도
모든 fitted 값을 동시에 쓰기보다 thermodynamic identity로 하나를 유도하도록 권고한다.
이번 직접 구현에서 원래 P와 n_B²d(E/A)/dn_B의 상대 차이는약−8.26%~+12.07%였다.
원전 approximate representation의 잔차를 숨기거나 해당 논문의 미시적 ETFSI 전체가 틀렸다고 일반화하지 않는다.

선택적 재구성은 원전 P(n_B)를 유지하고 outer 경계에서
e(n)=e(n0)+integral[P(u)/u²]du를 적분한다.
경계의 압력·Gibbs 에너지가 일치하고 밀도 점프는+0.745634%이며,
고정밀 적분과 interval 계산으로 선언 구간의 양의 압축률·인과성을 확인했다.
그러나 fitted 자유중성자 분율은 경계에서0.00491945이고,
composition-fit mu_n과 새 mu_B의 차이는−0.021495~+0.025630MeV이다.
이는 명시된 thermodynamic representation이지 새로 푼 ETF/neutron-drip 평형이 아니다.

원전 SPG surface-fit 동등성은 아직 확정하지 않았다.
Carreau의 명시적9핵종 중48Ni는 AME2016 추정값이며,
고정 symmetric saturation 면적과 density-dependent 면적 처방은 서로 다르다.
5모형×9/8핵종×2면적의20개 fit을 실제 계산한 결과,
격자 오차보다 처방 차이가 크고 기존 후보와 저Yp surface tension 차이는최대33.07%였다.
후보 calibration을 사용한 생성은 가능하지만 원 저자의 calibration 회복과 구별해야 한다.

추가로 종별 rest-energy를 복원한 SPG 후보는 beta 평형을 다시 풀었고,
5종 모두 명시된 BSk22 energy-fit branch와 공통접선을 갖는 것을 확인했다.
이 실험은 다른 candidate construction의 실행 가능성을 보여주며,
기존 common-mass 원전 모형이나 boundary-anchored BSk22와 자동 동일시하지 않는다.

## 적대적 검토에서 수정한 결함

|발견|수정·보존|
|---|---|
|M3 optimizer 실패 및 CLD gas proton 누락|반경 제거, analytic gradient polish, gas proton 자유도와 KKT를 추가|
|Outer screening b2 분자 누락|60자리 독립 에너지 미분 반례로 검출·수정|
|얇은 중간 핵종 층 누락|제3최저 Gibbs 핵종을 발견하면 전이 bracket을 재귀 세분화|
|TF 작은 volume cell의 전자 평형 실패|screened-electron 문제를 직접 제거하고 stale warm-start를 수정|
|전자 root bracket의 ulp 경계 실패|독립 adversarial warm-start 반례를 보존하고 outward padding|
|충분히 작은 에너지 변화에서 optimizer line-search 실패|물리 residual 기준은 유지하며 실제 profile restart·경계 초기값 수정|
|BSk22 재구성 후 이전 음속 metadata가 남음|새 에너지에 맞는 derivative metadata 재계산 및 회귀 추가|
|Core 좁은 요청구간 밖으로 derivative audit가 나감|실제 stencil을 구간 내부로 제한하고 step 기록|
|기준 commit CI618 checksum 불일치|실패 로그 보존; 실제 파일·Git blob 바이트 기준 manifest 갱신|

## 검증·게시와 남은 작업

Fortran build가 성공했고, 최종 frontend/pipeline 회귀는 **773passed,9skipped**,
generator 회귀는 **142passed**이다.
실행 당시 전체 계산 코드·테스트·script SHA를 기록했으며 검사 중 변경은 없었다.
Generator에는 기존 Matplotlib Axes3D 경고1개가 남는다.
독립 검토의 PASS는 각 명시적 범위에만 해당한다.

새 게시 HEAD의 remote CI 결과는 PR15 본문과 GitHub Actions 결과에 별도로 기록한다.
CI618의 실패는 MODEL_COVERAGE의 오래된 checksum1건이었으며,
물리 코드 실패로 재분류하지 않았다.

남은 실제 작업은 원 SPG calibration/저Yp 불확실성, 정밀 density·geometry 전이 재생성,
공통 lower-crust와 각 inner functional의 물리적 에너지·조성 기준 정합이다.
원본 압력 차이를 줄이기 위한 무단 refit이나 일괄 rejection은 하지 않는다.
Python 연구 구현에서 production Fortran으로의 이관은 별도 engineering 과제다.
Core–crust 조합과 새로운 stellar/mode 계산은 이번 루프에서 수행하지 않았다.

논문·입력부터 식·구현·실행·반례·수정까지의 두10단계 기록은
RESEARCH_LOOP.json에 있고, 후속 실행 지시는 HANDOFF_KO.md에 있다.
과거 source 비교, 실패, v6 proxy receipt는 별개의 역사적 근거로 보존된다.
