# Stage B1/B4: core source 및 variant ledger

기준 commit: 4069b258c906790d315e1f4760008bea9e41b598. Branch: research/core-support-20261117. 이 기록은 B의 source identity와 variant 결정이다. 실제 안정성·causality domain은 B2, whole-crust 상단 결합은 B3에서 판정한다.

## 결정

- **implementation-verified:** A와 extended SLy4의 계수·fitted kinetic coefficient·rest mass·lepton·에너지 기준이 같다. 서로 다른 grid와 support이므로 table byte identity는 주장하지 않는다.
- **결정:** 공통 K=20.73553 MeV fm²를 유지한다. 이번 지원 확장에 입자별 K variant는 필요하지 않으며 구현·채택하지 않는다.
- **numerically checked, 기존 결과 재사용:** GM1 exact-public-node 최대 pressure residual 0.4863857320%를 UNRESOLVED_SOURCE_DIFFERENCE로 보존한다. 이번 한정 탐색에서 exact provider input은 찾지 못했고 새로운 fitting은 하지 않았다.

## A와 extended SLy4의 연결

A core는 data/research_matching_20261117/integration/real_pairs_v2/core/RG_SLY4/generated, 생성 receipt는 data/research_matching_20261117/integration/real_pairs_v2/core/RG_SLY4/DISPATCH.json다. data/matched_eos_reference_20261117/campaign/RECEIPT.json의 core:RG_SLY4 seal로 연결했다.

| 항목 | A의 원래 core | 기존 extended core | 판정 |
|---|---|---|---|
| 실제 행 / 밀도 |129 / 0.08–0.5 fm⁻³|385 / 0.02–1.0 fm⁻³|실제 표 확인|
| metadata upper bound |1.4 fm⁻³|1.4 fm⁻³|1.0–1.4는 생성 support가 아님|
| 입자 / 온도 |npemu / 0 MeV|npemu / 0 MeV|동일|
| 공통 kinetic K |20.73553 MeV fm²|20.73553 MeV fm²|동일|
| n/p rest energy |939.5654 / 938.272 MeV|939.5654 / 938.272 MeV|동일|
| e/μ rest energy |0.51099895 / 105.6583755 MeV|0.51099895 / 105.6583755 MeV|동일|
| ħc |197.3269804 MeV fm|197.3269804 MeV fm|동일|
| epsilon |rest energy 포함|rest energy 포함|원래 물리 epsilon 보존|

| 계수 | 값 | 단위 |
|---|---:|---|
|t₀|−2488.91|MeV fm³|
|t₁|486.82|MeV fm⁵|
|t₂|−546.39|MeV fm⁵|
|t₃|13777.0|MeV fm^(3+3α)|
|x₀,x₁,x₂,x₃|0.834,−0.344,−1,1.354|무차원|
|α|1/6|무차원|

두 receipt의 source dictionary가 정확히 동일하다. A의 실제 마지막 n=0.5 fm⁻³ 상태는 P=68.33246940869708 MeV fm⁻³, epsilon=507.9711119469372 MeV fm⁻³다. A REPRESENTATION.json 마지막 core 상태와 부동소수점 값이 정확히 일치한다. 따라서 같은 anchor를 다시 생성할 필요가 없다.

Fortran core source hash는 A a5fd20940c286e1390ad95a5542cc1792ad1a9bb5a7b6a3ebea12b99a881d3c6, extended 396f24f0566a9eb03bbe7362cb959dd019a240c92945d533ac10cdbab818f2aa로 다르다. 실제 diff는 GM1/KDE0v1/SLy4 computational lower bound .08→.02와 설명 주석뿐이다. 이 파일의 functional·kinetic·beta solver 수식은 동일하다. 다른 historical dependency 차이는 JSON에 그대로 기록했다. 전체 historical binary identity는 NOT_INDEPENDENTLY_OBSERVED이며 binary 동일성을 주장하지 않는다.

### 에너지 기준

cold, charge-neutral, neutrino-free beta equilibrium은 μₙ−μₚ=μₑ, muon이 존재하면 μₑ=μ_μ, nₚ=nₑ+n_μ다. epsilon은 nucleon rest energy, fitted functional, rest energy를 포함한 cold free-lepton 기여의 합이다. CompOSE 복원은 각 파일 header를 사용하여 epsilon=n_B*m_n_ref*(1+Q7), P=n_B*Q1로 한다. mu_B=m_n_ref*(1+Q3)와 (epsilon+P)/n_B는 별도로 검토한다. header 재인코딩은 물리 epsilon을 보존한다. epsilon→epsilon+C*n_B shift는 적용하지 않았다.

## 공개본과 kinetic control

공개 [RG(SLY4), CompOSE134](https://compose.obspm.fr/eos/134)는 Gulminelli–Raduta unified EOS에 SLy4 interaction을 사용한다. 보존 submission PDF는 n,p,e,μ와 homogeneous phase3을 명시한다. native phase6/public phase3은 같은 homogeneous 영역을 고르는 서로 다른 label이며 crust construction 동일성을 뜻하지 않는다.

계수 출처로 source contract에 기록된 자료는 [Douchin–Haensel2000](https://arxiv.org/abs/astro-ph/0006135) Table1이다. 공통 K의 보존 conventional reference는 [HFODD output](https://www.fuw.edu.pl/~dobaczew/hfodd/pb208-a.out.txt)이며 이번 pass에서 새로 확인하지 않았다. 이 자료를 RG(SLY4) provider의 exact implementation 증거로 승격하지 않는다.

기존 diagnosis는 interaction·rest masses·lepton prescription을 고정하고 nuclear kinetic K만 Kₙ=20.72124899117293, Kₚ=20.74981305729148 MeV fm²로 바꾼 뒤 beta equilibrium을 각각 풀었다.0.08317638–0.4786301 fm⁻³의20개 exact public nodes에서 최대 상대 pressure 차이는 common K 0.1827747585%→split K 0.0112669185%, 최대 |ΔYp|는8.94044e-5→5.75593e-6이다. 이것은 **controlled difference**의 수치 증거다. provider가 split K를 사용했다는 attribution 및 모든 pressure/energy residual의 설명은 **unresolved**다.

고정 composition에서 τ_q=(3/5)(3π²)^(2/3)n_q^(5/3), Δepsilon=ΣΔK_q τ_q이다. 따라서 Δμ_q=ΔK_q(3π²n_q)^(2/3), ΔP=Σn_qΔμ_q−Δepsilon=(2/3)Δepsilon가 성립한다. K의 MeV fm²와 τ의 fm⁻⁵가 epsilon의 MeV fm⁻³를 준다. 이 항등식은 composition이 재조정된 두 beta-equilibrium trajectory의 전체 차이를 대체하지 않는다. 기존 analytic control residual3.40854e-11을 재사용했으며 새 production variant 검증으로 세지 않았다.

## GM1: 미해결 provider 차이

native GM1은 m_N c²=939 MeV, (gσ/mσ)²=11.79 fm², (gω/mω)²=7.149 fm², (gρ/mρ)²=4.411 fm², b=0.002947, c=−0.001070, isovector vertex gρ/2의 **independent rounded-parameter reconstruction**이다. native와 [CompOSE54](https://compose.obspm.fr/eos/54)는 모두 npe다. provider는 n_B<10⁻³ fm⁻³에 Douchin–Haensel lower EOS를 붙인다고 명시한다.

[APS 원 논문](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.67.2414)과 [저자기관 공개본](https://escholarship.org/uc/item/6sw424k1)의 bibliographic identity를 확인했다. 원 논문과 provider citation의 일치는 exact provider generator input 확보를 뜻하지 않는다. 보존 archive와 이번 유한 탐색에서 exact generator input/code는 찾지 못했다. 다른 곳에도 존재하지 않는다는 no-go 주장은 하지 않는다.

public GM1 header939.56533/938.272 MeV는 표 복원의 reference mass다. header만으로 provider kernel mass를 단정할 수 없다. 이전 header-mass control과 rounding-box diagnostic은 모든 차이를 설명하지 못했고 production에 채택하지 않았다. B에서 이들을 다시 실행하거나 refit하지 않았다.

132개 exact public nodes,0.0829066667–0.4986133333 fm⁻³에서 최대 |ΔP|=0.2413784708 MeV fm⁻³, 최대 |ΔP/P_ref|=0.004863857320이다. 서로 다른 극값 위치의 통계이므로 같은 density 값처럼 읽지 않는다. interpolation 포함 comparison의0.004863865296과 exact-node diagnostic 값을 JSON에서 구분했다. B 지원 확장으로 이 잔차가 해결되었다고 기록하지 않는다.

## Identity와 검증 범위

A campaign seal, extended receipt output hash 및 공개 reference intake hash를 실제 파일 SHA-256과 대조했다. 전체 path/hash/authority, exact anchor, 실제 support, kernel diff 및 prior comparison 수치는 data/core_support_20261117/SOURCE_LEDGER.json에 있다. 기존 raw reference, source epsilon, request/receipt/table은 수정하지 않았다.

source verification은 Exa search1/fetch1, SciSpace search1, Parallel search1/fetch1 총5회로 종료했다. 문헌 검색은 exact provider convention의 증거로 사용하지 않았다. publication fidelity, global ground state, stellar central-density sufficiency, TOV/Cowling은 주장하지 않는다.

## 사용한 skill context

- PhysMath Research Loop1.3.0: 상위 실행자가 전달한 selected-B/finite-evidence/preserve-source excerpt.
- Exa Search: SKILL.md, references/searching.md, patterns-papers.md.
- Parallel web-search: SKILL.md. SciSpace는 live tool contract에 따라 bibliographic discovery에만 사용했다.
