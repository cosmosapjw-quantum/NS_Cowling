# Stage C: TF continuation에 필요한 source-local 기록

기준 worktree: `/home/cosmosapjw/ns_cowling_wt/tf_pair_expansion_20261117`; supplied base `894cbf1c8315f07ee926c29b83766ffdab2ac7b7`.
범위는 기존 finite-range RMF TF의 실제 입력·stationarity·source 의미를 확인하여 Stage C continuation을 지원하는 것이다. 새 full audit, parameter fit, 전체 density campaign, provider-code 재현 주장은 하지 않는다. 아래 source 상태와 새 수치 실행 상태는 별개다.

## 직접 확인한 구현과 보존 convention

**repo-inspected:** `fortran_eos/src/compose_eos_spatial_tf.f90`의 `evaluate_spatial_tf_profile`, `screened_electrons_native`, `solve_spatial_tf`를 읽었다. TF는 local nucleon densities와 finite-range scalar/vector fields, screened local electrons를 쓰는 cold nucleonic npe 모형이다. DD vertices는 local n=n_n+n_p에 의존하고 local rearrangement를 포함한다. native crust의 Y_mu는 0이며 muon threshold flag는 local kinetic electron Fermi energy 최대값과 m_mu c²를 비교한다. 공개 표의 Y_mu=0만으로 이 local 검사까지 통과했다고 주장할 수 없다.

단위는 H=ħc (MeV fm), r,R (fm), n_i (fm⁻³), masses/rest energies와 raw fields (MeV), epsilon,P (MeV fm⁻³)다. M_i*=M_i−g_s phi, n_3=n_p−n_n, rho vertex g_r tau_3/2, Coulomb V=e A_0 convention을 보존한다. 총 epsilon에는 rest energy가 들어간다. public header reference mass는 CompOSE 복원 기준이며 kernel mass의 증거가 아니다. 새로운 epsilon+C n_B shift나 parameter refit은 하지 않는다.

local functional derivative에 들어가는 항은
Σ_R=g_w′ w n+g_r′ b n_3/2−g_s′ phi n_s,
mu_n=E_Fn*+g_w w−g_r b/2+Σ_R,
mu_p=E_Fp*+g_w w+g_r b/2+Σ_R+V,
mu_e=E_Fe−V.
이 식과 기존 TF_THEORY Eq.(1)의 부호는 source ledger의 Avancini2009 Eqs.(53)–(61), Avancini2010 Eqs.(9)–(11)에 결부되어 있다. 이 pass에서는 native 식을 직접 확인하고 해당 문헌 identity를 재조회했다. Exa 0812 fetch는 본문 수식 대신 bibliographic 요약을 반환했으므로 그 결과를 수식 재유도·재검증으로 승격하지 않는다. 기존 유효 source/functional evidence를 재사용한다.

## DD2 public seed가 주는 것

**repo-inspected + derived:** SEED_CATALOG의 actual_n_fm3는 mean baryon density이며 Ye=<n_e>/n_B, Ymu=<n_mu>/n_B다. globally neutral npe(mu)N matter에서 total proton charge fraction은 Yp,total=Ye+Ymu다. 이는 free proton fraction이나 cluster의 Z/A가 아니다. nuclei 안의 proton도 total charge에 포함된다. charge_fraction_axis=0은 beta-equilibrium 1D table의 축 표시이며 물질이 proton-free라는 뜻이 아니다.

GPPVA_DD2 seed rows217,225,245,265,282는 각각 actual n_B=0.0020000001,0.0099999998,0.029999999,0.050000001,0.067000002 fm⁻³이고 Ye=0.06988,0.03148,0.02341,0.026,0.03045, Ymu=0이다. 이들은 initial composition 좌표다. n_n(r),n_p(r),n_e(r), phi(r),w(r),b(r),V(r), WS R, cluster radius, dripped neutron 분해를 공급하지 않는다. catalog가 제공하는 P,epsilon,mu_B는 비교용 thermodynamics이며 초기 field 해답이 아니다.

**literature-supported:** official DD2 submission PDF는 eos.compo phase7을 inner crust로 정의한다. 별도의 eos.thermo 마지막 추가 column이 geometry code이며 1 droplet,2 rod,3 slab,4 tube,5 bubble,6 homogeneous다. phase7을 droplet 또는 하나의 geometry로 읽으면 안 된다. PDF는 eos.micro가 없다고 명시한다. eos.mr의 radius는 stellar radius(km)이며 WS radius(fm)가 아니다.

공개 DD2 crust는 BsK22 outer crust를 n_B<0.002 fm⁻³에 사용하고 TF inner crust는 0.067 fm⁻³까지다. unified eos217은 그 위에 homogeneous npeμ core를 쓴다. 이 수치는 공개 construction의 설명이지 independently generated branch의 transition을 강제하는 matching law가 아니다. 기존 catalog는 outer row216 n=.000269021에서 inner row217 n=.0020000001 사이를 explicit join-law 필요로 기록한다. 없는 public rows를 interpolation만으로 source TF 해답으로 만들 수 없다.

**derived constructive initialization:** nonnegative dense/hollow shape s(x), weighted mean a=<s>>0를 택하고 Yp=Ye+Ymu일 때
n_p(x)=n_B Yp s(x)/a,
n_n(x)=n_g+[n_B(1−Yp)−n_g]s(x)/a,
0≤n_g≤n_B(1−Yp)
로 두면 mean baryon number와 total initial proton fraction을 정확히 맞출 수 있다. shape amplitude·width·R은 공개 seed가 정하지 않는 native ansatz 선택이다. e는 native screened-electron/neutrality solve로 다시 구한다. fraction은 beta-equilibrium minimization 중 자유롭게 변해야 하며 public Ye를 최종 constraint로 고정하면 다른 문제다. 이 문단은 초기화 구성의 유도이며 실행·성공 주장도 아니다.

## WS geometry와 정확한 cell stationarity

**literature-supported:** Grill2014 Sec.II/III와 Avancini2010은 same n_B에서 competing WS geometries 및 homogeneous state의 자유에너지를 비교한다. T=0에는 total energy 비교로 줄어든다. five idealized topologies d=3 droplet/bubble,d=2 rod/tube,d=1 slab을 지원한다. source paper가 예시에서 droplet/rod/slab만 favor했다고 해서 나머지 family를 numerical failure로 지워도 된다는 뜻은 아니다. source TF는 semiclassical이고 shell/pairing 보정이나 calibrated nuclear masses를 제공하지 않는다.

**repo-inspected + derived:** radial mean은 <f>=d/R^d ∫_0^R r^(d−1)f(r)dr다. centre regularity와 outer zero normal derivative를 가진 finite volume grid를 쓴다. native `g%lap`는 positive −Δ_x operator다. grid weights w_i=(i/N)^d−((i−1)/N)^d, face flux가 centre/edge에서 zero이며 weighted self-adjointness와 constant null mode를 갖도록 구성되어 있다. Poisson은 charge의 mean을 빼고 inverse를 계산하므로 raw charge residual을 별도로 확인해야 한다. <V>=0은 gauge 고정이다. 실제 Neumann solvability는 <n_p−n_e>=0이다.

Off-shell meson gradient 부호는 scalar +G_s, vector −G_w−G_r이고 eliminated Coulomb E_C=0.5<qV>≥0이다. x=r/R와 density profiles를 고정해 field stationary solution을 따라 미분하면 meson-gradient magnitudes는 R⁻², neutral Poisson inverse는 R²로 scale한다. 따라서 envelope derivative는
d epsilon/d logR=−2G_s+2G_w+2G_r+2E_C.
interior radius stationarity는 G_s−G_w−G_r=E_C다. 이는 diffuse finite-range functional에서 유도한 식으로 sharp-interface surface-energy 규칙을 그대로 대입하지 않는다. derivative는 MeV fm⁻³, native `radius_stationarity_MeV=abs(derivative/n_B)`는 MeV다. epsilon을 미분한 gradient에는 stationary-field envelope가 필요하므로 failed PDE solve에는 사용하지 않는다. uniform branch는 R이 redundant하다.

기존 native initial clustered family는 평균 proton fraction .10, uniform guess는 .035를 쓴다. 이들은 public Ye reconstruction이 아니다. profile을 다른 mesh로 옮길 때 weighted means와 nonnegativity를 보존하거나 복원하고 실제 재최적화해야 한다. 단순 interpolation·reevaluation은 mesh convergence가 아니다.

## 현재 native residual 의미와 Stage C 해석

아래 값은 현재 읽은 코드의 stopping thresholds이며 새 실행 PASS 기록이 아니다.

| diagnostic | 현재 코드 의미 |
|---|---|
| scalar/vector PDE | max absolute field equation residual / H³; inner target2e−10, accepted≤2e−8 (fm⁻³ normalization) |
| screened e PDE | max abs[−ΔV−4πe²(n_p−n_e)]; final≤2e−10 MeV fm⁻² |
| global neutrality | abs< n_p−n_e >≤2e−11 fm⁻³ in screened electron solve |
| chemical KKT | active n:abs(mu_n−mu_B), active p:abs(mu_p+mu_e−mu_B), absent species:one-sided inequality; e electrochemical constancy included |
| STATIONARY | optimizer success, max KKT<.02 MeV, radius residual<2e−4 MeV, neither R nor density at search bound |
| radius/density bounds | default R∈[5,180]fm; n_n,n_p≤.4fm⁻³; these are numerical bounds, not physical phase boundaries |
| UNIFORM_STATIONARY | baryon contrast<1e−4, KKT<.02MeV, no density upper-bound hit; radius redundant |

STATIONARY is first-order finite-grid evidence. A finite multi-start winner does not prove a local Hessian minimum, continuum/global minimum, or every density support. Existing selector appropriately retains unresolved families and feasible lower unconverged attempts. Report density/charge constraints and individual field residuals alongside a status label.

For continuation, compare reoptimized same-model/same-particle-content candidates at same n_B; retain dense/hollow starts and uniform competitor. On representative states and close energy crossings, reoptimize N→2N and compare epsilon/n_B,P,R,Ye, density contrast/profiles, all residuals and the energy ordering. A positive energy gap only resolves a family ordering when larger than the measured discretization/optimization uncertainty. No new tolerance widening follows from this source note.

For an optimized smooth cold neutral beta-equilibrium branch, P=n_B mu_B−epsilon and independently P=n_B² d(epsilon/n_B)/dn_B. The latter requires reoptimized nearby densities and one-sided treatment at crossings. Radius convergence alone does not verify this density derivative. Source residual, numerical nonconvergence, implementation defect and model applicability remain separate classifications.

## 확인한 source와 증거 경계

- [CompOSE DD2 crust](https://compose.obspm.fr/eos/210), [provider PDF](https://compose.obspm.fr/download/1D/Crust/DD2/eos.pdf), [DD2 unified](https://compose.obspm.fr/eos/217): official source identity, contents, particle/construction and label semantics. Parallel search `search_31e54ab9a4f3f6acbf4188686fdc671a`, extract `extract_e0e2604a08f9693b6a44557891f48feb`. PDF short quotation: “eos.micro : not available”. Extract acquired_at dates2026-11-10(PDF),2025-12-21(eos217) are provider/cache metadata, not proof of live byte freshness.
- [Grill et al.2014](https://arxiv.org/html/1404.2753), Phys.Rev.C90,045803: TF/WS framework and geometry comparison. This paper's own listed model sample is not DD2; the DD2 provider's citation supports framework provenance but does not expose exact DD2 generator input/code. Exa search requested5results; relevant URL/bibliography verified. Its original outer-crust example uses BPS, whereas the current public DD2 provider uses BsK22; do not silently merge them.
- [Avancini et al.2009](https://arxiv.org/html/0812.3170), Phys.Rev.C79,035804 Eqs.(53)–(61), [Avancini et al.2010v3](https://arxiv.org/html/1010.3644v3), Phys.Rev.C82,055807 Eqs.(9)–(11): reuse existing source-local functional ledger `data/research_completion_20261116/tf_literature.json` and `docs/research_completion_20261116/TF_THEORY.md`. Exa fetch refreshed bibliographic identity and partial text; not all equations were re-extracted. Prior raw HTML paths under build/tf_literature are absent in this worktree, so historical raw SHA256 is ledger-recorded and not freshly byte-verified.
- [Maruyama et al.2005](https://arxiv.org/html/nucl-th/0503027), Phys.Rev.C72,015802: WS zero-gradient boundaries and diffuse profile/radius discussion, reused prior ledger without new source fetch.
- SciSpace search returned [Warm pasta phase](https://scispace.com/papers/warm-pasta-phase-in-the-thomas-fermi-approximation-l6jn1qcgw2xm?utm_source=chatgpt&utm_medium=tool) with2022 repository metadata rather than original2010 bibliographic date. It is discovery-only; exact date/equations use primary source and prior ledger. Other returned papers were outside the narrow model target and not adopted.

External source checks stop at5calls: Exa search1/fetch1, SciSpace search1, Parallel search1/fetch1. There were no source fetch failures; Exa's equation excerpt coverage was limited and SciSpace's metadata date was unsuitable for primary attribution. No provider equality, numerical campaign convergence, production matching or stellar result is asserted by this note.

## 사용한 skill context

PhysMath Research Loop1.3.0 SKILL.md and references/POLICY_INTEGRATION.md were actually loaded. Applied excerpts: “Retrieve only the sources and prerequisite evidence needed for selected work”; “Reuse unchanged valid work”; “Numerical/scientific output” requires justified residuals, invariants and convergence rather than byte equality; “Complete exactly the selected scope.” Parent assigned Stage C source-local continuation support and no recursive delegation; no child agents were spawned. Exa Search SKILL.md + references/searching.md + patterns-papers.md and Parallel web-search SKILL.md were loaded; SciSpace used its discovered live tool contract. The companion JSON preserves file identities, source URLs, tool IDs, evidence states and limitations.
