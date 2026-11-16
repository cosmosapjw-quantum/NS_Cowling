# SPG 표면 입력의 공개 원전 회수 — 2026-11-16

기준 commit: a49aaf9f82ddfda484e0f965aefdfe73aa7918f2. 선택 범위 R2,R3,R4,R8,R9. 생산 물리코드·기존 계산·fit 입력은 변경하지 않았다. 기존 AME2016 2,408핵종 일치는 선행 결과로 재사용했으며 새 발견으로 세지 않았다.

## 판정과 최소 unknown

이번 추가 공개자료 조사에서 M1–M5 원저자 sigma0,b 수치는 회수하지 못했다. 실제 arXiv source archive까지 열었고, 원고 뒤에 숨은 calibration 파일이나 주석 처리된 coefficient table 경로는 닫았다. 이것은 획득한 archive에 대한 관찰이며 공개 공간 전체의 부재 증명이 아니다.

가장 작은 원저자 입력은 각 M1–M5 bulk vector와 대응되는 (sigma0 [MeV fm^-2], b [dimensionless]) 열 10개 값 및 원본 file identity다. 이 입력이면 명시된 crust functional을 재평가할 수 있다. 2,408핵 mass fit 자체의 재현을 별도로 주장할 때에만 잔차·가중치·밀도·곡률·optimizer 계약 전체가 필요하다. 최종 EOS/AME 표를 다시 요청할 이유는 없다. 조건부 출력 역복원은 다른 provenance로 진행할 수 있다.

## 실제 원전·자료 읽기

| 자료 | 확인한 내용 | 근거 한계 |
|---|---|---|
| [arXiv 2402.15548v1 source](https://arxiv.org/src/2402.15548v1) | 전체 member 목록, main.tex, bibliography. TeX/bib/bbl/out와 EPS 그림만 존재. main.tex 199–212행은 Eqs.22–25; 304–312행은 Eq.31 및 2% 문장. | 이 source archive에는 계수표·코드·mass data 파일 없음. |
| [NSEoS repository](https://github.com/thomascarreau/NSEoS) | surface file commit history 26개와 현재 및 곡률 추가 전 2019 파일. | Carreau antecedent의 실제 구현. SPG 저자 코드라는 근거 없음. |
| [Zenodo 저자 query](https://zenodo.org/api/records?q=creators.name:%22Scurto,%20Luigi%22&size=25) | query 결과는 11084025 및 15024143. 전자는 기존 네 파일; 후자는 magnetized NL3wr EOS 출력 두 파일. | query-bound 결과; 저자의 모든 repository가 없다는 뜻 아님. |
| [QNP2024 author slides](https://indico.icc.ub.edu/event/180/contributions/2736/attachments/1274/2555/QNP2024.pdf) | 실제 PDF와 text, content slide8은 기존 surface form 및 AME2016 fit 반복. | 새 숫자·density·weight 계약 없음. |
| [Scurto2023 v2](https://arxiv.org/pdf/2212.09355v2) | 실제 PDF/text. CLD phase density와 R_d 정의. | CompOSE A/Z export map 회수하지 못함. |
| [Scurto2025 delta paper](https://www.aanda.org/articles/aa/pdf/2025/11/aa54729-25.pdf) | Parallel primary-paper excerpt의 AppendixA Eqs.A22–A26는 같은 CLD·두 계수식을 반복하고 fit는 2024 논문으로 환류. | 직접 PDF 획득 HTTP403; byte-verified 원본으로 주장하지 않음. |

재사용 thesis: [Coimbra2025](https://hdl.handle.net/10316/120319), SHA256 cf9ab5e4fc36d569975e3328c3d257e9e677e86a8b18f91ee0ef4072efd90929. Printed pp.22–24 Eqs.2.27–2.33, Eq.3.6, AppendixA.2 pp.109–111 및 A/Z 관련 전체-text 검색을 읽었다. 원문에서 명시하는 것은 dense/gas phase density, volume fraction f, R_d, sigma(yI)다. inspected source에서 exported A의 geometric/excess 정의는 찾지 못했다.

## 신규 회수: 곡률 추가 전 실제 두 계수 code

현재 NSEoS tree SHA e7657f8b0eb28eb015c93a7f7e1d134db8bca20d는 4계수 curvature fit다. [2019-03-22 commit 6366a09f4ea90df852b338c7e17390836911f141](https://github.com/thomascarreau/NSEoS/blob/6366a09f4ea90df852b338c7e17390836911f141/source/nseos/nuclear_surface_en.c)는 2계수 fit를 포함한다.

해당 nuclear_surface_en.c 213–300행:
- N=9, spherical_nuclei.data; 2,408개 cohort 아님.
- I=1−2Z/A, n0=nsat[1−3 Lsym I²/(Ksat+Ksym I²)].
- bulk는 이 n0,I에서 meta-model enpernuc. Surface r0=(4πn0/3)^(-1/3).
- residual Be−be[i], Be=e_bulk+e_surface+e_Coulomb. 음의 binding energy per baryon 관례, weights[i]=1.
- nb_of_params=2; initial guess (1.08,23.0), GSL lmsder, xtol=gtol=1e-8, ftol=0.
- 해당 fit surface에는 curvature 없음. 2019-05-22 history message가 curvature 추가를 명시.

따라서 antecedent를 인용했다는 이유만으로 SPG의 2% corresponding mass score를 unit B/A weights로 치환할 수 없다. 이 코드는 재현 가능한 별도 candidate 계약을 제공하지만 SPG fit-equivalence를 성립시키지 않는다. 이 lane에서 SPG mass refit나 EOS 생성은 수행하지 않았다.

## A/Z 의미의 조건부 복원 계약

원 SPG Eqs.21–25·thesis Eqs.2.27–2.33은 nI,ng,f,R_d를 정의한다. inspected provider PDF는 heavy-nucleus quadruple을 열거하나 저장된 A가 전체 dense volume인지 gas excess인지 식으로 고정하지 않는다. raw neutron/proton slot의0은 gas 부재의 증거가 아니다.

Vd=4πR_d³/3, proton gas=0라는 별도 가정 아래:

| Export map | 직접 유도한 관계 |
|---|---|
| Geometric A_G=nI Vd, Z=npI Vd | yI=Z/A_G; t_G=Ye/(Z/A_G)=fnI/nB; f=t_G nB/nI; ng=(1−t_G)nB/(1−f) |
| Excess A_X=(nI−ng)Vd, Z=npI Vd | t_X=Ye/(Z/A_X)=f(nI−ng)/nB; ng=nB(1−t_X); f=nB Ye Vd/Z; nI=ng+A_X/Vd; yI=Z/(A_X+ng Vd) |

Excess map에서 Z/A_X를 sigma에 들어가는 yI로 대입하면 틀린다. 이는 정의·baryon balance·neutrality에서 직접 유도한 조건부 식이며 실제 export 관례의 입증은 아니다. inverse-recovery lane에 두 map을 전달했다. 다중 density의 stationary equations 및 held-out thermodynamics로 조건부 map을 판별할 수 있다.

원본 two-parameter functional의 선형 재매개변수화:
D(y)=y^(-3)+(1−y)^(-3);
1/sigma=alpha+beta D;
b=alpha/beta;
sigma0=1/(alpha+16 beta).

가능한 map마다 상태·sigma를 조건부 복원하고 서로 다른 y의 여러 행에서 이 항등식을 검사하는 것이 직접적인 reconstruction route다. 독립 held-out density에서 variational·energy·pressure 잔차와 map identity를 동시에 검사해야 한다. Output-derived identification은 원본 coefficient-file acquisition과 다른 provenance다.

## 접근 실패와 한정 종료

- APS landing page HTTP403, official Harvest metadata/fulltext HTTP401: access blocked. Supplemental material 부재라고 판단하지 않았다.
- A&A PDF direct HTTP403. Parallel primary text excerpt와 original PDF byte acquisition은 다르다.
- NSEoS LICENSE404는 실제 tree의 LICENSE.txt 경로로 수정하여 회수 완료.
- Exa author/code/exact-DOI supplement 검색, SciSpace exact-paper search, Parallel primary-source search를 하나의 bounded pass로 수행했다. 새 유용한 자료 없는 반복검색은 종료. Parallel search_id=search_6c264d63b7e5484e6e3668fd4e611a06.
- 검색에서 못 찾음, 읽은 archive에 없음, endpoint access blocked를 분리했다. 현재 source gap은 물리적 no-go가 아니다.

필요하다면 사용자에게 요청할 자료는 M1–M5 원본 sigma0,b 파일과 A/Z export mapping으로 한정한다. 다만 진행 중인 conditional inversion 결과를 먼저 통합하여 unknown을 더 줄일 수 있다.

## 보존과 실행 경계

URL/hash/member·HTTP receipts와 MIT-licensed NSEoS evidence: data/research_blockers_20261116/spg_sources/. Full external paper/source/text payloads는 Git-ignored build/blocker_spg_sources/로 byte-identical 이동했다. SOURCE_CACHE_RELOCATION.json과 MANIFEST.json이 현재 위치 및 identity를 기록한다. 이 lane은 thesis payload를 새로 복사하지 않았다. Source-code mutation, author contact, paid access, full EOS refitting 없음.

skill_context_used: PhysMath Research Loop1.3.0, skill://plugins_6acdfbba921881919e0c17843706aef4/physmath-research-loop/SKILL.md; references/THEORY_LOOP.md의 R2,R3,R4,R8,R9; POLICY_INTEGRATION.md; SUBAGENT_MODE.md. 실제 source reuse, export maps, claim decision와 computational contract를 산출하고 one primary pass에서 종료. Recursive delegation 없음.

## 후속: beta 잔차에 대한 bounded formula check

2026-11-16 parent 요청으로 source equations와 production/inverse formulas를 읽고 760개 행의 잔차를 한 차례 대조했다. 새 광범위 문헌검색·production mutation·refit는 수행하지 않았다.

원문 main.tex109행은 M_n=M_p=938.9 MeV를 명시한다. native gdfm_mass=938.9와 일치한다. Native five-model n0·meson masses·12 coupling parameters는 보존된 Zenodo model_par 행을 옮긴 spg_author_candidates.json의 각 parameter_line과 수치 일치한다. Isovector vertex는 archived Gamma_rho=paper g_rho/2라는 기존 normalization 계약을 유지했다. Density-dependent rearrangement는 native neutron/proton chemical potentials에 동일하게 들어가므로 fixed-state mu_p−mu_n에서는 상쇄한다. 그 항 누락으로 charged residual만 설명할 수 없다.

proton-free gas, geometric A, original sigma(yI), T=0 아래 C=epsilon_C/f, S=epsilon_surface=2fC, y=yI라 두면 직접 변분은 다음을 준다:

- mechanical: P_I−P_g−C[1+f Phi'/Phi]=0;
- neutron: N=mu_nI−mu_ng−3y sigma'/(R nI)=0;
- beta: B=mu_pI−mu_nI+mu_e+2C/(nI y)+3sigma'/(R nI)=0.

여기서 Phi=(2−3f^(1/3)+f)/5이다. Beta는 nI,f에서 y를 바꾸는 허용 방향이므로 Phi는 변하지 않는다. 따라서 f Phi' term을 beta에 추가할 근거가 없다. 원 thesis final A46–A48과 production compose_eos_spg_native.f90의 g(3),g(4)는 위 조건으로 환원된다. Source의 bulk homogeneous beta 관계를 곧바로 gas boundary에 대입하면 surface/Coulomb constrained derivative를 놓칠 수 있다.

inverse probe의 chemical residual K는 항등적으로 K=N+yB다. 따라서 관측 K≈yB는 neutron surface condition은 근사적으로 만족하고 charged condition에 잔차가 집중됐다는 뜻이다. 원래 functional에서는 surface derivative의 계수는 전부 1이다.

source −938.9 Q4와 free-electron mu_e를 별도 재계산하여 비교했다. 전체760개 행의 최대 차이는0.00445MeV이며 beta 잔차의0.60–0.78MeV 규모를 설명하지 못한다. Electron mass를 source의 common nucleon mass와 섞거나 raw neutron/proton zero를 gas absence로 승격하지 않았다.

inverse lane이 발견한 가설 B≈f S_beta, S_beta=3sigma'/(R nI)를 같은 결과파일에서 독립 재집계했다. 정확한 input byte SHA 및 요약은 BETA_FORMULA_CHECK.json이다.

| Model | beta RMS [MeV] | neutron RMS [MeV] | RMS(B−f S_beta) [MeV] |
|---|---:|---:|---:|
| M1 | .266849 | .000516 | .004278 |
| M2 | .384866 | .000561 | .004631 |
| M3 | .273109 | .000631 | .004523 |
| M4 | .343167 | .000556 | .004654 |
| M5 | .313139 | .000481 | .004361 |

이는 표에 담긴 상태가 B_alt=mu_pI−mu_nI+mu_e+2C/(nI y)+(1−f)3sigma'/(R nI)≈0인 관계에 훨씬 가깝다는 numerical observation이다. **원 author code에서 (1−f) 오계수가 들어갔다고 입증한 것은 아니다.** 잔여 최대0.0073–0.0081MeV, 출력 반올림·추정 상태 오차·미기록 source convention이 남는다. B_alt는 원 energy functional의 beta derivative가 아니므로 production physics로 채택하지 않는다.

현재 판단은 (i) 표면 두 계수의 매우 정밀한 conditional output identification, (ii) 원 에너지 functional에 대한 source states의 charged stationary residual, (iii) 원저자 코드·export convention 미확정을 분리하는 것이다. Original functional 유지가 올바르다. 다음 최소 원저자 자료는 coefficient 수치 전체보다 **inhomogeneous beta-equilibrium를 구현한 식/코드와 A/Z export mapping**이 될 수 있다. 이것은 원저자 구현의 오류 단정이나 원 EOS 전체에 대한 no-go가 아니다.
