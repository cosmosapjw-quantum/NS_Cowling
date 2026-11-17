# G: TF crust source conventions and bounded precision lookup

작성일: 2026-11-17 UTC. 기준은 F commit `d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a`와 그때 동결한 source/parameter pack이다. 이 문서는 G의 native crust 계산을 위한 R2/R5 지원 결과다. G 전체 수치 캠페인의 통과 여부나 Maxwell 후보 승격을 판정하지 않는다.

## 1. 결론과 적용 범위

**기존 native spatial TF에는 finite-range meson fields, 셀 평균 중성 조건, 공간적으로 screening되는 전자, DD rearrangement가 이미 구현되어 있다.** TM1/NL3/DDME-X를 논의할 때 필요한 연속체 수준의 항은 Xia et al.의 TF 명세와 대응시킬 수 있다. 다만 공개 XMLSLZ 표와 native 계산은 mesh, 최적화, 실제 입력의 숫자 정밀도와 mass convention까지 같다고 확인되지 않았다. 특히 새 canonical TM1 crust는 독립적으로 재생성한 모델이며, 공개 TM1 표를 그대로 정밀화한 자료라고 부를 수 없다. [S1–S3, C1–C2]

**점검한 공식 경로에서는 더 높은 정밀도의 DDME-X 또는 NL3 crust 행을 확보하지 못했다.** 현재 CompOSE 직접 다운로드는 F의 행과 byte-identical이었다. arXiv source archive에는 논문에서 언급한 `RDF_beta.dat`가 없었다. APS/Zenodo 추가 경로는 401/403이어서 접근 가능한 자료가 아니었다. 이는 점검한 경로에 관한 유한한 결과이고, 저자에게 원자료가 없거나 더 정밀한 자료가 존재하지 않는다는 증거는 아니다. 데이터의 Q3를 덮어쓰거나 허용오차를 바꾸는 근거도 아니다. [D1–D3]

**원 DDME-X supplement의 \(J,L\) 정체성은 unresolved다.** 고정된 네 입력을 실제 native driver로 계산해 nucleon rest-energy와 symmetry-energy 정의의 두 가지 가능성을 검사했다. 어느 것도 supplement의 \(J=32.865829\), \(L=47.806151\) MeV를 재현하지 않았다. 이 검사는 매칭 입력을 수정하지 않은 별도 진단이며 fitting을 하지 않았다. [N1–N3]

## 2. source–native convention 대응

F의 원문과 parameter pack을 재사용했다. Xia 2022 Eq. (13)–(20), Niu 2025 Sec. II와 Table IV를 읽었으며, 수치 코드와 아래처럼 연결된다. [S1–S2]

| 항목 | 원문 또는 frozen source의 명세 | native 구현 및 해석 |
|---|---|---|
| 상태 | \(T=0\), neutrino-free β equilibrium | neutron/proton/electron TF. native crust selector는 muon threshold가 넘는 상태를 제외하므로 그 조건 안에서만 동등한 입자 내용이다. |
| 핵자 effective mass | source \(\sigma\) convention: \(M_i^*=M_i+g_\sigma\sigma\) | \(\phi=-\sigma\), \(M_i^*=M_i-g_\sigma\phi\). \(M_i=m_ic^2\)를 MeV로 저장한다. |
| isovector vertex | \(\Gamma_\rho\tau_3\) | \(g_\rho\tau_3/2\), \(g_\rho=2\Gamma_\rho\), \(\tau_{3p}=+1\). |
| DD coupling | \(x=n_b/n_0\)의 rational \(g_\sigma,g_\omega\), exponential \(g_\rho\) | `tf_couplings`는 셀의 평균이 아닌 **각 mesh 위치의 \(n_n+n_p\)** 를 대입하고 분석적 \(dg/dn_b\)를 계산한다. |
| 전하 중성 | 셀 적분으로 \(N_p=N_e+N_\mu\); 공간별 중성 조건은 아님 | 전자만 있는 admissible 영역에서 \(\langle n_p-n_e\rangle=0\). `screened_electrons_native`가 공간적인 \(n_e(r)\)를 계산한다. |
| Coulomb | 공간 전하 screening을 포함한 Poisson 문제 | \(V=eA_0\) convention. \(\mu_p\)에는 \(+V\), \(\mu_e\)에는 \(-V\); mean-zero gauge. Coulomb energy는 \(\tfrac12\langle(n_p-n_e)V\rangle\). |
| meson/Coulomb BC | 셀 중심과 경계의 reflective/Neumann 조건 | cell-centered finite-volume Laplacian, 경계 밖 face 없음. 원문의 cosine-transform solver와 이산화 방법은 다르다. |
| geometry | droplet/rod/slab/tube/bubble/uniform, \(D=3,2,1\) | \(D\) 및 droplet/hole/uniform seed를 통해 표현. 각 geometry의 radius 및 density profile에 대한 최소화가 필요하다. 초기 seed 자체는 최종 phase 판정이 아니다. |
| 표의 phase label | \(3,2,1,-2,-3,0\) | 각기 droplet, rod, slab, tube, bubble, uniform. 물리적 비교에는 label뿐 아니라 최종 density contrast/profile과 에너지가 필요하다. |
| numerical state | continuum Euler–Lagrange equations | scalar/vector/electron PDE residual, chemical KKT residual, radius stationarity, search-bound contact가 별도로 기록된다. `optimizer_success`만으로 stationary를 선언할 수 없다. |

Native electron chemical potential은 electrochemical 값이다. muon guard는 \(\max_r[\mu_e(r)+V(r)]>m_\mu c^2\)를 검사한다. G가 다루는 저밀도 범위라는 이유만으로 이 guard를 생략할 수 없다. 원문 Niu 표는 muon을 실제 포함한 더 넓은 밀도 범위를 다루므로, native 전자-only branch를 전체 공개 표의 입자 내용과 동일시하지 않는다. [S2, C1]

Xia의 극저밀도 rod/slab treatment는 유한 transverse size를 따로 다룬다. 현재 문제가 된 DDME-X \(n_-\simeq0.00857468\ {\rm fm}^{-3}\), NL3 \(n_-\simeq0.03146277\ {\rm fm}^{-3}\)는 그 극저밀도 특례의 영역이 아니다. 이 사실만으로 TF shell correction, lattice correction 등 원문에 없는 항을 새로 넣거나 빼는 변경을 정당화하지는 않는다. [S1, G support diagnosis]

### Frozen model identity

| 모델 | G에 관련된 frozen 입력 | 취급 |
|---|---|---|
| canonical TM1 | \(M_n=M_p=938\) MeV, \(m_\sigma c^2=511.198\), \(m_\omega c^2=783\), \(m_\rho c^2=770\) MeV; \(g_\sigma=10.0289,\ g_\omega=12.6139,\ \Gamma_\rho=4.6322,\ g_2=-7.2325\ {\rm fm}^{-1},g_3=0.6183,c_3=71.3075\) | F의 `TM1_CANONICAL` 그대로. 공개 TM1 header의 reference masses와 physical RMF masses를 혼동하지 않는다. F가 확인한 TM1/PK1 label ambiguity도 유지한다. |
| NL3 | F의 `NL3_PROVIDER_PRECISION` | original NL3의 rounded vector와 별개 숫자 입력이다. 공개 provider와 비교하려면 frozen 고정밀 벡터를 사용한다. |
| DDME-X | F에서 provider 비교에 사용한 \(n_0^{\rm ref}=0.151832\ {\rm fm}^{-3}\), \(M_n=M_p=938.5\) MeV 등의 전체 벡터 | 원 supplement의 NMP 표와 정체성 문제가 남아 있다. 아래 mass/definition 진단으로 이 입력을 교체하지 않았다. |

여기서 \(n_0^{\rm ref}\)는 density-dependent vertex에 넣는 기준 밀도이고, 실제 \(P_{\rm SNM}=0\)으로 계산되는 \(n_{\rm sat}\)와 구분한다. 수치상 근접하다는 이유로 서로 덮어쓰면 functional 자체가 바뀐다. 전체 값과 출처는 F의 `PARAMETER_INPUTS.json`이 SSOT다. [S3, N1]

## 3. functional, rearrangement와 압력의 조건

에너지 단위를 MeV, 길이를 fm로 두고 \(\hbar c=197.3269804\ {\rm MeV\,fm}\)를 명시적으로 사용한다. 아래 \(M,m_\sigma c^2,m_\omega c^2,m_\rho c^2,\phi,\omega,\rho,V\)는 코드에서 에너지 단위로 저장되는 양이다. 코드의 scalar/vector potential convention은

\[
U(\phi)=\frac{(m_\sigma c^2)^2\phi^2}{2}
+\frac{a_3\phi^3}{3}+\frac{a_4\phi^4}{4},
\qquad a_3=-g_2\hbar c,\quad a_4=g_3,
\]

\[
W(\omega,\rho)=\frac{(m_\omega c^2)^2\omega^2}{2}
+\frac{(m_\rho c^2)^2\rho^2}{2}
+\frac{c_3\omega^4}{4}
\]

이다. 현재 frozen TM1/NL3/DDME-X에는 \(\omega^2\rho^2\) cross term을 더하지 않는다. \(U,W\)를 \((\hbar c)^3\)으로 나누면 에너지 밀도 MeV fm\(^{-3}\)가 된다. [S3, C2]

Native off-shell energy expression은 local Fermi energy, \(U\), vector source terms와 \(-W\)를 포함하며 scalar gradient는 \(+\), vector gradient는 \(-\), Coulomb energy는 \(+\frac12\langle(n_p-n_e)V\rangle\)로 쓴다. 이것은 solved vector field를 먼저 대입한 on-shell positive-energy expression과 쓰는 방식이 다르다. 따라서 vector gradient의 부호만 떼어 보고 오류라 판정하거나 on-shell 식을 중복해서 더하면 안 된다. 실제 구현의 field equation과 같은 variational expression을 비교해야 한다. [C1]

DD interaction에서 필요한 local rearrangement는 위 convention으로

\[
\Sigma_R=-g_\sigma'(n_b)\phi(n_{s,n}+n_{s,p})
+g_\omega'(n_b)\omega n_b
+\frac{g_\rho'(n_b)\rho}{2}(n_p-n_n).
\]

세 항의 차원은 모두 MeV다. 이를 핵자 chemical potential 양쪽에 넣으면

\[
\mu_n=E^*_{F,n}+g_\omega\omega-\frac{g_\rho\rho}{2}+\Sigma_R,\qquad
\mu_p=E^*_{F,p}+g_\omega\omega+\frac{g_\rho\rho}{2}+\Sigma_R+V
\]

이고, \(\mu_e=E_{F,e}-V\)다. Native `evaluate_spatial_tf_profile`에 이 식이 직접 들어 있다. TM1/NL3의 constant coupling에서는 \(g'_i=0\)으로 돌아간다. [S1, C1–C2]

### Pressure identity의 전제

\(\epsilon\)을 셀 평균 에너지 밀도라 하고, profile과 radius를 포함한 자유도를 \(q,R\)라 하자. 평균 baryon density와 charge constraint를 고정한 variational functional은

\[
\mathcal L=\epsilon-\mu_B(\langle n_n+n_p\rangle-\bar n)
-\lambda_Q\langle n_p-n_e\rangle .
\]

내부 stationary point에서는 \(\mu_n=\mu_B\), \(\mu_p=\mu_B+\lambda_Q\), \(\mu_e=-\lambda_Q\)이므로 \(\mu_n=\mu_p+\mu_e\)가 된다. 경계에 붙은 0-density degrees of freedom은 equality 대신 KKT inequality를 만족해야 한다.

제약 조건을 제거한 좌표에서 \(\epsilon=\epsilon(\bar n,q,R)\)를 쓰면, stationary \(q_*(\bar n),R_*(\bar n)\)에 대한 envelope theorem으로

\[
\frac{d\epsilon_*}{d\bar n}=\mu_B,\qquad
P=\bar n^2\frac{d(\epsilon_*/\bar n)}{d\bar n}
=\bar n\mu_B-\epsilon_*.
\]

미수렴 경로 \(q(\bar n),R(\bar n)\)에서는 정확하게

\[
\bar n^2\frac{d(\epsilon/\bar n)}{d\bar n}
-\left(\bar n\frac{\partial\epsilon}{\partial\bar n}-\epsilon\right)
=\bar n\frac{\partial\epsilon}{\partial q}\frac{dq}{d\bar n}
+\bar n\frac{\partial\epsilon}{\partial R}\frac{dR}{d\bar n}
\]

가 남는다. 오른쪽이 stationarity에서 0이 되는 것을 Wolfram evaluator로 symbolic check했다. 따라서 `pressure=nbar*mu_B-epsilon`로 정의한 결과의 Euler identity가 성립해도, profile/radius stationarity와 독립적인 density derivative first-law 검증이 자동 완료되는 것은 아니다. 비미분적인 geometry switch에서는 각각의 branch에 대해 먼저 적용하고 transition 처리를 별도로 해야 한다. [W1]

이 symbolic check는 수치 수렴이나 Maxwell matching의 증명이 아니다. Wolfram context 호출은 internal error였고, 뒤이어 실행한 evaluator는 실제 algebraic result를 반환했다. 두 상태를 구별하여 receipt에 보존했다.

## 4. 고정밀 crust source에 대한 유한 조회 결과

Niu arXiv v2는 `RDF_beta.dat`와 `RDF_2D.dat`를 언급한다. 그러나 조회 결과를 실제 파일 확보와 구분해야 한다. [S2]

| 확인한 공식 경로 | 실제 응답 | G에서 가능한 결론 |
|---|---|---|
| CompOSE 현재 DDME-X `eos.thermo` | HTTP 200, 148965 bytes, SHA256 `373c8ba02b983231c5ec63f284901d67c5e8d96d14fcbd74ee643b8c7a301459` | F와 byte-identical. 원 C2 실패를 고칠 새 숫자를 제공하지 않는다. |
| CompOSE 현재 NL3 `eos.thermo` | HTTP 200, 148966 bytes, SHA256 `d2b373ea55402e8496e07860fea5716a18e2bb610239e3556aeeacd15d2a508b` | F와 byte-identical. 이전 interval의 interpolant-shape 문제를 제거할 새 행이 아니다. |
| arXiv `2506.11492v2` source | HTTP 200, 1081594 bytes, SHA256 `4c84cde3c09045d512deabe5faaee5f00402bd9a088230b9a106e669edbcaa15` | 16 members는 README, TeX, figure PDF. 광고된 두 numerical .dat 파일은 없었다. 이름에 “beta”가 있는 figure PDF를 데이터로 취급하지 않았다. |
| APS supplement URL, Exa fetch | `CRAWL_HTTP_401` | 접근 실패. supplementary table 내용 미확인. |
| DOI direct fetch | HTTP 403 | 접근 실패. |
| Crossref DOI metadata | HTTP 200 | DOI `10.1103/x9fk-ykm1`, PRC publication date 2026-02-11 및 official harvest link를 확인했다. 논문 본문/첨부 자료 확보와 별개다. |
| APS official harvest fulltext | HTTP 401 | 접근 실패; 이후 우회하지 않았다. |
| Zenodo record 14809151 metadata | HTTP 200 | 과거 DDME-X ZIP의 파일명, 173141 bytes, advertised MD5를 확인했다. |
| Zenodo archive API와 일반 web download | 각각 HTTP 403 | ZIP bytes를 얻지 못했으므로 수치 정밀도나 F와의 semantic identity는 비교하지 않았다. |

Exa/SciSpace는 source discovery에, Parallel은 DDME-X 원문 정의/erratum의 좁은 조회에 사용했다. Crossref의 실제 official metadata를 기준으로 게시일을 기록하며 검색 서비스의 날짜를 논문 출판일로 대체하지 않는다. Parallel 결과에서 본 문제를 해결하는 공식 correction은 확인되지 않았다. 추가 문헌 survey나 parameter search를 계속하지 않고 이 범위에서 동결한다. 원 tool results는 외부 evidence directory에, repository에는 URL·상태·hash와 필요한 결론만 남겼다. [D1–D4]

따라서 DDME-X의 C2 discrepancy 및 NL3의 부족한 retained crust support는 현재 지원 담당자의 진단 그대로 유지한다. 이 작업으로 더 정밀한 public crust를 새로 확보했다고 주장할 수 없다. 추후 실제 supplementary table이 확보되면 그 bytes와 provenance를 동결한 뒤 동일한 admission gate로 검사하거나, 별도 native TF model로 생성하고 그 차이를 명시해야 한다.

## 5. DDME-X NMP: mass 및 symmetry-energy 정의 진단

원 supplemental Table II는 \(n_0=0.151832\ {\rm fm}^{-3}\), \(E/A-M=-16.097256\), \(K=267.046210\), \(J=32.865829\), \(L=47.806151\) MeV를 싣는다. 기존 source identity assessment는 unresolved 상태다. 이번에는 frozen driver와 두 \(n_0^{\rm ref}\), 두 equal nucleon rest energy만 쓰는 사전에 정한 네 control을 계산했다. parameter를 목표 NMP에 맞추지 않았다. [N1–N3]

정의는

\[
\delta=1-2Y_p,\quad
J_2=\frac12\left.\frac{\partial^2(E/A)}{\partial\delta^2}\right|_{\delta=0,n_{\rm sat}},
\qquad L_2=3n_{\rm sat}\left.\frac{dJ_2(n)}{dn}\right|_{n_{\rm sat}}.
\]

대안 \(J_{\rm PNM-SNM}=E(n,Y_p=0)/A-E(n,Y_p=1/2)/A\)는 일반적으로 \(J_2\)와 다르다. \(E/A=E_0+S_2\delta^2+S_4\delta^4+S_6\delta^6\)이면 각각 \(S_2\)와 \(S_2+S_4+S_6\)임을 같은 symbolic check에서 확인했다. 고차 항이 있는 경우 정의 변경은 의미를 바꾸며 상호 교환할 수 없다. [W1]

Native `mode=baryons`에서 \(P_{\rm SNM}=0\)의 root를 찾고, \(\delta=0.01,0.005\), density step \(h=10^{-5},5\times10^{-6}\ {\rm fm}^{-3}\)를 사용했다. 네 control마다 37회의 native state 호출, 총 148회를 수행했고 process exit 0을 확인했다. 개별 native request/output/receipt는 `nmp_definition_diagnostic/`에 있다. 표의 \(J,L\)는 \(\delta=0.005,\ h=10^{-5}\) 결과를 반올림했다. [N2–N3]

| \(n_0^{\rm ref}\) (fm\(^{-3}\)) | \(M_n=M_p\) (MeV) | \(n_{\rm sat}\) (fm\(^{-3}\)) | binding (MeV) | \(K\) (MeV) | \(J_2\) (MeV) | \(L_2\) (MeV) | \(J_{\rm PNM-SNM}\) (MeV) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.151832 | 938.5 | 0.1516597953 | −16.062449 | 266.3934 | 32.2430 | 49.6524 | 33.3189 |
| 0.151832 | 939.0 | 0.1519176082 | −16.095094 | 267.3482 | 32.2577 | 49.6564 | 33.3346 |
| 0.152000 | 938.5 | 0.1515798747 | −16.064762 | 266.1020 | 32.2550 | 49.6888 | 33.3305 |
| 0.152000 | 939.0 | 0.1518376638 | −16.097391 | 267.0550 | 32.2698 | 49.6939 | 33.3462 |

마지막 control의 binding, \(K\), saturation density가 원 Table II에 더 가까워지는 사실은 **numerically checked**다. 그러나 이것만으로 저자의 실제 mass/reference convention을 확정할 수 없고, \(J,L\)는 여전히 목표와 다르다. 수치 미분 step의 변동보다 discrepancy가 훨씬 크며, PNM−SNM 정의를 써도 \(J\)는 33.32–33.35 MeV여서 원 값에 맞지 않는다. 그 정의의 slope 역시 52.66–52.71 MeV이다.

판정은 “mass-only explanation과 이 두 symmetry-energy 정의 사이의 선택만으로 original supplement의 모든 NMP가 해결되지 않았다”이다. 실제 원 저자의 입력/정의/오타 중 무엇이 원인인지는 **unresolved**로 유지한다. F/G matching parameter pack이나 source Q3는 이 진단 때문에 바뀌지 않았다.

## 6. 증거와 재현 경로

Repository-relative paths:

- **S1**: `data/new_model_maxwell_20261117/sources/papers/PROVENANCE.json`, Xia et al. [arXiv:2208.12893v1](https://arxiv.org/pdf/2208.12893v1), Eq. (13)–(20), Tables I–II. 원 PDF/추출문은 F worktree의 ignored evidence.
- **S2**: 같은 provenance, Niu et al. [arXiv:2506.11492v2](https://arxiv.org/pdf/2506.11492v2), Sec. II, Table IV; 출판 식별자 [DOI:10.1103/x9fk-ykm1](https://doi.org/10.1103/x9fk-ykm1).
- **S3 / N1**: `data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json`, SHA256 `f76284c2d291719767209a9530480cd82fd4f0f527eee7538d5706863ed04451`; `uniform_validation/SOURCE_IDENTITY_ASSESSMENT.json`. Original DDME-X [arXiv:1910.13007](https://arxiv.org/abs/1910.13007), supplement source provenance `sources/papers/ddmex_source_provenance.json`. TM1 coefficients are also independently tabulated in Long et al. [nucl-th/0311031](https://arxiv.org/abs/nucl-th/0311031), Table I.
- **C1**: `fortran_eos/src/compose_eos_spatial_tf.f90`: `make_tf_grid`, `tf_poisson_native`, `screened_electrons_native`, `evaluate_spatial_tf_profile`, `solve_spatial_tf`.
- **C2**: `fortran_eos/src/compose_eos_tf_models.f90`: `tf_couplings`, `tf_scalar_potential`, `tf_vector_potential`.
- **D1–D3**: `data/crust_refinement_20261117/source_physics/PRIMARY_SOURCE_FETCH.json`, `ADDITIONAL_PRIMARY_FETCH.json`, `FINAL_PRIMARY_FETCH.json`: exact URL, retrieval UTC, success bytes/hash or actual HTTP error.
- **D4**: `source_physics/SOURCE_PHYSICS_MANIFEST.json`: external raw connector result paths and hashes, source reuse, finite lookup scope.
- **W1**: `source_physics/SYMBOLIC_CHECK.json`: exact Wolfram code/result and context failure receipt.
- **N2**: `source_physics/nmp_definition_diagnostic.py`, `source_physics/nmp_definition_diagnostic/NMP_DIAGNOSTIC.json`.
- **N3**: `source_physics/NMP_EXECUTION_RECEIPT.json`; all 148 native calls retain `RUN.json`. Frozen F driver SHA256 `5d65389f8cb29902e051864f8040f4df0b4ef48ba504da9e979ed9ceaea5f2d4`.

새로 받은 full paper archive, figures, full connector excerpts는 `/home/cosmosapjw/ns_cowling_wt/scratch/G_source_physics_20261117/`에 둔다. Repository에는 논문 전문을 복사하지 않았다. 이 subtask는 commit/push를 수행하지 않았고, 재생성한 crust의 실제 convergence/Maxwell 판정은 G의 해당 결과를 따른다.
