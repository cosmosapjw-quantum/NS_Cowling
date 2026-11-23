# NS_Cowling — 원전 기반 EOS 생성·접합 완료 보고서

기준일: 2026-11-23 (Asia/Seoul)

연구 branch: research/eos-original-papers-20261123

출발 commit: 0fbeeda62f078209cad03798c435933c34e2e1a1

대상: 사용자 제공 BBP1971·NV1973·Myers1966 원전, 32개 기존 core, 새로 명명한 crust 표현 6개, matching method 7개.

## 1. 결과와 채택 범위

1,344개 조합을 실제로 실행했다. 입력·실행파일 해시가 실행 전후 동일한 최종 native campaign에서 816개가 접합·직렬화·재독출 gate를 통과했다. 이어 실제 출력에 남은 crust 행을 독립 Fortran 감사기에 넣어 138개를 제외했다. 최종 수는 **678개 표현별 matched candidate**다.

| 구분 | 수 |
|---|---:|
| 실행한 core × crust 표현 × method | 32 × 6 × 7 = 1,344 |
| Native 접합·재독출 통과 | 816 |
| 실제 retained-source 열역학 모순으로 제외 | 138 |
| 최종 표현별 candidate | 678 |
| 공개 후보 / 비공개 후보 | 650 / 28 |
| 공개 후보의 정확히 같은 native CSV/NML을 묶은 그룹 | 620 |
| 기존 NV1974 fit의 의미상 재계산 | 138 |
| NV 중복 제외 + 모든 crust 상 보존 + 소스 endpoint chain인 공개 후보 | **147** |
| 위 147개의 정확한 native 표 byte identity 그룹 | **142** |

678은 독립적인 미시물리 EOS 678종이라는 뜻이 아니다. 원전 근사·finite source view·matching 처방을 고정했을 때 통과한 조합의 수다. CSV/NML의 byte identity가 다르다고 물리적으로 서로 다른 함수라는 뜻도 아니다.

후속 계산의 첫 입력 목록은 selection/PHASE_COMPLETE_NEW_PUBLIC_CANDIDATES.csv다. 이 목록은 source의 모든 선언된 crust 상과 기존 core의 선택된 마지막 행을 보존하고, NV의 기존 fit 재계산을 제외한다. 이 coverage 조건도 원전 미시물리의 전역 최소 또는 완전한 surface-to-star 검증을 대신하지 않는다. 이번 source들은 outer-only 또는 inner-only이며, 전체 outer+inner crust를 새로 하나의 통일 모형으로 생성한 것은 아니다.

| 모델 | 실행 | Native 통과 | source 감사 후 | 공개 / 비공개 | 모든 crust 상 + endpoint chain | 그중 공개 |
|---|---:|---:|---:|---:|---:|---:|
| BBP1971_original_published | 224 | 138 | **0** | 0 / 0 | 0 | 0 |
| BBP1971_original_variational | 224 | 138 | 138 | 132 / 6 | **0** | 0 |
| BBP1971_variational_maxwell | 224 | 138 | 138 | 132 / 6 | 110 | **105** |
| BPS_Myers1966_magic14 | 224 | 132 | 132 | 127 / 5 | 25 | **21** |
| BPS_Myers1966_magic20 | 224 | 132 | 132 | 127 / 5 | 25 | **21** |
| NV1973_author_fit | 224 | 138 | 138 | 132 / 6 | 132 | 126 |

Raw variational BBP의 138개는 gas-free 또는 gas-bearing의 한 상만 남긴 접합이다. 비볼록한 두 상 전체가 채택된 것이 아니다. 원래 source 표 전체는 pressure decrease와 first-law interval 위반으로 refuted다. 이 distinction은 각 record의 retained_source_extent와 별도 coverage CSV에 보존했다.

| Matching method | 최종 후보 |
|---|---:|
| first_law_constrained_hermite_pressure_bridge | 103 |
| cubic_F | 115 |
| quintic_F | 115 |
| cubic_P | 115 |
| quintic_P | 115 |
| css | 115 |
| maxwell | 0 |
| 합계 | 678 |

여기서 matching method maxwell과 crust 내부를 구성하는 BBP1971_variational_maxwell은 다른 작업이다. 후자는 neutron-drip 부근의 두 crust branch를 공통 접선으로 연결한 source이고, 전자는 그 source와 core 사이의 직접 Maxwell 접합 방법이다.

528개 NO_ADMITTED_MATCH는 등록된 finite view와 현재 search protocol에서 채택 결과가 없었다는 실행 결과다. 모든 연속적 보간이나 다른 표현이 불가능하다는 정리가 아니다. 추가 source 감사에서 뒤집힌 138개는 모두 published BBP의 첫 source interval을 보존한다. 모든 허용 bridge가 최소 3개, native Maxwell이 최소 2개의 원본 행을 보존하므로, 이 경우에만 고정된 source-view protocol의 모든 prefix를 배제할 수 있다. 다른 후보를 추가로 재탐색하지 않은 source-specific overturn은 이번 최종 결과에 0개다.

## 2. 원전 identity와 구현

사용자 제공 PDF는 원전 근거로 읽고 페이지와 표를 시각 대조했다. PDF 본문은 공개 코드 저장소에 재게시하지 않았다.

| 원전 | DOI | 업로드 PDF SHA256 |
|---|---|---|
| Baym, Bethe & Pethick 1971 | [10.1016/0375-9474(71)90281-8](https://doi.org/10.1016/0375-9474(71)90281-8) | f5d15dbc97c8dae24557b55844c2bfc550984aaa50e514e8d9412e11a05ca495 |
| Negele & Vautherin 1973 | [10.1016/0375-9474(73)90349-7](https://doi.org/10.1016/0375-9474(73)90349-7) | 5611ae66700d983e7f66881839802c2204e16f91cafbf19c670fbd178042bb5b |
| Myers & Swiatecki 1966 | [10.1016/S0029-5582(66)80001-9](https://doi.org/10.1016/S0029-5582(66)80001-9) | 23d5031429ac8d635fca19db993972d1527b3be03a05f6647042faea01078561 |

생산 물리는 새 Fortran 모듈 세 개에 구현했다.

- fortran_eos/src/compose_eos_bbp_original.f90
- fortran_eos/src/compose_eos_nv_original.f90
- fortran_eos/src/compose_eos_myers1966.f90

source-prefix 열역학 감사도 fortran_eos/tests/source_prefix_audit.f90에서 실행한다. Python 변경은 요청 계약, native 호출, 행·파일 identity 검사, 조합 스케줄링, 결과 집계와 포장이다. 기존 matching kernel과 기존 네 primary 모델의 ID·원자료·수치 허용치는 유지했다.

### BBP: 원 논문의 처방과 열역학적 완성

원전 식 (3.8), (3.9), (3.18), (3.19), (4.30), (5.3), (5.9), (5.13), (6.4)를 따라 spherical compressible liquid-drop energy를 구현했다. 원전 p.251은 surface term에 포함된 neutron-gas energy의 특정 미분항을 생략했다고 명시한다. 이에 따라 세 계산을 별도의 model ID로 보존했다.

1. Published closure: 저자들의 미분항 생략을 유지한다.
2. Full variational: 같은 energy functional에서 모든 active-coordinate 미분을 유지한다.
3. Variational Maxwell: 선택된 gas-free/gas-bearing branch에 native common tangent를 구성한다.

Published closure에서 다시 풀어 확인한 같은 상 first-law 상대 residual 최대는 **0.06507544**다. 최대 true-energy gradient는 약 0.164934 MeV다. 원 논문 수치 처방을 구현했다고 해서 정확한 cold thermodynamic barotrope가 되지는 않는다.

Raw variational branch는 neutron drip 근처에서 압력이 감소한다. 변분 Maxwell 표현에서는

\[
n_a=2.7622175767948561\times10^{-4}\ {\rm fm^{-3}},\qquad
n_b=3.1566706457246754\times10^{-4}\ {\rm fm^{-3}},
\]

\[
P_t=5.4672144575670\times10^{-4}\ {\rm MeV\,fm^{-3}},\qquad
\mu_t=939.8883473106335\ {\rm MeV}.
\]

공통 상태로 직렬화하기 전 두 endpoint의 차이는
\(|\Delta P|=4.21\times10^{-15}\ {\rm MeV\,fm^{-3}}\),
\(|\Delta\mu|=9.21\times10^{-12}\ {\rm MeV}\)였다.
같은 상에서 독립적으로 재해결한 first-law residual 최대는 \(8.995\times10^{-9}\)다.

257개 밀도 요청 중 coexistence 내부 6개를 제거하고 양쪽 endpoint 2개를 넣어 253행을 생성했다. Gas-free face Hessian 검사를 추가했고, active face 밖의 미분 placeholder도 명시했다. 탐색은 gas-free seed 1개와 gas-bearing seed 최대 5개의 유한 국소 탐색이다. 이 결과를 전체 spherical branch, 전역 convex hull 또는 pasta의 전역 최소로 확대하지 않는다.

### BBP 고밀도 표의 내부 불일치

원전 p.252 마지막 비교 행의 인쇄값 \(n_b=0.0789\ {\rm fm^{-3}}\)를 그대로 사용하면

\[
P_{\rm native}=0.40107766,\qquad P_{\rm table}=0.384
\]

로 약 +4.4473% 차이가 난다. 원전 값이나 model constant를 조정해서 이 차이를 숨기지 않았다.

같은 상태를 나타낸다고 제시된 Table 1의 \(A,k,k_g,n_N\)에는 baryon identity

\[
n_b=A n_N+\left(1-\frac{A n_N}{n_i}\right)n_g,\quad
n_i=\frac{k^3}{1.5\pi^2},\quad n_g=\frac{k_g^3}{1.5\pi^2}
\]

를 적용할 수 있다. 명목값은 \(0.0769676427403\ {\rm fm^{-3}}\)다. 독립 Fortran의 16-corner 계산에서 인쇄 반단위 오차의 범위는

\[
[0.0759399063569,\;0.0780056334808]
\]

이며 Table 2의 \([0.07885,0.07895]\)와 겹치지 않는다. A의 허용폭을 ±5로 넓혀도 상한은 0.0780391289262다.

원전의 같은 행 \(\rho=1.30\times10^{14}\ {\rm g\,cm^{-3}}\)를 제어값으로 직접 역산하면 \(n_b=0.0768969977622\), \(P=0.384202233121\)이고 압력은 인쇄 반올림 범위 안에 들어간다. A=987.889858은 990보다 0.213% 작다.

따라서 **해당 표의 여러 열이 주어진 rounding 범위에서 동시에 일치할 수 없다는 결론은 derived/numerically checked**다. \(0.0789\to0.0769\)가 오식의 정답이라는 결론은 inference-only다. 원래 0.0789는 보존했다. 증거와 실행 probe는 research_myers/bbp_cross_table/에 있다.

### Myers–Swiatecki 질량식과 BPS 확장

질량식·shell correction·pairing·변형 항을 Fortran으로 계산했다. Deformation energy의 3차식이 큰 변형에서 내려가는 runaway를 bound nucleus의 ground state로 잘못 선택하지 않도록 안정한 국소 shape minima를 분류했다.

\(2\le N\le258,\ 2\le Z\le184\)의 유한 영역에서 원전의 입자 안정성 조건을 적용한 후보는 두 magic 선택 각각 8,528개다. Magic14는 13상, magic20은 11상을 선택한다. 최종 campaign의 요청 257점에서는 phase endpoint와 내부점을 보존하여 각각 294행, 288행을 출력했다.

두 모형 모두 순수 질량공식의 저밀도 출발 핵종은 \(^{52}{\rm Cr}\)이다. 기존 표·실험 입력을 사용한 \(^{56}{\rm Fe}\) 출발 BPS와 구별한다. Magic14 결과는 BPS Table 2 마지막 9개 핵종의 binding energy를 인쇄 반단위 \(5\times10^{-5}\) MeV/핵자 안에서 재현하지만 앞의 세 핵종에는 차이가 남는다.

드립은 두 선택 모두 \(n_b=2.587036333092080\times10^{-4}\ {\rm fm^{-3}}\), \(P=4.909067899495889\times10^{-4}\ {\rm MeV\,fm^{-3}}\)다. 기존 12핵종 BPS1971_sequence를 대체하거나 UCRL의 원래 질량표와 완전히 같다고 표시하지 않았다.

### NV1973

1973 p.317의 Table 4와 식은 1974 published fit과 같은 함수다. 4,097개 밀도점에서 float64 값의 정확한 일치를 확인했다. 1973 p.311 Table 3의 1.996 MeV와 1974의 1.966 MeV는 판본별 값으로 보존했다.

추가로 1973 Table 1의 bulk nuclear potential과 해석적 gradient/Hessian을 구현했다. 상대론적 전자와 proton Slater exchange를 포함한 균일 β평형 재구성은 별도 diagnostic이다. 원전 11점의 uniform-gas energy와 최대 0.0168109599 MeV/핵자 차이가 남는다. β평형 residual은 \(1.93\times10^{-12}\) MeV 이하, Euler residual은 \(2.27\times10^{-14}\ {\rm MeV\,fm^{-3}}\) 이하다. 이 계산을 비균일 Hartree–Fock crust 재계산으로 내보내지 않는다.

## 3. 실제 retained rows의 Fortran 판정

단위는 \(n\)에 fm\(^{-3}\), 압력과 에너지밀도에 MeV fm\(^{-3}\), chemical potential에 MeV를 쓴다. 파일의 mn은 \(m_nc^2\)의 MeV 값이다. \(T=0\), neutrino-free β평형, local charge neutrality를 source contract로 사용한다. \(c_s^2/c^2=dP/d\epsilon\)이다.

Cold first law와 비감소 압력으로부터

\[
\frac{d(\epsilon/n)}{dn}=\frac{P}{n^2},\qquad
P_1\left(\frac1{n_1}-\frac1{n_2}\right)
\le \frac{\epsilon_2}{n_2}-\frac{\epsilon_1}{n_1}
\le P_2\left(\frac1{n_1}-\frac1{n_2}\right)
\]

를 얻는다. CompOSE의 \(P=nQ_1\), \(\epsilon/n=mn(1+Q_7)\)를 대입한 상관된 형태로 계산한다. Fortran real128을 사용하며 원래 binary64의 직렬화 오차에 대한 \(64\epsilon_{64}\) 상자를 명시한다. 기존 matching 허용치를 변경하지 않았다.

Source 전체뿐 아니라 native survivor 816개 각각의 실제 archive에서 crust prefix를 읽어 검사했다. Core와 crust의 mass convention이 달라 재인코딩된 Q7도 실제 출력값으로 검사한다. Maxwell의 retained count는 native chemical-potential cutoff 결과를 따르며 별도의 density cutoff를 재정의하지 않는다.

| 전체 source | 압력 감소 | first-law interval 위반 | 첫 위반 |
|---|---:|---:|---:|
| Published BBP | 0 | 221 | 1 |
| Raw variational BBP | 1 | 1 | 4 |
| Variational Maxwell BBP | 0 | 0 | 없음 |
| NV1973 fit | 0 | 0 | 없음 |
| Myers magic14 | 0 | 0 | 없음 |
| Myers magic20 | 0 | 0 | 없음 |

Positive 결과명은 FINITE_TABLE_NECESSARY_BOUNDS_SATISFIED다. 유한 표의 필요조건 통과를 continuous thermodynamic consistency의 충분조건이나 전체 미시물리의 증명으로 부르지 않는다. Source별 analytic derivative·stationarity·원전 재현 검증은 별도로 제공한다.

## 4. 고친 구현·provenance 문제와 검증

독립 검토에서 실제 결함을 발견하여 수정했다.

- BBP의 CompOSE Q5에 μe를 넣던 오류를 고쳐 β평형의 μL=0을 정확히 기록했다. Q4=-μe/mn은 유지하고 write/read 왕복으로 검사했다.
- Gas-free 후보에도 face Hessian gate를 적용했다.
- Myers DOI의 S 누락을 수정했고 수치 CSV가 수정 전후 byte-identical임을 확인했다.
- Adjudication catalog와 component 파일을 실제 campaign SHA에 연결했다.
- Maxwell prefix 추출에서 native retention count를 사용하게 했다.
- 공개 failure trace 602개를 각 projection에 실제 복사·해시 검증했다.
- 공개 trace 경로 재지정 후 RESULTS와 SURVIVORS가 달라진 packaging 실패를 보존하고, 두 inventory를 함께 갱신하여 해결했다.

검증 환경은 GNU Fortran 13.3.0, Git 2.43.0이다.

| 실행 | 관측 결과 |
|---|---|
| BBP native 전용 suite | 27 passed |
| NV native 전용 suite | 43 passed |
| Myers native 전용 suite | 19 passed |
| 최종 source-prefix 감사기 suite | 20 passed |
| 최종 고정 생산 소스의 통합·회귀 대상 suite | 195 passed in 21.52 s |
| 광범위 frontend suite | 1,969 passed, 35 skipped, 4 subtests passed; source 변경 중 build guard exit75로 1개 중단 |
| 위 중단 case | 최종 고정 build 대상 195개 suite 안에서 재실행하여 통과 |
| 최종 native matching | 1,344/1,344, 41.4867 s, input/binary hashes unchanged |
| 공개 shard readback | 650개 모든 raw tar의 크기·SHA256 일치 |
| 공개 prefix proof archive readback | 782 candidate + 6 source, 2,364개 member SHA256 일치 |
| 실제 cell_000449 추출 후 native verify-csv | PASS |

광범위 suite를 무조건 전체 PASS라고 요약하지 않는다. Exit75는 compilation 도중 source가 달라지면 stale build를 게시하지 않도록 한 기존 guard이며, 물리식 회귀 실패와 구분했다. 초기 import-path 수집 실패와 packaging inventory 실패도 실행 기록에 보존했다.

Independent Myers 검산에서는 real128 Eq.(9) scan+golden-section과 production Eq.(10)의 분류를 비교했다. 24개 사례의 최대 질량 차이는 \(4.12\times10^{-13}\) MeV다. 34개 fixed-pressure query에서 8,528후보 brute-force와 cached search의 우승 핵종이 모두 같았다. 요청 resolution과 OpenMP thread 수 변경에 따른 전이 표의 byte identity도 검사했다.

## 5. 병렬화와 경량화

Generation은 최대 4 worker, matching은 8 worker로 실행했고 각 native child의 OMP/BLAS thread는 1로 고정하여 중첩 과다 병렬화를 피했다. Myers에서는 고정압력의 전자·격자 역산을 Z별로 공유하고 유한 catalogue를 검색한다.

34개 fixed-pressure query의 단일 측정에서 cached Newton은 0.010650215 s, 핵종별 bisection은 2.980707841 s였다. 약 280배는 이 조회 작업에만 해당하며 mass catalogue 생성·I/O·전체 matching의 speedup을 뜻하지 않는다.

공개 650개 gzip 후보의 50,178,788 bytes를 원래 tar byte 그대로 6개 solid XZ shard, **5,639,068 bytes**로 포장했다. 크기는 약 8.90배 줄었고 packaging은 15.24 s였다. 출력 손실이나 부동소수점 quantization은 없다. 공개 retained-prefix proof bundle은 1,013,148 bytes이며 비공개 수치 결과를 포함하지 않는다.

## 6. 산출물과 재현

모든 상대경로의 기준은 data/eos_original_papers_20261123/이다.

| 경로 | 내용 |
|---|---|
| selection/PHASE_COMPLETE_NEW_PUBLIC_CANDIDATES.csv | 우선 검토할 공개 147개 coverage-qualified 목록 |
| selection/PUBLIC_CANDIDATES_WITH_COVERAGE.csv | 공개 650개와 retained phase·행 수·정확한 표 identity |
| selection/EXACT_TABLE_IDENTITY_GROUPS.json | byte-identical native 표 grouping |
| selection/MODEL_COUNTS.csv | 모델별 최종 집계 |
| public/qualified_campaign/MODEL_METHOD_COUNTS.csv | 42개 모델×method 집계 |
| delivery/DELIVERY_INDEX.json | 650개 후보의 shard/member/hash index |
| delivery/shards/ | 원본 candidate member를 보존한 6개 공개 shard |
| prefix_evidence/ | 공개 source와 실제 retained-prefix 증거 archive 및 index |
| generated_sources/ | 6개 native 생성·재독출 결과와 diagnostics |
| research_bbp/, research_nv/, research_myers/ | 원전 추적, derivation, independent probes, 테스트 기록 |
| evidence/DELIVERY_VALIDATION.json | inventory·privacy·trace·전체 shard readback 결과 |

예시 추출은 실제 검증했다.

~~~sh
python3 scripts/extract_eos_reaudit_candidate_20261123.py \
  --index data/eos_original_papers_20261123/delivery/DELIVERY_INDEX.json \
  --cell cell_000449 \
  --out /tmp/ns-cowling-selected
~~~

새 요청을 생성하는 public frontend 예시는 다음과 같다.

~~~sh
PYTHONPATH=python_frontend python3 -c \
  'from compose_eos.generation_dispatch import run_generation; raise SystemExit(run_generation("BBP1971_variational_maxwell", "/tmp/ns-cowling-bbp-new"))'
~~~

출력 디렉터리는 새 경로를 사용한다. Full private replay에는 기존 private core 입력이 필요하다. Public projection은 42개 비공개 조합의 상태·identity만 제공하며 endpoint·수치표·candidate 파일은 게시하지 않는다.

Matching executable SHA256는 기존 frozen binary와 같다.

~~~text
d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588
~~~

최종 source-prefix auditor binary SHA256:

~~~text
2fea55d0b4a078fd4a090e4a151af33c244361c276a428c051ce6def3ed130d9
~~~

## 7. 아직 정확히 필요한 원전

### NV의 full nonuniform Hartree–Fock

1973 원전 p.306은 density-matrix expansion의 \(B_{IJ},C_{IJ},D_{IJ}\)와 spin-orbit W를 선행 논문에 의존한다. 이번 업로드 세 논문과 확보한 자료만으로 그 수치 입력을 완전히 채울 수 없었다.

필요 원전:

**J. W. Negele & D. Vautherin, “Density-Matrix Expansion for an Effective Nuclear Hamiltonian”, Physical Review C 5, 1472–1493 (1972).**

**DOI: [10.1103/PhysRevC.5.1472](https://doi.org/10.1103/PhysRevC.5.1472)**

APS full text는 접근 제한, 확인한 HAL 항목은 No file, 다른 경로는 metadata 또는 필요한 계수 미제공이었다. 원전 본문과 부록을 제공받아야 정확한 full HF 구현을 이어갈 수 있다. 현재 상태는 published fit과 source-supported uniform reconstruction은 구현·검증, full nonuniform HF는 missing primary coefficients로 blocked다.

### 원래 Myers 질량표의 판본까지 동일하게 재현

**UCRL-11980, Appendix D, 1965-05-27**의 전체 질량표가 필요하다. 1966 논문의 reference와 설명이 지목하는 보고서다. 식으로 8,528개 후보를 생성한 작업은 완료했지만, 그 보고서의 모든 행·판본과 같은지 확인하려면 Appendix D가 필요하다. 1966 논문의 DOI는 위에 제공한 S-prefixed DOI다. 보고서 자체의 확인되지 않은 DOI는 만들지 않았다.

## 8. Claim 상태와 재개 조건

- Original-paper formula 구현, CompOSE I/O 수정, native 실행·archive identity: implementation-verified.
- BBP selected common tangent, Myers finite catalogue, table-density inconsistency: derived 및 numerically checked, 명시한 유한 범위에 한정.
- NV1973/1974 fit 동일성과 판본별 Table 3 보존: source-supported 및 implementation-verified.
- BBP 0.0769 오식 후보: inference-only.
- Full nonuniform NV HF, UCRL Appendix D와 전체 표 동일성: blocked 또는 unresolved.
- 678개 후보 전체의 full-star/TOV/Cowling 또는 전역 microphysical admission: 이번 결과의 claim에 포함하지 않음.

실행 산출물과 source identity는 보존하며 이후 원전이 추가되면 새 model/representation identity로 확장한다. 기존 published 값, source catalog, 기존 네 primary model, 이전 결과를 덮어쓰지 않는 것이 재개의 기준이다.
