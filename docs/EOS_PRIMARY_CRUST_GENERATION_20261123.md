# NS_Cowling: primary-source crust generation and matching
## 2026-11-23 — 실행 결과와 과학적 범위

이전 감사에서 자료가 없어 실행하지 못했던 BPS outer crust, BBP inner crust, Negele–Vautherin 계열에 대해 원전 수식·계수·수치표를 확보하고 기존 Fortran EOS generator를 확장했다. 기존 32개 core source identity와 7개 matching method를 사용해 원전 기반 대표 3종의 672개 조합을 모두 실행했다. BBP의 압력 적분형을 별도의 파생 표현으로 추가하여 총 896개 조합을 계산했다.

처음 native matching과 유한 표 검사를 통과한 조합은 547개였다. 이 가운데 BBP 원문표를 사용한 137개는 원문 첫 구간의 연속 cold first-law 필요조건 위반을 그대로 보존하므로 제외했다. 최종 목록은 **표현·인터페이스 수준의 후보 410개**다. 공개 가능한 후보는 393개, 비공개 core 입력에서 나온 후보는 17개다.

이 숫자는 서로 다른 source/method 조합의 출력 수다. 미시물리적으로 서로 다른 EOS의 개수, 전체 역사적 모델의 정확한 재현 수, 완전한 별의 안정성 검증 수로 해석하지 않는다. 특히 BBP 파생형과 NV fit의 278개 후보에는 저밀도 에너지 기준과 neutron-gas 해석에 관한 제한이 남아 있다.

## 1. 실제로 구현한 모델

| Model ID | 원전과 구현 | 생성 행 수 | 범위 |
|---|---|---:|---|
| BPS1971_sequence | BPS1971 수식 + Table 2의 12개 평형핵·결합에너지. 고정 압력에서 Gibbs 최소화, 11개 전이의 양쪽 상태 보존 | 291 | outer crust; 전체 UCRL 질량 후보 탐색은 미재현 |
| BBP1971_table5 | BPS1971 Table 5가 직접 수록한 BBP의 n, rho, P. 전체 35행 보존 후 명시된 n≤0.08 fm^-3 영역의 26행 출력 | 26 | 원문 유한 표의 불일치를 그대로 보존하는 감사 대상 |
| BBP1971_pressure_integral | 원문 n/P knot를 구간별 polytrope로 연결하고 첫 원문 에너지 anchor에서 first law를 해석적으로 적분 | 281 | 별도의 pressure-primary 파생 EOS; 전체 CLDM·원문 표의 정확한 복원은 아님 |
| NV1974_author_fit | Negele1974 Table III의 저자 해석식을 Fortran Horner 방식으로 평가하고 P, mu, sound speed를 해석적으로 계산 | 257 | 원래 HF 계산을 다시 푼 결과가 아닌 저자 에너지 표현 |

Python은 요청 생성, 입출력, checksum, 실행 병렬화, archive 포장과 검증을 담당한다. 생산 EOS 수치 계산은 Fortran에 있다.

주요 생산 코드:
- fortran_eos/src/compose_eos_bps_primary.f90
- fortran_eos/src/compose_eos_bbp_primary.f90
- fortran_eos/src/compose_eos_nv_primary.f90
- fortran_eos/app/compose_eos_generate.f90
- python_frontend/compose_eos/primary_model_contracts.py
- python_frontend/compose_eos/generation_dispatch.py
- python_frontend/compose_eos/native_generation.py

공개 generator catalog는 18개에서 22개 모델로 확장됐다. 새로운 family는 primary_crust다. 기존 모델과 matching 허용오차는 그대로 유지했다.

## 2. 전체 조합 결과

각 칸의 분모는 32개 core source identity다. 수치 후보에 source-first-law 제외를 적용한 결과다.

| 표현 | first-law bridge | cubic_F | quintic_F | cubic_P | quintic_P | CSS | Maxwell | 합계 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BPS1971_sequence | 17 | 23 | 23 | 23 | 23 | 23 | 0 | 132 |
| BBP1971_table5 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| NV1974_author_fit | 23 | 23 | 23 | 23 | 23 | 23 | 0 | 138 |
| BBP1971_pressure_integral | 23 | 23 | 23 | 23 | 23 | 23 | 2 | 140 |
| 합계 | 63 | 69 | 69 | 69 | 69 | 69 | 2 | 410 |

first-law bridge의 정확한 구현 이름은 first_law_constrained_hermite_pressure_bridge다.

- 기존 자료 공백 3종에 대응한 672개 조합: 270개 표현별 후보, BBP 원표 224개 source-level 제외, 그 밖에 178개 미채택.
- 별도 BBP 파생형의 224개 조합: 140개 표현별 후보, 84개 미채택.
- 전체 896개: 후보 410개, BBP 원표 source-level 제외 224개, 그 밖의 미채택 262개.
- 실행·프로토콜 실패 0개. 입력 및 matching binary hash는 전체 실행 동안 불변.
- 모든 조합을 실제 native kernel에 전달했다. Python에서 물리적 성공 가능성을 예측해 실행을 생략하지 않았다.
- 후속 endpoint 탐색이 미완료인 실패 trial 25개는 23개 BPS CSS 조합에서 발견됐지만, 그 조합들은 후속 phase view에서 이미 양성 결과를 확보했다. 음성 조합을 이 미완료 탐색만으로 확정하지 않았다.
- Maxwell 12개 조합은 4개 core의 P(mu) interpolant가 기계적 안정성 조건을 만족하지 못해 미채택됐다. 해당 core는 2020-MVCD-D1MS-full, XMLSLZ_DDME-X, XMLSLZ_NL3, XMLSLZ_PK1이다. 이는 현재 표현·보간 경로의 거부 결과이며 물리적 해의 부재 증명은 아니다.

모든 7개 방법에 후보가 있는 source pair는 2개다. 둘 다 BBP pressure-integral 표현을 사용한다.

| Core source | Maxwell P_t [MeV fm^-3] | mu_t [MeV] | n_low [fm^-3] | n_high [fm^-3] |
|---|---:|---:|---:|---:|
| native_reconstructed:GPPVA_FSU2:core | 0.00651032911910 | 941.645977656 | 0.005275221563 | 0.012992151602 |
| native_reconstructed:PK1_CANONICAL:core | 0.0232806486897 | 943.566552599 | 0.012694333925 | 0.022908413490 |

이 두 결과에도 아래의 BBP 에너지 anchor 제한이 적용된다.

## 3. BBP 원문표에서 드러난 실제 오류

BPS1971은 해당 구간의 질량밀도를 rho=epsilon/c^2로 정의하고 BBP 데이터를 Table 5에 수록한다. 따라서 baryon density를 rho/m_n으로 바꿔 쓰지 않았다. 원문이 제공하는 n을 독립적인 열로 보존했다.

영온·단일 barotrope의 first law는

\[
d(\epsilon/n)=\frac{P}{n^2}\,dn
\]

이다. 안정한 구간에서 압력이 단조 증가하면 두 점 사이에

\[
P_1\left(\frac1{n_1}-\frac1{n_2}\right)
\leq
\frac{\epsilon_2}{n_2}-\frac{\epsilon_1}{n_1}
\leq
P_2\left(\frac1{n_1}-\frac1{n_2}\right)
\]

가 필요하다.

첫 두 원문 점의 rho, n, P 모두에 마지막 인쇄 자리 반 칸의 독립적인 닫힌 오차 범위를 허용한 결과:

| 필요조건 | 허용 범위 [MeV] |
|---|---:|
| Delta(epsilon/n) | [0.608800457685, 1.649410968788] |
| 압력 적분 | [0.268434540784, 0.285452553577] |

두 범위 사이의 최소 간격은 약 0.323347904108 MeV다. 인쇄 반올림만으로 일치시킬 수 없다. 이 계산은 Fortran 생성기의 진단과 별도의 Fortran probe에서 수행했고 원전 이미지를 독립적으로 대조했다.

기존 유한 표 verifier는 큰 격자 간격의 미분 stencil masking과 비거부 C7–C9 진단을 사용하므로 이 연속 first-law 모순을 잡지 못했다. 따라서 native 유한 표 PASS를 원전의 연속 열역학적 일관성으로 취급하면 실제 false positive가 생긴다.

이번에는 수치 생존으로 기록된 BBP 원표 137개 archive를 모두 열어 첫 두 원문 점의 보존을 검사했다. n과 P는 일치했고 epsilon의 최대 상대 직렬화 차이는 2.220446049250313e-16이었다. 모든 후보가 모순 구간을 보존하므로 제외했다. 원래 kernel journal과 EOS bytes는 변경하지 않았고 별도의 판정 단계와 증거를 남겼다.

이 결과의 범위는 **해당 원문 유한 표와 명시한 인쇄 반올림 가정**이다. 원래 BBP CLDM 자체의 불가능성을 뜻하지 않는다.

BBP 압력 적분형은 P(n)=P_i(n/n_i)^gamma_i를 사용하며, epsilon/n은 첫 원문 anchor에서 해석적으로 적분한다. gamma→1 극한에는 안정한 exprel 전개를 사용한다. 원문의 모든 n/P knot를 보존하지만 에너지는 달라진다. 최대 상대 에너지밀도 차이는 약 0.1414%이고, 원문 인쇄 반 칸 대비 최대 차이는 약 18.78배다. 이를 반올림 범위 안의 복원으로 부르지 않는다.

## 4. BPS와 NV의 중요한 해석 제한

### BPS

원전의 자유 전자기체와 bcc 격자식을 사용한다. 격자식의 전하 상수는 e_G^2=alpha hbar c이며, 전자 screening·exchange 등 현대 보정은 추가하지 않았다. 프로젝트의 명시된 기본 상수와 원전의 반올림된 binding energy를 사용했다.

- 고정 isotope에서 원전 Table 5의 43개 행과 비교: rho 최대 상대차 5.54e-4, P 4.21e-4, Gamma 3.28e-4.
- 원전 Table 5의 4개 고정밀도 순수상 표점은 정확한 Maxwell 밀도 간격 안에 있다. 전이 양끝을 재현해 이 차이를 분류했고 조성 값을 강제로 맞추지 않았다.
- 계산된 neutron drip: n=2.58699792268e-4 fm^-3, rho=4.32428292293e11 g cm^-3.
- 원문 Table 5 마지막 rho와 약 0.59% 차이가 있다. 경계를 원문 표 숫자에 강제로 고정하지 않았다.
- 16점과 513점 요청의 전이 파일이 byte-identical이었다.
- 실제 132개 matched 후보에서 보존된 Gibbs 전이의 최대 |Delta mu|는 3.41e-13 MeV였다.
- 132개 중 109개는 전체 sequence view에서, 23개는 첫 phase view에서 출발했다. 전체 11개 원래 전이를 끝까지 보존하는 후보는 26개다. 나머지는 matching interface 위치에 따라 뒤쪽 crust 구간을 bridge가 대체한다.
- outer-only 입력을 이용한 출력이므로, source endpoint-chain 분류를 원래 미시적 전체 crust가 보존됐다는 뜻으로 해석하지 않는다.

### NV

Negele1974의 Table III, p.24 해석식을 사용했다. x=ln[n/(1e-4 fm^-3)], Q(x)=sum(c_i x^(i-1)), F=exp(Q) MeV라고 두면

\[
\epsilon=n[m_nc^2+c_0+F],\quad
P=nFQ',\quad
\mu_b=m_nc^2+c_0+F(1+Q'),
\]

\[
\frac{c_s^2}{c^2}
=\frac{F[Q'+(Q')^2+Q'']}{\mu_b}.
\]

생성 범위는 원전 표의 2.79e-4≤n≤0.0789 fm^-3다. 전자 에너지는 원전 총에너지에 포함되어 있으므로 다시 더하지 않았다.

- Fortran real128 5점 에너지 미분과 31개 내부점 비교: P 최대 상대차 1.35e-14, sound-speed squared 최대 상대차 1.79e-14.
- 원전 11개 표점에 대한 에너지 fit 차이는 최대 0.115411 MeV/baryon이다. 구현 roundoff와 원전 fit 오차를 구분해야 한다.
- 기존 NUDA의 1973 라벨 표는 n=0.00577 fm^-3에서 E/A-m_n=1.996 MeV, 직접 확인한 1974 원전은 1.966 MeV다. 두 출처를 별도로 보존했다. 미확보 1973 본문의 값을 추정하지 않았다.
- 표의 neutron chemical potential은 마지막 점유 HF 궤도 고유값이다. 이를 정확한 fit 제약으로 강제하지 않았다. 표로부터 만든 Euler-pressure 진단과 fit 압력의 차이는 최대 약 27.20%다.
- Negele는 원문 p.24에서 에너지 fit의 두 번째 미분이 noise를 증폭하므로 adiabatic index의 세부 모양에 과도한 의미를 부여하지 말라고 설명한다.

특히 낮은 밀도에서 NV fit은 mu_b-m_n c^2=-0.2168231 MeV이며, 257점 중 7점에서 이 값이 음수다. BBP 파생형 역시 첫 에너지 anchor에서 mu_b=938.876804478 MeV로, 명시된 현대 neutron rest energy 939.56542052 MeV보다 낮다.

따라서 이 두 inner-crust 표현의 278개 양성 결과는 **안정하고 causal한 수학적 barotrope/interface 후보**로 보존하되, 원래 neutron-gas crust의 물리적 실현까지 검증된 것으로 취급하지 않는다. 에너지 anchor를 임의로 옮기거나 불리한 저밀도 행을 잘라내지 않았다. 이 제한은 모든 개별 archive의 SOURCE_ADMISSION.json에도 포함되어 있다.

## 5. 검증·성능·재현성

전체 회귀 검사:

~~~bash
PYTHONPATH=python_frontend \
COMPOSE_EOS_NATIVE_BUILD_DIR=/tmp/ns_eos_primary_crust_20261123/build \
COMPOSE_EOS_FORTRAN_EXE=/tmp/ns_eos_primary_crust_20261123/build/compose_eos_fortran \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m pytest -q python_frontend/tests pipeline/tests
~~~

실제 결과: **1985 passed, 35 skipped, 4 subtests passed in 233.50 s**.

주요 표적 검증은 잘못된 행/열 수, trailing garbage, NaN, overflow exponent, 잘못된 isotope/계수, 외삽, 출력 덮어쓰기, 원문 knot 보존, phase boundary, 독립 미분, native CompOSE readback, 진단 checksum 변조를 포함한다. NV 입력의 1e999가 숫자 변환 단계에서 FPE trap을 일으키는 결함은 변환 전 lexical exponent 검사로 수정했다. 조성이 없는 BBP/NV 출력은 가짜 abundance 파일을 만들지 않고 명시적인 thermo-only 경로로 처리한다.

- 매칭 병렬 worker 8개, native/BLAS thread는 worker마다 1개.
- 896개 매칭 wall time: 30.800277 s. 통제된 serial 대비 speedup 측정값은 아니다.
- Matching 실행파일 SHA256은 이전 감사의 binary와 정확히 같다:
  d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588.
- NV는 Q,Q',Q''를 한 번의 결합 Horner 평가로 계산한다.
- BBP 적분형은 구간별 해석 적분을 사용한다. 생산 경로에서 Python solver를 호출하지 않는다.
- BPS 전이근은 고정된 원전 입력과 압력 조건으로 구하고 샘플링 밀도와 분리했다.
- 공개 후보 393개의 gzip 합계 31,884,948 bytes를 원본 member bytes를 보존하는 solid-XZ 두 shard, 3,126,576 bytes로 포장했다. shard 기준 약 90.19% 감소다. index와 문서 크기는 별도다.
- 410개 qualified archive 전체의 기존 member bytes, 모든 sidecar, 두 shard 안의 393개 tar payload를 독립적으로 검사했다.
- Wolfram evaluator는 연결/내부 오류로 실행하지 못했다. Wolfram 검산 성공을 주장하지 않는다. 독립 Fortran real128·유한차분 검증과 직접 유도를 근거로 사용했다.

## 6. 생성·실행·추출 방법

아래 블록은 이번 작업에서 저장소 루트로부터 실행한 명령의 기록이다. 같은 경로로 다시 실행하면 생성기·판정기는 이미 존재하는 결과를 덮어쓰지 않고 거부한다. 게시된 저장소에도 generated_sources/와 sources/가 포함되므로, fresh clone에서 이 블록을 그대로 실행하는 방식 역시 재실행 절차가 아니다. 재실행할 때는 prepare의 --generated-root, --views-root, --catalog를 모두 새 경로로 지정하고, 뒤의 matcher·adjudicator에도 그 catalog와 새 --out, --private-artifacts, --work 경로를 연결해야 한다.

~~~bash
bash fortran_eos/build.sh /tmp/ns_eos_primary_crust_20261123/build

PYTHONPATH=python_frontend python3 scripts/prepare_eos_primary_crust_sources_20261123.py \
  --native-build /tmp/ns_eos_primary_crust_20261123/build

PYTHONPATH=python_frontend python3 scripts/run_eos_all_pairs_20261123.py \
  --catalog data/eos_primary_crust_generation_20261123/SOURCE_CATALOG.json \
  --native-exe /tmp/ns_eos_primary_crust_20261123/build/compose_eos_fortran \
  --out data/eos_primary_crust_generation_20261123/campaign_v1 \
  --private-artifacts /tmp/ns_eos_primary_crust_20261123/private_campaign_v1 \
  --work /tmp/ns_eos_primary_crust_20261123/campaign_work \
  --workers 8 --n-bridge 80 --timeout 180

python3 scripts/adjudicate_eos_primary_crust_20261123.py \
  --campaign data/eos_primary_crust_generation_20261123/campaign_v1 \
  --catalog data/eos_primary_crust_generation_20261123/SOURCE_CATALOG.json \
  --out data/eos_primary_crust_generation_20261123/adjudicated_v2 \
  --private-artifacts /tmp/ns_eos_primary_crust_20261123/private_adjudicated_v2
~~~

완전한 32개 core 재실행은 비공개 core 입력도 필요하다. 공개 source catalog에는 공개 source 35개와 비공개 source 1개의 제한된 메타데이터만 포함되어 있다. 공개 숫자만으로 재실행 가능한 조합은 868개다. 이미 비공개로 실행한 나머지 28개의 상태는 공개 ledger에서 메타데이터로 확인할 수 있다.

공개 후보 한 개 추출 예:

~~~bash
python3 scripts/extract_eos_reaudit_candidate_20261123.py \
  --index data/eos_primary_crust_generation_20261123/delivery/DELIVERY_INDEX.json \
  --cell cell_000056 \
  --out /tmp/ns_primary_candidate_000056
~~~

각 후보에는 CompOSE EOS, 내부 표, native matching/readback 검증 기록, 원래 provenance와 SOURCE_ADMISSION.json이 들어 있다. 원래 생성물의 모든 member는 byte identity를 유지하며, 추가 sidecar가 정확한 채택 범위를 설명한다.

## 7. 핵심 산출물과 공개 경계

data/eos_primary_crust_generation_20261123 아래:

- generated_sources/: 실제 Fortran 생성물과 생성 영수증.
- research_bps/, research_bbp/, research_nv/: 수치 전사, 원전 유도, native 검산 코드와 결과.
- sources/: 새 표현의 lossless matching views.
- public/SOURCE_CATALOG.json: 공개 source catalog와 비공개 source의 제한된 메타데이터.
- public/native_campaign/: 원래 native 실행 결과의 공개 projection.
- public/qualified_campaign/MODEL_METHOD_COUNTS.csv: 모델·방법별 최종 집계.
- public/qualified_campaign/SURVIVOR_CATALOG.csv: 410개 후보의 식별·방법·범위 메타데이터.
- public/qualified_campaign/SOURCE_ADJUDICATION.json: BBP 원표 제외 근거와 137개 실제 prefix 검사.
- public/qualified_campaign/INNER_SOURCE_QUALIFICATIONS.json: 278개 inner-source 후보의 해석 제한.
- delivery/DELIVERY_INDEX.json 및 두 shard: 공개 393개 실제 EOS.
- review/: 독립 원전·코드·archive·공개 경계 검토 근거.

공개 ledger의 896개 기록은 공개 숫자를 포함한 868개와 비공개 metadata-only 28개로 구성된다. 비공개 candidate의 native endpoint 수치, 표, archive 경로는 공개 projection에서 제외했다. 혼합된 원래 journal과 비공개 17개 qualified EOS는 별도 백업에 보존한다.

이번 실행의 원래 3개 역사적 placeholder를 새 모형과 동일하다고 바꾸지 않았다. 새 source ID가 어떤 자료 공백에 대응하는지는 catalog에 명시했다.

## 8. 원전 및 추가로 필요한 문헌

### 현재 계산에 사용했고 확보한 자료

1. Baym, Pethick & Sutherland (1971), The ground state of matter at high densities: equation of state and stellar models, Astrophysical Journal 170, 299–317. DOI: https://doi.org/10.1086/151216
   - 공개 PDF: https://adsabs.harvard.edu/pdf/1971ApJ...170..299B
   - 사용 부분: BPS Eqs. 1, 3–5, 14–19; Table 2; Table 5의 outer/BBP 구간.
2. J. W. Negele (1974), The Equation of State of Matter at Sub-Nuclear Density, IAU Symposium 53, 1–25. DOI: https://doi.org/10.1017/S0074180900099885
   - 사용 부분: Table II p.20, Table III p.23, 에너지 fit와 압력식 p.24.
   - 기존 NV1973/NUDA 파일은 변경하지 않았다.

### 사용자 제공 자료가 필요한 남은 원전 재현

1. **BBP 원래 CLDM 함수·계수 재현**
   - Baym, Bethe & Pethick (1971), Neutron star matter, Nuclear Physics A 175, 225–271.
   - **DOI: 10.1016/0375-9474(71)90281-8.**
   - 필요한 자료: 전체 PDF의 nuclear/neutron-gas functional, surface prescription, 계수와 평형식. 본문을 확보하지 못했으므로 식 번호·세부 페이지를 추정하지 않았다.
   - 저자 유래 반올림 전 n, epsilon, P 표와 단위·에너지 기준·provenance도 중요한 대안이다.
2. **1973 NV 원문과 1974 재록의 직접 대조·전체 HF 재현**
   - Negele & Vautherin (1973), Neutron star matter at sub-nuclear densities, Nuclear Physics A 207, 298–320.
   - **DOI: 10.1016/0375-9474(73)90349-7.**
   - 필요한 자료: 원문 PDF, 특히 원래 표·functional·계수·경계조건. 1974 재록과 기존 NUDA 표의 한 행 차이를 1973 본문으로 확인해야 한다.
3. **BPS 전체 원래 mass-candidate 집합 탐색**
   - Myers & Swiatecki, UCRL-11980, 27 May 1965의 질량표.
   - 이 보고서 자체의 DOI는 확인하지 못했다. 관련 1966 논문의 DOI를 이 보고서의 DOI로 대신 쓰지 않는다.
   - 현재 구현한 Table 2 제한 모델의 실행에는 추가 업로드가 필요하지 않다.

## 9. 게시·복구 범위

작업 baseline은 7ad3040d93215db67b23041987e17ea19b6c6715이며 연구 브랜치는 research/eos-primary-crust-generation-20261123이다. 게시 대상은 이 연구 브랜치와 Draft PR이며, 기본 브랜치 병합이나 release tag 생성은 수행하지 않는다.

원래 dirty Dropbox checkout과 이전 source/catalog/journal은 보존했다. 새 source 사실, 코드, 공개 결과는 Git으로 추적하고, 전체 실행 기록·비공개 입력 유래 결과·원전 작업 PDF·runtime 증거 및 incremental Git bundle은 별도 백업으로 남긴다. incremental Git bundle을 적용하려면 위 baseline commit을 보유한 기존 저장소가 필요하다. 공개 연구 브랜치와 별도 작업 파일 백업도 함께 보존한다. 원전 정확도와 전체 stellar validation에 대한 남은 한계는 위에 명시한 범위대로 이어서 검토해야 한다.
