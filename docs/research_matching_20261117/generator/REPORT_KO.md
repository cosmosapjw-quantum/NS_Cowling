# 공개 generator 계약과 18모델 재사용 — 2026-11-17

기준 소스는 `478c7b5de1514ac512f91156f36551478f36f688`이다. 이 lane은 공개 generator 입출력 계약 C2–C8/R3의 병렬 검증으로, 기존 18모델 물리 캠페인을 다시 실행하지 않았다. 기존 자료의 446개 등록 output hash(공개 cohort 278 + TF source-grid 168)를 현재 파일에 대조했고 모두 일치했다. 18모델 native 생성 성공은 유지되며, 전체 원저자 EOS 동일성은 여전히 부분 완료다. 전체 microscopic BHF/ETFSI 재계산이나 CompOSE 숫자 일치를 독립 생성의 공통 필요조건으로 추가하지 않는다.

## 실제 결함과 수정

`native_generation.verify_native_generation`은 `SUPPORT.json`을 hash 대상에 넣지 않았다. `RESULT.json.output_sha256`에 SUPPORT hash를 명시한 뒤 `cross_phase_interpolation=UNAPPROVED`를 `APPROVED`로 바꾸거나 파일을 삭제해도 검증이 성공했다. 원인은 optional 진단 loop가 `sound_speed.csv`와 `attachment_diagnostics.json`만 다루고, 비교가 실제 수집한 hash의 키만 순회하기 때문이다. 이는 물리 해 실패가 아니라 지원구간 계약의 provenance 누락이다.

SPG 구현 lane이 기존 optional 증거 처리 두 곳에 `SUPPORT.json`을 추가했다. 생성된 파일은 hash를 얻고, 이미 봉인된 파일의 변경·삭제는 거부한다. 과거의 미봉인 receipt는 읽을 수 있다. 여기서 checksum은 내용 보존만 확인하며 지원구간의 물리적 승인은 별도다.

새 테스트 파일은 `python_frontend/tests/test_generation_closure_20261117.py`다. 수정 전 실제 결과는 **3 failed, 15 passed**, 수정 후 관련 네 파일 합계는 **83 passed in 0.79s**였다. 세 실패는 SUPPORT 누락·변경·삭제 반례다. 실제 Fortran reader로 Q1–Q7 각각의 NaN/Infinity 14개를 거부했고, 과거 DD2 nonuniform 출력의 미정의 Yn/Yp/A/Z/Yh mask와 C3/C4 비적용 상태도 보존됨을 확인했다. 음수/0 neutron mass 및 음수 proton mass는 기존 reader가 이미 거부하므로 새 결함으로 올리지 않았다.

## 18모델의 현재 closure와 남은 범위

표의 행 수는 기존 공개 CLI 실행 출력이다. TF 추가 열은 별도 기존 source-grid 실행 행 수다. 모든 행의 생성·readback은 원 실행 증거를 재사용했고, 이번에는 등록 파일의 byte identity를 직접 확인했다. 이 표는 새로운 물리 PASS나 연속 지원구간 인증이 아니다.

| 모델 | 공개 native 행 | 추가 TF 행 | 생성 표현과 남은 범위 |
|---|---:|---:|---|
| APR_APR | 163 | — | 출판 21계수 functional + bulk two-charge Gibbs; 원 기본구간 고밀도 C10 실패. 인과 경계 약 0.86895315 fm⁻³를 넘는 사용은 별도 제한 |
| BBB_BHF2 | 161 | — | 출판 SNM/PNM 에너지의 parabolic/PCHIP 재구성; 저자 smooth-fit/미분 동일성 미해결, 원 구현 .7<n≤.8도 .9 knot에 의존. 별도 source-only9행 표현과 구분 |
| GM_GM1 | 161 | — | 출판 rounded nonlinear RMF, 공통 정지에너지 939 MeV, npe; 조사 기본구간에서 기존 consistency 통과 |
| RG_KDE0v1 | 161 | — | 출판 비대칭 Skyrme functional; 원 기본구간 고밀도 C10 실패, 인과 경계 약 1.41096589 fm⁻³ |
| RG_SLY4 | 161 | — | 출판 비대칭 Skyrme functional; 원 기본구간 고밀도 C10 실패, 인과 경계 약 1.26500567 fm⁻³ |
| SPG_M1 | 17 | — | GDFM spherical CLD, 명시적 candidate surface; 원저자 calibration/변분·export 동일성 미확정 |
| SPG_M2 | 17 | — | 같은 표현; 원저자 calibration/변분·export 동일성 미확정 |
| SPG_M3 | 17 | — | 같은 표현; 원저자 calibration/변분·export 동일성 미확정 |
| SPG_M4 | 17 | — | 같은 표현; 원저자 calibration/변분·export 동일성 미확정 |
| SPG_M5 | 17 | — | 같은 표현; 원저자 calibration/변분·export 동일성 미확정 |
| GPPVA_DD2 | 3 | 66 | RMF spatial TF 유한 geometry/seed/mesh 최소화; 연속 mesh·global minimum·미표본 보간 미인증 |
| GPPVA_DDME2 | 3 | 71 | 같은 수치 범위; 연속 mesh·global minimum·미표본 보간 미인증 |
| GPPVA_TW | 3 | 72 | 같은 수치 범위; 연속 mesh·global minimum·미표본 보간 미인증 |
| GPPVA_FSU2 | 3 | 52 | 같은 수치 범위; 연속 mesh·global minimum·미표본 보간 미인증 |
| GPPVA_FSU2R | 3 | 79 | 기존 압력감소→slab/tube 공존 후보; 최신 기준 N256 내부 tube 표본 HOLD는 TF lane의 현재 후속 결과와 합쳐 판정 |
| GPPVA_FSU2H | 3 | 83 | 기존 N16 압력감소와 정련 결과 보존; 유한 표본을 모든 밀도 승인으로 확장하지 않음 |
| GPPVA_NL3wrL55 | 3 | 81 | 공개 n=.075,N16은 noncanonical_radial_profile, Iphase=9; canonical/uniform으로 바꾸지 않음 |
| GPPVA_TM1e | 3 | 86 | slab/tube 공존 후보의 N256 근거 존재; global minimum 또는 전체 연속 crust 동일성 증명은 아님 |

공개 cohort의 main 행은 916개다. readback segment 키를 모두 세면 52개 표이고, main과 byte-identical인 TF `mesh_16` 8개를 중복 제외하면 기존 인계의 44개 표다. 합계 정의의 차이를 생성 성공 수의 불일치로 해석하지 않는다. TF source-grid는 8모델 590행이다. 각 모델의 실제 밀도 구간, phase ID, 원 receipt 경로·hash, 모든 세그먼트 행 수는 `data/research_matching_20261117/generator/RETAINED_COHORT_MATRIX.json`에 있다.

## 새 조건부 SPG 출력

병렬 SPG lane의 repaired native 실행으로 `public_composition_inferred_20261116` 5모델×129행, 총 645행이 생성되었다. 이 lane에서도 RUN_RECEIPT 및 독립 검산 SUMMARY의 SHA256, 각 states.csv identity를 대조하고, **실제 공개 request/verify API로 다섯 결과를 다시 읽어 모두 129행을 확인**했다. 재생성은 반복하지 않았다. 확인 결과는 `INFERRED_SPG_READBACK.json`과 모델별 matrix 항목에 추가했다.

| 모델 | 새 조건부 행 | 시작 n (fm⁻³) | 끝 n (fm⁻³) |
|---|---:|---:|---:|
| SPG_M1 | 129 | 0.000233173690936 | 0.0765391794443 |
| SPG_M2 | 129 | 0.000241621932934 | 0.0892191511393 |
| SPG_M3 | 129 | 0.000235541029948 | 0.0811501611471 |
| SPG_M4 | 129 | 0.000247560285015 | 0.0988255585432 |
| SPG_M5 | 129 | 0.000236112742687 | 0.0730367964506 |

조건부 조성 역추정 surface input을 명시한 생성·readback은 완료되었으나, 원저자 calibration과 변분/export 동일성은 미확정이다. 새 시작점들도 lower-source-fit의 .002 끝과 겹친다. 마지막 상태의 `fixed_density_source` 교차는 에너지 동일 밀도 교차이며, 실제 transition.json에 남아 있는 ΔP와 Δμ가 0이 아니므로 common-tangent 접합으로 승격하지 않는다. 독립 검산의 scope는 별도 명시적 CLD 에너지와 공유 native bulk/lepton 함수, 전체 645행의 에너지·고정좌표 열역학 및 선택한 constrained derivative/reoptimized 검사다. 원저자 코드의 독립 구현으로 부르지 않는다.

## Matching으로 넘길 때 보존할 의미

공개 TF smoke의 밀도는 .003, .030, .075 fm⁻³이고 lower-source-fit 끝은 .002이므로 실제 gap이 있다. SPG 시작점은 약 .000193862–.000205001이고 같은 lower 끝 .002와 겹친다. 따라서 두 계열 모두 명목 .002 경계만으로 단순 연결할 수 없다. `cross_phase_interpolation`과 `cross_segment_interpolation`은 기존 출력에서 UNAPPROVED다. 별도 common-tangent/branch 지원 근거를 적용할 때만 그 구체적인 인터페이스가 추가로 승인된다.

TF 비균일 cell은 total profile 및 lepton fraction을 제공하지만, free/cluster 분해 mask가 없다. 이에 대한 0은 관측된 0이 아니다. 실제 readback에서 C3/C4가 적용되지 않는 상태를 테스트했다. NL3wrL55의 phase9도 같은 제한에 속한다. 원 source-energy reference와 CompOSE header normalization은 구분해야 하며, 이 lane은 에너지 상수를 조정하지 않았다. 접합 방법의 수리·기준은 matching lane 소유다.

공개 17행 SPG와 3행 TF의 기존 C10/C11 실패는 미분 stencil 지원 부족을 포함하며, 물리적 불안정의 일반 증명으로 쓰지 않는다. 반대로 APR/KDE0v1/SLy4 고밀도 비인과성은 별도 native 해석/차분 연구로 확인된 모델 제한이다. 두 원인 분류를 유지한다.

## 재현과 provenance

실제 읽기 binary는 부모가 확인한 `/tmp/ns_blockers_a49aaf9_20261116/build/compose_eos_fortran/compose_eos_fortran`을 `COMPOSE_EOS_FORTRAN_EXE`로 지정했다. 1 worker, OMP/OPENBLAS 1 thread를 사용했다. 새 build, 물리 campaign, git 조작, 원 자료 변경은 하지 않았다.

```bash
PYTHONPATH=python_frontend \
COMPOSE_EOS_FORTRAN_EXE=/tmp/ns_blockers_a49aaf9_20261116/build/compose_eos_fortran/compose_eos_fortran \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 -m pytest -q \
  python_frontend/tests/test_generation_closure_20261117.py \
  python_frontend/tests/test_native_generation_receipts.py \
  python_frontend/tests/test_native_generation_architecture.py \
  python_frontend/tests/test_generation_dispatch.py
python3 data/research_matching_20261117/generator/summarize_retained_cohort.py
```

최종 재실행도 **83 passed in 0.88s, exit 0**이고 대상 Python 소스·테스트 identity는 실행 전후 같았다. binary 및 소스 hash, 실제 명령은 `EXECUTION.json`에 있다. 근거는 `red_support.log`, `green_support.log`, `final_support_tests.log`, `INFERRED_SPG_READBACK.json`, `RETAINED_COHORT_MATRIX.json`이며 기존 물리 해석의 authority는 `docs/research_blockers_20261116/{REPORT_KO,HANDOFF_KO}.md`, `docs/research_coexistence_20261116/{HANDOFF_KO,TF,SPG}.md`와 각 원 receipt다. 이번 전체 task의 이후 SPG/TF 결과는 부모 통합 보고서가 별도로 결합한다.

사용한 연구 지침: 실제 cloud PhysMath Research Loop 1.3.0의 SKILL.md, POLICY_INTEGRATION.md, SUBAGENT_MODE.md, CODING_LOOP.md 및 systematic-debugging/verification-before-completion. 한 병렬 primary 검증 pass 안에서 실제 반례→최소 수정→해당 회귀로 종료했으며 추가 reviewer tree를 만들지 않았다.
