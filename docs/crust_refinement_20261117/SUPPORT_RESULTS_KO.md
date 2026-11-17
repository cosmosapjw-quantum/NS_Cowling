# G: F에서 남은 source support와 정밀도 실패의 실제 재계산

기준은 F commit `d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a`이다. 네 개의 구체적 실패에 대해 462개 상태점을 지정하고 native 공간 Thomas–Fermi 계산을 실제 수행했다. 674회 native 호출 중 453개 격자 상태가 기존 stationary gate를 통과했다. 연결된 유효 source component에 대한 29회 native Maxwell 호출에서 **16개 mesh/anchor 변형 출력, 네 가지 물리적 pair 유형**을 얻었다. 네 유형 모두 N64 출력이 있다. 16개를 독립된 물리 모델 조합 16개로 세지 않는다.

모든 성공 출력은 CompOSE 직렬화, 실제 readback, 독립 interface admission, 직렬화된 출력에 대한 native consistency 검증을 통과했다. F의 순수 cubic-envelope 검사 함수를 재사용한 독립 계산도 16/16 PASS다. 이는 지정한 두 connected source component의 공통 support에 대한 결과다. 원래 공개 NL3와 DDME-X 파일의 실패 판정은 그대로 보존했고, 새로 계산한 native 대응 source에는 별도 identity를 부여했다.

수치 결과의 SSOT는 [SUMMARY.json](../../data/crust_refinement_20261117/support_repair/SUMMARY.json), 전체 전이값은 [SUCCESS_INDEX.csv](../../data/crust_refinement_20261117/support_repair/SUCCESS_INDEX.csv), 원시 matching 결과는 [MATCHING.json](../../data/crust_refinement_20261117/support_repair/MATCHING.json)이다.

## 1. 무엇이 부족했는가

F의 동일 source bytes와 동일 component에서 독립적인 cubic root를 계산했다. native retained-knot 조건은 `mu_tol=1e-10*max(1,abs(mu_t))`를 사용하며, 전이점보다 낮은 low-source knot와 높은 high-source knot가 각각 두 개 이상이어야 한다. 아래 수치는 matching을 승인한 결과가 아니라 기존 실패의 진단값이다.

| F 대상 | 원래 실패 | 저밀도/고밀도 retained knots | 독립 진단 μt [MeV] | Pt [MeV fm⁻³] |
|---|---|---:|---:|---:|
| GPPVA_NL3wrL55–native PK1, candidate 1 | source knots 부족 | 1 / 752 | 943.293695726 | 0.017463803032 |
| GPPVA_TM1e–native PK1, candidate 1 | source knots 부족 | 1 / 751 | 943.317620997 | 0.017937500864 |
| public XMLSLZ_NL3–native GM1, candidate 8 | source knots 부족 | 1 / 317 | 943.462434084 | 0.048519111580 |
| public XMLSLZ_DDME-X–native GM1, candidate 4 | C2 열역학 정합성 실패 | 14 / 384 | 942.297750162 | 0.011545715907 |

첫 세 사례는 모두 **crust 쪽 support 부족**이다. core를 더 촘촘히 계산해도 이 원인은 없어지지 않으므로 F의 PK1/GM1 core를 그대로 재사용했다. 공개 DDME-X에서는 해당 component의 C2 상대오차 최대값이 약 `4.40239e-6`으로 기존 `2e-6` gate를 넘었고, `|μ_Q−μ_E|` 최대값은 약 `0.00413732 MeV`였다. Q3를 덮어쓰거나 gate를 넓히지 않았다.

근거: [DIAGNOSIS.json](../../data/crust_refinement_20261117/support_repair/DIAGNOSIS.json), [F REJECTION_DIAGNOSIS.json](../../data/new_model_maxwell_20261117/matching_validation/REJECTION_DIAGNOSIS.json). 원래 공개 파일의 고정밀도 대체 자료를 얻지 못한 유한 범위의 탐색과 source convention은 [SOURCE_PHYSICS_KO.md](SOURCE_PHYSICS_KO.md)에 따로 기록했다. 접근하지 못한 보충자료의 부재를 주장하지 않는다.

## 2. 실제 계산한 source

밀도 구간과 mesh를 field 계산 전에 고정했다. 실패점은 연결 구간을 끊으며, 실패점 사이를 보간하여 source knot를 만들어 넣지 않았다. 각 연속 accepted run을 내보낸 뒤 기존 native phase/shape 진단으로 최대 admissible component를 구했다. core 교차점에 맞춘 source trimming은 하지 않았다.

| 새 native source | nB [fm⁻³] | 실제 격자 간격 | 초기 profile/radius의 출처 | mesh |
|---|---|---|---|---|
| GPPVA_NL3wrL55 | 0.010–0.015 | 0.00025, 추가 0.01201/0.01202/0.01204/0.01208 | 기존 selected node 11 attempt 1, node 10 attempt 5 | 16, 32, 64 |
| GPPVA_TM1e | 0.006–0.012 | 0.00025 | 기존 selected node 7 attempt 1, node 6 attempt 5 | 16, 32, 64 |
| NL3_PROVIDER_PRECISION | 0.026–0.038 | 0.0005 | 공개 NL3의 WS radius를 초기값으로만 사용 | 16, 32, 64 |
| DDMEX_ORIGINAL_TABLE_RHO0 | 0.005–0.012 | 0.00025 | 공개 DDME-X의 WS radius를 초기값으로만 사용 | 16, 32, 64 |

기존 native anchor는 `selected.csv`에서 선택된 정확한 attempt 번호를 `attempts.csv`와 `profiles.csv`에 연결했다. 두 모델의 낮은 밀도 anchor는 초기 seed 명칭이 hole이었으나 최종 선택 geometry는 droplet이다. 초기 seed 명칭을 최종 geometry와 동일시하지 않았다. 공개 두 모델에서는 radius만 초기 추정에 사용했으며, 에너지·압력·화학퍼텐셜은 모두 새 field solution에서 계산했다.

native source는 T=0, 구형 WS cell의 전자 screening을 포함한 공간 TF droplet branch다. 실제 profile, field residual, KKT, radius stationarity, baryon/charge residual과 optimizer 상태를 모두 저장했다. `maxiter=1800`, `ftol=1e-14`, 기존 KKT/radius/bound gate를 유지했다. optimizer mode 8/9로 종료하되 물리 residual 조건과 bound 조건을 이미 만족하는 경우에만 반환된 유한 profile에서 **한 차례** optimizer 재시작을 허용했다. 재시작의 tolerance와 admission gate는 같다. 재시작 후에도 optimizer가 성공하지 않은 9개 상태는 불합격으로 남겼다.

질량 convention은 NL3wr와 NL3에서 mn=mp=939 MeV, TM1e에서 938 MeV, DDME-X에서 938.5 MeV다. DDME-X 대응 계산의 `n0=0.151832 fm^-3`는 명시적으로 고정한 모델 convention이다. 이를 공개 supplemental source identity가 완전히 확인되었다는 주장으로 사용하지 않는다. 평균 charge neutrality와 beta equilibrium은 field 계산의 근거이며, CompOSE 수송은 thermo-only다. 따라서 수송된 composition에 대한 C3/C4 검증을 수행했다고 주장하지 않는다.

계획과 각 source/state 해시는 [PLAN.json](../../data/crust_refinement_20261117/support_repair/PLAN.json), [REGENERATION.json](../../data/crust_refinement_20261117/support_repair/REGENERATION.json), [PREPARED_MANIFEST.json](../../data/crust_refinement_20261117/support_repair/PREPARED_MANIFEST.json)에 있다. 각 `states/<model>/<anchor>/N*/point_*/`에 실행 input, stdout/stderr, STATE.json, POINT.json을 보존했다.

## 3. N64에서 검증된 출력

아래는 각 pair 유형의 N64 성공 출력이다. 압력과 에너지 밀도의 단위는 MeV fm⁻³, 밀도는 fm⁻³, 화학퍼텐셜은 MeV다. 전체 정밀도는 CSV/JSON을 따른다.

| Native crust → 고정 F core | N64 anchor | μt | Pt | nlow | nhigh | Δε |
|---|---|---:|---:|---:|---:|---:|
| GPPVA_NL3wrL55 → PK1 | lower_adjacent_branch | 943.278256791 | 0.017162040741 | 0.011608021369 | 0.019445364164 | 7.392795049 |
| GPPVA_TM1e → PK1 | upper_F_component | 943.275885737 | 0.017115971331 | 0.008175489943 | 0.019414485077 | 10.601473090 |
| NL3_PROVIDER_PRECISION → GM1 | provider_radius_seed | 943.455020372 | 0.048225847921 | 0.031359895559 | 0.039518915270 | 7.697668107 |
| DDMEX_ORIGINAL_TABLE_RHO0 → GM1 | provider_radius_seed | 942.297837384 | 0.011547470119 | 0.008595666217 | 0.020113542449 | 10.853269864 |

전이는 P와 μ가 같은 두 상 사이의 밀도 점프를 유지한다. 네 N64 출력의 연결된 저밀도 source support 시작점은 각각 0.01075, 0.006, 0.026, 0.007 fm⁻³다. 고밀도 끝은 모두 1 fm⁻³다. 따라서 이 표들을 표면까지 완성된 neutron-star EOS로 바로 해석하지 않는다. outer-crust/surface support는 별도 연결 문제다.

출력은 `data/crust_refinement_20261117/support_repair/matching/<model>/<anchor>/N64/` 아래에 있고, 각 출력의 `RESULT.json`에는 source component, output hashes, readback, admission, consistency 결과가 함께 들어 있다. 정확한 경로는 SUCCESS_INDEX.csv에 있다.

## 4. mesh 및 anchor 의존성

N32→N64 상대 변화는 `abs(X64-X32)/abs(X64)`다. N64 성공 anchor를 그대로 사용하여 비교했다.

| Pair 유형 | Pt 변화 | nlow 변화 | nhigh 변화 | Δε 변화 |
|---|---:|---:|---:|---:|
| NL3wr–PK1 | 0.35241% | 0.61869% | 0.20777% | 1.43219% |
| TM1e–PK1 | 0.82399% | 1.09653% | 0.48545% | 1.63699% |
| regenerated NL3–GM1 | 0.05194% | 0.00208% | 0.01646% | 0.07183% |
| regenerated DDME-X–GM1 | 0.36711% | 0.22227% | 0.34788% | 0.44185% |

이 값들은 측정된 민감도이며, 새 임의 tolerance를 붙여 “완전 수렴”으로 바꾸지 않았다. N128 계산은 이 bounded pass에서 수행하지 않았다. 특히 NL3wr/TM1e의 전이 latent energy에는 N32→N64에서 약 1.4–1.6% 변화가 남는다.

두 기존 native anchor의 같은 밀도에서 에너지/바리온을 비교하면 N16에서는 NL3wr 최대 절대차 0.03498 MeV, TM1e 0.14513 MeV로 서로 다른 local branch 또는 underresolution 영향이 나타난다. 반면 N32/N64의 공통 accepted 상태점에서는 두 anchor의 E/A 차이가 모두 `9.3e-11 MeV`보다 작다. N32에서 전이점도 사실상 일치한다. 따라서 N16의 anchor 두 개를 독립된 안정 물리상 두 개로 세거나, N64에서 다른 anchor가 실패한 사실을 물리상 부재로 해석하지 않는다.

N64 NL3wr upper anchor는 nB=0.01125, TM1e lower anchor는 nB=0.00775에서 optimizer mode 8이 남아 support가 끊겼다. 이후 성분의 첫 source knot는 각각 0.0115와 0.008 fm⁻³이고, 전이점 아래 retained knot가 하나뿐이다. 정확한 독립 재진단에서 두 경우 모두 low/high=1/753이다. 반대 anchor의 실제 accepted source는 필요한 support를 제공하여 같은 고해상도 물리 branch의 matching에 성공했다. 서로 다른 anchor의 점들을 합쳐 gap을 메우지는 않았다.

전체 29회 Maxwell 호출 중 나머지 13회는 교차 없음 8회, 공통 μ support 없음 3회, retained source knots 부족 2회다. 후자의 기존 wrapper category `UNKNOWN`은 원시 결과에 그대로 보존하고, SUMMARY의 해석 필드에서 정확한 native 오류 문구와 독립 root 계산에 따라 분류했다.

## 5. 검증과 representation 오류의 분리

실제 production에서는 native 계산이 모두 끝난 뒤 첫 export에서 phase-only CompOSE metadata 오류가 났다. 중간 API는 thermo-only view를 쓰기 전에 적용했으나 writer는 non-default phase label의 eos.compo sidecar를 위해 `composition_file_present=True`를 요구했다. 최종 API는 phase sidecar를 직렬화한 다음 native 분석에 thermo-only view를 돌려준다. field solver, 에너지/압력/Q값, 초기 profile 보간, compiler build는 이 수정에서 바뀌지 않았다.

이 실패를 성공 실행으로 덮어쓰지 않았다.

| 실행 단계 | 실제 child exit | 근거와 결과 |
|---|---:|---|
| 최초 production | 1 | native 674회/462격자 상태 완료 후 export metadata 오류. 원래 stderr/계산 결과 보존 |
| 최종 API로 cached matching | 0 | TF 재호출 0회, native Maxwell 29회, verified output 16개 |
| 독립 summary/envelope 계산 | 0 | native 호출 0회, 16/16 envelope PASS, original F source와 output hashes 재검사 |

[TRANSPORT_REPAIR.json](../../data/crust_refinement_20261117/support_repair/TRANSPORT_REPAIR.json)은 중간/최종 API SHA와 `run_state`, `interpolate_initial`, `build_driver`의 동일 AST를 기록한다. 중간 API snapshot은 보존했고 기존 PLAN의 SHA를 최종 SHA로 바꾸지 않았다. 실행 근거는 [PRODUCTION_RUN.json](../../data/crust_refinement_20261117/support_repair/PRODUCTION_RUN.json), [MATCHING_RUN.json](../../data/crust_refinement_20261117/support_repair/MATCHING_RUN.json), [SUMMARY_RUN.json](../../data/crust_refinement_20261117/support_repair/SUMMARY_RUN.json)이다.

독립 envelope 계산은 각 common-μ union interval에서 cubic polynomial의 모든 실근을 검사하고, 하나의 forward switch, 양의 density jump, 양 끝의 압력 차 부호 및 native root 일치를 확인했다. 16개 출력의 최대 root 차는 `3.83125e-11 MeV`, 최대 압력 상대 residual은 `2.60664e-11`이다. 이는 동일 connected component 밖의 다른 geometry나 전역 microscopic 최소에 대한 검증이 아니다. 근거: [ENVELOPE_AUDIT.json](../../data/crust_refinement_20261117/support_repair/ENVELOPE_AUDIT.json).

## 6. 재현 및 identity

실행 runner는 [repair_maxwell_support_20261117.py](../../scripts/repair_maxwell_support_20261117.py)다. production 당시 정확한 snapshot은 `support_repair/PRODUCTION_RUNNER.py.txt`에 있다. 독립 계산은 [summarize_support_repair.py](../../data/crust_refinement_20261117/support_repair/summarize_support_repair.py)로 다시 수행할 수 있다.

| 대상 | SHA256 |
|---|---|
| Production runner | b47a18d1f114e0c29d7bd8be74f45511e991c7184c4b26a751a05ada6302d7d8 |
| Production native TF executable | 21f00ecd9cd492e37329e5690d342842c19025d4b649f97f99e84e592522c189 |
| Maxwell/consistency native executable | 3893fbe061c2564ed4138211f7ea379a973851eb6e785fad22edcec0875152b1 |
| Production 당시 중간 Python API | 9e056c9f68cada2f0561134960240b8c84f61c408b7468632cd5bfa37425c36c |
| 최종 Python API | 61ad299bfeab5ee5588758b6367f7b69e3911317257dcecea23665f08a3c4c18 |
| 독립 summary script | 1c44a5f51dd7a8ee3d4841030261b79d8e4879fa4da65e1da8d7bfaed013b912 |

binary를 저장소에 넣지 않아도 재빌드 근거가 남도록 compiler/version, 명령, 각 compile exit, 전체 Fortran source SHA를 `TF_BUILD_PRODUCTION.json`에 복사해 두었다. pilot의 v1 build receipt는 `TF_BUILD_PILOT.json`이며 main production의 674회에 pilot 4회를 합산하지 않았다.

기존 결과의 산술·해시 검증은 G worktree에서 다음처럼 수행한다.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONPATH="$PWD/python_frontend:/tmp/ns_cowling_audit_20261117_pydeps" \
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python \
data/crust_refinement_20261117/support_repair/summarize_support_repair.py
```

최종 adapter로 기존 state의 export/matching을 진행한 정확한 command는 MATCHING_RUN.json에 있다. 역사적 PLAN.json과 production runner의 bytes는 그대로 유지했다. 새 재현에는 별도의 [PLAN_FINAL_API.json](../../data/crust_refinement_20261117/support_repair/PLAN_FINAL_API.json)을 사용한다. 이 계획은 **앞으로 실행할 재현 계획이며, 원래 production의 실행 증거가 아니다.** 기존 계획과 비교하여 provenance 설명 이외에 바뀐 항목은 Python transport API의 SHA 하나뿐이다. 물리 모델, density grid, meshes, gates, Fortran source와 executable identity는 동일하다.

[launch_final_api_reproduction.py](../../data/crust_refinement_20261117/support_repair/launch_final_api_reproduction.py)는 모든 source hash와 물리 계획의 동일성을 확인한 뒤, 원래 runner의 코드를 바꾸지 않고 별도 output BASE에서 실행한다. 출력 위치는 새 디렉터리이거나 동일 prospective PLAN.json만 들어 있는 준비된 디렉터리여야 하므로 원래 production을 덮어쓰지 않는다. 준비 단계인 `prepare`는 실제 exit 0으로 실행했으며 native 호출은 0회다. 근거는 [FINAL_API_PLAN_CHECK.json](../../data/crust_refinement_20261117/support_repair/FINAL_API_PLAN_CHECK.json)이다. 새 field/Maxwell 재계산은 수행하지 않았다.

최종 API로 처음부터 재현하는 command는 다음과 같다. `run`을 `prepare`로 바꾸면 계산 없이 source/plan을 검사하고 독립 출력 디렉터리만 준비한다.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONPATH="$PWD/python_frontend:/tmp/ns_cowling_audit_20261117_pydeps" \
COMPOSE_EOS_FORTRAN_ROOT="$PWD" \
COMPOSE_EOS_FORTRAN_EXE="/home/cosmosapjw/ns_cowling_wt/stellar_consumers_20261117/build/compose_eos_fortran/compose_eos_fortran" \
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python \
data/crust_refinement_20261117/support_repair/launch_final_api_reproduction.py run \
--plan data/crust_refinement_20261117/support_repair/PLAN_FINAL_API.json \
--output data/crust_refinement_20261117/support_repair/reproduction_final_api
```

PLAN_FINAL_API.json SHA256은 `919ac0d3aaff303be0e76d28abda763dd2ea812234bff67e5fe83716a831968a`, launcher SHA256은 `7692378604cb02d12833714a8ff2f606f1c8b58eb3ab90f4944e97ab04d2b444`다. 원래 PLAN.json의 SHA256 `0eef613ca56cefc9015c4c3d4fa2015b0720893b475f30a01e11e9b88d1cdca0`는 바뀌지 않았다.

## 7. 결론의 범위와 다음 생산 계산

G support 작업의 완료된 산출은 네 가지 pair 유형의 실제 sharp Maxwell 출력과 그 source/실행/직렬화 근거다. 기존 native NL3wr/TM1e는 실제 source 재계산으로 support 부족을 해결했고, 공개 NL3/DDME-X 문제는 새로운 native 대응 모델로 별도 성공 사례를 얻었다. 원래 공개 파일의 retained-knot/C2 실패는 해결됐다고 표시하지 않는다.

다음 정확도 개선에서는 네 N64 대표 branch의 N128 비교, 남은 optimizer mode 8 상태점의 원인과 유효 support, 저밀도 outer-crust/surface 연결을 따로 다루면 된다. 현재 데이터만으로 전 geometry의 global minimum, 공개 보충자료와의 완전 source identity, stellar surface까지의 완성 EOS를 주장할 수 없다.
