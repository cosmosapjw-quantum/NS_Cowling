# G: native crust 재생성, mesh 및 sharp Maxwell 결과

작성일: 2026-11-17 UTC. 기준 F commit: `d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a`.
수치 SSOT: `data/crust_refinement_20261117/tf_regeneration/SUMMARY.json`.

## 1. 결과와 주장 범위

**DDME2 native crust–PK1 native core와 새 canonical TM1 native crust–GM1 native core에서 실제 재생성한 stationary crust로 sharp Maxwell matching이 성립했다.** DDME2는 N=16,32,64,128, TM1은 N=16,32,64의 7개 유한 mesh variant가 각각 native interval 검사, Maxwell 구성, 독립 admission, 물리량·metadata write/readback 및 직렬화 후 native consistency 검사를 통과했다. 두 모델의 가장 미세한 source는 각각 9개 요청 밀도 전체가 stationary이며, source interval rejection이 없다. [R1–R3]

이 결과는 계산한 TF branch와 밀도 구간에 관한 **numerically checked / implementation-verified** 결과다. 모든 geometry와 초기 basin을 망라한 global minimum, 연속체 극한의 증명, 공개 CompOSE crust와의 numerical identity, 완전한 stellar EOS를 뜻하지 않는다. 원래 F가 고른 DDME2 N16 state와 G의 첫 초기 basin이 달랐음을 확인하고, 정확한 F-selected profile을 이용한 네 개 anchor control을 추가했다. N16에서는 서로 다른 stationary droplet이지만, N32/64/128의 같은 anchor에서 두 초기 basin은 측정된 수치 정밀도 안에서 같은 상태로 모인다. 이 확인을 전체 밀도 구간의 유일성으로 확장하지 않았다. [R4]

가장 미세한 admitted 결과는 다음 디렉터리다. 각 `RESULT.json`에 native report, 독립 admission과 readback 결과가 있고, 그 아래 `compose/`에 실제 표가 있다.

- `tf_regeneration/tables_v2/DDME2/N128/block_00/match_00_00/`
- `tf_regeneration/tables_v2/TM1_CANONICAL/N064/block_00/match_00_00/`

두 표는 유한한 lower density에서 시작한다. Surface/outer closure와 sharp-interface stellar equations는 별도 검증 대상이며, 이 문서의 성공 수에 포함하지 않는다.

## 2. 정의에서 계산까지

길이는 fm, baryon density는 fm⁻³, 에너지·chemical potential은 MeV, 압력·에너지 밀도는 MeV fm⁻³를 쓴다. 기존 코드의 `tf_hbarc`를 그대로 사용하고, nucleon rest energies와 CompOSE reference masses를 구별한다. 표의 reference mass를 맞추기 위해 물리적 energy density를 이동시키지 않았다.

기존 `compose_eos_spatial_tf.f90`의 screened-electron TF functional과 field equations를 재사용했다. Cell-centered finite-volume grid의 정규화된 cell weight는 차원 \(d\)에 대해

\[
w_i=(i/N)^d-[(i-1)/N]^d,\qquad \sum_i w_i=1
\]

이다. 제약은 \(\langle n_n+n_p\rangle=\bar n\), 셀 평균 charge neutrality \(\langle n_p-n_e\rangle=0\)이다. 각 위치의 전하밀도를 0으로 강제하지 않는다. Meson/Coulomb field는 중심·cell boundary의 Neumann 조건, Coulomb mean-zero gauge, 공간적으로 screening되는 electron profile을 유지한다. Density-dependent coupling에는 local density와 해당 rearrangement 항을 사용한다. 문헌 convention 대응과 field energy 부호는 [SOURCE_PHYSICS_KO.md](SOURCE_PHYSICS_KO.md)의 S1–S3/C1–C2를 따른다.

고정 평균 밀도에서 nucleon profile과 radius를 변화시킨다. 내부 active degree of freedom에서 beta equilibrium은 \(\mu_n=\mu_p+\mu_e=\mu_B\)이며, 0-density boundary에서는 native KKT 조건을 쓴다. 사용한 압력은

\[
P_E=\bar n\mu_B-\epsilon.
\]

Stationary branch에서는 envelope theorem으로 \(d\epsilon_*/d\bar n=\mu_B\)이므로 \(P_E=\bar n^2d(\epsilon_*/\bar n)/d\bar n\)이다. 이 관계를 코드에 입력한 Euler identity만으로 확인했다고 하지 않고, 별도로 재최적화한 인접 밀도 상태의 에너지 미분과 비교했다. [R5]

Native muon guard \(\max_r[\mu_e(r)+V(r)]>m_\mu c^2\)는 유지했다. 가장 미세한 DDME2/TM1 grid에서 이 최대값은 각각 53.100234/62.030326 MeV이고 threshold-crossed flag는 모두 0이다. Electron-only native branch의 적용 범위에 관한 확인이며, 공개 표 전체에서 muon을 무시해도 된다는 주장이 아니다. [R1]

| 항목 | 구현 또는 입력 | 확인 |
|---|---|---|
| Native functional/fields | 기존 `compose_eos_spatial_tf.f90` | 파일 불변; derivative, PDE, constraint 기록 |
| 기존 DDME2 coupling | 기존 TF registry case | 변경하지 않음; homogeneous regression 유지 |
| Canonical TM1 | F의 `PARAMETER_INPUTS.json` | \(M_n=M_p=938\) MeV, frozen coefficients 그대로 |
| NL3 / DDME-X 지원 case | `NL3_PROVIDER_PRECISION`, `DDMEX_ORIGINAL_TABLE_RHO0` | Additive case; 원 DDME-X NMP identity는 unresolved 유지 |
| Thin adapter | `fortran_eos/tools/tf_candidate_driver.f90` | Native field/optimizer 호출 및 scalar/profile 출력 |
| 실행·증거·직렬화 | `python_frontend/compose_eos/tf_candidate.py` | Input/stdout/stderr/STATE, SHA, 실제 child exit 기록 |
| 생산 grid | `scripts/run_crust_refinement_20261117.py` | Continuation, connected support, 독립 matching |
| 정확한 F-selected control | `scripts/check_tf_selected_anchor_20261117.py` | 같은 역사적 profile에서 네 mesh를 독립 초기화 |

TM1/NL3/DDME-X는 기존 registry case를 덮어쓰지 않고 세 case만 추가했다. Insert를 제거하면 원 registry source가 byte-identical하게 복원됨을 `REGISTRY_ADDITION.json`에 기록했다. Scalar convention은 F의 \(a_3=-g_2\hbar c,\ a_4=g_3\), isovector vertex는 \(g_\rho=2\Gamma_\rho\)를 그대로 매핑했다. 새 case의 homogeneous fixed-profile limit은 F parameter-driven native RMF와 \(\epsilon,P,\mu_B,Y_p\)를 실제 비교했다. 자유 변분으로 uniform seed가 droplet로 변하는 경우와 fixed-profile 평가를 혼동하지 않았다. [R6–R7]

## 3. 사전에 정한 계산과 실제 실행

| 모델 | 밀도 grid (fm⁻³) | 첫 mesh | 관측 후 추가 | 초기 profile 역할 |
|---|---|---|---|---|
| DDME2 | .006–.014, 간격 .001, 9점 | 16,32,64 | N128의 동일 9점 | 기존 DDME2 droplet-seeded checkpoint |
| TM1_CANONICAL | .036–.052, 간격 .002, 9점 | 16,32,64 | 추가 full grid 없음 | TM1e profile/radius만 초기값; 실제 functional은 canonical TM1 |

창은 F의 유망 transition을 둘러싸도록 고정했다. 이후 core에 맞추어 energy, pressure, coupling, 허용오차를 조정하지 않았다. Source validity와 connected support를 core와의 matching 전에 독립 검사했다.

SLSQP `maxiter=1800`, `ftol=1e-14`를 유지했다. Mode 8/9가 반환한 실제 profile이 있을 때만 optimizer state를 초기화하는 **한 번의 warm restart**를 허용했다. 원 실패는 별도 directory에 남고, 같은 physical gates와 `ftol`을 다시 적용했다. Native status, 실제 optimizer success, baryon/charge errors, muon guard 모두를 검사한다. Radius/density search boundary에 닿은 상태를 admissible data로 바꾸지 않는다.

전체 요청은 grid 63개, rod 2개, density derivative 8개, F-selected anchor 4개로 77개다. Warm restart 25회를 합쳐 실제 native invocation은 102회였고, 최종 76개 요청이 stationary다. 남은 하나는 DDME2 N64의 \(n=.007\): 한 번 재시작 뒤에도 mode 8, `REQUIRES_REFINEMENT`여서 버렸다. Native process exit 0과 physical acceptance를 구별한다. 102회 중 26개 초기/재시작 상태는 refinement-needed였으며 이 기록을 보존했다. [R1, R8]

N64의 실패 지점을 보간하지 않았다. 해당 mesh의 연결된 export support는 .008–.014이고, 단독 .006 state를 이 구간에 연결하지 않았다. N128에서 .007을 포함한 9개를 별도로 실제 풀어 .006–.014 support를 회복했다. N16에서는 \([.008,.009]\) interval 하나가 `INTERPOLANT_SHAPE_FAIL`이므로 두 source component로 나뉜다. 낮은 component는 core와 common \(\mu\) support가 없어 거절됐고, 높은 component의 matching만 admitted다.

## 4. Mesh sensitivity와 stationarity

Anchor는 DDME2 \(n=.01\), TM1 \(n=.044\) fm⁻³이다. 아래는 production branch 값이다.

| 모델 | N | P (MeV fm⁻³) | \(\mu_B\) (MeV) | R (fm) | Chemical KKT (MeV) |
|---|---:|---:|---:|---:|---:|
| DDME2 | 16 | .013998089877 | 943.101559457 | 30.51944906 | 6.90e−6 |
| DDME2 | 32 | .014627914677 | 943.166293853 | 33.11796641 | 1.51e−5 |
| DDME2 | 64 | .014472607153 | 943.152238240 | 32.41041569 | 3.84e−5 |
| DDME2 | 128 | .014445908694 | 943.150006648 | 32.31410479 | 5.04e−5 |
| TM1_CANONICAL | 16 | .132714529507 | 945.249090792 | 18.54814977 | 6.91e−6 |
| TM1_CANONICAL | 32 | .132525742609 | 945.245200526 | 18.49862358 | 5.57e−4 |
| TM1_CANONICAL | 64 | .132515973150 | 945.245068296 | 18.51109748 | 1.08e−4 |

DDME2 N32→64의 anchor pressure 차이가 약 1.06%였으므로 N128을 선택적으로 계산했다. N64→128에서도 약 .184%의 차이가 남아 있어 continuum convergence를 선언하지 않는다. TM1 N32→64 anchor pressure 차이는 약 \(7.37\times10^{-5}\) relative이며, 추가 전체 N128 grid는 수행하지 않았다.

| 가장 미세한 source grid의 최대값 | DDME2 N128 | TM1 N64 |
|---|---:|---:|
| Chemical KKT (MeV) | 7.30e−5 | 4.13e−4 |
| \(\lvert d\epsilon/d\log R\rvert/\bar n\) (MeV) | 9.49e−8 | 9.26e−8 |
| Native scalar PDE residual | 8.44e−11 | 1.80e−10 |
| Native vector PDE residual | 1.04e−10 | 1.03e−10 |
| Native electron PDE residual | 4.54e−14 | 1.61e−14 |
| Baryon constraint error (fm⁻³) | 8.28e−16 | 1.81e−16 |
| Charge constraint error (fm⁻³) | 1.66e−19 | 3.94e−19 |
| Radius/density bound 및 muon threshold flag | 모두 0 | 모두 0 |

PDE residual은 native field equation의 정규화 convention을 따른다. 이를 압력 오차나 continuum error bar로 바꾸지 않았다.

### 에너지 미분으로 확인한 압력

각 anchor의 \(\bar n\pm h,\bar n\pm2h\), \(h=10^{-5}\) fm⁻³를 실제 다시 최적화했다. 8개 모두 같은 droplet geometry에서 stationary이며,

\[
P_h=\bar n^2\frac{(\epsilon/n)(\bar n+h)-(\epsilon/n)(\bar n-h)}{2h}
\]

와 \(P_{2h}\)를 계산했다.

| 모델/mesh | \(P_E\) | \(\lvert P_h-P_E\rvert\) (MeV fm⁻³) | Relative 차이 | \(\lvert P_h-P_{2h}\rvert\) |
|---|---:|---:|---:|---:|
| DDME2 N128 | .0144459086937 | 4.711e−9 | 3.261e−7 | 1.166e−8 |
| TM1 N64 | .132515973150 | 2.600e−10 | 1.962e−9 | 1.843e−9 |

이 비교는 정의한 Euler identity와 별개의 local derivative check다. 전체 density interval의 고차 regularity나 geometry switch의 미분 가능성을 증명하지 않는다. [R5]

## 5. 실제 sharp matching과 mesh 이동

모든 값은 native Maxwell report에서 가져왔다. \(\mu_*,P_*\)의 equality, density jump와 \(\Delta\epsilon=\mu_*\Delta n\)를 검사하며, reference mass 차이는 physical-energy-preserving Q 재인코딩으로 처리한다.

| Crust→core | N | \(\mu_*\) (MeV) | \(P_*\) (MeV fm⁻³) | \(n_-\) (fm⁻³) | \(n_+\) (fm⁻³) | \(\Delta\epsilon\) (MeV fm⁻³) |
|---|---:|---:|---:|---:|---:|---:|
| DDME2→PK1 | 16 | 943.110335315 | .014086040601 | .01004379865 | .01715672903 | 6.708278162 |
| DDME2→PK1 | 32 | 943.108976643 | .014062743373 | .00972160751 | .01713739549 | 6.993896212 |
| DDME2→PK1 | 64 | 943.106509515 | .014020506557 | .00977341248 | .01710225763 | 6.911881567 |
| DDME2→PK1 | 128 | 943.105867281 | .014009525845 | .00977355895 | .01709310405 | 6.903105930 |
| TM1→GM1 | 16 | 945.252555232 | .132866995372 | .04401760900 | .05352328539 | 8.985264900 |
| TM1→GM1 | 32 | 945.250709261 | .132768203652 | .04402783884 | .05351171506 | 8.964640725 |
| TM1→GM1 | 64 | 945.250293008 | .132745929778 | .04402652029 | .05350910552 | 8.963416462 |

DDME2 N64→128에서는 \(P_*\)가 약 \(7.83\times10^{-4}\) relative 이동하고, \(n_-\)는 약 \(1.50\times10^{-5}\) relative 이동했다. TM1 N32→64의 \(P_*\) 이동은 약 \(1.68\times10^{-4}\) relative이다. 숫자상 transition이 유지되는 것은 확인했지만, 이 차이를 엄밀한 extrapolation error bar로 해석하지 않았다.

가장 미세한 DDME2 matching의 low/high output row 수는 5/762, TM1은 6/289다. Maxwell output은 원 source의 물리적 thermodynamic rows를 보존하고 새 interface row를 생성한다. Abundance 없이 phase label만 직렬화하므로 C3/C4는 `NOT_AVAILABLE`이고, 실제 native profile에서 구한 cell neutrality와 별개다. C7–C9는 finite-difference diagnostics이며 hard gate로 바꾸지 않았다. [R2–R3]

## 6. Basin/geometry control이 실제로 확인한 것

DDME2의 역사적 `selected.csv`는 \(n\simeq.01\)에서 attempt 5를 고른다. Initial seed는 `hole`이지만 최종 geometry는 `droplet`이다. G 생산의 첫 profile은 dimension 3, initial seed `droplet`인 attempt 1이었다. Initial seed와 최종 geometry, source selection을 구분해야 한다. `SEED_IDENTITY.json`이 실제 row와 hashes를 보존한다.

정확한 F-selected profile을 각각 N16/32/64/128에 독립 보간해 같은 density, model, executable, gate로 다시 풀었다. Production state는 교체하지 않았다. 아래 차이는 selected-control minus production이다.

| N | \(\Delta(E/A)\) (MeV) | \(\Delta P\) (MeV fm⁻³) | \(\Delta R\) (fm) | Nucleon profile weighted RMS 차이 (fm⁻³) |
|---:|---:|---:|---:|---:|
| 16 | −4.947804e−3 | +1.199692e−3 | +4.940425 | 3.485315e−3 |
| 32 | −4.511946e−11 | +1.160972e−9 | +1.031930e−5 | 6.670371e−9 |
| 64 | −5.027090e−11 | +7.037180e−10 | −8.021870e−7 | 7.156779e−9 |
| 128 | −6.981082e−11 | −2.920295e−10 | −8.454235e−6 | 4.301451e−9 |

따라서 N16 production을 F-selected minimum의 exact refinement 시작점이라고 부르지 않는다. N32 이상에서 **이 anchor의 두 basin은 관측 가능한 수치 차이 안에서 같은 stationary state와 양립한다.** 전체 density grid의 basin uniqueness나 모든 local minimum의 소멸은 미확인이다. [R4]

별도의 N64 rod control은 두 모델 모두 stationary였다. 같은 anchor에서 rod minus droplet의 \(E/A\)는 DDME2 +.03737935 MeV, TM1 +.00441493 MeV다. 이 두 geometry 비교는 droplet 선택의 제한된 지지이며, slab/tube/bubble 및 추가 profile basin까지 포함한 global minimum 판정이 아니다. TM1e에서 가져온 초기 profile 역시 canonical TM1 physical parameters를 대체하지 않는다.

## 7. 보존한 실패와 생산적 수정

1. **초기 regression 설계 오류.** 자유 변분에 uniform seed를 넣고 uniform 상태가 그대로 남는다고 가정한 테스트가 실패했다. 실제로 droplet로 가는 것이 가능한 문제이므로, homogeneous functional 비교는 `evaluation_only=True`의 fixed-profile native evaluation으로 분리했다. 이 mode는 `PROFILE_EVALUATED`이며 stationary data로 export되지 않는다. Physical optimizer를 바꾸지 않았다. 수정 뒤 네 물리 regression이 실제 통과했다.
2. **Phase-only transport 표현 오류.** 최초 export의 abundance presence/applicability masks는 이미 false였다. 문제는 phase sidecar를 위한 `composition_file_present=True`가 native thermo-only diagnosis에는 unsupported composition으로 보였다는 점이다. 직렬화 전에 이 flag를 끈 중간 수정은 nonzero phase labels의 writer contract를 위반해 실패했다.
3. **최종 표현 수정.** Phase-only sidecar가 있는 원 data를 먼저 직렬화하고, native 분석 호출에만 기존 `_thermo_only_view`를 반환한다. Q1–Q7, source state, profile, phase label, tolerance는 그대로다. 다섯 번째 실제 export/diagnosis regression까지 통과한 뒤 저장된 native states를 재사용해 matching을 다시 실행했다.
4. **문서 직렬화 수정.** 첫 Markdown write에서 JavaScript string이 TeX escape를 소비한 것을 확인했다. Raw literal로 다시 쓰고 C0 control byte와 math delimiter 보존을 확인했다. 수치 파일은 변경하지 않았다.

Original/interim/final API snapshot과 실패 로그를 모두 보존했다. API AST 비교에서 `run_state`, `interpolate_initial`, `build_driver`가 불변이고 `export_compose`만 바뀌었음을 `PROJECTION_REPAIR.json`에 기록했다. 원 `tables/`의 composition rejection도 남아 있으며, 유효한 최종 출력은 `tables_v2/`다. 실패한 writer 중간 수정 때문에 실제 state를 다시 계산하지 않았다. [R6–R8]

| 실제 실행 receipt | Child exit | Wall time | 결과 |
|---|---:|---:|---|
| Primary production | 0 | 14.02 s | 54 grid + 2 rod 요청, 원 표현 rejection 보존 |
| Final repaired unit regressions | 0 | 7.20 s | 5 passed; pytest 6.62 s |
| Cached-state reproject | 0 | 2.40 s | N16/32/64의 두 pair 실제 admission |
| DDME2 N128 | 0 | 37.32 s | 동일 9점 실제 재생성 및 matching |
| Density derivative | 0 | 12.09 s | 8개 reoptimized state 및 두 step 비교 |
| Selected-F anchors | 0 | 5.19 s | 4개 anchor, 2회 single warm restart |

Wall time은 process wrapper의 실제 monotonic 시간이다. 초기/중간 실패의 exit code와 상세 메시지는 별도 JSON/log에 있고, 이후 성공 receipt로 덮어쓰지 않았다.

## 8. 증거·재현 경로

아래 경로는 `data/crust_refinement_20261117/`에 상대적이다. Test source 경로만 repository root에 상대적이다.

- **R1**: `tf_regeneration/SUMMARY.json`, `NATIVE_STATE_MANIFEST.json`: 102 native 호출의 exact state/input/executable hash, acceptance와 final aggregate.
- **R2**: `tf_regeneration/tables_v2/DDME2/N*/SUMMARY.json`, 각 `RESULT.json`, `CONSISTENCY.json`, `NATIVE_INTERVAL_DIAGNOSTIC.json`.
- **R3**: 같은 구조의 `tables_v2/TM1_CANONICAL/N*/`.
- **R4**: `tf_regeneration/SEED_IDENTITY.json`, `selected_anchor_control/SEED.json`, `selected_anchor_control/SUMMARY.json`.
- **R5**: `tf_regeneration/DDME2_THERMODYNAMIC_CHECK.json`, `TM1_CANONICAL_THERMODYNAMIC_CHECK.json` 및 참조하는 실제 `states/*/STATE.json`.
- **R6**: `tf_regeneration/REGISTRY_ADDITION.json`, `BUILD_PROVENANCE.json`, `PROJECTION_REPAIR.json`, 세 API snapshot 및 production runner snapshot.
- **R7**: `python_frontend/tests/test_tf_candidate.py`와 `tf_regeneration/TEST_RED*`, `TEST_GREEN_INITIAL*`, `TEST_GREEN*`, `TEST_EXPORT_RED*`, `TEST_FINAL*`, `TEST_FINAL_REPAIRED*`.
- **R8**: `tf_regeneration/INPUT.json`, `PRODUCTION_RUN.json`, `REPROJECT.json`, `REFINE128.json`, `DERIVATIVE.json`, `SELECTED_ANCHOR.json` 및 대응 stdout/stderr.

Native TF executable SHA256는 `21f00ecd9cd492e37329e5690d342842c19025d4b649f97f99e84e592522c189`, 재사용한 D native Maxwell executable은 `3893fbe061c2564ed4138211f7ea379a973851eb6e785fad22edcec0875152b1`다. Final Python API SHA256는 `61ad299bfeab5ee5588758b6367f7b69e3911317257dcecea23665f08a3c4c18`이다. 전체 source/driver/compiler 명령은 build receipt와 summary에 있다.

실행 환경은 `/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python`, OMP/OpenBLAS/MKL thread 수 1, G의 `python_frontend`와 `pipeline`, 기존 audit dependency path다. Runner의 `production`, `refine128 --models DDME2`, `derivative`, `reproject` phase가 역할을 분리한다. 기존 결과가 있으면 source/executable identity를 확인하는 state checkpoint를 사용한다. Driver build도 다른 source identity를 기존 build directory에 덮어쓰지 않는다.

이 subtask는 commit/push하지 않았다. G 전체의 독립 review, stellar-admission 결과와 publication decision은 상위 작업의 최종 기록을 따른다.
