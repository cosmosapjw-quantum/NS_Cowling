# H: native 저밀도 crust 연장과 N128 전이창 계산

## 결과와 범위

G commit `905b845ad0f14d2f1e7f67f9eb802691144850cd`의 물리 모형·에너지 기준·native acceptance gate를 유지하면서 DDME2, DDMEX_ORIGINAL_TABLE_RHO0, NL3_PROVIDER_PRECISION의 구형 droplet branch를 저밀도로 실제 연장했다. DDME2는 G의 N128 상태 9개를 byte identity로 재사용했다. 나머지 두 모델은 G 전이창 전체를 N128에서 다시 계산했다. 세 모델 모두 `n_b=0.0003 fm^-3`의 stationary 상태에 도달했다.

요청한 102개 밀도점 가운데 최종 accepted 상태는 99개다. 새 native 호출은 145회이며, 그중 STATIONARY 90회, REQUIRES_REFINEMENT 55회다. 재사용한 G 상태 9개는 새 호출에 포함되지 않는다. 모든 native process의 실제 종료 코드는 0이다. 종료 코드 0인 REQUIRES_REFINEMENT는 물리 상태의 승인으로 해석하지 않았다.

| 모델 | 요청 / 최종 accepted | G N128 재사용 | 새 native 호출: production + repair | 남은 실패 밀도 / fm^-3 | strict Maxwell에 연결된 source 구간 / fm^-3 |
|---|---:|---:|---:|---:|---|
| DDME2 | 19 / 18 | 9 | 16 + 1 | 0.005 | 0.006–0.014 |
| DDMEX_ORIGINAL_TABLE_RHO0 | 38 / 37 | 0 | 68 + 2 | 0.002 | 0.003–0.012 |
| NL3_PROVIDER_PRECISION | 45 / 44 | 0 | 57 + 1 | 0.004 | 0.005–0.038 |

원래 낮은 쪽 source block도 보존했다. 실패 밀도를 가로질러 두 block을 합치거나 다른 spatial mesh의 상태를 끼워 넣지 않았다. 실제 outer-crust closure와 stellar admission은 별도 H full-EOS 결과에서 다룬다. 이 문서는 finite native source와 strict Maxwell 결과까지 기록한다.

## 정의와 계산 규약

- `n_b`: cell-averaged baryon number density, fm^-3.
- `epsilon`, `pressure`: MeV fm^-3. `mu_B`, `mu_e`: MeV. Cell radius: fm.
- 기존 frozen RMF/Thomas–Fermi functional과 구형 Wigner–Seitz droplet convention을 그대로 사용했다. N128은 radial cell discretization이다. 정적 matter 계산이므로 metric을 새로 선택하거나 상대론적 stellar 방정식을 수정하지 않았다; 프로젝트 signature는 (-,+,+,+)다.
- 저온 thermodynamic representation은 `mu_B = (epsilon + P)/n_b`, `dP/dmu_B = n_b`를 사용한다. Native interval 진단은 source knot 사이의 열역학적 admissibility를 검사한다. 실패 knot 사이의 미계산 구간을 source로 복원하지 않는다.
- CompOSE Q encoding의 기준질량은 DDME2와 NL3에 mn=mp=939 MeV, DDME-X에 mn=mp=938.5 MeV이다. 이는 기존 모형의 frozen convention이다. Matcher의 physical-energy-preserving Q reencoding과 `physical_energy_shift_mev=0`을 확인했다.
- 물리적 c 또는 hbar를 새로 1로 설정하거나 복원하지 않았다. Native code의 frozen MeV–fm conversion을 재사용했다.
- exported source는 thermo-only table이며 phase label을 포함한다. Nuclear abundance를 복원하거나 composition validation을 새로 주장하지 않는다.

## 동결한 입력과 실제 실행

첫 expensive native 계산 전에 `data/low_density_crust_20261117/low_extension/PLAN.json`을 저장했다.

- PLAN SHA256: `02e93d047518c32ce12a577bc86b3fea4c364ac30423bc87602a6221f818be4a`.
- TF driver SHA256: `21f00ecd9cd492e37329e5690d342842c19025d4b649f97f99e84e592522c189`.
- General native matcher SHA256: `533d55ba2b5264e3e5248395baefc9463b4732b818e3e8af14beb774779f5356`.
- 상세 model registry, driver source, Python API, G state index, core manifest 및 모든 G seed state의 SHA256는 PLAN에 있다.
- 실행 entry point: `scripts/extend_low_density_crust_20261117.py`.
- N128, SLSQP ftol=1e-14, maxiter=1800, 기존 native stationarity / optimizer-success / radius / field / conservation gate를 유지했다.
- 모든 OMP / OpenBLAS / MKL thread 수는 1로 설정했다.

Production은 2026-11-17 14:26:08–14:29:12 UTC에 실제 실행되었고 process receipt의 wall time은 184.126530694 s이다. Repair는 이 실행의 독립된 완료 source에서 시작했으므로 일부 시간이 겹친다. 실제 PID와 exit code는 다음 receipt에 있다.

`low_extension/PROCESS_DDME2_DDMEX_ORIGINAL_TABLE_RHO0_NL3_PROVIDER_PRECISION.json`

`low_extension/repaired/{MODEL}/PROCESS.json`

Production 원본 index / POINT receipts / input.txt / stdout.txt / stderr.txt / STATE.json은 변경하지 않았다. PLAN에서 동결한 source bytes가 production 동안 바뀌지 않았다는 기록과 최종 hash 재확인을 모두 남겼다.

## 저밀도 연장과 제한된 수치 repair

하강 밀도 계획은 0.005, 0.004, 0.003, 0.002, 0.0015, 0.001, 0.00075, 0.0005, 0.0004, 0.0003 fm^-3이다. NL3에는 G lower edge 0.026부터 이 구간까지 연결하는 0.024, 0.022, 0.020, 0.018, 0.016, 0.014, 0.012, 0.010, 0.008, 0.006도 실제 계산했다. 모든 정확한 density values와 ordering은 PLAN에 있다.

첫 실패가 optimizer mode 8/9이고 finite returned profile, KKT < 0.02, radius residual < 2e-4, inactive bounds 조건을 만족하는 경우에만, 반환 profile로 optimizer-state restart를 한 번 허용했다. 이는 gate나 tolerance의 완화가 아니다.

그 뒤 connection gap에 남은 실패점에는 한 번의 서로 다른 productive repair를 시행했다. 최초 진행 방향의 반대편에서 이미 accepted된 이웃 profile을 초기값으로 사용하고, 같은 N128 / ftol / native gate로 다시 풀었다. 반복적인 warm restart는 추가하지 않았다.

| target | 반대편 accepted 시작점 | 추가 실제 native solve 결과 |
|---|---:|---|
| DDME2 0.005 | 0.004 | REQUIRES_REFINEMENT, mode 8, KKT 2.80540098e-5 |
| DDME-X 0.002 | 0.0015 | REQUIRES_REFINEMENT, mode 8, KKT 1.40793708e-5 |
| DDME-X 0.0004 | 0.0003 | STATIONARY, mode 0, KKT 1.55649086e-5 |
| NL3 0.004 | 0.003 | REQUIRES_REFINEMENT, mode 8, KKT 1.09774229e-4 |

위 실패는 현재 optimizer/acceptance의 수치적 미해결 상태다. 이 유한 시도만으로 물리적 평형해가 없거나 해당 모델이 불가능하다고 결론내리지 않는다. DDME-X의 G N64 실패점 0.00675와 0.01025는 이번 N128 production에서 모두 stationary로 계산되었다. N256의 추가 계산은 시행하지 않았다.

Repair plan / source identity / 실패 근거와 실행 코드는 `low_extension/repaired/{MODEL}/PLAN.json`, `REPAIRS.json`, `low_extension/repair_failed_knots.py`에 있다.

## Strict Maxwell 결과

독립된 각 connected accepted source block을 export한 후 native thermodynamic interval diagnosis를 시행했다. 세 모델의 모든 exported connected block에서 interval rejection은 0개다. 각 최종 upper block은 아래 frozen core와 strict Maxwell construction, source-support admission, CompOSE serialization/readback, metadata readback 및 serialized native verify를 통과했다.

| crust → core | mu_t / MeV | P_t / MeV fm^-3 | n_low / fm^-3 | n_high / fm^-3 | Delta epsilon / MeV fm^-3 |
|---|---:|---:|---:|---:|---:|
| DDME2 → PK1_CANONICAL fine | 943.1058672809809 | 0.014009525845495105 | 0.00977355894875663 | 0.017093104048714065 | 6.903105929597611 |
| DDME-X → GM_GM1 | 942.2975616250997 | 0.011541924902351135 | 0.008596606764460788 | 0.020104352588030134 | 10.84372082935082 |
| NL3 → GM_GM1 | 943.4548545156440 | 0.048219293613915755 | 0.03136030129260822 | 0.03951721251204436 | 7.695677487830146 |

`Delta epsilon = mu_t * Delta n`의 native relative residual은 각기 1.2866e-16, 1.6381e-16, 1.1541e-16이다. 새로운 상수 energy shift는 0이다.

N64에서 N128로 바꾸었을 때 P_t의 상대 변화 `abs(P128-P64)/abs(P64)`는 DDME-X 0.0004802105376, NL3 0.0001359085879이다. 각각 약 0.0480%, 0.0136%의 유한 mesh sensitivity이다. DDME2의 전이창은 exact G N128 reuse이므로 기존 N128 값과 동일하다. 이러한 비교는 continuum extrapolation이나 오차 상한의 증명이 아니다.

최종 세 strict-match table 경로는 다음과 같다. 각 `compose`의 부모 폴더에 `RESULT.json`이 있다.

- `data/low_density_crust_20261117/low_extension/repaired/tables_v2/DDME2/N128/block_01/match_00_00/compose`
- `data/low_density_crust_20261117/low_extension/repaired/tables_v2/DDMEX_ORIGINAL_TABLE_RHO0/N128/block_01/match_00_00/compose`
- `data/low_density_crust_20261117/low_extension/repaired/tables_v2/NL3_PROVIDER_PRECISION/N128/block_01/match_00_00/compose`

## 최종 source / run / summary manifest

`data/low_density_crust_20261117/low_extension/MANIFEST.json`

SHA256: `e98b442907450ec21b10adc77cf79ba686947ae453beee0d1703ad559df13329`.

이 파일은 145개 실제 native invocation을 중복 없이 집계하고 다음 사항을 확인한다.

1. PLAN에 고정한 source / executable / G source-state bytes와 현재 bytes의 SHA256 일치.
2. Production 및 모델별 repair process가 모두 COMPLETED이며 actual_exit_code=0.
3. 모든 native child의 actual_exit_code=0 및 원본 실패 status의 보존.
4. 최종 accepted state의 전부가 N128이며, path / SHA256 / n_b / epsilon / P / mu_B를 제공함.
5. 최종 strict-match table의 실제 bytes SHA256 및 readback / metadata / independent admission / native verification 통과.

최종 원본 state index는 모델별로 다음 위치에 있다.

`low_extension/repaired/{MODEL}/FINAL_N128_STATES.json`

이 index는 실패 source도 보존한다. Accepted source를 이용할 때에는 순서대로 실패점에서 block을 끊어야 한다. 실패점을 제거한 뒤 인접한 것처럼 연결하면 이 결과의 source-support 계약을 위반한다.

재확인 명령:

```bash
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python data/low_density_crust_20261117/low_extension/verify_and_summarize.py
```

재실행 entry points:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/extend_low_density_crust_20261117.py prepare
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/extend_low_density_crust_20261117.py run
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python data/low_density_crust_20261117/low_extension/repair_failed_knots.py run --models DDME2 DDMEX_ORIGINAL_TABLE_RHO0 NL3_PROVIDER_PRECISION
```

현재 runner는 명시된 sibling G worktree의 frozen executable와 source states를 사용한다. 다른 기계에서 replay하려면 PLAN의 경로와 SHA256에 해당하는 bytes를 먼저 준비해야 한다. 이 실행에서 source 파일이나 code identity가 달라진 것을 허용하는 자동 fallback은 없다.

## 남은 결론의 경계

이 기록은 finite 구형 droplet branch의 stationary source와 strict Maxwell matching에 대해 implementation-verified / numerically checked evidence를 제공한다. 다른 geometry보다 낮은 global ground state, 연속 mesh limit, 모든 중간 밀도에서의 실제 microscopic stationarity, outer-crust 전체 closure 및 stellar 안정성을 이 결과 하나로 주장하지 않는다.

세 실패 knot와 original provider source의 식별 문제는 보존한다. 공개 원표의 실패 결과를 이 native counterpart 성공으로 덮어쓰지 않는다. 다음 H 통합 결과는 위 strict tables를 입력으로 frozen outer endpoint와의 full-EOS admission 및 TOV/Cowling 비교를 수행한다.
