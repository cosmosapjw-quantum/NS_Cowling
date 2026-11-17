# H: NL3wrL55·TM1e 실제 N128 공간 격자 및 저밀도 연장

## 결과

G commit `905b845ad0f14d2f1e7f67f9eb802691144850cd`의 선택된 N64 상태를 초기 guess로 사용하여 두 모델을 실제 N128 native TF solver에서 다시 최적화했다. 원래 70개 밀도점에 대한 126회 native 호출 중 67개 source state가 optimizer/physical gate를 통과했다. TM1e의 전이 support를 막는 n=0.008 fm^-3 한 점에만 추가한 opposite-neighbor control은 거절됐다. 총 native 호출은 127회이며, 원실패·warm restart·bounded control을 모두 보존했다.

NL3wrL55→PK1_CANONICAL은 n=0.0003–1.0 fm^-3의 sharp Maxwell 출력이 수치 readback, metadata readback, serialized independent admission, native verify를 통과했다. TM1e는 저밀도 N128 상태 자체는 n=0.0003 fm^-3까지 얻었지만, 전이 근처 optimizer failure가 source support를 끊어 N128 Maxwell 출력은 미성립이다. 이 문서는 whole-EOS 표면 연결 또는 TOV/Cowling 승인을 주장하지 않는다.

## 고정한 물리·수치 정의

- 밀도 n: fm^-3; 에너지 밀도 ε와 압력 P: MeV fm^-3; μ_B, μ_e 및 핵자 질량 기준: MeV; Wigner–Seitz 반경: fm. 이 계산은 기존 cold beta-equilibrated native TF functional, Wigner–Seitz cell 전체 전하 중성 조건, 그 ħc 상수를 그대로 사용한다. GR signature는 프로젝트 기본 (-,+,+,+)이고 이 국소 TF 계산에는 직접 등장하지 않는다.
- NL3wrL55는 G의 lower_adjacent_branch, TM1e는 G의 upper_F_component를 선택했다. N64 profile을 보간한 값은 N128 solver의 초기 guess일 뿐 EOS source row로 내보내지 않았다.
- spherical droplet, 3차원 Wigner–Seitz cell을 유지했다. N128 결과의 geometry는 전부 droplet이다. 다른 geometry와의 전역 minimum 또는 continuum limit은 검증 범위에 포함되지 않는다.
- 실제 최적화: maxiter=1800, ftol=1e-14. 원래 stationary/optimizer success, neutrality, baryon conservation, KKT, radius stationarity, radius/density bound gate를 유지했다.
- 첫 실패가 G에서 정한 유한·작은 잔차의 optimizer mode8/9 조건을 충족할 때만 동일 tolerance의 warm restart를 한 번 수행했다. 최종 실패점은 accepted source의 연결을 끊는다. 실패점을 가로지르는 interpolation, N64–N128 splice, ε/n 또는 μ_B의 임의 이동을 하지 않았다.
- 실제 source에서 계산한 μ_B=(ε+P)/n의 단위와 물리 에너지를 보존하는 CompOSE Q 재표현을 사용했다. output physical_energy_shift_mev는 0이다. abundance 파일은 제공하지 않으며 phase label만 포함한다. native source neutrality와 export abundance availability를 구분한다.

## 밀도 계획과 source identity

계산 전에 `data/low_density_crust_20261117/support_refinement/PLAN.json`에 70점 계획, solver/API/code SHA, 선택 G source state SHA를 고정했다. NL3wrL55는 0.010–0.015의 기존 전이 창을 0.00025 간격으로 유지하고 0.01201/0.01202/0.01204/0.01208을 포함했다. 여기에 0.0003/0.0005/0.00075/0.001/0.0015/0.002/0.003/0.004/0.006/0.008/0.009를 추가하여 36점을 계산했다. TM1e는 0.006–0.012를 0.00025 간격으로 유지하고 0.0003/0.0005/0.00075/0.001/0.0015/0.002/0.003/0.004/0.005를 추가하여 34점을 계산했다.

| 고정 대상 | SHA256 |
|---|---|
| H PLAN.json | `8a1d4084dc0a6d9d155ff929dcfe2bbcba8b0b5fd837cfd35565478282549707` |
| H main runner | `f65410a92b9cd7960e46dabe737cc08cf4c1412044d63f91f0db5705e84f11df` |
| frozen tf_candidate.py | `61ad299bfeab5ee5588758b6367f7b69e3911317257dcecea23665f08a3c4c18` |
| TF executable | `21f00ecd9cd492e37329e5690d342842c19025d4b649f97f99e84e592522c189` |
| Maxwell executable | `533d55ba2b5264e3e5248395baefc9463b4732b818e3e8af14beb774779f5356` |
| G NL3wrL55 N64 n=.012 source | `541dc8d52d8bf3c477f589d344a28b1e087c1afe30ec63b6ebbc7f73c85061b8` |
| G TM1e N64 n=.008 source | `9adb00fae0428769978944b9aff8ddbedc602a36d5643bbf7cf4ba25ce626838` |

첫 prepare는 seeds directory 미생성으로 actual exit1이었고 계산 시작 전에 directory 생성만 추가했다. 실패 stdout/stderr/receipt와 이전 PLAN을 보존했다. 수정 후 prepare는 actual exit0이며 production 동안 위 source bytes는 바뀌지 않았다. shared solver/API/model 파일은 수정하지 않았다.

## N128 실제 상태 및 Maxwell 결과

| 항목 | NL3wrL55 | TM1e |
|---|---:|---:|
| 요청 source state | 36 | 34 |
| accepted source state | 36 | 31 |
| 최저 accepted n [fm^-3] | 0.0003 | 0.0003 |
| 전이 쪽으로 연결된 저밀도 component 상한 [fm^-3] | 0.015 | 0.00775 |
| accepted states 최대 KKT [MeV] | 6.25483e-4 | 2.03909e-4 |
| 최대 radius stationarity [MeV] | 4.29238e-7 | 1.21743e-6 |
| 최대 baryon error [fm^-3] | 1.06165e-15 | 8.74301e-16 |
| 최대 charge error [fm^-3] | 2.72745e-19 | 2.82909e-19 |
| N128 strict Maxwell | 검증된 출력 1개 | source gap 때문에 미성립 |

NL3wrL55→PK1의 N128 전이는 μ_t=943.2775231078147 MeV, P_t=0.017147777503688004 MeV fm^-3, n_-=0.011619307065360653 fm^-3, n_+=0.019435813564131422 fm^-3이다. ε_-=10.94308341133884, ε_+=18.316218300855468 MeV fm^-3이며 Δε=7.373134889516628 MeV fm^-3이다. Δε=μ_t Δn의 상대 잔차는 1.2046143641871617e-16이다. 이 잔차는 endpoint construction identity를 검사하며 독립적인 functional accuracy 증거로 사용하지 않는다.

| NL3wrL55 selected branch | N32 (G) | N64 (G) | N128 (H) |
|---|---:|---:|---:|
| P_t [MeV fm^-3] | 0.017222521540651185 | 0.01716204074067177 | 0.017147777503688004 |
| n_- [fm^-3] | 0.011536203530662753 | 0.011608021369174801 | 0.011619307065360653 |
| Δε [MeV fm^-3] | 7.498674069657492 | 7.392795049231822 | 7.373134889516628 |

N64→N128의 변화는 P_t −0.0831092%, n_- +0.0972233%, n_+ −0.0491150%, Δε −0.265937%다. G N32→N64보다 P_t 변화가 줄었다. 이것은 실제 mesh refinement의 수치 결과이며 임의 오차 목표의 달성이나 continuum 한계 인증은 아니다. 이번 finite-domain 판정에서 N256을 필요로 하는 별도 위험은 제기하지 않았으므로 N256 계산은 하지 않았다.

## TM1e의 구체적 blocker와 한 번의 진단적 repair

n=0.008, 0.0085, 0.0095 fm^-3가 primary 및 한 번의 warm restart 뒤에도 optimizer mode8로 남았다. native process exit는 모두 0이지만 optimizer_success=0이므로 accepted_stationary=False이다. 작게 남은 physical residual을 근거로 gate를 완화하지 않았다. source connected component는 0.0003–0.00775, 0.00875–0.00925, 0.00975–0.012로 갈라진다. 각 component의 native Maxwell 시도는 NO_COEXISTENCE_CROSSING이었다. 고립된 n=0.00825 점은 최소 component 크기를 충족하지 못한다.

G에서 얻은 TM1e 전이가 n_-=0.00817549 부근에 있어 n=0.008의 실패가 전이 bracket을 실제로 막는다. 따라서 선택된 동일 branch의 accepted N128 n=0.00825 profile에서 n=0.008로 역방향 초기화를 한 번만 추가했다. 이 control은 G N64 profile을 N128로 옮긴 첫 seed와 구분되며, solver/physical gate 및 maxiter/ftol은 동일하다. 결과는 다시 mode8이었다. 원 failed iterate와 E/A 차이는 −2.3067059373715892e-10 MeV였으나 이 사실은 optimizer 승인을 대체하지 않는다. 추가 control 이후에는 같은 문제를 반복 실행하지 않았고 .0085/.0095는 전이 밖이어서 별도로 재시도하지 않았다.

따라서 TM1e에 관해서 확인된 것은 accepted 저밀도 N128 native 상태 및 구체적 optimizer/source-support blocker다. N128 matching 성공이나 물리적 비존재를 뜻하지 않는다. N64 결과와 N128 저밀도 상태를 합친 출력도 만들지 않았다.

## 실행·출력 검증 및 경로

원 production wrapper actual exit0, 103.263247 s; bounded control wrapper actual exit0, 1.520536 s. bounded control의 물리 상태는 거절됐다는 별도 결과를 유지한다. 원 production은 native Maxwell을 4회 호출했고 NL3wrL55 1개만 출력했다. 그 출력 773행의 CompOSE numerical 및 metadata readback 차이는 모두 0, serialized independent admission은 ADMITTED, native verify는 PASS다. hard gate C1–C6/C10/C11을 통과했고 C7–C9는 기존 finite-difference diagnostics로 유지했다.

기준 root는 `/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117`이다.

- 요약: `data/low_density_crust_20261117/support_refinement/SUMMARY.json`
- 실제 accepted source별 값·경로·SHA: `data/low_density_crust_20261117/support_refinement/ACCEPTED_SOURCE_STATES.json`
- 실제 all-state 시도·실패: `data/low_density_crust_20261117/support_refinement/REGENERATION.json` 및 `states/`
- CompOSE source component: `data/low_density_crust_20261117/support_refinement/PREPARED_MANIFEST.json` 및 `prepared/`
- strict match 판정: `data/low_density_crust_20261117/support_refinement/MATCHING.json`
- NL3wrL55 N128 출력: `data/low_density_crust_20261117/support_refinement/matching/GPPVA_NL3wrL55/continued_G_lower_adjacent_branch/N128/run0_component_000__component_000/compose/`
- bounded control: `data/low_density_crust_20261117/support_refinement/bounded_anchor_control/`

NL3wrL55 출력 eos.thermo SHA256는 `b96bbf5918163e01c1a23444cc0e94aa4844161d90f0138aed43291c54029ce9`, eos.nb SHA256는 `6d2b041259eeaf346461c3c2c3e647e0e44eb88cd71bb230f913efe53ec89f8e`다. 전체 파일별 SHA는 SUMMARY.json과 MATCHING.json에 보존했다.

## 재현

아래 runner는 동일한 worktree에서 기존 state의 exact executable/model/mesh/density identity를 확인해 cache를 재사용한다. 최초 실행은 위 PLAN/source identity를 고정한 뒤 수행됐다.

```bash
cd /home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python \
 scripts/refine_supported_crust_20261117.py all
```

bounded control의 실행 소스와 PLAN은 해당 하위 폴더에 있다. 그 script는 이미 존재하는 native work_dir를 덮어쓰지 않으므로 별도 지시 없이 재실행하지 않는다. 검증 결과 요약은 `support_refinement/summarize.py`로 실제 산출물에서 다시 작성할 수 있다. 이 작업은 commit/push를 수행하지 않았다.

## 후속 whole-EOS 검사와 범위

상위 H 작업의 `data/low_density_crust_20261117/stellar/h_initial/FULL_CONTRACT_ADMISSION.json`은 위 NL3wrL55 N128 output 파일들의 동일한 SHA를 입력으로 사용했다. 고정된 E outer/surface reference와 direct Hermite closure에서는 n=0.0003–0.0115 fm^-3의 모든 기존 source knot에서 lower join이 거절됐다. 따라서 이 문서의 finite-domain sharp Maxwell PASS를 surface-complete stellar EOS 승인으로 읽으면 안 된다. 이 판정은 고정 reference/closure의 결과이며 일반적인 causal connection의 비존재를 뜻하지 않는다. 해당 receipt의 이 문서 작성 시 SHA256는 `723332f6c5d7d4dcde413ca33229f98e3506457393c4883eb31f8d9e99d02516`다. TOV/Cowling 입장 여부와 후속 bridge 선택은 상위 H 문서에서 다룬다.
