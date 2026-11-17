# H: opt-in causal lower closure와 1.4 M_sun TOV inventory

## 유한 실행 결과

동결된 H strict Maxwell source 네 개에 명시적인 causal exponential lower closure를 적용했다. 각 source의 밀도 순서에서 첫 번째와 다음 admissible 실제 source knot를 관측량 계산 전에 선택했다. 총 8개 full-EOS contract가 WholeCrustBarotrope validation을 통과했고, 모두 1.4 M_sun target-mass TOV background를 계산했다. Primary 4개에는 같은 중심밀도에서 r0/2 및 더 엄격한 tolerance를 적용한 TOV sensitivity를 추가했다.

실제 완료된 TOV solver 호출은 154회다. 48개 고정 density scan 배경, target mass root 평가, 최종 8개 배경 및 4개 tighter control이 포함되며 실패한 scan point는 없다. Admission과 TOV process는 각각 shell exit code 0으로 완료되었다. 원본 source, outer/surface reference, contract 및 실행 code의 before/after SHA256가 일치한다.

이는 선택한 phenomenological closure의 정적 TOV 계산이다. 추가 Thomas–Fermi source 생성, 실패 native knot의 복구, 새로운 Cowling spectrum 계산은 수행하지 않았다.

## Closure의 정의와 근거

Cold Euler chemical potential을 mu=(epsilon+P)/n으로 정의하고 x=(mu-mu_L)/(mu_R-mu_L), y=n/mu로 둔다. 이번 opt-in 연결식은

```text
y(mu) = y_L + (y_R-y_L) Phi_kappa(x),
Phi_kappa(x) = expm1(kappa*x)/expm1(kappa),   Phi_0(x)=x,
n(mu)=mu*y(mu),
P(mu)=P_L + integral(mu_L..mu) n(u) du,
epsilon(mu)=mu*n(mu)-P(mu).
```

이다. 증가하는 y를 사용하므로 dn/dmu=y+mu*y' > 0이고, c_s^2/c^2=n/(mu*dn/dmu)=y/(y+mu*y')는 (0,1]에 있다. Endpoint pressure integral은 kappa를 유일하게 결정한다. kappa에 독립적인 자유 fitting을 추가하거나 stellar observable을 입력하지 않았다.

사용한 구현은 `CausalExponentialInterval`이며 contract kind는 `causal_exponential_n_over_mu`다. 실제 pressure integral, endpoint preservation, positivity, causality 및 saved-kappa 재검증을 통과한 contract만 선택했다. 기존 Hermite default를 대체하지 않는 opt-in 경로다.

단위는 n: fm^-3, mu: MeV, P/epsilon: MeV fm^-3이다. Stellar solver는 m=GM/c^2와 r을 km로 적분하며 metric signature (-,+,+,+), 기존 c와 G 변환 및 nominal solar-mass convention을 유지한다. Sharp Maxwell 양쪽 source row와 유한 density/energy jump는 그대로 보존한다.

Lower closure의 전 구간은 명시적인 phenomenological barotrope다. DDME2 0.005, DDME-X 0.002, NL3 0.004 fm^-3의 미해결 native knot가 이 구간에 들어가며, contract의 `known_failed_native_knots_in_closure`에 기록했다. 이들은 microscopic TF support로 승격되지 않는다. Disconnected lower source state는 upper table에 가져오지 않았다.

## 관측량 이전에 선택한 endpoint

입력은 `stellar/h_final/INPUTS.json`의 strict matched source 4개다. 고정 E outer/surface reference는 `data/stellar_extension_20261117/eos/fine/REPRESENTATION.json`이다. 아래 선택과 contract SHA256를 FULL_CONTRACT_ADMISSION.json, TOV_PLAN.json에 저장한 뒤 TOV process를 시작했다.

| 모델 | primary n_R / fm^-3 | next n_R / fm^-3 | primary kappa | next kappa |
|---|---:|---:|---:|---:|
| DDME2 | 0.006 | 0.007 | 2.4902682205225344 | 2.4918654790434274 |
| DDME-X | 0.003 | 0.004 | 1.2752551528538176 | 1.4185597969522754 |
| NL3 | 0.005 | 0.006 | 3.3096869584088755 | 3.4819008235094207 |
| NL3wrL55 | 0.001 | 0.0015 | 7.656727519722228 | 3.801794872440535 |

## 1.4 M_sun primary 배경

질량 bracket을 찾기 위한 고정 중심밀도는 [0.08, 0.16, 0.32, 0.5, 0.8, 1.0] fm^-3이다. 각 contract에서 첫 sampled ascending target-mass bracket [0.16,0.32]를 사용했다. 이 선택은 radial stability 판정이 아니다.

| 모델 / core | 중심 n_c / fm^-3 | R / km | R-r_t / km | lower closure radial width / km | outer+surface width / km |
|---|---:|---:|---:|---:|---:|
| DDME2 / PK1_CANONICAL | 0.3038516328930928 | 14.4204166641948 | 0.942862141621395 | 0.223743437860577 | 0.663763230614588 |
| DDME-X / GM1 | 0.3192437118952588 | 13.8751739805798 | 0.813926977331622 | 0.113539931722807 | 0.606289448769589 |
| NL3 / GM1 | 0.3192417336380075 | 13.8751978011773 | 0.882922185285990 | 0.140355790197395 | 0.606291908220531 |
| NL3wrL55 / PK1_CANONICAL | 0.3038512982908763 | 14.4204215817763 | 0.953962047040502 | 0.0808378450886931 | 0.663763759578131 |

r_t는 retained sharp Maxwell interface의 radius다. R-r_t는 선택한 whole-EOS closure에서의 interface-to-surface thickness이며, bridge 내부를 새로 계산한 microscopic crust thickness로 해석하지 않는다.

8개 target-mass 배경의 최대 |M-1.4 M_sun|은 4.4290435541e-8 M_sun으로 사전에 고정한 1e-7 M_sun acceptance 이내다. 각 interface에서 실제 두 phase segment의 dense output을 양쪽에서 읽어 얻은 r_jump와 m_jump는 모두 0이다. 원래 density jump는 유한하게 유지되어 있다.

## Join과 수치 sensitivity

Primary 설정은 r0=0.001 km, rtol=1e-10, atol=1e-12, max_step_h=0.002다. Tighter control은 같은 중심밀도에서 r0=0.0005 km, rtol=2e-11, atol=2e-13, max_step_h=0.001로 계산했다.

| 모델 | next-primary Delta R / km | tighter-primary Delta R / km | next-primary lower closure 폭 변화 / km |
|---|---:|---:|---:|
| DDME2 | 9.47732628e-7 | -5.49221824e-9 | 0.0164231234594 |
| DDME-X | 1.09790230e-7 | 3.83496381e-8 | 0.0223887380378 |
| NL3 | 1.79293131e-8 | 2.09883417e-8 | 0.00789041049068 |
| NL3wrL55 | 3.02394310e-7 | 1.20018129e-8 | 0.0236908870469 |

Next-primary 비교는 각각 target mass를 다시 맞춘 서로 다른 closure의 유한 sensitivity다. Tighter-primary 비교는 같은 중심밀도에서 solver controls를 바꾼 유한 수치 sensitivity다. 어느 값도 엄밀한 오차 상한은 아니다. 특히 NL3의 next-primary radius 차이는 tighter-control radius 차이와 같은 크기이므로 분해된 물리 차이로 주장할 수 없다.

Lower closure width는 선택한 endpoint가 바뀌면서 정의 구간도 달라진다. 폭 변화 자체를 새 물리 shell의 생성으로 해석하지 않는다.

정적 radius의 작은 변화는 compressibility 또는 mode spectrum의 동등성을 뜻하지 않는다. NL3wrL55 primary bridge의 max(c_s^2/c^2)는 0.1146466493이고 next bridge에서는 0.00414096408이다. 두 값 모두 causal gate 안에 있지만 lower closure의 국소 압축률은 다르다. 이 inventory에서 새 Cowling eigenmode를 계산하지 않았으므로 mode 영향은 미해결이다.

## 실행·파일 identity

실행 전에 h_thermo_construction의 code freeze와 targeted test 71 PASS를 확인했다. Production code identity:

`python_frontend/compose_eos/whole_crust.py`

SHA256 `7bb94e914019045a636f2008d5ed4ce2b6e346c8f79e19f0100416c5e924b513`.

소유 실행 코드는 `scripts/prepare_causal_low_density_stellar_20261117.py`, `scripts/run_causal_low_density_backgrounds_20261117.py`다. 새 native TF 호출은 없다.

출력 기준 폴더:

`data/low_density_crust_20261117/stellar/causal_extension/`

| 파일 | SHA256 |
|---|---|
| ADMISSION_PLAN.json | 7eea8891121048de5f49002c00113a4d4400c663e0017d90838db61cea78c832 |
| FULL_CONTRACT_ADMISSION.json | 38551d8be9066908e14a3862f9ec584f11b100a730d72381e8f27dcee6d63974 |
| TOV_PLAN.json | 42eef4c3d4cb774ae2c88cb9410ab5deb834af79b94e6394bdd6cacd8d0ab7c4 |
| TOV_RESULTS.json | f8d2ab4398a8333f19410605b7f4d4739e0c7af5bbddfd5749ec96c650eafe1b |

모든 contract는 `contracts/{source_id}/{primary|next_admissible_join_control}/REPRESENTATION.json`에 있고, 각 실행 세부 기록은 `calculations/`에 있다. `ADMISSION_PROCESS.json`, `TOV_PROCESS.json`에 내부 receipt가 있으며 `REMOTE_PROCESS_EXIT_RECEIPTS.json`은 remote process 도구가 독립 확인한 종료 결과다. Admission shell PID 313130은 exit0/runtime21.24s, TOV shell PID314880은 exit0/runtime30.89s로 완료되었다.

최종 요약과 파일 identity는 `FROZEN_SUMMARY.json`에 있다. 기존 결과가 있으면 runner는 덮어쓰지 않는다.

## 근거 상태

- **derived / numerically checked:** endpoint-determined causal lower connector의 first-law 및 causality 조건.
- **implementation-verified:** 8 whole-EOS contract admission, exact endpoint/source retention, 154 completed TOV solver calls, 8 target-mass backgrounds, 4 tighter controls, two-sided interface continuity, source/code hash preservation.
- **unresolved:** bridge의 microscopic nuclear physics, native failed-knot stationarity, 다른 geometry와 비교한 global minimum, radial stability, 새 4모델의 mode spectrum과 closure dependence.

이 결과는 직접 cubic Hermite 연결이 거절된 경우에도 명시적으로 선택한 다른 causal phenomenological closure로 whole-EOS TOV 계산이 가능함을 해당 네 source에서 보인 것이다. 원래 Hermite rejection이나 미해결 native state를 삭제하거나 그 판정을 뒤집지 않는다.
