# H 연구 결과: 저밀도 source와 전체 EOS

상태: 계획한 native source, causal whole-EOS/TOV, DDME-X Cowling 비교를 완료했다. 잔차 기준을 넘은 low-frequency 후보 네 기록은 거절 상태로 보존한다. 최종 commit과 게시 상태는 RESUME 및 외부 delivery receipt를 따른다.

## G 게시

G commit 905b845ad0f14d2f1e7f67f9eb802691144850cd를 research/crust-refinement-20261117에 push했고 git ls-remote의 정확한 SHA 일치를 확인했다. H는 research/low-density-crust-20261117에서 진행했다. H push는 이 단계에 포함하지 않는다.

## 실제 저밀도 계산과 sharp matching

N128의 다섯 crust model에서 요청한 밀도점은 172개다. 실제 새 native 호출 272회 중 STATIONARY 157회, REQUIRES_REFINEMENT 115회였고, exact G N128 상태 9개를 재사용해 최종 accepted source 상태는 166개다. Native process exit 0과 optimizer/physical state admission을 구분했다. 실패 6개는 삭제하지 않았으며 source를 그 점에서 끊었다.

| Crust → core | 요청/accepted | 남은 실패 n / fm^-3 | Maxwell에 연결된 lower edge / fm^-3 | H N128 P_t / MeV fm^-3 |
|---|---:|---|---:|---:|
| DDME2 → PK1 | 19/18 | .005 | .006 | .0140095258454951 |
| DDME-X original rho0 → GM1 | 38/37 | .002 | .003 | .0115419249023511 |
| NL3 provider precision → GM1 | 45/44 | .004 | .005 | .0482192936139158 |
| NL3wrL55 → PK1 | 36/36 | 없음 | .0003 | .0171477775036880 |
| TM1e → PK1 | 34/31 | .008, .0085, .0095 | 전이를 포함하는 N128 component 없음 | N128 미생성 |

DDME2의 전이창은 exact G N128 source reuse다. 나머지 성공한 세 모델의 P_t는 G N64 대비 각각 DDME-X -0.0480211%, NL3 -0.0135909%, NL3wrL55 -0.0831092% 변했다. 이는 finite mesh sensitivity이며 continuum 오차 상한이 아니다. 네 strict matched output은 원 source 값, Maxwell identity, 독립 source-support admission 및 serialized native readback을 통과했다.

다섯 모델 모두 .0003 fm^-3의 actual stationary state는 얻었다. 그러나 DDME2, DDME-X, NL3에서는 중간 실패점이 저밀도 상태와 전이창을 분리한다. TM1e의 .008 실패는 N128 전이 bracket을 막는다. 실패점을 건너뛴 native interpolation이나 N64/N128 혼합 EOS는 만들지 않았다. 세부 원본은 LOW_EXTENSION_KO.md와 SUPPORT_REFINEMENT_KO.md 및 대응 source manifest에 있다.

## 하부 연결식의 제약과 새 구현

기존 고정 E outer/surface와 direct Hermite lower closure로는 H에서 DDME-X의 .003/.004 연결만 수용됐다. G의 DDME-X primary lower join은 .007이므로 실제 연결 가능한 microscopic source 범위를 아래로 확장했다.

다른 세 모델의 direct cubic 거절을 일반적인 causal barotrope 불가능으로 해석할 수는 없다. y=n/mu를 쓰면 c_s^2/c^2=n/(mu dn/dmu)≤1은 y'≥0과 동치다. 서로 다른 positive y endpoint의 smooth causal connection에는

n_L/(2mu_L) (mu_R²-mu_L²) < P_R-P_L < n_R/(2mu_R) (mu_R²-mu_L²)

가 필요충분하다. 고정 endpoint의 단일 cubic에는 이보다 좁은 제약이 추가된다. DDME2 N128 n_R=.006에서 실제 평균밀도 .00205096598은 cubic 하한 .00217950741보다 작지만 일반 causal 허용구간 안에 있다.

이에 CausalExponentialInterval를 명시적 opt-in으로 구현했다. n/mu=y_L+(y_R-y_L)expm1(kappa*x)/expm1(kappa), x=(mu-mu_L)/(mu_R-mu_L)로 두고 압력 적분이 endpoint P를 보존하도록 유일한 kappa를 결정한다. Stellar observable을 fitting에 사용하지 않는다. 기존 Hermite 기본값과 본문, stellar solver는 그대로다. 원 source energy shift는 0이다.

이 family는 명시한 현상론적 lower closure다. 해당 구간의 핵종·neutron drip·에너지 최소를 미시적으로 계산한 것은 아니다. 실패 native knot가 closure 안에 있으면 그 위치를 contract에 기록했고 TF 수렴 성공으로 세지 않았다. 유도·안정적 수치 primitive·독립 quadrature/CAS 검산은 LOWER_THERMODYNAMICS_KO.md와 LOWER_EXPONENTIAL_IMPLEMENTATION_KO.md에 있다.

## 조건부 1.4 M_sun TOV 결과

최초 및 다음 causal admissible source knot를 관측량 전에 고정하여 8 contracts를 계산했다. Primary 4개에는 같은 중심밀도의 tighter TOV control도 적용했다. 이 inventory의 actual TOV solve는 154회이며, 모든 target-mass background가 사전 기준을 통과했다.

| Crust → core | primary/next lower join n / fm^-3 | primary R / km | sharp interface r_t / km |
|---|---:|---:|---:|
| DDME2 → PK1 | .006/.007 | 14.420417 | 13.477555 |
| DDME-X → GM1 | .003/.004 | 13.875174 | 13.061247 |
| NL3 → GM1 | .005/.006 | 13.875198 | 12.992276 |
| NL3wrL55 → PK1 | .001/.0015 | 14.420422 | 13.466460 |

M=1.4 M_sun에서 최대 target-mass 차이는 4.43e-8 M_sun이고, 실제 two-sided dense solution의 r/m jump는 0이다. 유한 에너지밀도 epsilon과 수밀도 n의 jump를 보존했다. 여기서는 첫 sampled ascending mass bracket을 선택했으며 radial stability나 maximum mass를 결론내리지 않는다.

Next-join의 반지름 차이는 1.79e-8–9.48e-7 km다. Tighter setting의 반지름 차이는 최대 3.83e-8 km였고, 이는 엄밀한 오차 상한이 아니다. 특히 NL3 변화는 수치 sensitivity와 비슷하여 분해된 물리 효과로 해석할 수 없다. 작은 R 변화가 compressibility나 Cowling spectrum의 동등성을 뜻하지 않는다. 예를 들어 NL3wrL55 두 lower closure의 max(c_s²/c²)는 .114647과 .00414096으로 다르다. 자세한 결과는 CAUSAL_STELLAR_KO.md와 TOV_RESULTS.json을 따른다.

## DDME-X Cowling 비교

별도의 고정된 계획에서 기존 Hermite 연결만 사용해 N64/.007, N128/.007, N128/.003, N128/.004를 비교한다. 실제 frozen G stellar code를 사용하므로 spatial/source refinement와 lower-join 변경의 효과를 분리할 수 있다. l=2, equilibrium bulk Gamma1, perfect fluid Cowling이며 slow/rapid conversion은 명시적 이상화한 interface 경계조건이다.

50–3500 Hz의 고정 scan에서 처음 검출된 slow sign-change root 두 개와 rapid root 하나를 선택했다. 모든 spectrum의 완전성을 주장하지 않는다. Root-finder convergence만으로 수용하지 않으며 determinant/join<1e-6, surface/junction<1e-10을 검사한다. 거절된 root는 원본을 보존하고 한 번의 tighter control만 허용했다.

선택한 ordinary 12개 root record와 numerical refinement 3개를 모두 계산했다. 최종 잔차 gate는 11 PASS, 4 FAIL이다. 기존 raw runner는 하나라도 거절되면 exit1을 반환하므로 actual exit1, raw status REJECTED_OR_EMPTY를 그대로 보존했다. 계산이 중단되거나 결과가 비어 있었던 것은 아니다. 별도 MODE_ACCEPTANCE.json의 지위는 COMPLETED_WITH_REJECTED_CANDIDATES다.

| Case | R / km | Slow 약1.905 kHz branch / Hz | Rapid branch / Hz | 낮은 slow candidate |
|---|---:|---:|---:|---|
| H_DDMEX_N128_primary | 13.875173909 | 1905.456530662 | 1905.438805040 | 211.115932 Hz, FAIL |
| G_DDMEX_N64_primary | 13.875173405 | 1905.457980401 | 1905.440175219 | 211.209128 Hz, FAIL |
| H_DDMEX_N128_G_join | 13.875173392 | 1905.458132454 | 1905.440453666 | 211.103771 Hz, FAIL |
| H_DDMEX_N128_next_join | 13.875173849 | 1905.456726576 | 1905.439021280 | 211.114515 Hz, PASS |

같은 .007 join에서 N64→128의 주파수 변화는 slow 0.00015205293 Hz, rapid 0.000278446402 Hz다. N128의 .003→.004 next-join 변화는 각각 0.000195914876, 0.000216239801 Hz다. 반면 primary의 tighter numerical control은 slow 0.000196191269, rapid 0.000166469104 Hz 변했다. 따라서 mesh-only와 next-join의 작은 차이는 수치 sensitivity와 같은 규모이고, 이를 정밀하게 분리된 물리적 shift로 주장하지 않는다. .007→.003의 더 넓은 join 선택 변화는 두 branch에서 약0.0016 Hz의 유한 차이를 보였다. 이 역시 엄밀한 오차 상한이나 일반적인 closure 독립성의 증명은 아니다.

Primary .003의 약211.116 Hz slow candidate는 tighter 계산 후에도 determinant/join=3.21422e-6로 거절됐다. G N64 .007, H N128 .007 및 primary numerical refinement의 낮은 후보도 거절됐다. H N128 .004의 약211.1145 Hz 후보 하나는 잔차 gate를 통과했다. 이 단독 통과를 다른 미승인 후보와의 검증된 모드 차이나 전체 저주파 spectrum 확정으로 확대하지 않는다. 근 찾기 convergence와 실제 경계/접합 residual을 분리했고, 물리적 모드 부재 주장은 하지 않는다.

Fresh review는 numerical refinement가 거절될 때 ordinary surface_h로 되돌아갈 수 있는 runner 설정 오류도 확인했다. 완료된 원 runner의 byte snapshot과 원 raw output을 보존한 뒤, 재사용 runner가 요청된 surface_h를 유지하고 effective settings를 저장하도록 수정했다. 실제 발생한 failed low-root fallback 1개는 planned paired refinement에서 제외했다. 승인된 high-frequency refinement 두 개는 원래 계획된 surface_h=3e-8에서 계산되었다. 세 surface 설정의 별도 configuration regression은 모두 통과했으며 과거 결과를 덮어쓰지 않았다.

## 구현 검증

새 connector 26 case와 기존 관련 회귀를 포함한 targeted 71 tests가 실제 통과했다. 이후 전체 frontend/pipeline suite는 1453 passed, 9 skipped, 4 subtests passed, actual exit 0이었다. 해당 실행의 frontend/pipeline/native 213개 source bytes는 전후 동일했고 git diff --check도 통과했다. Native 실행 파일은 exact G build를 재사용하며 H에서 변경한 native source는 없다.

독립 H 검토의 범위와 최종 판정은 INDEPENDENT_REVIEW_KO.md 및 review/REVIEW_RECEIPT.json에 기록한다. 한 번의 fresh review 안에서 source identity, 새 causal core, 실제 TOV/Cowling 영수증 및 원 실패 보존을 확인한다. 개별 수치 실패는 전체 suite 통과로 지우지 않는다.

## 다음 I의 우선순위

1. DDME-X low-frequency candidate의 접합 잔차가 기준 아래로 내려오지 않는 원인을 finite controlled test로 분리한다. 알려진 EOS knot/음속 불연속에서의 적분 분할과 match 위치·수치 조건을 검토하되 잔차 기준은 유지한다. 승인된 해를 확보한 뒤 새 causal lower closure의 mode sensitivity를 비교한다.
2. Native mode8의 여섯 밀도점을 원 functional·parameter·acceptance를 고정한 채 변수/constraint conditioning과 gradient/line-search 진단으로 다룬다. 미시적 저밀도 support를 복구하면 명시적 phenomenological lower 구간을 줄일 수 있다.
3. 그 결과로 whole-EOS 후보를 다시 고르고, 필요한 geometry/mesh 비교를 특정 위험에 한정해 수행한다. 이미 얻은 finite TOV 해를 global equilibrium이나 continuum convergence로 승격하지 않는다.
