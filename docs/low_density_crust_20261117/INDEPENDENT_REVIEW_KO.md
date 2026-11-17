# H 독립 검토

**판정: ACCEPT_SCOPED. Blocking review findings: 0.**

G `905b845ad0f14d2f1e7f67f9eb802691144850cd`를 기준으로 H의 실제 N128 source, opt-in lower closure 구현, whole-EOS/TOV와 DDME-X Cowling 산출물을 검토했다. 이 판정은 아래에 한정한 구현·유한 수치 결과 및 실패 보존의 승인이다. Native 실패 6개와 Cowling 미승인 후보 4개가 해결됐다는 뜻은 아니다. 과학적 미해결 결과를 지우지 않고 보고한 현재 H 범위에는 게시 전 수정을 요구할 blocker가 남지 않았다.

검토자는 생산 코드, native state, EOS, stellar 계산 결과를 수정하거나 새 TF/TOV/Cowling solve를 실행하지 않았다. 원 input/stdout/STATE, CompOSE 수치, 실행 receipt, 정확한 source hash와 코드·수학을 독립적으로 읽었다. 최종 사실 집계는 [H_SUMMARY.json](../../data/low_density_crust_20261117/H_SUMMARY.json), 종합 서술은 [RESULTS_KO.md](RESULTS_KO.md)의 아래 고정 bytes에 한정한다.

## Source와 sharp Maxwell

272개 실제 native invocation의 input hash와 raw stdout/STATE를 비교했다. STATIONARY 157개와 REQUIRES_REFINEMENT 115개, 정확한 G N128 재사용 9개를 확인했다. 요청 172개 밀도점 중 최종 accepted는 166개다. Native process exit 0을 stationary admission으로 간주하지 않았다.

| 모델 | 요청/accepted | 원실패가 유지된 밀도 / fm^-3 | 전이를 포함하는 N128 source |
|---|---:|---|---|
| DDME2 | 19/18 | .005 | .006–.014 |
| DDME-X original rho0 | 38/37 | .002 | .003–.012 |
| NL3 provider precision | 45/44 | .004 | .005–.038 |
| NL3wrL55 | 36/36 | 없음 | .0003–.015 |
| TM1e | 34/31 | .008, .0085, .0095 | 전이 bracket 미성립 |

각 source block은 미리 요청한 순서에서 연속하는 accepted 상태로만 구성된다. 실패점 삭제 후 gap을 보간하거나 서로 다른 spatial mesh를 합친 출력은 없다. TM1e의 고립된 .00825 accepted state도 보존하며, 최소 component 크기를 충족하지 않아 export하지 않았음을 확인했다.

Raw CompOSE에서 `n,P,epsilon,mu`를 별도로 계산해 native 상태와 비교한 최대 상대 차이는 2.10e-16이었다. 네 strict matched output의 원 source 행, endpoint, 파일 SHA 및 Maxwell identity를 확인했다. Physical energy shift는 0이고 protected TF functional/model parameter bytes는 G와 동일하다. DDME2 전이창은 새 N128 계산으로 세지 않고 exact G reuse로 분류했다.

이 검사는 정해진 gate에서 stationary인 finite droplet branch와 explicit source support의 증거다. Global geometry minimum, continuum 수렴, 원 public EOS의 모든 물리적 convention까지 일치한다는 결론은 포함하지 않는다. 세부 독립 기록: [SOURCE_REVIEW.json](../../data/low_density_crust_20261117/review/SOURCE_REVIEW.json).

## Lower connector의 수학과 구현

고정 양끝 상태의 cold Euler chemical potential을 mu=(epsilon+P)/n으로 두면, 양의 n,mu에서 causality는 y=n/mu의 비감소 조건 y'≥0로 표현된다. 서로 다른 y endpoint에는 압력 적분의 strict causal bounds가 필요하다. Exponential family의 적분은 kappa에 대해 연속·엄격히 감소하고 두 극한 causal bounds에 접근하므로 admissible pair마다 유일한 유한 kappa를 결정한다. 이는 stellar 관측량에 맞춘 parameter fitting이 아니다.

다음 구현 사항을 직접 검토했다.

- 양·음의 큰 kappa와 kappa≈0에서의 안정적 primitive가 직접 적분식과 일치한다.
- Susceptibility와 causal margin의 extrema는 양끝과 2+kappa*mu/d=0으로 구할 수 있다.
- Sound-speed stationary 함수 F=y+mu*kappa*A/d는 F'=kappa*y/d이므로 내부 root가 최대 하나다. 구현의 양끝 및 유일한 부호변화 root 검사가 충분하다.
- Equal-y luminal domain, strict endpoint bounds, 저장 kappa와 kind/support 변조 거절을 명시적으로 처리한다.
- Positive endpoint energy와 d epsilon/dmu=mu dn/dmu>0가 연속구간 energy positivity를 보장한다.

G와의 AST 비교에서 추가된 public class는 `CausalExponentialInterval` 하나이고, 기존 symbol 변경은 `WholeCrustBarotrope.from_contract`의 opt-in dispatch와 pressure-inverse type annotation뿐이었다. 기존 Hermite 본문과 stellar consumer files는 변하지 않았다. Targeted 71 tests는 실제 exit 0이며 독립 quadrature/first-law, 극한 kappa, endpoint identity, contract/inverse, sharp TOV의 독립 pressure-coordinate 비교를 포함한다. [CORE_REVIEW.json](../../data/low_density_crust_20261117/review/CORE_REVIEW.json)에 exact source/test hash와 판정을 남겼다.

새 lower interval은 현상론적 barotrope다. DDME2 .005, DDME-X .002, NL3 .004를 포함하는 closure는 native 실패를 복구한 미시적 source가 아니다. 이 의미가 모든 해당 contract에 명시되어 있고, 연결되지 않은 낮은 밀도 상태는 upper source에 합쳐지지 않았다. 기존 cubic failure를 일반 causal connection의 불가능성으로 확장하지 않은 판단은 정당하다.

## 조건부 TOV inventory

8개의 first/next admissible 계약을 관측량 전에 고정했다. 고정 E surface/outer를 그대로 유지했고, upper n/P/epsilon과 원 sharp interface를 보존했다. 실제 matched source와 contract upper array의 최대 상대 차이는 3.21e-16이었다.

8개의 1.4 M_sun background 및 4개의 primary tighter control, 총 실제 TOV 호출 154회를 확인했다. 최대 target-mass residual은 4.43e-8 M_sun으로 사전 기준 이내다. 실제 양쪽 dense output의 r/m jump는 모두 0이다. Central-density 선택은 첫 sampled ascending bracket이며, 이것만으로 radial stability 또는 maximum mass를 주장하지 않는다.

Next-join Delta R은 1.79e-8–9.48e-7 km, tighter Delta R 최대값은 3.83e-8 km였다. 특히 NL3의 작은 차이는 수치·mass-root sensitivity에서 분리됐다고 볼 수 없다. 이 작은 반지름 변화가 compressibility 또는 mode spectrum의 동등성을 뜻하지 않는다. NL3wrL55의 두 lower closure는 max(c_s²/c²)가 약 .11465와 .00414로 다르다. Causal inventory의 검증 범위는 TOV까지이고 새 causal 모델의 Cowling spectrum은 계산되지 않았다. [CAUSAL_TOV_REVIEW.json](../../data/low_density_crust_20261117/review/CAUSAL_TOV_REVIEW.json)을 참조한다.

## DDME-X Cowling과 남은 수치 후보

이 lane은 기존 Hermite lower closure와 frozen G stellar physics를 사용했다. 4개 ordinary case의 12 selected record와 primary numerical refinement 3 record가 모두 계산됐다. Raw 잔차를 독립적으로 다시 판정한 결과는 **11 PASS, 4 FAIL**이며, 각 profile에서 실제 interface의 high/low 두 행을 확인했다. 총 profile point 18,015개의 유한성도 확인했다.

실패한 낮은 slow 후보는 다음과 같다. 이 주파수는 미승인 수치 후보이며 검증된 모델 간 shift로 사용하지 않는다.

| Case | 주파수 / Hz | determinant 절댓값 |
|---|---:|---:|
| H N128, join .003 | 211.115932444 | 3.21422e-6 |
| G N64, join .007 | 211.209128318 | 1.45219e-6 |
| H N128, join .007 | 211.103771492 | 5.54877e-6 |
| Primary numerical refinement의 fallback | 211.115931158 | 3.25532e-6 |

Gate는 determinant/join <1e-6, surface/junction <1e-10으로 유지됐다. H N128 join .004의 낮은 후보 하나가 gate를 통과했지만, 이를 미승인 후보와의 mode 차이 또는 전체 저주파 spectrum 확정으로 해석하지 않았다. 모든 ordinary case의 약1.905 kHz slow/rapid branch 및 계획된 high-frequency refinement 두 개는 수용됐다.

Native source refinement 및 next-join의 높은 branch 주파수 차이는 약1.5e-4–2.8e-4 Hz이고 tighter numerical 변화도 약1.7e-4–2.0e-4 Hz다. 따라서 현 결과는 작은 효과를 정밀한 물리적 shift로 분리한 증거가 아니다. 보고서는 이를 수치 sensitivity와 함께 제시한다.

Raw runner의 actual exit 1과 `REJECTED_OR_EMPTY`는 모든 계획 계산을 마친 뒤 일부 residual rejection 때문에 나온 결과로 보존됐다. 별도 최종 집계 `COMPLETED_WITH_REJECTED_CANDIDATES`와 모순되지 않는다. Suite 통과가 이 네 root의 수용을 대체하지 않는다. [COWLING_REVIEW.json](../../data/low_density_crust_20261117/review/COWLING_REVIEW.json)에 raw files/hash/acceptance readback을 남겼다.

## 검토 중 수정된 사항

1. TF 문서의 neutrality를 pointwise local neutrality로 오해할 수 있는 문구를 Wigner–Seitz cell 전체 전하 중성으로 정정했다. 계산 데이터는 바꾸지 않았다.
2. Numerical refinement의 bounded retry가 ordinary surface_h로 돌아갈 수 있는 runner 오류를 발견했다. 원 실행 runner와 raw outputs를 보존한 뒤 재사용 runner만 요청된 surface_h를 유지하고 effective settings를 기록하도록 수정했다. 실제 해당 fallback이 일어난 failed low-root 한 기록은 paired refinement에서 제외했다. 세 surface 설정의 configuration regression과 수정된 solve callsite를 확인했다.
3. 에너지밀도와 수밀도의 jump를 뜻하는 문구를 epsilon/n 비율의 jump와 혼동하지 않도록 정정했다.

원 runner snapshot SHA는 `780f8df344a93af3805fca3702187376f85b44b40ba0cd8fedf1f829b2cdf731`이고 수정된 runner SHA는 `723bb8085eb863fe8e07378dc25bb9112458b2117657a2eecefd6facfec239a3`다. Snapshot은 byte 재현 자료이며, `__file__`에 의존하는 ROOT 때문에 원 `scripts/run_low_density_stellar_20261117.py` 경로에 복원하여 실행해야 한다. 수정된 재사용 runner로 과거 수치 결과를 소급해서 계산했다고 표기하지 않았다.

전체 frontend/pipeline suite는 실제 **1453 passed, 9 skipped, 4 subtests passed, exit 0**이었다. 실행 전후 및 검토 시점의 213개 source SHA가 동일했다. Runner의 이후 configuration 수정은 이 core source 범위 밖이며 별도 세 설정 검증과 원 실행 identity를 제공한다. [FULL_SUITE_REVIEW.json](../../data/low_density_crust_20261117/review/FULL_SUITE_REVIEW.json)을 참조한다.

## 최종 범위와 다음 계산

현재 H는 finite native matching, 명시적 causal whole-EOS 구성, 조건부 TOV, 수용·거절을 구분한 finite Cowling 결과를 재현 가능하게 제공한다. 남은 여섯 native source gap, 낮은 Cowling 후보의 접합 잔차, causal lower closure의 spectrum sensitivity는 다음 계산의 대상이다. 전역 평형, continuum limit, 물리적 conversion rate, radial stability, complete spectrum 또는 lower closure의 미시적 타당성은 이번 승인 범위에 없다.

I에서 낮은 mode의 residual 원인을 구체적 integration/matching control로 분리하고, source mode8의 conditioning/gradient/line-search를 고정 functional 안에서 진단하겠다는 우선순위는 실제 남은 문제에 대응한다. 실패를 물리적 no-go로 승격하거나 잔차 기준을 느슨하게 만드는 단계는 제안하지 않는다.

최종 review receipt는 [REVIEW_RECEIPT.json](../../data/low_density_crust_20261117/review/REVIEW_RECEIPT.json)이다. 이 검토 뒤 추가될 RESUME 인계 정보, 파일 manifest 및 commit/package receipt는 과학 source/results를 변경하지 않는 metadata로 취급한다.
