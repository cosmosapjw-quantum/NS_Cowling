# TF 나머지 압력 감소 구간: 실제 정련과 두 상 공존 후보

기준 `a49aaf9f82ddfda484e0f965aefdfe73aa7918f2`. 나머지 원래 N16 감소 **50구간을 모두 계산**했다. 기존 최악 3구간의 정련은 재실행하지 않았다. 47구간은 N64의 원 endpoint 두 점에서 감소가 사라졌고, TM1e 한 구간은 N128의 endpoint·중간점 세 점에서 해소되었다. 남은 FSU2R/TM1e slab–tube 두 구간은 N128에서도 원시 압력이 감소한다. 그 원시 실패를 보존한 채 실제 branch 에너지로 두 개의 common-tangent 연구 후보를 계산하고 별도 열역학 CSV를 native 생성·재읽기했다. **48/50의 표본 감소 해소와 2개의 유한 branch 공존 후보**이며 전체 연속 EOS·전역 최소 승인이 아니다.

## 계산 및 실패의 정확한 범위

| 모델 | 남은 원 구간 | endpoint 합집합 | 신규 N32/N64 mesh-density 점 | N64 음의 ΔP 구간 |
|---|---:|---:|---:|---:|
| FSU2R | 17 | 33 | 64 | 1 |
| FSU2H | 17 | 32 | 62 | 0 |
| TM1e | 16 | 32 | 64 | 2 |
| 합계 | 50 | 97 | 190 | 3 |

기존 정련에서 FSU2R 0.048, FSU2H 0.050999999의 N32/N64 결과, 총 4 endpoint-mesh 값을 재사용했다. N64의 세 미해소 구간에만 N128 endpoint와 실제 중간점 9개를 추가했다. 따라서 **신규 generator 고유 model–mesh–density 점은 199개**다. N64 97개 endpoint의 중심 재최소화와 ±h, ±h/2는 별도의 485회 branch solve이며, 중심은 기존 density와 겹친다. 공존 후보 두 실행 76회와 TM1e 정밀화 26회, 총 102회 branch solve를 별도로 기록한다. 서로 다른 종류의 실행 수를 고유 density 수로 더하지 않는다.

새 generator 시도 1,412개는 STATIONARY 1,010개, UNIFORM_STATIONARY 226개, **EVALUATION_FAILED 160개**, **REQUIRES_REFINEMENT 16개**다. 160개 오류는 모두 `vector field Newton line search failed`이고, 16개는 optimizer mode 0이어도 native 물리 stationarity gate를 충족하지 못했다. 후속 family 탐색·정련 시도와 함께 attempts/profiles에 남겼다. 485+102 branch solve는 모두 STATIONARY다. 과거 87개 실패 기록과 직전 정련 8개 실패/정련요구는 기존 경로·archive 그대로 유지했다.

Fortran 물리 producer는 TF 모델·공간 functional·screened electron·native SLSQP를 그대로 사용한다. Python은 요청, profile 전달, 하위 실행, 파일 진단·요약만 한다. 최대 4개 계산 process, 각 BLAS/OMP/MKL thread 1이었다. generator는 완성 node마다 attempts/profiles/selected를 flush하고 각 실행 receipt를 보존했다. 물리 stopping tolerance, source density, 상 선택 guard, 생산표를 변경하지 않았다.

## 같은 실제 branch에서의 압력–에너지 검증

TF 내부 convention은 질량에너지 MeV, 길이 fm, `ħc=197.3269804 MeV fm`, T=0다. density는 fm⁻³, ε와 P는 MeV fm⁻³, μB는 MeV다. 정의는 P=nμB−ε이고, β 평형의 중성 branch에서 dε/dn=μB다.

`blocker_tf_branch_derivative.f90`는 선택된 실제 nn/np profile과 반경을 읽어 중심을 재최소화하고, 각 ±h, ±h/2에서 중심 profile을 density 비율로 rescale한 뒤 동일 dimension/seed의 실제 native TF를 다시 최소화한다. h=min(10⁻⁵,10⁻³n) fm⁻³다. 독립 에너지 차분으로
`μ_fd=[ε(n+h)−ε(n−h)]/(2h)`, `P_fd=n μ_fd−ε(n)`
를 계산한다. 이는 P=nμ−ε를 다시 대입하는 Euler 항등식 검사가 아니다.

97 stencil, 총 485 상태 모두 STATIONARY이고 stencil 내 geometry 변경은 없었다. h/2의 최대 |P_fd−P_native|는 **8.68944×10⁻⁷ MeV fm⁻³**다. 정규화한 density profile의 최대 weighted 거리 0.002043, 최대 상대 반경 변화 0.0006141도 저장했다. phase=9라는 같은 정수만으로 branch를 연결하지 않았다. 일반 checker의 phase 처리도 변경하지 않았다.

N32→N64 최대 endpoint 압력 상대 차이는 FSU2R 1.7338%, FSU2H 2.6475%, TM1e 2.2317%다. 독립 최소화와 branch 순위 변경이 섞인 관측값이며 순수 공간 truncation 오차나 continuum bound는 아니다. N32의 비정형 profile은 N64에서 canonical rod 등으로 바뀐다.

## 잔여 원 구간과 분류

| 모델·원 bracket [fm⁻³] | N64 ΔP | N128 세 점의 ΔP | 판정 |
|---|---:|---|---|
| FSU2R [0.073000006, 0.074000001] | −0.00120285 | +0.00105239, −0.00225933 | slab–tube branch 교차가 남음 |
| TM1e [0.061999999, 0.063000001] | −0.000103121 | +0.00205552, +0.00203415 | N128 세 점 slab로 해소 |
| TM1e [0.081, 0.082000002] | −0.000545146 | +0.00166969, −0.00222468 | slab–tube branch 교차가 남음 |

N128 FSU2R의 tube−slab 에너지/핵자 차이는 좌끝 +3.39462×10⁻⁴, 중간 +2.65567×10⁻⁵, 우끝 −2.78111×10⁻⁴ MeV다. TM1e의 두 번째 구간은 각각 +4.74649×10⁻⁴, +1.76768×10⁻⁴, −1.14534×10⁻⁴ MeV다. 모든 여섯 native family의 원시 에너지를 보존했다.

따라서 48구간의 개선은 **수치 mesh 및 그에 따른 유한 branch 순위의 민감성**으로 분류한다. 남은 두 구간은 독립 에너지 미분 오차보다 큰 압력 감소와 실제 canonical branch 에너지 교차가 함께 있어 **탐색된 branch 최소-envelope의 비볼록성**을 지지한다. 아직 탐색하지 않은 더 낮은 구조를 배제하는 전역 최소 증명은 아니다. vector-field 평가 실패는 별도 수치/최적화 경로 오류이며 물리 평형의 부재를 뜻하지 않는다.

## 공존 후보의 유도와 실제 계산

각 상은 screened charge profile을 가진 Wigner–Seitz cell이며 **cell 평균으로 각각 중성**이다. 공간의 모든 점에서 np=ne를 강제한 것은 아니다. 두 중성 상의 거시적 혼합에서 표면·계면 에너지의 체적당 기여를 무시하고 β 평형 T=0의 유효 보존량을 baryon 하나로 둔다. 고밀도 상 체적분율 χ에 대해
n=(1−χ)nA+χnB, εmix=(1−χ)εA(nA)+χεB(nB)다.
고정 n에서 이 식을 최소화하면 εA′(nA)=εB′(nB)=μt와
εB−εA=μt(nB−nA), 즉 PA=PB=Pt를 얻는다. 혼합선은 εmix(n)=μt n−Pt다.

이는 **별도로 중성인 두 상의 Maxwell형 연구 구성**이다. 전역 중성만 요구하여 두 상의 charge density도 변분하는 two-charge Gibbs 문제를 풀지 않았다. 전자의 전기화학 퍼텐셜은 전자 전하 −e에서 μ̃e=μe−eΦ이므로 각 상의 전위 기준·Galvani 차이를 정하지 않고 순수한 μe 값의 같음을 추가하지 않는다. 계면 전위나 계면 에너지를 계산했다고 주장하지 않는다.

`blocker_tf_common_tangent.f90`는 모든 함수값을 실제 warm-start native 최소화로 계산한다. energy derivative와 2차 차분으로 Newton update를 수행하고 branch geometry 및 native stationarity를 검사한다. 원 bracket 밖으로 이동한 값도 실제 solve로 계산했으며 표를 외삽하지 않았다. FSU2R 실제 평가 domain은 [0.072738561,0.074390446] fm⁻³, TM1e 최종 정밀화 domain은 [0.081213880,0.082384568] fm⁻³다.

| 최종 후보 | nA [fm⁻³] | nB [fm⁻³] | Pt [MeV fm⁻³] | μt [MeV] | h/2 ΔP | h/2 Δμ |
|---|---:|---:|---:|---:|---:|---:|
| FSU2R slab/tube | 0.07275856095 | 0.07431969379 | 0.2342521711 | 952.0117462 | +5.53141×10⁻⁸ | +7.43228×10⁻⁷ |
| TM1e slab/tube | 0.08122389021 | 0.08237448523 | 0.3275351091 | 953.1378791 | −8.63041×10⁻⁸ | −1.05018×10⁻⁶ |

두 후보 모두 사전에 둔 root 목표 |Δμ|<2×10⁻⁶ MeV, |ΔP|<2×10⁻⁷ MeV fm⁻³를 최종 h/2 검사에서도 만족했다. 이는 계산 목표이며 기존 물리 tolerance를 완화한 것이 아니다. TM1e 최초 h/2 ΔP=−3.27336×10⁻⁷은 목표를 넘었으므로 그 기록을 보존하고 h=10⁻⁵, SLSQP ftol=10⁻¹⁴의 한 차례 정밀화로 해결했다. FSU2R 추가 실행은 하지 않았다.

각 공존구간 내부의 1/4, 1/2, 3/4 지점에서 양 branch를 실제 재최소화했다. 공존선보다 위에 있는 최소 에너지 차이는 각각 2.20274×10⁻⁶ 및 1.72545×10⁻⁶ MeV fm⁻³다. 유한 내부 표본에서 두-branch 혼합이 에너지를 낮춘다는 검산이며 모든 구조·모든 density에 대한 convex-hull 증명은 아니다.

`candidate_chords/{FSU2R,TM1e}/candidate_thermo.csv`는 Fortran이 생성하고 다시 읽은 endpoint·중간점 3개씩의 **별도 연구 후보**다. 정의된 n, ε, P, μB, χ만 출력한다. μQ와 핵/lepton 조성을 임의의 0으로 채우지 않았으며 `SCOPE.json`에서 모든 조성 presence를 false로 표시했다. 원시 압력 감소 표는 그대로 남아 있고 생산 registry나 EOS를 후보로 대체하지 않았다.

## Provenance 및 전달

초기 producer 42개 의존의 byte identity는 `SOURCE_IDENTITY.json`에 PASS로 남아 있다. N128 두 실행은 native child가 각각 exit 0, 9개 선택 점과 표를 생성했지만 병렬 작업 중 **연결되지 않은** `compose_eos_core_domain.f90`, `compose_eos_crust_attachment.f90` 내용이 바뀌어 상위 source-set binding이 실패했다. 원 DISPATCH의 EVIDENCE_OR_EXECUTION_FAILED는 보존했다. `ambiguous_N128/READBACK_RECOVERY.json`은 실제 isolated binary 생성·실행에 관계된 물리·입출력·실행 의존 15개가 변하지 않았음을 비교한 후 native table readback을 별도로 통과시킨 기록이다. 물리 재계산으로 덮지 않았다. 이후 부모의 production diagnostic 통합은 이 과거 실행의 source identity와 구별한다.

실행·재현 근거는 `PLAN.json`, 모델별 `RECEIPT.json`, 원 DISPATCH, attempts/profiles/selected, 97 derivative receipt, 3 tangent receipt, `COMMON_TANGENT_SUMMARY.json`에 있다. `MANIFEST.json`과 `SHA256SUMS.txt`는 전달 시점 파일 identity다. Python replay wrapper는 최초 실행 뒤 올바른 build-directory 환경변수와 덮어쓰기 방지 guard를 추가했으므로 그 최종 hash를 최초 wrapper의 byte identity로 주장하지 않는다. 실제 native producer 소스·binary와 명령·입력·출력 identity는 별도로 기록되어 있다.

다음 행동은 두 후보를 **명시적인 local-neutral 연구 분기**로 소비하여 이후 observables 계산에 연결하는 것이다. 원 TF raw support를 이 분기로 자동 승인하지 않는다. 전역 phase 선택을 주장하려면 공존 domain의 경쟁 구조와 mesh 수렴을 추가해야 하고, two-charge/계면 물리가 목표이면 그 물리 문제를 별도로 정식화해야 한다. 이미 계산된 48개 표본 개선과 두 후보를 다시 동일 입력으로 반복하는 것은 필요하지 않다.
