# A whole-crust 최종 독립 검토

검토 기준: `8c917baf3f5975a6109c3b6e59177079c285c4e0`, 2026-11-17.
권위 문서는 `docs/matched_eos_reference_20261117/PLAN_KO.md`의 A와
`docs/whole_crust_20261117/PLAN_KO.md`다. 생산 구현과 분리된 검토이며 구현 파일은 수정하지 않았다.

## 판정 범위

물리 구성식과 출처 분리는 A의 허용 범위에 맞는다. pressure-primary outer는
원본 n/P를 유지하고 첫 epsilon/n에만 적분상수를 고정하는 명시적인 새로운 표현이다.
원본 에너지·chemical potential의 모든 knot를 보존한 미시적 BSk22 복원으로
해석하면 안 된다. lower bridge와 surface는 그 한계를 문서·contract에 기록하고 있다.
상단은 기존 생성·matching 결과의 knot를 재사용하지만 그 사이의 연속 law는
별도로 검증한 cubic P(mu_E)이다. 그러므로 원래 quintic/CSS 생성 law 자체의
모든 점을 정확히 재사용한다는 주장은 성립하지 않으며 생산 RESULT의 limitations는
이 차이를 명시한다.

**초기 판정: 조건부 적합.** 공개된 consumer 계약에서 허용한 CSS endpoint 압력
불일치에 대한 구체적인 실패를 발견했다. 아래 R1의 수정·회귀 확인과 생산 파일의
실제 native write/readback 확인 전에는 전체 A 완료를 이 문서만으로 인증하지 않는다.

## 물리·수학 검토 결과

- **derived:** F=epsilon/n에 대해 dF/dn=P/n^2를 적분하는 Python/Fortran 식을
  각각 대조했다. Gamma=1 극한, mu=F+P/n, cs2/c2=Gamma P/(epsilon+P)가 일관된다.
  양의 첫 F와 증가하는 P,n에서는 F가 양수이며 각 구간 sound-speed 극값은 양끝에 있다.
- **derived:** cubic P(mu)의 n=P_mu, epsilon=mu*n-P는 양끝 n/P/epsilon과 first law를
  동시에 보존한다. susceptibility와 causal margin mu*n_mu-n의 극값 검사 계수는 맞다.
  양의 mu,n_mu이면 epsilon_mu=mu*n_mu>0이므로 양의 왼쪽 energy는 구간 양성을 보장한다.
  Python은 별도로 epsilon 극값도 검사한다.
- **derived:** surface P=K*n^gamma, epsilon=a*n+P/(gamma-1)는 a>0,
  1<gamma<=2에서 안정·인과적이다. n=0에서 P=epsilon=0, mu=a>0이며
  h=log(mu/a)는 surface를 기준으로 한 무차원 enthalpy다. 유한 cutoff 변화는
  생략 enthalpy의 수렴을 검증할 수 있지만 stellar M–R 수렴을 뜻하지 않는다.
- **implementation-reviewed:** header mass 변경은 Q3/Q6/Q7 등을 재표현하여 차원 있는
  물리 상태를 보존한다. 상대 source energy에 임의의 C*n을 더하는 경로는 발견하지 않았다.
  첫 outer energy와 source native upper energy를 조건으로 bridge를 구성한다.
- **implementation-reviewed:** 17개 isotope gap은 explicit regularization으로 표시되고
  그 내부 composition은 미정이다. EOSPiece의 row provenance, native 미사용 prefix,
  원본 snapshot, 출력 hash가 서로 다른 출처를 구분한다. 새 native kernel이 계산했다고
  해서 public outer를 native microscopic EOS로 이름 바꾸는 경로는 발견하지 않았다.

## 독립적으로 실행한 consumer 확인

기존 테스트의 단순 analytic fixture를 이용해 생산 구현을 직접 호출했다.
실행 환경은 `CODEX_PRIMARY_RUNTIME_PYTHON`,
`PYTHONPATH=python_frontend:/tmp/ns_cowling_comparison_pydeps`다.

세 smooth seam P=0.1, 0.4, 0.9에서 nextafter 양쪽의 n/epsilon/mu 상대차 최대는
각각 1.7763568394002506e-16, 0, 0이었다. 정확한 sharp fixture의 P=1.6에서
low/high density 차는 1이고 latent identity 잔차는 -5.329070518200751e-15였다.
이 확인은 서로 다른 region law를 사용하면서 sharp density gap을 보간하지 않는
정상 경로를 실제로 통과했다는 증거다. 생산 전 구간을 대신하는 검증은 아니다.

### R1: 허용된 CSS pressure mismatch 내부에서 consumer가 실패함

**initially implementation-verified failure; 아래 후속 검토에서 해결 확인.**
`from_contract`는 doublet의 압력 상대차를 2e-12까지 허용하고,
`_upper_search_pressure`에서는 high row의 search pressure를 low pressure로 바꾼다.
그러나 `_state_one`의 transition 판단 폭은 8*machine epsilon뿐이며 high interval의
실제 왼쪽 압력은 그대로다. 두 tolerance 사이의 pressure에서는 high interval을
선택하면서 root bracket 밖의 목표를 전달한다.

재현은 existing `consumer_contract(sharp=True)`에서 `upper.p[2] *= 1+1e-12`만 바꾼다.
이때 admitted doublet는 [1.6, 1.6000000000016001]이며 midpoint
1.6000000000008001에서 `state_at_pressure`가
`ValueError: f(a) and f(b) must have different signs`를 냈다.
정확한 pressure doublet 정상 테스트로는 이 실패를 발견하지 못한다.

수정은 admitted endpoint pressure band와 실제 소비자의 경계 규약을 일치시켜야 한다.
원본 물리 행을 조용히 바꾸거나 밀도 jump 내부를 보간하는 것으로 해결하면 안 된다.
positive/negative pressure mismatch 양쪽, band 내부와 양쪽 smooth 영역을 검증할 필요가 있다.
생산 doublet가 이보다 훨씬 정밀하더라도 공개 API가 허용한 contract의 문제는 남는다.

## 생산·출처 관련 확인 조건

root의 실제 생산 진단에 따르면 최초 native writer 실패 원인은 raw public outer의
비표준 group abundance였다. provider tuple 값 1을 native reader가 표준 nuclear
number fraction Yi로 해석하여 Yh=A*Yi=A가 되고 row_status=3이 된다. 이는 파생
row_status=0 문제와 구별된다. phase-only derived row를 status=1로 canonicalize하는
처리도 필요하지만 최초 실패의 원인으로 혼동하면 안 된다.

**derived / source-semantic correction:** 217개 확인된 pure-isotope outer knot에만 한정해
Yn=Yp=0, 단일 nuclear group이라는 inventory가 성립하면 baryon conservation이
Yh=1을 고정한다. Ymu=0, Ye=Z/A는 charge closure의 독립 교차검사다. 따라서 표준
number fraction은 Yi=1/A이다. raw tuple=1을 유지해 Yh=A로 읽는 것은 그 inventory와
모순된다. A>0, 0<=Z<=A, single-group 및 field presence 확인, raw tuple=1 조건,
source identity를 요구하는 명시적인 adapter로 표준 표현을 구성하는 것은 타당하다.
필드 부재를 zero로 간주하거나 inner crust까지 확장해서는 안 된다. raw 파일을 보존하고
행별 변환량과 semantic assumption을 기록해야 하며, 이 수정 이후에는 composition이
모두 unchanged라고 기록해서도 안 된다. 열역학 energy gauge는 이 adapter로 바뀌지 않는다.
이 문서 작성 시점에는 최종 생산 artifact가 아직 제공되지 않았다. 기존 집중 검사 19 Python/29 native 통과는 root 보고이며,
본 검토자가 같은 48개를 재실행해 얻은 결과로 혼동하지 않는다. 전체 과거 1301개 suite는
이번 변경 검토에 재실행하지 않았다.

runner의 `code_commit=HEAD`는 새 파일이 untracked인 생산 시점에는 그 파일의 identity를
보장하지 않는다. native executable hash 외에 runner/consumer/kernel/test의 source hash,
실제 patch 또는 최종 commit에 연결된 기록이 필요하다. 입력 reference와 기존 upper의
hash 검증, output hashes, 원본 snapshot은 구현되어 있다. source-provenance의 129행
상수는 현재 129-knot source campaign에 한정된 로직이며 범용 source loader로 볼 수 없다.

초기 검토 파일 SHA256:

| 파일 | SHA256 |
|---|---|
| whole_crust.py | 052a669c3941ae747bebc7383188fe9fd8c2293fd623fb7f79c6c8fad4876839 |
| whole_crust_kernel.f90 | bc0cc951da37b7b736a80c3c6fe69cabd8e890c45ee0de919c04e58476e71dec |
| build_whole_crust_20261117.py | b2dba772a9c89784e25c9ffe617a06c887118b6aa48ad102f858bf62b292722d |
| test_whole_crust.py | 60e820c6816ee6e77c684f4b47372fd87fafbbf86029118d531680f2da3edcc3 |
| test_whole_crust_native.py | f3bf6d4450e4b78f6300d0794356f2afef9e18b762e8b6728b08c5272b91fb58 |

## 실제 사용 범위

A의 usable baseline은 `[0,0.5] fm^-3`의 명시적 equilibrium barotrope를 뜻한다.
압력 입력 consumer는 P=0 vacuum, analytic surface, outer, lower bridge, native upper
연속 구간 및 CSS one-sided 상태를 제공한다. CoMPOSE의 finite-positive-density 표만
단독 전역 spline해서는 같은 EOS가 되지 않는다. 반드시 REPRESENTATION.json과
WholeCrustBarotrope를 함께 사용하거나 같은 piece laws를 구현해야 한다.

이 결과로 목표 별의 central density가 지원된다고 결론내리거나 composition-dependent
buoyancy/frozen composition closure, Seidov flag로부터 실제 mode 안정성, TOV/Cowling
수렴을 주장할 수 없다. A에서 surface와 하단 연결을 만들고 B에서 core support를
확장하며 D에서 stellar consumer를 검증하는 기존 순서와 일치한다.


## 후속 독립 검토: R1과 composition adapter

**R1 resolved / implementation-verified.** 수정본은 admission을 통과한 doublet의 실제
압력 encoding band `[min(P0,P1),max(P0,P1)]`와 8 epsilon 여유를 명시적으로 보존한다.
band 안에서는 source one-sided n/epsilon를 선택하고 smooth root solver를 호출하지 않는다.
원본 배열은 바꾸지 않는다. 이는 유한 폭의 물리적 mixed phase를 추가하는 규약이 아니라
동일 전이압의 수치 encoding 차이를 처리하는 규약이다. 이 band 안에서 반환한 값으로
finite-difference first law를 검사하거나 물리적 compressibility를 추론해서는 안 된다.

검토자가 수정본에 대해 다음 선별 실행을 직접 수행했다.

```text
pytest python_frontend/tests/test_whole_crust.py -q \
  -k 'encoding_band or provider_outer or phase_only'
14 passed, 19 deselected
```

이는 signed band 2개, phase-only canonical metadata 1개, adapter 정상 1개,
모호하거나 혼합된 source 거절 10개다. 추가로 band 경계의 바로 바깥 `nextafter`
압력을 직접 평가했다. 양의/음의 1e-12 상대 discrepancy 모두 아래쪽 smooth branch는
n<4, 위쪽 smooth branch는 n>5였고, band 경계·중간에서는 low/high가 정확히 n=4/5를
반환했다. 최초 midpoint bracket 실패는 재현되지 않았다.

adapter는 모든 일곱 composition presence mask, phase0, 정확한 free-nucleon/muon
zero, 양의 정수 A/Z, Z<=A, 기존 decoded Yh=A, printed charge closure를 검사한다.
runner는 앞서 raw NP=4, species IDs={10,11,0,1}, 하나의 group, raw group value=1,
실제 source hash 및 217개 행을 확인한다. 따라서 public whole table 전체를 임의로
재해석하는 경로가 아니다. accepted correction ledger는 원본 row_status=3과
새 status=2를 함께 보존한다. thermo-only derived 행은 PARTIAL=1이며 species masks는
false여서 missing composition을 관측된 zero abundance로 바꾸지 않는다.

runner는 `implementation_sha256`도 추가하여 생산 source identity를 보강했다.
후속 검토 파일 SHA256은 다음과 같다. 이 값은 이후 production receipt에 있는 값과
직접 비교해야 하며, commit 또는 실행 성공 자체를 대신하지 않는다.

| 파일 | SHA256 |
|---|---|
| whole_crust.py | f1bfdfbd5615733e7f2636b2b2ad826805d1ce8b5d84e37ac4f0ee868ca046be |
| build_whole_crust_20261117.py | cfa6988b7cbd23c21a423dbc45531af0308671fbd2a18d5ade2d8bfb62803e11 |
| test_whole_crust.py | 0d982d97aba3741e6f4c8207b07cd22242214ada5fd83db614e1c9111f34d31b |

**이 시점 판정:** 발견된 R1은 해결되었고 adapter와 canonical metadata는 확인했다.
최종 native production/write/readback artifact 확인만 남아 있으므로 전체 A의 최종 PASS는
아직 선언하지 않는다.

## 최종 생산 artifact 독립 확인과 A 판정

**최종 판정: A 범위 적합, implementation-verified.**
`data/whole_crust_20261117/production_v2/`의 실제 파일을 확인했다.
생산 source commit은 `5c0e94649cd586d63323a9ec3fe908fd91586d6c`, artifact checkpoint는
`f7a4037bb089a167098265f4c7f3820c734b9dec`다. 위 초기·후속 검토의 pending 상태는
실패를 감추지 않기 위한 이력이며 이 최종 판정으로 갱신한다.

검토자가 생산 receipt를 읽는 것 외에 독립적으로 수행한 검증은 다음과 같다.

- 구현 파일, public source snapshot, 기존 두 upper table snapshot과 PAIR receipt,
  8개 output table의 SHA256 **68개를 실제 파일에서 다시 계산하여 일치**를 확인했다.
  production의 public source는 재다운로드 snapshot이며 앞선 PHYSICS_REVIEW의
  `eos_generator/eos_data/SPG_M1` 진단 파일과 byte identity가 같다는 뜻은 아니다.
  생산의 권위 있는 입력 hash는 CAMPAIGN 및 복사된 PROVENANCE.json이다.
- native output `eos.nb`와 `eos.thermo`를 직접 읽어 n/P/epsilon을 복원하고, 각
  REPRESENTATION.json에서 새 consumer를 생성하여 **8개, 총 5638행**을 재평가했다.
  source doublet의 low/high를 구분했으며 n과 epsilon의 최대 상대 오차는 각각
  **3.1086244689504383e-15**였다. 이 검사는 생산 runner의 저장된 PASS를 반복해서
  읽는 것이 아니라 저장된 실제 수치와 연속 consumer를 다시 대조한 것이다.
- 8개 consumer 모두 P=0에서 n=P=epsilon=h=0을 반환했다. 모든 RESULT의
  native write/readback, upper exact-interval certification 및 native hard gate가
  PASS인 기록과 row provenance 길이를 확인했다.
- CSS 4개 variant의 실제 transition에서 delta n은
  **0.00012336473282301141 fm^-3**이며
  `abs(delta epsilon-mu_low*delta n)`은 모두
  **1.942890293094024e-16 MeV fm^-3**였다. sharp gap 보간 없이 양쪽 상태를 재평가했다.
- `validation/FOCUSED_FINAL.json`은 실제 원격 pytest의 returncode=0,
  **62 passed** stdout을 보존한다. 이 62개는 root 실행 기록이고, 앞서 검토자가 직접
  수행한 14개 회귀 및 위 artifact 수치 대조와 근거의 출처를 구분한다.

| construction | baseline | lower row36 | lower row37 | surface gamma4/3 |
|---|---:|---:|---:|---:|
| quintic_P source knots | 704행 | 708행 | 707행 | 704행 |
| CSS source knots | 702행 | 706행 | 705행 | 702행 |

과학적 baseline은 quintic_P/row40/gamma5/3 하나다. 나머지 일곱 표는 endpoint,
upper construction, surface closure 민감도 조합이다. 8개 독립 microscopic EOS를
생성한 결과로 세면 안 된다. 모든 표의 실제 core upper support는 0.5 fm^-3다.

생산 outer의 최대 energy-per-baryon 보정은 0.004915606957638374 MeV,
outer cs2/c2 범위는 6.7299963573426004e-6–0.0026519594679298157이며,
17 isotope gap label이 보존되어 있다. lower bridge의 최대 cs2/c2는 row40에서
0.006974341835140935, row36에서 0.0338087582965788, row37에서
0.016125555528160033이다. 이는 임의 에너지 이동 없이 선택한 endpoint 민감도다.
표면 export cutoff n/n0=1e-8에서 생략 h는 gamma5/3의
2.400330294743814e-10, gamma4/3의 1.7826705611459798e-7로 기록되어 있다.

별도의 `OUTER_COMPOSITION_NORMALIZATION.json`은 217개 source-specific 변환을
보존하며 실제 printed-charge closure 최대 잔차는 4.7619047593983055e-9다.
`OUTER_RECONSTRUCTION.json`의 composition unchanged는 **이미 명시적으로
normalize한 입력에서 pressure reconstruction이 추가로 composition을 바꾸지 않는다**는
뜻으로 읽어야 한다. raw provider tuple 전체를 그대로 보존했다는 뜻은 아니다.

native overall PASS는 C1–C6/C10/C11 hard gate의 의미다. C7/C8 finite-difference
최대 상대 잔차 약 **5.27475% / 1.88479%**는 그대로 남아 있으며 이를 zero로
인용하거나 모든 수치 미분 항등식이 동일 정밀도로 만족한다고 주장하면 안 된다.
생산 interval law의 first law는 analytic construction 및 exact interval certification으로
검증했고, 여러 piece/phase를 가진 유한 표의 별도 finite-difference diagnostic과
구분한다. generic tabulated consumer를 사용할 때의 재보간 오차까지 인증한 것은 아니다.

이 증거는 원래 A의 최소 종료조건인 source-labeled surface/outer/inner/core 후보,
명시적인 boundary law, continuous stability/causality, source-preserving energy policy,
실제 native 생산과 CoMPOSE 왕복을 충족한다. B의 필요한 core support 확장 및
D의 TOV/Cowling·composition closure 검증은 완료 범위에 포함되지 않는다.
