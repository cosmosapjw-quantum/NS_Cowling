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

**initially implementation-verified failure; 수정 확인 대기.**
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
