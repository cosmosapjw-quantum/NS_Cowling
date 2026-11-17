# E — SLy4 고밀도 연장, 최대질량과 선택된 Cowling 모드

D 커밋 `300e5a0532dd137cb09738568fca0f35cfcc9282`의 push와 원격 SHA 확인을 먼저 마쳤다.
현재 사용자 요청은 E의 실행과 push도 승인한다. 이 문서는 게시 전 검증된 E 산출물을 설명하며,
게시 후 원격 SHA receipt는 다음 단계의 입력으로 보존한다.

## 결론

보존된 B hybrid crust와 동일 native SLy4 core를 사용하는 평형 barotrope에서
최대질량 중심밀도는 **\(n_c\in[1.2060,1.2061]\;\mathrm{fm}^{-3}\)**다.
대표 계산은 **\(M_{\max}\simeq2.0503813\,M_{\odot,\mathrm{nom}}\)**,
**\(R\simeq10.0064\;\mathrm{km}\)**를 준다.
기존 상한 \(n_c=1\;\mathrm{fm}^{-3}\) 이후에 남아 있던 최대점 위치를 이번에 결정했다.
기존 D 자료는 그 상한에서 질량이 아직 증가한다는 사실을 올바르게 기록하고 있었다.

최대점의 native 평형 응답은 \(\chi_{\rm eq}=dP/d\epsilon=0.96807415\)다.
\(\chi_{\rm eq}=1\)은 \(n_B=1.2650056715\;\mathrm{fm}^{-3}\)에서 나타나므로
이번 최대점은 해당 barotropic 인과 경계보다 낮다.
이는 full microscopic causality나 모든 조성을 고정한 음속에 대한 증명이 아니다.

## 실제 확장과 보존한 부분

기존 native generator로 \(n_B>1\)의 상태를 새로 계산했다. source 매개변수,
rest energy 기준, 전자·뮤온 처리와 beta-equilibrium 조건은 유지했다.
외삽과 에너지 이동은 사용하지 않았다.

기존 773개 행을 보존하고 64/128/256개 native tail interval을 각각 추가한
coarse/medium/fine contract를 만들었다. 공통 실제 상한은
\(n_B=1.2649856715193053\;\mathrm{fm}^{-3}\)다.
native \(\chi_{\rm eq}=1\)보다 \(2\times10^{-5}\;\mathrm{fm}^{-3}\) 낮게 정해,
source knots 사이의 연속 보간도 \(\chi\le1\)을 만족시켰다.
이는 지원구간을 제한한 것이며 인공적인 causal completion을 붙인 것이 아니다.

기존 \(n,P,\epsilon,\mu,h\), 조성·mask, 하부 contract와 기존 upper 상태 배열은 그대로다.
기존 \(n=1\) 끝점에서 왼쪽과 새 오른쪽의 일방적인 \(\chi\)는 구분한다.
열역학 상태의 보존과 2차 미분의 좌우 차이를 동일한 주장으로 다루지 않았다.

| Contract | 새 tail 구간 | 연장부 연속 최대 \(\chi\) |
|---|---:|---:|
| coarse | 64 | 0.9999918439 |
| medium | 128 | 0.9999902103 |
| fine | 256 | 0.9999898003 |

주 산출물은 `data/stellar_extension_20261117/eos/fine/REPRESENTATION.json`이다.
hash와 원본 보존 비교는 `EOS_RESULT.json`, `SOURCE_BINDING.json`, 각 grid의 receipt에 있다.

## TOV 최대점과 수치 정확도

사용한 식은 D에서 검증한 동일 exact-contract enthalpy TOV solver다.
서명은 \((-+++)\), 반경은 km, 질량은 \(GM/c^2\)와 IAU nominal solar mass 기준으로 변환한다.
실제 진공 표면 \(h=0\)까지 적분하고 중심 정칙 전개를 사용한다.
새로운 TOV 물리식이나 consumer 구현을 추가하지 않았다.

최종 최대점의 권위 파일은 **`tov/MAXIMUM_REFINEMENT.json`**이다.
재현용 대표값은 다음과 같다.

| 물리량 | 대표 계산 |
|---|---:|
| \(n_c\) | 1.2060328050005111 fm\(^{-3}\) |
| \(P_c\) | 855.9246351636019 MeV fm\(^{-3}\) |
| \(M\) | 2.0503812906965106 nominal \(M_\odot\) |
| \(R\) | 10.006409438136822 km |
| \(GM/(Rc^2)\) | 0.3025705044461724 |

대표 좌표의 끝자리는 재현할 모델을 식별한다. 중심밀도의 과학적 정밀도 주장은
실제 확인한 \([1.2060,1.2061]\) bracket이다.
초기 adaptive DOP853 설정에서는 매우 평평한 최대점 부근의 한 질량값에
약 \(10^{-6}M_\odot\) excursion이 있었다. 원자료와 실행 snapshot을 보존했다.
더 작은 ODE 오차와 step으로 계산한 중심차분은 세 source grid 모두에서
bracket 왼쪽에서 양수, 오른쪽에서 음수였다.
강화된 source-grid/ODE/center 비교의 관측 mass 변화는 \(2.07\times10^{-9}M_\odot\) 이내다.
이 수치는 검사한 수치 설정에 대한 민감도이며 crust/matching 모형 불확실성을 포함하지 않는다.

초기 선택점 \(1.2060130400651856\)은 최종 bracket 안에 남아 있다.
Cowling 계산은 이 고정 밀도와 그 0.98배 밀도에 결박되어 있으며,
정밀화 대표 좌표로 입력을 사후 교체하지 않았다.
전체 radial stability spectrum은 계산하지 않았다.

## 응답과 인과 경계의 의미

찬 물질의 최소화된 에너지 \(\epsilon(n,y)\)에서 \(\epsilon_y=0\), \(\epsilon_{yy}>0\)이면

\[
\frac{d^2\epsilon_{\rm eq}}{dn^2}
=\epsilon_{nn}-\frac{\epsilon_{ny}^2}{\epsilon_{yy}},\qquad
\chi_{y}-\chi_{\rm eq}
=\frac{n\,\epsilon_{ny}^2}{\mu\,\epsilon_{yy}}\ge0.
\]

native routine의 `frozen` 출력은 \(Y_p\)를 고정하되 전자·뮤온을 공통 chemical potential로
재분배하는 응답이다. full fixed-species response와 구분했다.
이 응답은 \(n_B=1.2397422126\;\mathrm{fm}^{-3}\)에서 1을 지나며,
최대점에서는 0.97994079다. 두 검사된 응답의 경계 모두 최대점보다 높다.
정확한 정의와 native 식은 `eos/RESPONSE_SCOPE.json`에 있다.
위 Schur-complement 식과 scalar convention의 대수 검산은 `SYMBOLIC_CHECK.json`에 보존했다.

## Cowling 관측량

평형 perfect fluid, \(\Gamma_1=\Gamma_{\rm eq}\), \(\ell=2\), \(\delta g_{\mu\nu}=0\)을 사용했다.
주파수는 무한대 시간 기준이다. 중심 정칙 조건, 표면 \(\Delta P=0\),
연속 EOS knot의 \(W,V\) 연속 조건과 실제 mode profile의 절점 수를 검사했다.
sharp Maxwell interface를 이 E 별에 추가하지 않았다.

두 고정 중심밀도는 \(1.1818927792638818\)과 \(1.2060130400651856\;\mathrm{fm}^{-3}\)다.
반복 가능한 상세 주파수·배경·설정은 `cowling/SUMMARY.json`과
`cowling/TIGHT_BACKGROUND_COMPARISON.json`에 있다.
대표 보고값은 아래 machine-readable `CURRENT_OUTCOME.json` 및
`COWLING_RESULTS_KO.md`와 일치한다.

| 고정 중심밀도 [fm^-3] | 선택된 구성 | f [Hz] | p1 [Hz] | W/V 절점 |
|---|---|---:|---:|---|
| 1.1818927792638818 | tight mode ODE | 2886.216009 | 7493.942557 | f 0/0, p1 1/1 |
| 1.2060130400651856 | tight background + mode ODE | 2899.488673 | 7504.249128 | f 0/0, p1 1/1 |

실제 10개 주 계산 구성의 20개 mode와 추가 tighter-background 구성의 2개 mode를 계산했다.
모두 기존 \(2\times10^{-6}\) matching/join 기준을 만족했고 거절된 mode는 0개였다.
\(f\)의 \(W,V\) 절점은 0, \(p_1\)은 1이다.
source-grid와 ODE/표면/중심 설정을 바꾼 관측 주파수 span은 수 mHz 수준이다.
추가 tighter-background 비교의 변화는 \(f\)에서 \(3.84\times10^{-5}\) Hz,
\(p_1\)에서 \(-1.65\times10^{-4}\) Hz였다.
이는 fixed-density mode의 측정된 수치 민감도이며 엄밀한 전체 오차 상계는 아니다.
frozen response, crust elasticity, superfluidity, full-GR mode/damping은 이 계산의 범위 밖에 남는다.

## 구현·검증·다음 단계

새 runner·증거 회귀·결과 문서를 추가했고 기존 D consumer/native 코드는 유지했다.
전체 `python_frontend/tests pipeline/tests`는 **1404 passed, 9 skipped, 0 failed**,
실제 exit 0이었다. 7개 E 회귀는 이 1404개에 포함된다.
skip은 기존 raw fixture 미배치 항목이며 이번 물리 계산의 실패가 아니다.
테스트 전후 source/test hash가 같았다. 실제 명령·환경·stdout은
`FRONTEND_PIPELINE.json/log`에 있다.
한 차례의 새 독립 검토 결과와 선택적인 재계산은 `INDEPENDENT_REVIEW.json`에 기록한다.

Douchin–Haensel의 unified SLy4 문헌값 \(M_{\max}\simeq2.05M_\odot\)는
크기 수준의 독립 참고다. 이번에 보존한 hybrid lower EOS와 같지 않으므로
문헌 끝자리와의 일치를 검증 기준으로 삼지 않았다.
원전: https://arxiv.org/pdf/astro-ph/0111092v2 ,
DOI https://doi.org/10.1051/0004-6361:20011402 .

사용자가 지정한 다음 작업은 기존 후보 밖으로 넓힌 sharp Maxwell 탐색이다.
CompOSE의 NL3, PK1, TM1, DDME-X를 원전 매개변수로 재생성하고,
새 provider crust와 기존 admissible source를 함께 평가한다.
성공은 공통 \(P,\mu_B\), 양의 density jump, 올바른 envelope 방향,
구간별 안정성·인과성과 실제 CompOSE readback 근거로 판정한다.
이 후속 작업은 E와 별도 branch/receipt를 갖는다.
