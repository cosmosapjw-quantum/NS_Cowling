# EOS 생성 완료 범위와 core–crust matching 수정 — 2026-11-17

기준 HEAD는 478c7b5de1514ac512f91156f36551478f36f688, 기준 원격 CI는 #629 성공이다.
이 작업은 실제 EOS 생성의 남은 수치·입력 문제와 노출된 matching 방법의 이론/구현 결함을
병렬로 조사한 뒤, 실제 생성 테이블을 연결하는 계산까지 수행했다.
생산 물리 계산은 Fortran에 두었으며 Python은 frontend·transport·검증 보고에 한정했다.
독립 연구 계산과 legacy Python 참조 구현은 그 역할을 명시했다.

**18개 명시적 native 생성 경로의 실행·readback은 확보되어 있다.
18개 원저자 full EOS의 동일성이나 모든 연속 crust 구간이 완성되었다는 뜻은 아니다.**
남은 제약과 새 실행 결과는
[data/research_matching_20261117/CURRENT_OUTCOME.json](../../data/research_matching_20261117/CURRENT_OUTCOME.json)에,
모델별 표현·지원범위는 generator/REPORT_KO.md와
RETAINED_COHORT_MATRIX.json에 구분했다.

## 생성기: 실제 진전

- 기존 18모델 cohort와 TF source-grid의 등록 산출물 446개 hash가 일치했다.
  unchanged 물리 계산을 무조건 반복하지 않았다.
- 공개 조성에서 복원한 SPG surface input을
  public_composition_inferred_20261116이라는 명시적 variant로 추가했다.
  입력 JSON SHA256은 7d6974d11c5986eb0b10c04aa6c085b3f1508711735edd9f895abca0b7e0cfcd이며,
  기본 AME2408_candidate는 유지했다.
- 이 조건부 처방으로 M1–M5를 각129행, 총645행 실제 재최적화하고 조성 포함 native
  readback을 통과했다. 최대 KKT 잔차1.9967e-8MeV, 별도 CLD 에너지 평가와의 최대 차이
  8.53e-14MeV fm^-3, 고정 좌표 압력 차이6.75e-10MeV fm^-3였다.
  중요40상태의 다섯 차분 간격·네 이웃, 총800번 재최적화도 수행했다.
  상세 방법과 drip 경계의 느린 차분 수렴은 spg/REPORT.md에 있다.
- M5의 실제 실패는 부피분율 약1.01e-8인 균일상 극한을 인공적인 탐색 하한의
  정확한 좌표와 일치하지 않는다는 이유로 lower unconverged cluster로 취급한 문제였다.
  희박한 부피분율과 독립 beta-equilibrated bulk 에너지 일치를 함께 확인하도록 수정했다.
  유한 phase나 에너지가 다른 희박 상태는 계속 거부한다.
- SUPPORT.json의 변경·삭제가 봉인된 RESULT 검증을 통과하던 결함을 고쳤다.
  support의 hash 보존은 물리적 보간 승인을 뜻하지 않는다.

SPG 원저자 calibration과 변분/export 동일성은 미확정이다.
공개 source에서 관찰한 표면 미분 패턴과 완전 변분식 사이의 차이를 보존했다.
고정 밀도 에너지 교차의 ΔP 약0.00445–0.00807MeV fm^-3,
Δμ 약0.0606–0.0817MeV는 common tangent의 작은 잔차가 아니다.
원문·공개 코드·자료를 확인했지만 이 조건부 처방을 원저자 설정으로 승격하지 않았다.

## FSU2R의 N256 공존 후보

기존 실패 profile의 장 방정식·경계값 문제와 optimizer 종료 문제를 구분했다.
SLSQP mode8은 positive directional derivative 및 Hessian reset 소진으로 발생했다.
원 실패의 normalized equality residual은1.33670852165e-13으로,
변경하지 않은 내부 종료 기준10*ftol=1e-13보다 조금 컸다.
별도 seed의 동일 기전도 재현했다. IEEE_INVALID는 통과한 상태에도 나타났으며,
trap을 걸어 vendor의 사용하지 않는 bound multiplier를 위한 의도적인0/0에서
발생함을 확인했다. 이를 물리 상태의 장 실패로 해석하지 않았다.

반환된 유한 profile에서 새 최소화를 실제 수행하자 mode0으로 수렴했다.
실패를 PASS로 이름만 바꾸거나 ftol·KKT·field 허용오차를 넓히지 않았다.
생산 solve_tf_family_native에는 이미 실패 보존과 이런 warm restart가 있었고,
즉시 중단하던 연구용 공존선 driver가 이를 이용하지 못했던 것이었다.

새 N256의 같은 유한 공존 후보는
n_low=0.072763037238587971, n_high=0.074330203388966634 fm^-3,
mu=952.01186827331685MeV, P=0.23426065228209580MeV fm^-3다.
h/2 검산의 Δμ=3.21165317e-7MeV, ΔP=2.39340778e-8MeV fm^-3는 기존 목표 이내다.
복구 후 다섯 밀도에서 slab/tube 각 하나, 총10번의 실제 최소화가 모두 mode0으로 끝났고,
세 내부 점에서 후보 에너지는 chord보다 높았다. 최소 차이는
2.21068928852e-6MeV fm^-3였다.
이는 조사한 두 branch의 유한 공존 근거이며, 전역 최소나 연속 mesh 극한의 증명이 아니다.
전체 profile·실패·restart·readback은 tf/ 보고서와 data/tf/에 있다.

## Matching 이론과 실제 수정

T=0의 한 보존 baryon charge barotrope에 대해
F=epsilon/n, mu=F+P/n, P=n^2 dF/dn, dP=n dmu를 사용했다.
표의 Euler 관계가 맞는 것만으로 보간 구간의 미분형 first law가 보장되지는 않는다.
네 smooth bridge는 선택한 구간의 유효 barotrope이며,
두 원 microscopic phase의 평형 교차라는 주장은 하지 않는다.
CSS는 crust 측 density jump가 있는 별도 구성이고,
sharp Maxwell은 실제 두 P(mu) branch의 동일 압력·화학퍼텐셜 교차다.

| 재현한 문제 | 최소 수정 및 검증 |
|---|---|
| 허용된 source Q3 잔차가 일부 방법의 endpoint·결과를 바꿈 | 모든 bridge가 같은 Euler mu=F+P/n를 사용. source Q3 gate는 유지 |
| C1 spline을 한 번의 Gauss 적분으로 지나 내부 first law 불일치 | 실제 knot별 동일 quadrature로 적분; 미분 관계 반례로 검증 |
| 401노드 사이의 비인과적 구간이 통과 | 각 다항식 구간의 stationary roots와 endpoint에서 positivity·안정성·인과성 검사 |
| 안정한 감소/상수 P/n를 거부하거나0으로 나눔 | 실제 물리 조건으로 canonical polynomial을 검사하고0 나눗셈 제거 |
| 차가운 입력 열·operation/error receipt·phase 수·stellar pair report 변조를 놓침 | 실제 열/metadata와 serialized receipt를 다시 결합해 검증 |
| SUPPORT의 미승인 경계와 변경·삭제를 무시 | API와 pipeline에 선언된 source 지원범위와 봉인 검증 적용 |
| 같은 물리 EOS를 다른 header 질량 기준으로 쓰면 Maxwell 거부 | 각 source의 물리 P,epsilon,mu를 유지하며 공통 output 기준으로 재인코딩 |

감소하는 q=P/n의 직접 반례는 x=ln[n/(0.04fm^-3)],
q=(10-x)MeV, F=(940+10x-x^2/2)MeV다.
dq/dx<0이지만 s=c_s^2/c^2=(q-1MeV)/(F+q)는 해당 구간에서0<s<1이다.
Wolfram 검산, native 수정 전 실패와 수정 후 통과를 함께 보존했다.

다항식 극값 검사는 부동소수점 근 분리와 평가다. 엄밀한 interval arithmetic 증명으로
표현하지 않는다. fixed diagnostic grid에서 놓치던 구간을 검사하는 구현상의 개선이다.
Maxwell의 세 내부 교차점, 역방향 density jump, 공통 mu support 부재도 별도 반례로 검증했다.

Header 기준 M에서 M'로 바꿀 때 Q7'=(M/M')(1+Q7)-1로 쓰는 것은
epsilon=nM(1+Q7)를 보존한다. 반대로 F→F+C는 epsilon→epsilon+nC, mu→mu+C라는
물리적 변경이다. 이번 수정은 후자를 수행하지 않는다.
동등한 여섯 encoding의 출력 epsilon 차이는 최대2.84e-14MeV fm^-3였고,
물리적인+100MeV/핵자 offset을 준 경우는 계속 거부했다.
상세 유도·원전·counterexample은 theory/THEORY_AUDIT.md,
native_matching/REPORT.md 및 MAXWELL_REFERENCE.md에 있다.

## 실제 생성 EOS 연결과 남은 입력 범위

core5를 각각129점 요청으로 n≤0.5fm^-3까지 새 native CLI로 생성했다.
APR은 자체 phase endpoint가 추가되어131행이다. 이 상단은 이번 matching 실험 범위이며,
원 모델의 고밀도 인과성 문제를 해결한 값이 아니다.

조건부 SPG5×core5에 네 smooth 방법과 CSS 및 sharp Maxwell을 적용했다.
최초150개 시도 중 smooth/CSS100개가 construction·native verifier·readback을 통과했다.
APR의 전체131행은 원래 Gibbs phase 경계를 현재 matcher가 지울 수 있으므로 거부했다.
원자료 rows0–64, phase6, n=0.085–0.20439495489345874fm^-3를 그대로 선택한 APR은
SPG_M1과 다섯 방법 모두 통과했으며 retained F의 오차는0.0MeV였다.

이 제한을 물리적 no-transition으로 해석하지 않는다.
전체 APR과 다른 여러 phase source를 지원하려면 source boundary 보존, 허용 segment
interpolant, native export, CSS stellar intake를 함께 구현해야 한다.
phase label을 지우거나 boundary를 넘어 임의 보간하는 것으로 우회하지 않는다.
pipeline_matching/APR_SUPPORT_SCOPE.md가 재현 가능한 현재 범위를 제공한다.

초기 Maxwell 조합의 header 기준 제약은 실제 입력을 통해 드러났고 위와 같이 수정했다.
이어 pipeline이 문서와 달리 조성을 그대로 넘겨 거부되던 문제를 수정했다.
public API는 drop_composition=True를 명시한 때만 복사된 thermo-only 입력을 만들고,
원본·모든 물리 열·phase·extra quantities·지원구간을 유지한다. 기본 거부는 그대로다.
생산 pipeline도 그 명시적 정책과 변환 receipt를 보존한다.

수정 후25조합은 실제로 다시 실행했다. APR 전체 입력5개는 phase 보존 미지원으로,
나머지20개 중13개는 제공된 공통 mu 구간에서 교차가 없음으로,
6개는 공통 mu 지원구간 부재로,1개는 여러 교차의 현재 미지원으로 끝났다.
이 조합 중 sharp Maxwell로 승인된 것은0개다. 공통구간 부재와 구현 제한을
물리적인 no-transition으로 부르지 않는다. 분류 수정 전 raw 결과도 보존했고,
최종 분류는 integration/TAXONOMY_FINAL.json에 별도로 기록했다.

성공한100개와 APR 단일 branch5개, 총105개 출력은 마지막 native 실행 파일로
다시 읽어 consistency를 통과했다. 전체 frontend/pipeline 회귀1,140passed,9skipped,
generator142passed였고, 마지막 분류 수정은5개 red 반례를 포함한99개 집중 회귀로
검증했다. 추가5개 분류 반례와 native matching·Maxwell·SPG 반례를 원격 CI에도 포함했다.
실제 게시 commit의 CI 결과가 최종 전체 실행 결과다. 원자료는 integration/에 있다.

## 판정과 이어갈 작업

- 구현/수치적으로 수리한 항목: SPG 희박 극한 분류, support 무결성, TF 연구 driver의
  실제 재시작, matching의 endpoint·적분·구간 검사·잘못된 단조성 조건·receipt 계약,
  Maxwell encoding 의존성.
- 모델 제한: 기존 독립 분석에서 확인한 APR/KDE0v1/SLy4의 고밀도 비인과성을
  matching을 위해 임의 수리하지 않았다.
- source fidelity 미해결: SPG 원 calibration/변분 일치, BBB 저자 smooth-fit 미분 표현,
  BSk source-fit의 조성·열역학·원문/배포 코드 차이. 각 명시적 representation은 따로 유지한다.
- 수치/구간 미해결: lower/inner gap·overlap의 일반 결합, TF 미표본 연속 지원,
  여러 source phase를 보존하는 전체 EOS matching 및 downstream 연결.

재현·local 인계는 HANDOFF_KO.md를 따른다.
새 반례 회귀는 CI에도 추가했다. 실패 원자료와 최종 산출물의 hash를 함께 게시하며,
모든 모델과 모든 구간에 대한 무조건적인 PASS를 선언하지 않는다.
