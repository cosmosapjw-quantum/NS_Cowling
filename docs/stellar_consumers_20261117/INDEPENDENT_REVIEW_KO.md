# D 독립 검토

C commit `539e6bd8de73f8654d121c0bc30023fb868dc05b`에서 시작한 D의 APR 소비자, TOV 배경, 평형 완전유체 Cowling 구현·산출물에 대해 한 번의 fresh-context 검토를 수행했다. 검토 범위에 열린 blocking finding은 없다. 실제 발견 F1은 소스에 근거한 최소 metadata 수정과 selective readback으로 종료했다. 재귀 reviewer나 A–C campaign 재실행은 없었다.

## F1: native Cowling 보간기의 오기 — CLOSED

Native runner와 `cowling/NATIVE_COMPARISON.json`은 provider를 Hermite P(mu)로 적었다. 그러나 `_compile_solver`가 빌드하는 `oscillation_source/solutions/src/eos/ns_table_eos.f90`는 g=epsilon−P의 ln(P) PCHIP를 구성하고 first law로 n, mu를 복원한다. `solutions/app/cell_solve.f90`는 이 모델을 배경 적분과 mode matching에 전달한다.

Native TOV와 native Cowling은 같은 EOS provider를 공유한다. executable·적분 경로가 다르며 exact-contract Python과 비교한다는 설명은 유효하다. 별도의 Hermite provider라는 설명은 잘못이었다. 담당자가 runner와 receipt의 label을 수정하고 실제 provider SHA-256 `466d9a548f6b0282b0f9d938d41f1069d5e726d80e6de450de06ac95f518682c`를 붙였다. 원 실행 로그의 오기는 correction record로 보존했다. reviewer가 수정된 receipt와 source hash를 확인했다. 수치 결과는 변하지 않았으므로 계산 재실행은 요구하지 않았다.

## 직접 확인한 수학·물리

서명 (−,+,+,+), rest-energy 포함 epsilon, geometric km 변환과 nominal solar GM가 일관된다. TOV enthalpy 변수변환, 유한 중심 초기조건, 해당 surface cutoff의 Schwarzschild 시간 정규화가 구현에 반영된다. canonical profile에서 mu exp(Phi)의 상대 spread를 독립 재계산한 값은 4.556e−16이다. 이는 표현·시간 gauge의 내부 일관성 검사이며 EOS 정확도의 독립 증명은 아니다.

Cowling scaled 식과 원래 W,V 식이 일치한다. 정칙 중심은 V=−W/(ell r)이며 두 baseline profile의 중심 비율은 1.0이다. 표면 Q=0, 두 방향 적분의 normalized determinant, join mismatch, root bracket, node 수를 함께 사용한다. 선택한 Gamma1=Gammaeq에서 N²=0이며 frozen-species 부력을 추정하지 않는다. 연속 Gibbs 경계에서 W,V를 연속으로 전달하고 one-sided susceptibility를 유지하는 해석이 일관된다.

APR의 각 cubic epsilon(n) interval에서 epsilon의 2차 도함수는 선형이다. 양 끝의 양의 2차 도함수는 내부 압축률 양성을 보장한다. mu−n epsilon''의 도함수는 −n epsilon'''이므로 n>0에서 극값은 endpoint에 있어 causal inequality의 endpoint 검사가 유효하다. 이는 binary64 다항식 검사이며 엄밀한 outward-rounded interval 증명으로 격상하지 않는다.

## 실제 검증과 identity

독립 focused pytest 실행은 remote PID 1964180에서 exit 0, **19 passed in 1.19 s**였다. 정확한 Schwarzschild interior TOV control, 약중력 incompressible Cowling control, 식 변환, EOS domain, 경계·보간·mask tests를 포함한다. 정확한 실행 명령은 동반 JSON에 있다.

55개 APR/crust source 파일을 RUN_CONTEXT의 SHA-256과 다시 비교하여 모두 일치함을 확인했다. B primary REPRESENTATION은 C base commit의 bytes와 같다. 25개 APR assembled table을 실제 readback하여 retained core의 n,P,epsilon,mu,F 최대 절대차 0, composition masks와 phase 배열 exact retention, bridge의 unavailable masks를 확인했다. 물리량 동일성과 파일 byte identity는 따로 판정했다.

TOV sequence는 실제 37행, nc=0.1–1.0 fm⁻³다. 질량은 모든 인접 표본에서 증가하고 끝값은 2.0183528731300773 Msun이다. 이 결과는 최대질량을 bracket하지 않는다. 두 native lane의 입력 table·executable hash와 exit 0 receipt를 직접 대조했다. APR matching 검증과 B hybrid 별의 관측량을 혼합하지 않는다.

## Cowling 수렴과 최종 회귀의 구분

최종 `cowling/SUMMARY.json`의 **12개 설정**을 읽고 24개 mode profile을 독립 재계산했다. 모든 f profile의 W/V node 수는 0/0, p1은 1/1이다. `ode_loose`를 제외한 baseline 대비 최대 절대 주파수 변화는 f 0.000537006 Hz, p1 0.002611512 Hz다. 이는 시험한 설정의 sensitivity이며 엄밀한 오차 상계는 아니다.

`ode_loose` f-mode는 |matching residual|=3.8695e−6으로 고정 criterion 2e−6을 넘고 `residual_accepted=false`로 보존된다. 느슨한 진단을 PASS로 세거나 전체 12개 설정이 모두 승인됐다고 쓰면 안 된다. 이 제외는 사후 tolerance 완화가 아니다.

Native sub1→sub4의 exact finite-cutoff 대비 차이는 f 0.001013735→0.000205347 Hz, p1 0.084873189→0.007416953 Hz로 줄어든다. 동일 central pressure·finite surface cutoff 비교이며 원 vacuum 시간 gauge와 혼합하지 않는다. Native provider가 두 lane에서 공유된다는 한계는 F1 교정에 반영했다.

Parent가 수행한 최종 frontend/pipeline 회귀는 별도 증거다. `FRONTEND_PIPELINE.json`의 exit 0, source/tests unchanged와 실제 log의 **1397 passed, 9 skipped in 117.40s**를 확인했다. reviewer는 전체 suite를 다시 돌리지 않았다. 9 skips는 staging되지 않은 기존 raw fixture의 제한이다.

12개 설정의 진행 중 상태는 최종 receipt로 해소됐다. 새로운 reviewer cycle이나 추가 전체 campaign은 필요하지 않다. 검토자는 이 문서와 동반 JSON만 생성했으며 코드 수정·commit·push는 하지 않았다. Mmax, APR whole-star observable, frozen-species 부력, elasticity/superfluid/full-GR 효과에 대한 새 결론은 없다.
