# E 독립 검토: native SLy4 지원 범위, TOV 최대점, 선택 Cowling 모드

## 판정

**ACCEPT_WITH_SCOPE. 열린 BLOCK_NOW finding은 0개다.** D commit `300e5a0532dd137cb09738568fca0f35cfcc9282`에서 재사용한 물리·소비자 구현과 E의 실제 생성 결과를 대상으로 한 차례 독립 C9 검토를 마쳤다. EOS/TOV와 최종 Cowling 결과, 추가 tighter-background control까지 읽고 아래 대상 검산을 실제 실행했다. 이 판정은 아래에 고정한 수치·물리 범위에 대한 것으로, 전체 미시적 인과성이나 radial stability spectrum을 추가로 확정하지 않는다.

검토자는 새 하위 reviewer를 만들지 않았고, 기존 캠페인을 반복하지 않았다. 작성한 파일은 이 문서와 `data/stellar_extension_20261117/INDEPENDENT_REVIEW.json`, `INDEPENDENT_REVIEW.log`뿐이다. 소스·생성 결과·물리 수용 기준을 수정하거나 commit/push하지 않았다. 게시·통합은 부모 작업의 책임으로 남긴다.

## 검토 대상과 실제 실행

`extend_stellar_support_20261117.py`, `run_stellar_extension_modes_20261117.py`, `plot_stellar_extension_20261117.py`, E artifact tests 및 재사용된 core extension/TOV/Cowling 구현을 읽었다. native analytic Hessian의 실제 lepton susceptibility, thermodynamic derivative와 cold Euler identity를 직접 대조했다. 최종 수치는 `eos/EOS_RESULT.json`, `eos/RESPONSE_SCOPE.json`, `tov/MAXIMUM_REFINEMENT.json`, `cowling/SUMMARY.json`, `cowling/TIGHT_BACKGROUND_COMPARISON.json`에 결속했다.

| 검증 | 실제 결과 | 근거 상태 |
|---|---|---|
| E artifact pytest | `7 passed in 4.52s`, exit 0 | implementation-verified |
| D source/input 바이트 대조 | 50개 파일이 지정 D commit 및 기록 hash와 일치 | BYTE_IDENTITY |
| 실제 native generator/linked objects | 4개 hash 일치 | BYTE_IDENTITY |
| Cowling 재사용 D 경로 | 7개 파일, source/test/control receipt hash 일치 | BYTE_IDENTITY |
| native 직접 재호출 | n=1, 대표 최대점, chi_eq=1 경계의 3점; 모두 exit 0 | numerically checked |
| 대표 TOV 재실행 | strict 설정 1회, 기록 대비 M·R 차이 0.0 | NUMERICAL_IDENTITY |
| TOV 원자료 산술 | 6개 fit, 4개 derivative bracket, 6개 convergence row 재계산 일치 | numerically checked |
| Cowling raw profile | 22개 CSV, 각 1199행, hash·node·정칙성·BC 확인 | numerically checked |
| 부모 전체 suite | 보존 log/receipt를 읽음: 1404 passed, 9 skipped, exit 0 | reused execution evidence |

전체 suite는 검토자가 다시 실행하지 않았다. 9개 skip은 staged되지 않은 기존 raw CompOSE 입력 관련 항목이며, 이번 E artifact tests의 skip은 없다. 실제 대상 검증 출력은 `INDEPENDENT_REVIEW.log`에 보존했다. 이 검토는 현재 읽은 소스·실행 기록과 해당 D 재사용 경계에 대한 것이며, 도구와 과거 감사의 계보를 재귀적으로 검토하지 않았다.

## native 연장과 정확한 prefix 보존

세 export는 같은 native SLy4 generator를 실제 실행한 65/129/257개 상태에서 각각 n=1 anchor를 제외한 64/128/256개 새 상태를 붙인다. assembled tail의 n, P, epsilon을 native 생성 표와 직접 비교해 상대 2e-14 이내 일치를 확인했다. 이는 단지 support metadata를 늘린 결과가 아니다. 기존 Fortran source·에너지 기준을 바꾸거나 인공 causal completion을 추가하지 않았다.

기존 773개 행의 모든 row field와 composition mask, surface/outer/lower bridge 및 기존 upper state 배열은 보존된다. artifact tests에서 원래 구간의 n, P, epsilon, mu, h를 array equality로 검증했다. 기존 내부 구간의 chi도 같다. n=1의 원래 끝점에는 기존 마지막 구간의 왼쪽 chi와 새 첫 구간의 오른쪽 chi가 구분된다. 예를 들어 fine에서 왼쪽 값은 0.8318776264232889, 오른쪽 값은 0.8318620397476921이다. 이 일방 미분의 차이는 명시된 Hermite knot 효과이며, thermodynamic state의 변경이나 상전이로 해석하지 않는다.

3개 contract의 모든 upper interval을 직접 validate했다. coarse/medium/fine 각각 364/428/556개 구간이 통과했고 연속 최대 chi는 각각 0.9999918438563513, 0.9999902103496859, 0.999989800271408이다. 검사는 knot에서의 양수성만 확인한 것이 아니라 기존 polynomial interval validator의 구간 내 극값 조건을 사용한다. support 바깥 입력을 거절하는 artifact test도 통과했다.

## 응답 정의와 causal 경계

밀도 n의 단위는 fm^-3, epsilon과 P는 MeV fm^-3, mu는 MeV다. chi=(c_s/c)^2는 무차원이다. lepton 분배를 최소화한 epsilon(n,y), y=Y_p에 대해 epsilon_y=0, H_yy>0이면

\[
\frac{dy}{dn}=-\frac{H_{ny}}{H_{yy}},\qquad
\chi_{\rm eq}=\frac{n}{\mu_B}\left(H_{nn}-\frac{H_{ny}^{2}}{H_{yy}}\right),\qquad
\chi_{Y_p}-\chi_{\rm eq}=\frac{nH_{ny}^{2}}{\mu_B H_{yy}}\geq0.
\]

Fortran analytic routine의 `sus`는 공통 전자·뮤온 화학퍼텐셜에 대한 두 lepton susceptibility의 합이다. 따라서 native `cs2_frozen`은 고정 Y_p·평형 lepton 분배의 응답이며, 모든 입자종의 개별 fraction을 동결한 응답이 아니다. `RESPONSE_SCOPE.json`의 명시적 해석이 구현과 일치한다.

| 대상 | 확인값 |
|---|---:|
| native equilibrium chi_eq=1 경계 | 1.2650056715193052 fm^-3 |
| 안전한 연속 export 상한 | 1.2649856715193053 fm^-3 |
| fixed-Y_p/lepton-equilibrated 경계 | 1.2397422126050377 fm^-3 |
| 대표 최대점의 chi_eq | 0.9680741481811074 |
| 대표 최대점의 chi_Yp | 0.9799407903603997 |

세 native 재호출에서 Euler identity와 analytic/FD chi를 대조했다. 경계 재호출은 chi_eq=1.0000000000000007, FD=0.9999999999831508이었다. 보존된 다섯 FD 간격의 경계 차이는 최대 8.50e-12 fm^-3다. 이는 native branch의 유한 scan과 국소 root의 수치 근거다. 전역 단조성 정리, 전체 frozen-species causality 또는 다른 고밀도 상의 배제는 이 결과에서 나오지 않는다.

## TOV 최대점과 초기 실패의 처리

기존 D의 exact-contract enthalpy DOP853 구현을 재사용한다. metric signature는 (-,+,+,+), 기하학적 길이는 km, m=GM/c^2이며 h=log(mu/mu_vac)이다. 중심 정칙 전개와 h=0 진공 표면은 그대로 유지한다.

초기 optimizer의 n_c=1.2060130400651856 fm^-3와 원래 작은 one-sided derivative의 약 1e-6 M_sun excursion은 보존되어 있다. `TOV_RESULT_INITIAL.json` 및 `ENDPOINT_DERIVATIVE.json`을 최종 정밀 최대점의 근거로 오인하지 않도록 `TOV_RESULT.json`은 authoritative refinement를 가리킨다. initial runner hash가 `RUNNER_INITIAL_SNAPSHOT.py`와 일치하며, refinement snapshot과 metadata 정정 receipt도 분리되어 있다. 초기 결과의 numerical defect를 삭제하거나 물리적 특징으로 해석하지 않은 처리는 적절하다.

최종 수치 근거는 다음과 같다.

- 직접 확인된 부호 bracket은 1.2060 < n_c,max < 1.2061 fm^-3다.
- 세 source grid 및 fine step-limit 검사에서 왼쪽 중심차분은 +(3.72–3.85)e-5 M_sun fm^3, 오른쪽은 -(7.47–7.59)e-5 M_sun fm^3다. 원자료 37개 call의 질량에서 4개 bracket과 6개 quadratic fit을 재계산해 저장값과 일치함을 확인했다.
- 대표 입력 n_c=1.2060328050005111 fm^-3, M=2.0503812906965106 M_sun,nom, R=10.006409438136822 km를 strict 설정으로 한 번 재실행했다. 설정은 r0=0.001 km, rtol=1e-12, atol=1e-14, max_step_h=0.0001, surface_h=0이고, mass·radius가 기록값과 동일했다.
- 저장된 source/step/center controls의 최대 관측 mass 차이는 2.0658848e-9 M_sun이다. 대표점부터 bracket 가장자리까지의 국소 곡률 estimate는 2.5847980e-9 M_sun이다.

이 숫자는 명시된 equilibrium barotrope와 검사한 수치 controls에 대한 것이다. 많은 끝자리는 재현할 대표 입력과 출력의 식별에 쓰고, 중심밀도는 실제 bracket으로 보고해야 한다. 반지름의 대표값과 fixed-n 수치 spread는 최대점 위치 또는 crust/matching 선택의 총 불확실성을 대신하지 않는다. complete radial stability spectrum은 계산되지 않았다.

## 최종 Cowling 결과와 추가 배경 검증

두 별의 중심밀도는 1.1818927792638818 및 1.2060130400651856 fm^-3로 고정되어 있다. 두 번째 별은 최대점 bracket 안의 선택 모델이며, 위 TOV 대표 모델의 중심밀도와 같다고 주장하지 않는다. closure는 Gamma_1=Gamma_eq, ell=2의 평형 perfect-fluid Cowling이고, 주파수는 무한원 관측자 기준이다. 기존 중심 정칙 조건, 연속 EOS knot의 W,V 연속성, 양의 h_cut에서 Delta P=0 조건을 그대로 사용한다.

| 선택 모델 | baseline f [Hz] | baseline p1 [Hz] | 검사한 최대 절대 shift f / p1 [Hz] |
|---|---:|---:|---:|
| near_endpoint | 2886.2159749082025 | 7493.943366276016 | 0.000299304 / 0.000808912 |
| endpoint | 2899.489457960419 | 7504.249074578502 | 0.001200881 / 0.001117290 |

각 모델의 baseline, mode ODE tolerance, surface cutoff, center radius, source-grid control을 포함한 10개 구성의 20개 모드가 모두 수용되었고 rejected attempt는 0개다. 표의 shifts는 JSON에 명시된 각 reference에 대한 관측 sensitivity이며, 일반적인 물리 오차 한계가 아니다. 특히 coarse source-grid 비교는 fine의 tight mode 설정을 reference로 사용한다.

TOV의 초기 noise가 mode 관측량에 미치는 영향을 확인하는 추가 검사는 두 번째 별의 중심밀도를 그대로 유지하고 배경 ODE만 strict하게 계산했다. 결과는 f=2899.4886731199385 Hz, p1=7504.249127615556 Hz이고, 동일한 tight mode 설정의 원래 배경에 비해 각각 +3.83857596e-5 Hz, -1.65482778e-4 Hz 이동했다. 실제 실행 receipt는 exit 0이며 source·contract·reference hash가 일치한다.

총 22개 raw CSV 각각에 대해 byte hash, 1199개 유한 행, 증가하는 r, 감소하는 h, 중심 정칙비 -ell*r*V/W, 세 cutoff 1e-10/1e-9/1e-8에서의 node 수를 검산했다. 모든 f에서 W,V node 수는 0, 모든 p1에서는 1이다. 최대 matching residual은 6.32570074e-7, join mismatch는 6.32570074e-7로 기존 2e-6 수용 기준 이내다. 최대 surface BC 상대 residual은 9.60585788e-17, 중심 정칙비의 최대 1과의 차이는 2.22e-16이다. 추가 배경 control의 두 profile도 같은 검사에 통과했다.

이 결과는 지정된 평형 Cowling f/p1에 한정된다. frozen-composition modes, g-modes, elasticity, superfluidity, full-GR mode 또는 radial 안정성을 확정하지 않는다.

## Findings와 종료

| Finding | 대상·실제 증거 | 영향·최소 수정 | 최종 상태 |
|---|---|---|---|
| E-C9-01 | `EOS_TOV_RESULTS_KO.md`의 atol/Delta LaTeX 역슬래시 및 범위 끝의 중복 역슬래시 | 수식 렌더링 문제; 해당 문단 표기만 교정, 숫자 유지 | NONBLOCKING, CLOSED: 부모 수정 후 관련 두 문단 readback |

초기 adaptive derivative noise, equilibrium/fixed-Y_p response 구분, 초기/정밀 최대점의 권위 분리는 이미 이루어진 수정을 이번 검토에서 확인한 항목이다. 이를 새 blocker로 재개하지 않는다. 열린 과학·구현 finding은 없다. 위 효과 범위에서 C9를 종료하며, 추가 검토 사이클이나 전체 캠페인 반복을 요구하지 않는다.

사용한 skill 문맥은 PhysMath Research Loop 1.3.0의 main, `POLICY_INTEGRATION.md`, `CODING_LOOP.md` C9, `SUBAGENT_MODE.md`다. 생산 결과와 독립 검토의 근거를 구분했고, typed identity·유한 검토·실패 보존·최소 수정을 적용했다. artifact별 SHA256과 실제 검증 상태는 동반 JSON에 있다.
