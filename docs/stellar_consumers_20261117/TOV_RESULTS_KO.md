# D: exact-contract TOV background

2026-11-17. 기준 commit `539e6bd8de73f8654d121c0bc30023fb868dc05b`.
이 기록은 B의 reference-assisted hybrid whole-crust+SLy4 core baseline을 그대로
소비한 국소 D 계산이다. 새로운 core 생성, fitted input 변경, composition 보충을 수행하지 않았다.
37개 중심밀도 모형과 1.4 nominal solar mass 모형을 계산하고, 정확한 저밀도 law,
중심·표면·ODE 수렴 및 기존 native solver의 다른 보간법을 확인했다.

## 입력과 구현

- SSOT: `data/core_support_20261117/primary_v1/REPRESENTATION.json`.
  SHA-256 `591ac27aa99448512e291684c98f8273efc729d17c63be88e331c20e83bd2653`.
- 물리 CSV SHA-256 `48661daebaaba64277d96ba5ddd5ed16cf4639915010c9fa4ea94f61c7bdc79d`.
- 지원: n=0–1 fm^-3, P=0–506.42133443967515 MeV fm^-3.
  1 fm^-3는 실제 생성 끝이며 maximum-mass endpoint라는 뜻이 없다.
- 구현: `python_frontend/compose_eos/stellar_background.py`.
  `WholeCrustBarotrope`의 승인된 surface, outer power laws, lower P(mu) bridge,
  upper piecewise P(mu) laws를 enthalpy 좌표로 직접 평가한다.
  표본 사이를 새 spline으로 대체하지 않았다.
- 선택 closure: cold effective beta-equilibrium barotrope. B의 missing composition masks는 보존한다.
  이 모듈은 smooth baseline 전용이며 sharp interface contract를 명시적으로 거절한다.
- solver는 SciPy DOP853의 **Python reference backend**다.
  기존 Fortran production backend와의 별도 비교를 아래에 기록했다.

## 정의, 유도와 경계조건

Metric signature는 (-,+,+,+), metric은
ds²=-exp(2 Phi)c²dt²+(1-2m/r)^(-1)dr²+r²dOmega²다.
기하학적 mass m=GM/c²와 areal radius r는 km다.
epsilon과 P의 내부 기하학적 값은 G/c^4를 곱하고 m^-2에서 km^-2로 바꾼다.

c=299792458 m/s, G=6.67430e-11 SI,
1 MeV fm^-3=1.602176634e32 Pa,
변환계수는 1.3238333135663825e-6 km^-2/(MeV fm^-3)다.
질량 단위는 native와 같은 IAU nominal GM_sun=1.3271244e20 m³/s²,
M_sun=GM_sun/G=1.988409870698051e30 kg다.

Cold first law d epsilon=mu dn, P=n mu-epsilon으로 dP=n dmu이고,
따라서 h=int_0^P dP/(epsilon+P)=ln(mu/a),
a=929.7530114783499 MeV다. Surface polytrope는 expm1(h)로 안정적으로 역산하고,
outer power-law interval은 mu의 analytic inverse, upper/bridge는 원래 P(mu) 다항식을 평가한다.
dP/dh=epsilon+P의 차원과 first law를 회귀에서 검사했다.

기하학적 변수의 TOV 식은

```text
dm/dr = 4 pi r² epsilon
dh/dr = -(m+4 pi r³ P)/(r(r-2m))
dr/dh = -r(r-2m)/(m+4 pi r³ P)
dm/dh = 4 pi r² epsilon dr/dh
dPhi/dh = -1
```

중심의 regular expansion은 m(r0)=4 pi epsilon_c r0³/3,
h_c-h(r0)=2 pi(epsilon_c/3+P_c)r0²다.
생략항은 각각 O(r0^5), O(r0^4)다. 기준 r0=1 m이며 0.3 m, 10 m와 비교했다.
정확한 analytic surface h=0, P=epsilon=n=0까지 적분한다.
Phi(h)=0.5 ln(1-2M/R)-h이며 표면 Schwarzschild 시간과 일치한다.
양의 cutoff h_s를 요청한 비교에서는 Phi(h)=0.5 ln(1-2M_s/R_s)-h+h_s다.
이는 finite-cutoff Schwarzschild normalization이라는 별도 경계조건이다.

## 실제 관측량과 지원 제한

| 항목 | 계산값 |
|---|---:|
| canonical 중심밀도 | 0.5377339776363226 fm^-3 |
| canonical 중심압력 | 84.65082808258397 MeV fm^-3 |
| canonical h_c | 0.24177274388331796 |
| canonical M | 1.399999994656446 nominal M_sun |
| canonical R | 11.756410784055323 km |
| canonical GM/(Rc²) | 0.17584236238014914 |
| n_c=1 endpoint M | 2.0183528731300773 nominal M_sun |
| n_c=1 endpoint R | 10.536944116358379 km |

위의 많은 자릿수는 재현용 raw output이다. 아래 실측 수렴보다 높은 정확도를 뜻하지 않는다.
논문·요약에서는 canonical R≈11.75641 km, endpoint M≈2.01835 M_sun 정도로 읽는다.
n_c=0.1–1을 0.025 fm^-3 간격으로 계산했다. 표본 간 mass slope는 모두 양수이고,
마지막 구간 0.975→1의 slope는 0.3880442555501237 M_sun/(fm^-3)다.
**Mmax는 현재 지원에서 bracket되지 않았다.** 이를 정하려면 기존 EOS의 유효성과 causality를
검사하면서 별도의 upper-core 생성 확장이 필요하다. 이 D의 canonical 관측량에는 확장이 필요하지 않아
수행하지 않았다. 표본의 mass 증가만으로 전체 branch의 radial-mode stability를 증명하지 않는다.

## 수치 검증과 오차의 분리

기준 controls: rtol=2e-10, atol=1e-12 (기하학적 r,m 각각 km),
max_step_h=0.002, r0=1 m. 실제 `convergence.csv`, `endpoint_convergence.csv`를 보존했다.

| 비교 | ΔM (nominal M_sun) | ΔR (m) |
|---|---:|---:|
| canonical, rtol=3e-11 / step=.001 | +4.57292e-8 | -6.30914e-5 |
| canonical, r0=10 m | -1.28274e-8 | -7.93904e-6 |
| canonical, r0=0.3 m | -8.77559e-9 | -1.75884e-5 |
| canonical, h_s=1e-4 | -3.34466e-12 | -4.33245433 |
| canonical, h_s=1e-6 | -7.70495e-14 | -0.04334465 |
| canonical, h_s=1e-8 | -8.88178e-16 | -0.00043345 |
| n_c=.1, tighter ODE | -3.37768e-10 | -0.00032448 |
| n_c=.1, r0=10 m | -1.15719e-8 | +0.00467512 |
| n_c=1, tighter ODE | +6.89036e-7 | -0.00040024 |

rtol 변화의 결과는 piecewise-law 구간을 가로지르는 adaptive integration 때문에 단조롭지 않다.
local rtol 숫자를 최종 질량 상대오차라고 해석하지 않는다.
관측된 canonical ODE mass spread는 3e-7 M_sun 이하, endpoint의 tested spread는
7e-7 M_sun 이하이며, 수치 오차와 실제 EOS/model uncertainty는 별도다.
endpoint 체크의 h_s=1e-8 radius 효과는 저질량 n_c=.1에서 0.06747 m,
n_c=1에서 0.00016179 m다. 모든 기준 결과는 h_s=0를 사용한다.

구현 검증 5개는 원 pressure consumer와 state 일치, enthalpy first law,
vacuum와 지원 밖 거절, analytic incompressible Schwarzschild interior 해,
배경 미분방정식·surface lapse·radius inversion을 다룬다. Incompressible 모형은
ODE의 독립 analytic control이며 causal physical EOS로 채택한 결과가 아니다.
정확한 knot에서는 cs²의 서로 다른 one-sided 값이 존재하므로 thermodynamic state는
knot에서도 비교하고 derivative의 양 소비자 비교는 interval interior에서 수행한다.
5 passed, 실제 pytest exit=0이다.

## 기존 native Fortran과 실제 보간 민감도

원 contract의 모든 law join을 유지하여 514개 및 2053개 샘플을 만든 뒤,
기존 `case_driver real_table_star`를 실제 호출했다. 양쪽 모두 admission과 solve에 성공했다.
실제 출력에 기록된 interpolant는
`monotone_PCHIP_epsilon_minus_pressure_of_logP_with_firstlaw_density_v1`이다.
초기 예상의 P(mu) Hermite가 아니므로 이를 동일 EOS 표현이라고 표기하지 않았다.

P_c=1.3562557880266705e34 Pa,
P_cut=7.774846713458777e21 Pa (원 contract h_s=1e-4)를 동일하게 둔다.
Exact-contract reference의 cutoff radius는 11752.078329725505 m다.

| 샘플 수 | native R (m) | native−exact ΔR (m) | native−exact ΔM (M_sun) |
|---:|---:|---:|---:|
| 514 | 11749.923704798046 | -2.154624927458 | -1.965253126901e-7 |
| 2053 | 11751.844352381375 | -0.233977344129 | +9.469875106305e-9 |

서로 다른 적분 변수/solver를 통한 질량과 반경 cross-check이며 표본 세분화로
native radius 차이가 줄어드는 것을 확인했다. 잔여 0.234 m는 이 비교의 representation/ODE
오차 합이며, 전부 exact-contract solver의 오차라고 귀속하지 않는다.
기준 B 파일을 바꾸거나 q-linear/PAV repair를 하지 않았다.
두 native run의 exit=0, raw stdout/stderr, namelist, binary SHA-256와 명령은
`NATIVE_COMPARISON.json`, `native_sub*.nml`, `native_sub*.stdout.log`에 있다.

## 실행·실패 보존과 재현

작업 디렉터리는 `/home/cosmosapjw/ns_cowling_wt/stellar_consumers_20261117`다.

```bash
export PYTHONPATH=python_frontend:pipeline:/tmp/ns_cowling_audit_20261117_pydeps
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python -u scripts/run_stellar_background_20261117.py
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python data/stellar_consumers_20261117/tov/check_endpoints.py
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python -m pytest -q python_frontend/tests/test_stellar_background.py
```

첫 runner attempt는 CSV field schema에 surface_h가 빠진 구현 오류,
둘째는 scripts/compose_cowling.py가 pipeline package를 가리는 import 경로 오류로 종료됐다.
둘 다 원 traceback/exit를 `attempt1_*`, `attempt2_*`에 보존하고,
CSV schema 통일 및 명시적 pipeline 우선 경로로 수정했다.
첫 회귀의 knot cs² 비교 실패는 one-sided derivative 비교 정의 문제였으며
thermodynamic 연속량과 interior derivative 비교를 분리하여 수정했다.
원 실패 로그를 보존했다. 수정 후 전체 production runner, endpoint 검사, pytest의
실제 exit는 모두 0이며 결과를 실패 이력과 혼동하지 않는다.

정확한 입력·구현·runner·test hash와 Python/NumPy/SciPy 버전은 `PROVENANCE.json`,
배경 profile은 `canonical_profile.csv`, 전체 sequence는 `sequence.csv`다.
최종 생성 파일 목록 및 content hash는 `MANIFEST.json`이다.
이 lane은 local 산출물만 만들었으며 commit/push하지 않았다.
