# SPG 구형 CLD 조성 테이블 재생성

이 노드는 `c8ff2d6c97a93f1ccf621051861d3091bfe94d81`의 25개 후보 셀 계산을 M1–M5 각각의 독립적인 구형 inner-crust 조성 테이블 생성으로 확장한다. 표면 계수는 기존 후보 mass fit을 그대로 사용한다. 원 저자 fit과의 동등성, 전역 최솟값, 전체 형상 census를 주장하지 않는다. CompOSE와의 일치는 구현 통과 조건이 아니며 비교 오차를 별도로 보존한다.

## 문헌과 모델 계약

- [Scurto, Pais, Gulminelli 2024, arXiv:2402.15548v1](https://arxiv.org/html/2402.15548v1), Eqs. 21–25: 두 핵자 상의 bulk 에너지, 양성자 밀도 대비를 사용한 Coulomb 항, 구형 표면 항. Eq. 25의 모델별 표면 fit은 후보 계수와 원전 계수의 동일성을 자동으로 보장하지 않는다.
- [CompOSE SPG(M2) 제출 문서](https://compose.obspm.fr/download/1D/Crust/M2/eos.pdf)는 BSk22를 `n_B < 0.002 fm^-3`에 사용한다고 명시한다. 이 밀도는 실제 outer-crust neutron drip보다 높다. 이 노드의 CLD를 자체 neutron drip까지 연장한 부분은 명시적인 후보 확장이다. 원 제출물의 BSk22 저밀도 구간을 재생성했다고 해석하지 않는다.
- 원 논문은 같은 밀도에서 CLD와 균질물질의 에너지 교차를 비교한다. 이 값은 receipt에 보존한다. 추가로 구현한 국소 중성 두 상의 공통 접선은 같은 functional의 열역학적 지지구간을 정하는 별도 구성이다. 출력은 cluster 쪽 끝점에서 종료하며 core 테이블과 접합하지 않는다.

`gdfm_uniform.py`의 세 meson bulk 및 rearrangement를 사용하고 핵자 공통 정지에너지 `M = 938.9 MeV`를 유지한다. 공간 좌표는 fm, 수밀도는 fm^-3, 에너지밀도와 압력은 MeV fm^-3이다. 전자는 유한 질량 자유 Fermi gas이며 muon 삽입 비용을 검사한다. 이번 지지구간에서는 muon 삽입 불안정성이 발견되지 않았다. 물리적 neutron/proton 질량을 복원하면 beta 조성이 바뀌므로 기존 테이블에 상수만 더하여 BSk22와 연결할 수 없다.

## 유도와 구현

독립 변수는 `x=(f,R,n_i,y_i,y_g)`이며

```
n_g = (n_B - f*n_i)/(1-f)
n_e = f*n_i*y_i + (1-f)*n_g*y_g
Delta = n_i*y_i - n_g*y_g
E = f*E_i + (1-f)*E_g + E_e + C + S
C = 2*pi*e2*f*Phi*R**2*Delta**2
S = 3*f*sigma(y_i)/R
```

여기서 `e2=1.43996448 MeV fm`이다. 반경은 `R**3 = 3*sigma/(4*pi*e2*Phi*Delta**2)`로 제거한다. 따라서 각 평가에서 `S=2*C`를 만족한다. 최적화 변수는 `t=f*n_i/n_B`, `n_i`, `y_i`, `v=y_g/y_i`이다. 구현에는 원 좌표의 다섯 해석적 도함수와 축약된 네 도함수가 함께 있다. 유한차분 검증은 에너지 함수에서 독립적으로 수행한다.

정지 조건에서만 압력 공식을 간소화하면 neutron drip과 유한 residual에서 오류가 생긴다. 구현은 다음 envelope 표현을 그대로 사용한다.

```
B = mu_p_g - mu_n_g + mu_e - 2*C/(Delta*(1-f))
partial_E_partial_nB_at_fixed_x = mu_n_g + y_g*B
mu_B = partial_E_partial_nB_at_fixed_x + t*d(E/n_B)/dt
P = n_B*mu_B - E
```

`t=1`에서는 gas가 없고 KKT 조건은 `d(E/n_B)/dt <= 0`이다. 자체 neutron drip은 이 도함수가 0이 되는 경계로 계산한다. `y_g=0` 경계에서는 proton 삽입 비용이 음수가 아니어야 한다. 공통 정지질량 convention에서 gas 진공의 `mu_n=mu_p=M`을 사용한다.

L-BFGS-B로 후보를 찾은 뒤 해석적 gradient의 root를 least-squares로 polish한다. 반경 제거가 기존 M3 `ABNORMAL` 실패를 해소한다. 통과 gradient 기준은 이전 `1e-4`보다 엄격한 `1e-5`이며 실패 threshold를 완화하지 않았다. 최저 에너지 시도가 미해결이면 `HOLD`를 유지한다. 국소 Hessian, 계산 경계 활성 여부, muon/proton 삽입 조건을 별도로 검사한다. 고정 초기값 세 개와 밀도 continuation을 사용하며 대표 밀도 네 곳에서는 양성자 gas 초기값 세 개를 강제로 추가한다. 모든 초기 검색 결과와 polish 결과를 보존한다.

탐색은 `0<t<=1`, `n_i>n_B`, `0<y_i<0.5`, `0<=y_g<y_i`의 구형 cluster branch이다. 계산용 경계와 이 유한 초기값 탐색은 전역 최소나 모든 pasta 형상의 완전성을 증명하지 않는다. 구형 functional의 재생성 가능성과 원전 전체 상도 재현은 다른 주장이다.

## 반례와 수정

1. 기존 M3 `n_B=0.005011872336 fm^-3` 최저 시도는 optimizer `ABNORMAL` 때문에 `HOLD`였다. 반경 제거와 정지점 polish 이후 같은 bulk/surface 입력으로 엄격한 residual 조건을 만족한다. 과거 실패 row 전체와 수정 결과를 `legacy_m3_repair.json`에 보존한다.
2. 순수 neutron gas로만 제한하면 M1/M3의 전이 부근 proton 삽입 비용이 각각 약 `-0.111/-0.098 MeV`였다. 이는 모델 전체가 거절된 증거가 아니라 제한된 변분 공간의 결함이다. 원전 Coulomb 대비항에 따라 gas proton 자유도를 추가한 후 작은 양의 gas proton branch를 찾았다. 제한된 branch 결과는 `logs/transition_probe.log`에 남긴다.
3. 조성 export의 `Yh`는 cluster 수밀도 분율이 아니라 **cluster baryon fraction `t`**다. 핵 수분율 `t/A`는 receipt에 별도로 보존한다. `Yn+Yp+Yh=1`, `Yp+(Z/A)*Yh=Ye+Ymu`를 native readback에서 확인한다. 연속 CLD의 `A,Z`는 정수 핵종으로 오해하지 않는다.

## 재현

저장소 루트에서:

```sh
bash fortran_eos/build.sh build/compose_eos_fortran
PYTHONPATH=python_frontend python3 -m pytest -q python_frontend/tests/test_candidate_gdfm_cld.py python_frontend/tests/test_spg_crust_generator.py
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 scripts/generate_spg_crust.py build/spg_crust_reproduced --points 513
```

출력 디렉터리는 새 경로여야 한다. `--points`는 두 중첩 격자를 위해 `8k+1 >= 17`이어야 한다. 각 모델의 `compose/`는 실제 native writer로 생성하고 reader로 다시 읽은 CompOSE 파일이다. `all_attempts.json.gz`는 coarse/fine 모든 시도, drip/전이, 독립 proton probe를 포함한다. `receipt.json.gz`는 수치 검증, 입력/소스 SHA-256, 조성 의미론과 CompOSE 차이를 포함한다. 보존 경로는 `data/research_completion_20261116/spg/`다.

독립적인 표본 재최적화·국소 곡률·격자 세분화·native roundtrip은 구현 근거다. 차분 sound speed와 causal 구간은 표본 진단이며 연속 구간이나 임의 보간함수에 대한 정리가 아니다. 저 `Y_p` 표면 extrapolation, 후보 mass fit의 원 저자 fit 동등성, source의 BSk22 저밀도 segment, 다른 geometry는 각각 별도의 의무로 남는다.

## 보존한 최종 실행 결과

대상 테스트 14개가 통과했다. 각 모델에서 257개 coarse와 513개 fine 밀도점을 실제로 재최적화하여 총 2,565개의 native 조성 행을 생성하고 readback했다. 129개 진단 격자는 coarse의 매 두 번째 해를 사용한다. 아래 오차는 coarse PCHIP 값을 추가 fine 노드의 독립 해와 비교한 최대 상대 압력 차이다. 이는 측정된 interpolation 오차이며 연속 함수에 대한 엄밀한 오차 상계가 아니다.

| 모델 | 자체 drip [fm^-3] | cluster 끝점 [fm^-3] | P: 129→257 | P: 257→513 | 양성자 gas 행 |
|---|---:|---:|---:|---:|---:|
| SPG_M1 | 0.000194020641192 | 0.0744899608593 | 0.130718% | 0.0618844% | 23 |
| SPG_M2 | 0.000205001488779 | 0.0874137415362 | 0.129145% | 0.0600569% | 0 |
| SPG_M3 | 0.000193862106305 | 0.0789288768376 | 0.132456% | 0.0627911% | 20 |
| SPG_M4 | 0.000202108600111 | 0.0968089199567 | 0.134019% | 0.0627986% | 0 |
| SPG_M5 | 0.0002029589416 | 0.0715959858761 | 0.125462% | 0.0583437% | 0 |

모든 native C1–C6/C10/C11 검사가 통과했고 각 필드와 조성 presence mask가 readback 후 유지되었다. 최악 projected gradient는 `1.3291e-8`이며 허용치 `1e-5`보다 작다. 독립적인 양성자 gas 초기값 검사는 다섯 모델×네 밀도에서 수행했고, 각 점의 모든 여섯 시도를 보존했다.

CompOSE의 phase 7 중 겹치는 밀도에서 비교한 압력 차이는 약 20–22%, 전자분율 차이는 약 10–13% 수준이다. 이 차이는 후보 fit·bulk 및 모델 선택의 추가 연구 대상이며 구현 검사를 자동 거절하는 조건으로 사용하지 않았다. 정확한 모델별 수치는 `SUMMARY.json`과 receipt에 있다.
