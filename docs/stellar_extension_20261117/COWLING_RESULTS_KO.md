# E Cowling: causal 최대질량 근방의 평형 f/p1

B primary의 reference-assisted hybrid를 유지하고 native SLy4 고밀도 꼬리만 확장한 E fine 계약에서, 두 고정 중심밀도의 ell=2 평형 완전유체 Cowling 모드를 계산했다. 표시 주파수는 근접 별 f=2886.22 Hz, p1=7493.94 Hz, 최대질량 밀도 bracket 안의 선택 별 f=2899.49 Hz, p1=7504.25 Hz다. 모두 무한대 시간으로 정규화된 주파수다. 실제 11개 설정, 22개 모드가 고정 matching/join 기준 2e-6을 통과했고 거절된 모드 시도는 0개다.

## 입력과 별의 선택

| 선택 | n_c [fm^-3] | M [M_sun] | R [km] | f [Hz] | p1 [Hz] |
|---|---:|---:|---:|---:|---:|
| 근접 별 | 1.181892779263882 | 2.05003895160 | 10.0645678753 | 2886.22 | 7493.94 |
| 최대질량 밀도 bracket 안의 선택 별 | 1.206013040065186 | 2.05038127757 | 10.0064566024 | 2899.49 | 7504.25 |

두 번째 중심밀도는 초기 TOV 최적화에서 선택한 입력을 그대로 고정한 것이다. 이후 [MAXIMUM_REFINEMENT.json](../../data/stellar_extension_20261117/tov/MAXIMUM_REFINEMENT.json)이 세 tail grid와 더 엄격한 배경 적분에서 확인한 최대질량 밀도 bracket [1.2060,1.2061] fm^-3 안에 있다. 최종 대표 최적화 밀도 1.2060328050005111 fm^-3의 모드로 재지정하지 않는다. 첫 번째 별은 고정한 두 번째 밀도의 0.98배다. 원 선택 파일의 identity와 당시 이유는 [SELECTION_USED.json](../../data/stellar_extension_20261117/cowling/SELECTION_USED.json)에 보존했다.

두 입력은 native 평형 음속의 chi=1 crossing 1.2650056715193052 fm^-3 및 export support 1.2649856715193053 fm^-3보다 작다. 최대질량 수열의 정량 판정은 TOV 자료에 근거하며, 양의 ell=2 f/p1만으로 radial 안정성이나 전체 mode 안정성을 판정하지 않는다.

Fine 계약 SHA-256은 `f1bd5c3c41b45a07a7299e665594173473ac790af6f519bf53085db843bd936d`, coarse 계약은 `993e9e669f43e96e70bca55878c0473b1a51fe8a98c951e8fa572d38a08486b5`다. 각각 native tail 256/64 구간이며, 원래 D의 저밀도 법칙을 포함한 source binding은 EOS lane을 따른다. 이 모델을 native unified SLy4 crust로 재해석하지 않는다.

## 방정식·중심·표면·시간 정규화

방정식과 정칙성 유도는 [D FORMULATION_KO.md](../stellar_consumers_20261117/FORMULATION_KO.md)의 변경 없는 식을 사용한다. Signature는 (−,+,+,+), 길이는 km, `k=2*pi*f/c`, `c=299792.458 km/s`, `Phi`는 시간 metric 지수의 절반이다. 평형 응답 `Gamma1=Gamma_eq`를 선택한다. EOS 연속 knot에서 W,V를 연속으로 연결하며 sharp Maxwell jump는 입력에 없다.

중심은 `W=r^(ell+1)w`, `V=r^ell v`, `(w,v)=(1,-1/ell)`인 정칙 series로 r0=1 m에서 시작했다. r0=0.1 m 설정으로 중심 절단 민감도를 실제 비교했다. Background CSV의 r=0 행은 analytic center limit이며, 적분이 특이점 r=0에서 직접 시작했다는 뜻이 아니다.

TOV 배경은 analytic vacuum h=0까지 적분했다. 모드는 h=1e-9에서 Lagrangian pressure 조건 Q=0을 부과하고 내부 matching을 수행했으며 h=1e-7과 비교했다. 표면에서 부과한 작은 Q 잔차만을 고유해의 근거로 쓰지 않았다. 실수 주파수 bracket의 부호 변화, 내부 determinant, join mismatch, W/V node와 설정 민감도를 함께 확인했다. `Phi(h)=0.5*ln(1-2GM/(c^2 R))-h`인 vacuum 시간 정규화를 사용한다.

250–12000 Hz를 250 Hz 간격으로 실제 scan했다. 근접 별은 f=[2750,3000], p1=[7250,7500] Hz, 두 번째 별은 f=[2750,3000], p1=[7500,7750] Hz bracket에서 근을 구했다. 다음 acoustic root의 scan 정보는 보존했지만 이번 관측량으로 계산하지 않았다. 모든 f의 W/V node 수는 0/0, p1은 1/1이다.

## 실제 수치 민감도

Baseline mode 설정은 rtol=1e-8, atol=1e-11, h_cut=1e-9, h_match=0.15 h_c, max_step_h=0.002, root_xtol=1e-5 Hz다. 배경은 rtol=1e-10, atol=1e-12, max_step_h=0.005다. Root 좌표 xtol은 실제 주파수 오차 인증값이 아니다. 아래 많은 자릿수는 재현용 기록값이며, 선택한 모델의 물리적 정확도를 뜻하지 않는다.

| 별 | 설정 | f [Hz] | p1 [Hz] | baseline 대비 f [Hz] | baseline 대비 p1 [Hz] |
|---|---|---:|---:|---:|---:|
| near_endpoint | baseline | 2886.215974908 | 7493.943366276 | +0.000000000 | +0.000000000 |
| near_endpoint | ode_tight | 2886.216009046 | 7493.942557364 | +0.000034138 | -0.000808912 |
| near_endpoint | surface_h1e_7 | 2886.215866630 | 7493.944046294 | -0.000108278 | +0.000680018 |
| near_endpoint | center_0p1m | 2886.216274212 | 7493.943724834 | +0.000299304 | +0.000358558 |
| endpoint | baseline | 2899.489457960 | 7504.249074579 | +0.000000000 | +0.000000000 |
| endpoint | ode_tight | 2899.488634734 | 7504.249293098 | -0.000823226 | +0.000218520 |
| endpoint | surface_h1e_7 | 2899.488257080 | 7504.250191869 | -0.001200881 | +0.001117290 |
| endpoint | center_0p1m | 2899.488291754 | 7504.248103858 | -0.001166207 | -0.000970720 |
| endpoint | background_tight | 2899.488673120 | 7504.249127616 | -0.000784840 | +0.000053037 |

`ode_tight`는 mode rtol=1e-9, atol=1e-12다. `surface_h1e_7`과 `center_0p1m`은 각각 표면 cutoff와 중심 반지름만 바꾼다. TOV 최대질량 좌표 정밀화 중 발견된 배경 적분 잡음이 선택한 모드에 영향을 주는지 확인하기 위해, 두 번째 별에만 한 번의 `background_tight`를 추가했다. 이 설정은 배경 rtol=1e-12, atol=1e-14, max_step_h=1e-4이며 mode 설정은 `ode_tight`와 같다. 기존 `ode_tight` 대비 이동은 f=+0.000038386 Hz, p1=-0.000165483 Hz였다. 관측한 배경 민감도는 다른 설정의 주파수 spread보다 작다.

| 별 | fine 설정의 f 범위 [Hz] | fine 설정의 p1 범위 [Hz] |
|---|---:|---:|
| near_endpoint | 0.000407582 | 0.001488930 |
| endpoint | 0.001200881 | 0.002088010 |

동일 중심밀도와 동일 tight mode controls에서 coarse−fine 차이는 다음과 같다. 두 source-grid 결과 모두 같은 native EOS를 다른 tail 분할로 표현한 비교다.

| 별 | coarse−fine f [Hz] | coarse−fine p1 [Hz] |
|---|---:|---:|
| near_endpoint | +0.000020628 | +0.000006824 |
| endpoint | +0.000023470 | -0.000214190 |

이 source-grid 차이는 위 적분 설정 민감도보다 작다. 두 grid의 차이로 수렴 차수나 엄밀한 전역 오차 상계를 만들지 않는다. 0.01 Hz 표시는 관측한 수치 민감도를 고려한 반올림이며 EOS 모델·closure·Cowling 근사의 물리적 불확실성을 포함한 오차막대가 아니다.

## 실행·검증·산출물

주 campaign은 실제 exit 0, 587.837 s였고, 표적 tight-background 확인은 exit 0, 75.999 s였다. 원 runner는 D 소비자를 변경하지 않는다. D commit `300e5a0532dd137cb09738568fca0f35cfcc9282`의 두 소비자, analytic control tests, formulation과 validation/review 기록을 byte identity로 확인하고 재사용했다. D의 12설정 campaign이나 전체 회귀를 다시 실행하지 않았다.

Production readback은 exit 0이며 22개 eigenprofile의 각 1199행을 읽어 유한 값, 반지름 순서, 중심 비율, W/V node를 다시 계산했다. 두 baseline background는 analytic center와 vacuum surface를 포함하는 각 801행이다. CSV에서 재구성한 `mu*exp(Phi)` 상대 spread의 최댓값은 7.7484e-16였다. 전체 모드의 최대 absolute matching residual은 6.32570074e-07, 최대 relative join mismatch는 6.32570074e-07로 고정 2e-6보다 작다. 이것들은 주어진 수치 해와 provenance의 검증이며 별도 EOS 정확도 증명은 아니다.

[SUMMARY.json](../../data/stellar_extension_20261117/cowling/SUMMARY.json)은 주 10설정을, [TIGHT_BACKGROUND_COMPARISON.json](../../data/stellar_extension_20261117/cowling/TIGHT_BACKGROUND_COMPARISON.json)은 추가 1설정을 보존한다. [MODE_TABLE.csv](../../data/stellar_extension_20261117/cowling/MODE_TABLE.csv)는 11설정을 함께 모은 표이고, [VALIDATION.json](../../data/stellar_extension_20261117/cowling/VALIDATION.json)은 전체 22개 profile readback을 기록한다. 원 scan, 설정별 JSON, W/V eigenprofile, 배경 CSV, 실제 명령·환경·exit code, 로그는 같은 디렉터리에 있다.

재현 명령은 다음과 같다. 실제 실행 argv와 종료 상태는 EXECUTION.json, EXECUTION_TIGHT_BACKGROUND.json, VERIFICATION_EXECUTION.json이 기준이다.

```bash
cd /home/cosmosapjw/ns_cowling_wt/stellar_extension_20261117
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=python_frontend:/tmp/ns_cowling_audit_20261117_pydeps
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python -u scripts/run_stellar_extension_modes_20261117.py --selection data/stellar_extension_20261117/cowling/SELECTION_USED.json
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python -u data/stellar_extension_20261117/cowling/tight_background_check.py
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python data/stellar_extension_20261117/cowling/verify_results.py
```

이 결과의 범위는 비회전·비자성 평형 완전유체의 ell=2 Cowling f/p1이다. Frozen-composition response, buoyancy g-mode, 탄성, 초유체 상대운동, metric perturbation, GR radiation damping과 full-GR mode는 계산하지 않았다. 최대질량 bracket 내부의 고정 모델에 대한 모드이며, 미세하게 다른 최적화 중심밀도의 동일 관측량으로 취급하지 않는다.
