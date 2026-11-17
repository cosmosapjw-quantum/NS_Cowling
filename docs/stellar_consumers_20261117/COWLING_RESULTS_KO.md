# D Cowling: 평형 barotropic 모드와 실제 수렴

B primary의 연속 EOS 법칙에 묶인 nominal 1.4 M_sun 별에서 l=2 Cowling f 및 p1을 계산했다. 최종 표시값은 f=2401.38 Hz, p1=6721.96 Hz다. 이는 무한대에서 측정한 주파수이며 물리 모델 불확실성을 포함한 오차막대는 아니다. 정확 계약 Python 경로와 독립 native PCHIP 소비자 계산을 모두 실행했다.

## 물리·구현 범위

- Gamma1=Gamma_eq인 평형 완전유체 closure. 조성 마스크를 채우지 않았으며 frozen-species 부력, crust 전단 탄성, 초유체 상대운동, 회전·자기장·metric perturbation은 포함하지 않는다.
- baseline central density 0.5377339866073595 fm^-3, mass 1.40000006478 M_sun, radius 11.7564106743 km. 근접한 TOV canonical 밀도와의 차이는 약 9e-9 fm^-3다.
- 내부 길이는 km, k=2 pi f/c; Phi는 시간 metric의 절반 지수다. h=ln(mu/mu_vac)를 적분 변수로 쓰고 원래 연속 EOS 조각을 직접 평가한다.
- w=W/r^(l+1), v=V/r^l. 중심 (w,v)=(1,-1/l), 표면 Q=0에서 양방향으로 적분하여 내부에서 정규화 determinant를 0으로 맞춘다. W,V는 연속 EOS knot에서 연속이다. sharp Maxwell jump는 이 경로의 입력이 아니다.
- 250–10500 Hz의 실제 250 Hz scan에서 f=[2250,2500], p1=[6500,6750], 다음 root=[9250,9500] Hz bracket을 얻었다. 처음 두 root의 W,V node 수는 각각 0,1이다. 세 번째 모드는 이번 관측량 계산 대상이 아니다.
- 방정식·중심 유도·단위·closure 근거는 FORMULATION_KO.md, 구현은 python_frontend/compose_eos/stellar_cowling.py다.

## 12개 실제 수렴 설정

baseline: mode rtol=1e-8, atol=1e-11, surface h=1e-9, h_match=0.15 h_c, max step h=0.002, r0=1 m. 배경 rtol=1e-10, atol=1e-12. root coordinate xtol=1e-5 Hz보다 실제 ODE 잡음 바닥이 크므로 이 xtol을 주파수의 인증 오차로 해석하지 않는다.

| 설정 | f [Hz] | p1 [Hz] | f shift [Hz] | p1 shift [Hz] | matching 판정 |
|---|---:|---:|---:|---:|---|
| baseline | 2401.381789874 | 6721.957630556 | +0.000000000 | +0.000000000 | accepted |
| ode_loose | 2401.380070391 | 6721.933133142 | -0.001719483 | -0.024497415 | f rejected; diagnostic |
| ode_tight | 2401.381839494 | 6721.957071699 | +0.000049620 | -0.000558857 | accepted |
| surface_h1e-6 | 2401.381979377 | 6721.957396874 | +0.000189503 | -0.000233683 | accepted |
| surface_h1e-7 | 2401.382043925 | 6721.957402900 | +0.000254051 | -0.000227656 | accepted |
| surface_h1e-8 | 2401.381932110 | 6721.957605252 | +0.000142235 | -0.000025304 | accepted |
| match_h007 | 2401.382092050 | 6721.957812834 | +0.000302176 | +0.000182278 | accepted |
| match_h030 | 2401.381350721 | 6721.957386810 | -0.000439153 | -0.000243746 | accepted |
| step_h1e-3 | 2401.382325588 | 6721.957359276 | +0.000535714 | -0.000271280 | accepted |
| center_10m | 2401.382326881 | 6721.957070991 | +0.000537006 | -0.000559565 | accepted |
| center_01m | 2401.381585015 | 6721.957226756 | -0.000204859 | -0.000403800 | accepted |
| background_loose | 2401.381810272 | 6721.955019045 | +0.000020397 | -0.002611512 | accepted |

고정 matching/join criterion은 2e-6이다. loose rtol=1e-6의 f matching residual=−3.8695e-6은 거절된 진단값으로 보존했다. root 좌표 수렴을 mode residual 통과로 바꾸지 않았다. 초기 runner 중단 로그와 수정 후 전체 결과가 함께 있다. 나머지 11개 설정에서 관측한 주파수 범위는 f 0.000976159381 Hz, p1 0.00279378961 Hz다. 이 범위는 설정 민감도의 관측값이며 엄밀한 전역 오차 상계가 아니다.

표면 h=1e-6부터 1e-9, 중심 10/1/0.1 m, h_match/h_c=0.07/0.15/0.30, ODE tolerance, step cap 및 배경 tolerance를 실제로 바꿨다. 표면 조건의 작은 잔차는 초기조건으로 부과한 결과다. 고유해의 독립적인 수치 근거는 내부 matching, node 수, 설정 민감도와 다른 구현의 비교다.

## native 소비자와 표현 민감도

native cell_solve가 실제로 컴파일한 provider는 oscillation_source/solutions/src/eos/ns_table_eos.f90의 monotone PCHIP(g=epsilon−P)(ln P)이며, n과 mu는 first law quadrature로 복원한다. native TOV와 동일 provider를 쓰지만 적분 경로는 다르다. 초기 runner의 Hermite P(mu) 설명은 잘못된 추론이었다. 원 로그는 보존하고 NATIVE_COMPARISON.json과 runner의 설명·source binding을 교정했다. 수치 출력을 재계산하거나 바꾸지 않았다.

동일 central pressure, h_cut=1e-4와 해당 절단 별의 Schwarzschild 시간 정규화에서 비교했다. exact 기준은 원 EOS 법칙을 직접 평가한 절단 배경이며 vacuum 배경의 lapse와 혼용하지 않는다. exact cutoff의 f=2401.381874932 Hz, p1=6721.957056883 Hz다.

| 원 EOS 구간 세분 | sample rows | native f [Hz] | native p1 [Hz] | native−exact f [Hz] | native−exact p1 [Hz] |
|---|---:|---:|---:|---:|---:|
| 1 | 514 | 2401.382888667 | 6722.041930072 | +0.001013735 | +0.084873189 |
| 4 | 2053 | 2401.382080279 | 6721.964473836 | +0.000205347 | +0.007416953 |

두 native 계산은 exit 0이며 f/p1 normalized matching residual은 3.1e-10 이하, W 및 V node 수는 각각 0/1이다. sample refinement에서 exact와의 p1 차이가 0.0849에서 0.00742 Hz로 감소했다. 원 B EOS를 다시 fitting하거나 energy offset을 넣지 않은 소비자 표현 민감도다. 결과의 0.01 Hz 표시는 이 관측된 수치 차이를 고려한 반올림이며 실제 중성자별의 물리 예측 정확도를 뜻하지 않는다.

## 검증·실행·보존

- analytic control: 약한 중력의 비압축성 Schwarzschild interior에서 Cowling omega²=l GM/R³ 극한을 검사했다. 이 control의 무한 음속은 물리 baseline이 아닌 제한해 검산용이다.
- scaled RHS를 원 W,V 식으로 변환한 값, 잘못된 scan/domain의 거절을 포함한 3개 dedicated test가 실제 통과했다. 전체 회귀 결과는 상위 D validation receipt를 따른다.
- exact campaign은 exit 0, 실행시간 502.50 s. 12 설정의 JSON과 f/p1 eigenprofile CSV, 실제 scan, native namelist/stdout/stderr/background/eigenprofiles를 data/stellar_consumers_20261117/cowling/에 보존했다.
- source SHA, native executable SHA와 전체 compiled source digest는 COWLING_SOURCE_IDENTITY.json 및 NATIVE_COMPARISON.json에 있다. VALIDATION.json은 75개 source identity, 12 설정, node 판정, native exit code 및 source EOS hash의 실제 readback을 기록한다.

재실행은 repository root에서 다음 환경으로 scripts/run_stellar_cowling_20261117.py --central-density 0.5377339866073595 및 scripts/run_stellar_cowling_native_20261117.py를 실행한다:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
PYTHONPATH=python_frontend:pipeline:/tmp/ns_cowling_audit_20261117_pydeps \
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python -u scripts/run_stellar_cowling_20261117.py --central-density 0.5377339866073595
```

Frozen-species response의 주파수 이동은 미계산이다. 필요한 off-equilibrium derivatives가 thermo-only crust 계약에 없기 때문이다. 이 계산은 f/p1의 선택한 평형 Cowling 관측량을 검증하며, GR radiation damping, 전체 mode 안정성, radial 안정성이나 최대질량을 증명하지 않는다. APR Gibbs consumer 결과를 이 SLy4 hybrid 별의 APR 모드로 해석하지 않는다.
