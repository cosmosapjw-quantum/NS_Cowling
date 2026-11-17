# EOS generator 연구·구현 루프 — 2026-11-16

기준 커밋은 677585947e4d126bc27cc3d8f30c265c407a7773이다. 목표는 원래의 core 5종과 crust 13종을 원전 물리와 조성에 따라 독립적으로 재생성하는 것이다. 전체 목표의 상태는 PARTIALLY_COMPLETE다. 이번 실행은 밀도별 유효구간 진단, core 조성 재생성, SPG 원전 bulk 입력 복원, 전체 source component 확보를 진전시켰다.

## 구현과 실제 계산

- 원본 component 경로는 core 5종, crust 13종 모두 생성·native readback을 수행했다. Crust는 총4172행이다. TW의 첫217행 열역학과 조성이 DD2/DDME2와 일치함을 확인했다. 원자료의 BPS/BSk22 서지 충돌은 receipt에 기록된다.
- 새 physical_domain_report는 각 밀도의 열역학·조성 가용성, 인접 구간의 안정성·인과성, 최초 위반 구간을 기록한다. 고밀도 문제를 포함한 전체 원자료와 유효한 저밀도 구간을 함께 보존한다. 조밀하지 않은 격자는 연속 보간의 불확실성으로 기록한다.
- GM1은 원전 비선형 RMF와 npe 조성으로 독립 계산했다. 290개 원본 core 행과 비교한 최대 압력 차이는0.4864%, 에너지밀도 차이는0.06161%다.
- KDE0v1·SLy4는 완전한 비대칭 Skyrme functional과 해석적 화학퍼텐셜에서 npemu 베타평형을 풀었다. 274+54개 원본 격자점에서 조성 포함 CompOSE 출력과 native readback을 수행했다. 최대 압력 차이5.619%/13.089%는 저밀도 쪽에 집중되어 있고, 정확한 provider 상수와 압력 구성 방식의 대조가 남아 있다.
- APR은 공개된21개 계수로 low/high-density Hamiltonian, 해석적 화학퍼텐셜과 npemu 조성을 계산했다. 국소 중성 Maxwell과 전역 중성 two-charge Gibbs 공존을 각각 구현했다. 2025-03-14 공식 제공자 문서는 Gibbs 처방을 명시하며, 최신 eos.thermo는 저장소 원자료와 SHA-256이 일치한다. Gibbs 경계는 nB=0.2043949549–0.2373842489 fm^-3, 압력은5.839084–5.902855 MeV fm^-3이다. 원본과의 최대 압력 차이55.07%는 nB=0.2에서 남아 있어 상 처방 교정만으로 원본 일치가 완성되지는 않았다. 내려받은 APREOS 저자 Fortran 루틴을 그대로 컴파일한12개 독립 대조에서는 압력이1.23e-6 MeV fm^-3 이내로 일치했다. nB=0.2의 원본과 생성값은 각각 자기 에너지 미분과 일치하므로, 남은 차이는 bulk functional/data 또는 제공자 규약에 관한 문제다.
- BBB2는 원전의 Paris+Urbana SNM/PNM 에너지 표와 parabolic asymmetry에서 조성을 계산했다. 양성자 분율의 최대 절대차는0.000978이다. 압력은 원전 Table3에 대해 최대65.18% 차이를 보여 PRESSURE_MISMATCH_REQUIRES_ORIGINAL_FIT로 기록했다. 후속 검산에서 baryon beta-equilibrium 경로를 미분한 후 lepton 압력을 더하는 처방이 차이를 설명하는 유력한 원인으로 나타났다. CubicSpline 대조에서는 이 처방으로 실제 CompOSE의15개 압력을 최대0.6413% 차이로 재현했다. 원전 저자 코드의 확인은 아직 없으므로 이는 수치적 근거를 가진 추론이다. 기본 generator는 총 에너지에 대해 열역학적으로 일관된 압력을 유지한다. 원전의 정확한 미분 처방과 smooth-fit 계수의 확보가 남는다.
- Native DD2·DDME2·TW에서 총771개 npemu 상태를 생성하고 조성·열역학의 CompOSE readback을 검증했다.
- 저자 Zenodo 데이터의447346개 파라미터 행 중 SPG M1–M5에 대응하는 후보를 추출했다. 새3-meson GDFM 구현이 저자 데이터의 포화 특성과 대칭에너지를 재현했다. 기존4-meson surrogate와 구별되는 source bulk 경로에서1285개 npemu 상태를 생성했다. M2/M4는 논문의 인쇄된 Ksat와 아카이브 수치 사이의 작은 차이를 후보 identity 기록에 남겼다.

- SPG M1–M5 모두에서 AME2016 후보2408개 핵에 대한 표면계수 fit 실험을 실행했다. vacuum spherical droplet 가정과 고정된 bulk 입력에서 sigma0=1.0178–1.0659 MeV fm^-2, b=14.0848–34.8553을 얻었다. 질량 RMS 잔차는15.2184–22.4442 MeV이며, 원전 fit의 재현으로 판정하지 않는다. 질량/결합에너지 가중치, holdout, 밀도격자 정련과 multistart를 비교했다. 독립 연속밀도 최소화와 격자 근사는 핵 질량에서1.04e-4 MeV 이내로 일치한다. 학습핵의 최소 Yp는0.3076923으로, neutron-rich crust 조성까지의 surface 함수는 외삽이다.

## 물리적 유효구간

원본 APR 표의 최초 초광속 secant 구간은 nB=[0.88,0.90] fm^-3이고 KDE0v1은 [1.425315,1.442943] fm^-3이다. 이 값은 표의 두 행 사이 구간 진단이다.

독립적으로 구현한 매끄러운 Skyrme 균일상에서는 미분 간격을 정련해 c_eq^2/c^2=1 위치를 KDE0v1 1.4109658882 fm^-3, SLy4 1.2650056682 fm^-3로 구했다. 마지막 간격 정련에 따른 변화는 각각1.55e-8,9.87e-9 fm^-3이다. 압력의 원본 대비1% 이내 일치 구간은 원본 격자에서 KDE0v1 0.1137985–0.1314325 및0.2102486–1.767846, SLy4는0.1 단일점 및0.2089296–1.21 fm^-3이다. 이는 원본 일치도 구간이며 물리적 인과성 구간과 구별한다. SLy4 원본 비교 격자의 끝은1.21 fm^-3이다. 이 경계는 명시한 functional, 상수, 균일상 및 약상호작용 평형 가정에 적용된다. 스핀·유한파수 안정성은 별도 의무다.

BBB 원전은 약0.8 fm^-3까지 BHF 계산의 수치적 안정성을 설명하고 그 위에는 명시적 rational extrapolation을 사용한다. 표의 최대밀도1.5 fm^-3, microscopic 계산 범위, 현재 PCHIP 재구성 범위를 구분해 기록했다.

밀도 미분은 조성 평형 응답을 포함한다. THERMODYNAMIC_CONTRACT.md에 자유에너지 최소화와 Schur complement 유도, Wolfram 검산, 고정 조성과 평형 음속의 구분을 정리했다. 원본의 Euler 관계 불일치는 원전 모델의 고밀도 인과성 한계와 별도로 남아 있다.

## 지원 범위가 좁았던 원인과 교정

1. Source extraction, 독립 미시계산, 원본 수치 일치, 유효구간 판정이 이전 상태표에서 혼합되어 있었다. 실행 결과와 모델별 남은 의무를 분리해 MODEL_COVERAGE.json에 기록한다.
2. TW는 실제 바이트/물리필드 비교 없이 서지 충돌만으로 전체 component가 보류되었다. 직접 비교 근거에 따라 원본 추출을 허용하고 충돌 정보를 출력에 보존했다.
3. 기존 SPG는 원전과 다른4-meson surrogate였다. 원전은 sigma/omega/rho의3-meson GDFM이며 x=n/n_sat이고 delta meson이 없다. 저자 파라미터를 사용한 신규 구현과 NMP 대조로 이 경로를 교정했다. 과거 M3/M5의 optimizer·spline 실패는 그 surrogate와 당시 격자에 대한 결과로 유지한다.
4. GPPVA crust는 공간 의존 finite-range Thomas–Fermi 모델이다. 이를 재현하려면 해당 meson gradient·Coulomb·pasta geometry를 포함한 방정식을 풀어야 한다. SPG surface fit을 공통 선행조건으로 삼는 분류를 교정했다.

## 전체 재생성 완료에 필요한 작업

SPG5종의 source bulk는 계산 가능해졌다. 후보 surface sigma_0,b도 이번에 계산했으나, 원전 crust 재현에는 원전과 동등한 finite-nucleus fit 정의와 전체 CLD 연결 검증이 필요하다. Zenodo에는 bulk parameters/NMP와 nB,epsilon,P,Yp의 EOS 열이 있고 surface 계수는 수록되어 있지 않다. 확보한 AME2016 후보 집합을 사용한 명시적 fit 실험은 실행했으며, 해당 실험과 원전의 fit 결과 사이의 동일성은 아직 확립되지 않았다. candidate_surface_fit.md에 에너지식, 질량 기준, 잔차와 민감도를 기록했다.

GPPVA8종은 source-backed bulk 계수, 공간 의존 meson/Coulomb 장, 입자 밀도와 셀 크기·geometry를 함께 최적화하는 Thomas–Fermi 구현이 필요하다. 이번에는 FSU2 계열과 TM1e 등의 정확한 문헌 입력을 추가 수집했다. 벌크 균일상 계산으로 finite-range crust 검증을 대신하지 않는다.

BBB의 압력 미분 처방, Gibbs 처방 확인 후에도 남는 APR bulk/provider 차이, Skyrme 저밀도 압력 차이를 해소해야 원본 일치도를 완료로 판정할 수 있다. 새 Python 물리 구현은 연구용으로 실행·검증되었으며, 기존 native-only production 목표에 대한 Fortran 이관도 남은 의무로 기록한다.

## 재현·근거

data/research_generator_20261116 아래에 전체 모델 원본 진단, 실제 생성된 CompOSE 및 조성 행, 원본 비교, 소스 해시, 로그와 체크섬이 있다. core_literature.md와 crust_literature.md가 원전·표·식과 구현의 대응을 제공한다. 최종 local 통합검사는 frontend/pipeline 636 passed, 9 skipped; generator 142 passed였다. 기존 Matplotlib Axes3D 경고는 남아 있다. 검증 명령·버전·실행 코드 해시와 로그는 EXECUTION.json에 기록한다.

R1–R4/C1–C4에서 목표·문헌·정의·구현 경로를 확인했고, R5–R7/C5–C8에서 위 functional·조성·유효구간을 실제 계산했다. R6/C9 독립 검토에서 cold F=epsilon 검사 누락과 비중성 경로의 mu_B 단조성 판정 문제를 발견해 교정했다. R8–R10/C10의 전체 목표 판정은 PARTIALLY_COMPLETE이며 남은 작업은 위 원전 일치도와 crust 방정식 구현이다.
