# 신규 RMF 균일 핵물질 생성: 구현과 검증

## 범위와 보존한 입력

대상은 NL3_PROVIDER_PRECISION, PK1_CANONICAL, TM1_CANONICAL,
DDMEX_PROVIDER_ROUNDED_N0의 네 신규 canonical 모형과, 별도 진단인
DDMEX_ORIGINAL_TABLE_RHO0다. 입력은 원전 대조를 거친 PARAMETER_INPUTS.json에
고정되어 있다. DDME-X 두 입력은 n0=.152와 .151832 fm^-3을 구분하며
질량은 두 경우 모두 명시한 938.5 MeV convention을 유지한다.
어떤 provider 잔차도 parameter나 물리적 에너지의 조정에 사용하지 않는다.

물리 계산은 새 standalone Fortran driver가 기존
compose_eos_tf_models.f90의 tf_couplings, tf_scalar_potential,
tf_vector_potential, tf_fermion을 그대로 호출한다. 기존 GPPVA/TF 생산 경로와
parameter registry는 수정하지 않는다. Python은 frozen parameter namelist,
native process 실행, 검증 기록과 CompOSE 입출력만 담당한다.

## 식과 native 계산의 대응

표의 질량은 정지에너지 m_i c²로 MeV 단위를 사용하고, 길이와 밀도는 fm과
fm^-3이다. ħc=197.3269804 MeV fm을 명시적으로 유지한다.
T=0의 npeμ 물질이며 외부 전자기장, 중성미자 포획은 없다.

M_i*=M_i-g_sigma phi, rho vertex=g_rho tau3/2, tau3(p)=+1이다.
원전의 sigma=-phi와 Gamma_rho를 바꾸면 g_rho=2 Gamma_rho이며

U(phi)=m_sigma² phi²/2 + a3 phi³/3 + a4 phi⁴/4,
a3=-g2 ħc, a4=g3

이다. quartic 부호를 추가로 뒤집지 않는다.
W=m_omega² omega²/2+m_rho² rho²/2+c3 omega⁴/4를 사용한다.
이 driver는 cross coupling=0, c3≥0인 이번 cohort와 DD2 control을 처리하며
그 밖의 functional 입력을 명시적으로 거부한다.

각 (n,Yp)에서 scalar 방정식은

U'(phi)/(ħc)³ - g_sigma(n_s,n+n_s,p)=0

이다. 0≤phi<min(M_n,M_p)/g_sigma에서 256개 구간 scan에 검출되는 상향 부호 전환을 찾고
괄호를 보존하는 이분법으로 근을 계산한다. 양의 effective mass와
U''/(ħc)³+g_sigma² Σ_i ∂n_s,i/∂M_i*>0인 scalar 근 중
kinetic+scalar energy가 가장 작은 검출 근을 택한다.
유한 구간 scan은 임의로 가까운 root pair의 부재나 전체 함수의 global minimum을
증명하지 않는다. 이번 결과는 검출된 양의 질량을 가진 국소 안정 후보라는 범위다.
유효질량이 음수인 해와 해당 구간 밖의 비물리적 stationary branch를
균일물질 core로 채택하지 않는다. 이 기준은 scalar 방향의 국소 안정성
검사이며 완전한 물리적 영역 판정은 이후 native domain screen이 담당한다.

omega는 양의 선형해를 상한으로 이분법을 적용하고 rho는 선형식을 사용한다.
DD 재배열 항은

Sigma_R=g'_omega omega n + g'_rho rho (n_p-n_n)/2
        -g'_sigma phi (n_s,n+n_s,p)

이며 neutron/proton chemical potential에 동일하게 더하고 P에 n Sigma_R를
포함한다. epsilon에는 별도 재배열 에너지 항을 더하지 않는다.
각 종의 P_kin=n_i E_F,i-epsilon_kin 관계로 압력을 구성한다.

n_p=n_e+n_mu는 lepton chemical potential의 bracketed solve로,
mu_n-mu_p-mu_e=0은 0≤Yp≤1의 bracketed solve로 처리한다.
mu_e≤m_mu c²에서는 muon 밀도가 정확히 0이다. composition endpoint의
KKT 분기를 구현하되 생산 결과의 실제 분기는 states.csv에 따로 기록한다.
총 에너지와 압력에는 lepton 정지에너지를 포함한다.

## 산출물과 검증 계약

generated/<MODEL>/coarse 및 fine은 각각 385, 769개 실제 native 상태다.
.015–.16 fm^-3에는 선형 격자, .16–1 fm^-3에는 로그 격자를 사용한다.
각 폴더의 compose/는 기존 native writer로 기록하고 기존 native reader로
되읽은 CompOSE table이다. states.csv에는 field, scalar curvature,
effective masses, 재배열, beta/charge/Euler residual, scalar root 개수를
보존한다. RUN.json과 BUILD.json에는 실행 명령, exit, 입력/실행파일 hash,
compiler와 thread 환경이 담긴다.

검증은 다음의 서로 다른 증거를 구별한다.

- 구현 검증: 모든 생성 상태의 field/beta/charge/Euler residual과 native
  CompOSE readback.
- 과학 검증: 대칭핵물질의 P=0 saturation density, binding energy, K, J, L;
  원전 표의 control은 계산값과 별도로 기록.
- 수치 검증: .015,.03,.06,.1,.2,.4,.8,1 fm^-3에서 상대 density step
  2e-4와 1e-4로 d epsilon/dn=mu_B와 dP/dn=n dmu_B/dn 확인.
- 독립 구현 control: 기존 compose_eos_rmf.f90:evaluate_ddrmf의 DD2 결과를
  .1,.4 fm^-3에서 새 driver와 대조.
- source identity 진단: .1,.2,.4,.8 fm^-3의 exact common homogeneous
  source rows에서 다섯 생성 입력과 네 provider를 교차 비교.

provider raw epsilon은 각 원래 thermo header의 reference neutron mass로
복원한다. 따라서 canonical TM1의 fit mass와 provider header가 다르더라도
차이를 없애기 위한 n에 비례하는 energy shift를 허용하지 않는다.
mu_E=(epsilon+P)/n와 source Q3의 mu_Q는 각각 기록한다.
TM1/PK1 label 충돌은 직접 계산값으로 판단하고, 가까운 다른 모형이
있다는 사실만으로 table의 전체 물리적 identity를 확정하지 않는다.

## 실제 결과

E 게시 commit 00bb5e1adcdbb8dfa0a0b2d5d7dd6f9037a4ddc3을 기반으로 F 작업에서
실제 실행했다. 최초 새 테스트는 구현 전 import 실패(exit 2)를 기록했고,
구현 후 물리 검사 6개가 통과했다. 생산은 13.27초에 exit 0으로 종료했다.
실행 기록은 uniform_validation/PRODUCTION_RUN.json과 stdout/stderr다.

각 입력에 대해 385/769점 table을 생성하고 native write/read를 완료했다.
모든 생산점에서 scalar root는 하나, effective mass와 scalar curvature는 양수,
composition은 interior beta branch였다. coarse의 모든 상태는 fine의 짝수 index
상태와 수치적으로 정확히 같았다. 이는 같은 밀도에서 native 계산이 일관됨을
확인하며 Maxwell root의 격자 수렴성은 별도 matching 결과에서 판정한다.

전체 fine table의 최대 field residual은 1.736e-13 fm^-3, beta residual은 2.387e-12 MeV, charge residual은 3.331e-16 fm^-3, Euler residual은 9.095e-13 MeV fm^-3이다. native readback의 최대 scaled residual은 3.418e-16다.

두 finite-difference step 중 2e-4에서 first-law 최대 상대 잔차는 2.124e-9,
Gibbs-Duhem 최대 scaled 잔차는 2.832e-8이었다. 1e-4에서는 각각
5.316e-10, 8.367e-9로 감소했다. 독립 기존 DD2 native 결과와의
최대 scaled 차이는 7.031e-13이었다.

대칭핵물질의 직접 계산 결과는 다음과 같다. n_sat은 coupling 입력 n0가 아닌
P(n,Yp=.5)=0의 계산 근이다.

| 입력 | n_sat [fm^-3] | E/A-M [MeV] | K [MeV] | J [MeV] | L [MeV] |
|---|---:|---:|---:|---:|---:|
| NL3_PROVIDER_PRECISION | 0.148256510 | -16.249638 | 271.735432 | 37.416617 | 118.595719 |
| PK1_CANONICAL | 0.148200835 | -16.267938 | 282.697339 | 37.642147 | 115.882379 |
| TM1_CANONICAL | 0.145222939 | -16.263118 | 281.164522 | 36.892427 | 110.795492 |
| DDMEX_PROVIDER_ROUNDED_N0 | 0.151579875 | -16.064762 | 266.102053 | 32.255023 | 49.689391 |
| DDMEX_ORIGINAL_TABLE_RHO0 | 0.151659795 | -16.062449 | 266.393389 | 32.243000 | 49.651969 |

NL3는 원전의 근사 NMP control과 부합한다. PK1과 TM1의 원전 control과
작은 잔차는 그대로 보존한다. 특히 PK1 K의 계산값 282.69734 MeV와
표의 282.644 MeV는 정확한 숫자 동일성을 주장하지 않는다.

DDME-X의 원전 supplemental NMP와의 불일치는 남아 있다.
ORIGINAL_TABLE_RHO0 입력의 J=32.2430 MeV, L=49.6520 MeV는 원표의
32.865829, 47.806151 MeV와 다르다. 이 차이의 원인은 현재 입증되지 않았으며
원표의 모든 output까지 재현했다는 claim은 보류한다. 이 사실과 아래의
현재 provider 열역학 일치는 서로 다른 evidence scope다.

### provider 대조에서 해결한 구분

exact common uniform density .1,.2,.4,.8 fm^-3에서 NL3/PK1/TM1은 각각
같은 이름의 current raw provider와 max scaled P/epsilon/mu residual
대략 1e-6–5e-6 수준으로 일치한다. PK1↔TM1을 교환하면 P 잔차가 13–15%다.
따라서 paper parameter table의 label 충돌이 current provider 파일의
thermodynamic curve까지 바꾸었다는 추정은 이 대조로 지지되지 않는다.
TM1의 다른 reference mass header도 물리적 energy shift 없이 처리된다.

DDME-X에서는 원래 고정했던 n0=.151832 진단 입력이 max scaled P=3.956e-6,
epsilon=1.566e-6, mu_E=1.163e-6으로 raw provider에 가깝다.
반면 n0=.152 입력은 P=1.658e-3, epsilon=3.258e-5, mu_E=4.941e-5다.
이는 사전에 이름을 붙여 분리한 두 입력의 forward evaluation 결과이며 fitting이
아니다. 둘의 이름과 수치를 유지한다. 현재 provider를 대표할 입력을 고른다면
.151832 variant의 표본점 근거가 더 강하다. 전체 표의 exact model identity와
원전 supplemental NMP identity까지 확정한 것은 아니다.

PROVIDER_COMPARISON.json에는 5×4 전체 행렬과 원래 header로 복원한 physical
source 값, 별도의 raw mu_Q 차이, species fractions와 source row/hash가 있다.
SOURCE_IDENTITY_ASSESSMENT.json에는 위의 claim 범위를 기계 판독 가능하게 보존했다.
