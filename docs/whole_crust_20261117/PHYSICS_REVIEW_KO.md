# A. Whole-crust 연결의 독립 물리 검토

기준: `8c917baf3f5975a6109c3b6e59177079c285c4e0`, 2026-11-17.
이 문서는 생산 구현과 독립적으로 원자료를 읽고 계산한 진단·구성식을 기록한다.
최종 선택한 endpoint, 출력 개수, 실행 검증은 해당 생산 보고서를 따른다.

## 범위와 원래 계획

`docs/matched_eos_reference_20261117/PLAN_KO.md`의 A는 공개 BSk22 lower와
native SPG–SLy4 `quintic_P` 결과를 연결하는 `reference-assisted hybrid` baseline이다.
CSS는 construction sensitivity 대상으로 남는다. 이 A 후보는 이전 sharp Maxwell의
세 출처 계열을 바꾸거나 교차시키지 않는 별도 산출물이다. 다음 B는 필요한 core
지원 확장 및 source-variant 정책, C는 TF continuation과 나머지 방법 확대,
D는 APR 경계와 TOV/Cowling 소비자 연결이다.

조건은 찬 barotrope, 전체 baryon 수 한 conserved charge, local neutrality,
beta equilibrium을 전제로 한다. epsilon은 정지 에너지를 포함한 에너지밀도,
n은 baryon 수밀도, mu는 에너지 단위의 baryon chemical potential이다.
아래 식의 sound-speed 비는 `c_s^2/c^2`이며 `c=1`을 물리 단위로 가정하지 않는다.
접합 대상 source의 epsilon에 임의의 `C*n`을 더하지 않는다.

## 공개 lower 실제 자료와 숨은 경계

직접 읽은 자료는 `eos_generator/eos_data/SPG_M1/{eos.nb,eos.thermo,eos.compo,eos.pdf}`다.
공급자 PDF는 BSk22 lower를 n<0.002 fm^-3에서 쓴다고 설명한다. 실제 파일에는
첫 217행이 n=4.67962e-10에서 2.69021e-4 fm^-3까지 있고, 다음 점은
0.002511886432 fm^-3이다. 따라서 공급자의 문장만으로 그 사이의 표 값을
존재하는 것으로 취급할 수 없다. `eos.compo`의 Iphase는 첫 217행에서 0,
그 뒤 148행에서 7이다. `eos.thermo` 끝의 0은 추가 열 개수이며 phase 번호가 아니다.

첫 217행의 composition에는 18개 isotope 구간과 17개 전환이 있다. 각 전환을
가로지르는 연속 두 행은 같은 압력의 Maxwell doublet가 아니다. 압력은 약
0.3% 증가하며 Euler mu 차이는 2.03e-4에서 5.53e-3 MeV이다. 공급자의 Q3와
`(epsilon+P)/n`도 다르므로 어느 chemical-potential 표현을 사용하는지 명시해야 한다.

| 항목 | 직접 복원한 값 |
|---|---:|
| 첫 n [fm^-3] | 4.6796200000e-10 |
| 첫 P [MeV fm^-3] | 9.000000079422e-12 |
| 첫 epsilon/n [MeV] | 929.78185997485 |
| 마지막 lower n [fm^-3] | 2.6902100000e-4 |
| 마지막 lower P [MeV fm^-3] | 4.981910045986e-4 |
| 마지막 lower mu, 저장 Q3 복원 [MeV] | 938.8962300409328 |
| 마지막 lower mu, Euler [MeV] | 938.8975403205699 |

모든 n,P,epsilon knot를 그대로 고정한 cubic P(mu_Euler) 보간은 216개 interval 중
40개에서 양의 compressibility를 위반했다. 좁은 같은-isotope 구간의 인쇄 정밀도와
원 표의 열역학적 잔차가 영향을 준다. 몇몇 isotope gap은 다음 필요조건도 위반한다.

\[
n_L < \frac{P_R-P_L}{\mu_R-\mu_L} < n_R.
\]

이는 mu가 증가하고 n=dP/dmu가 증가하는 매끄러운 구간에 평균값 정리를 적용한
조건이다. 해당 endpoint를 고정한 안정한 구간의 부재를 말하며, BSk22 물리 모형의
일반적 불가능성을 주장하지 않는다. raw lower를 한 개의 매끄러운 phase로 spline하는
방법은 이 문제와 17개 isotope 경계를 함께 숨긴다.

## 압력 중심의 명시적인 lower 표현

실제 사용할 수 있는 한 구성은 공개 n,P knot를 보존하고, 첫 epsilon/n만 에너지
적분의 경계조건으로 고정하는 것이다. 이 선택은 이전 세 계열 작업의 energy-primary
`public_reconstructed`와 별개의 A용 표현이다. 원본 epsilon 및 Q3와의 차이를
수정량으로 기록하고 raw 파일은 그대로 보존한다.

각 tabulated 구간에서

\[
\Gamma_i=\frac{\ln(P_{i+1}/P_i)}{\ln(n_{i+1}/n_i)},\qquad
P(n)=P_i(n/n_i)^{\Gamma_i},\qquad F=\epsilon/n
\]

로 정하면 first law는

\[
\frac{dF}{dn}=\frac{P}{n^2},\qquad
F(n)=F_i+\frac{P_i/n_i}{\Gamma_i-1}
\left[(n/n_i)^{\Gamma_i-1}-1\right],\qquad
\mu=F+P/n.
\]

Gamma=1의 극한은 `F=F_i+(P_i/n_i)*ln(n/n_i)`이며, 수치 계산은
`expm1((Gamma-1)*ln(n/n_i))/(Gamma-1)` 또는 그 연속 극한을 사용한다.
이 표현의 압력과 chemical potential은 knot에서 연속이다. sound speed는

\[
\frac{c_s^2}{c^2}=\frac{\Gamma_i P}{\epsilon+P}.
\]

각 구간의 Gamma는 상수이므로 이 값의 극값은 양의 에너지 영역에서 endpoint로
검사할 수 있다. 같은-isotope 구간과 isotope 사이의 표본 공백은 서로 다른 origin
label을 갖는다. 후자는 **유효 barotropic gap 표현**이며 미시적 isotope 공존점의
복원으로 부르지 않는다. 공백의 composition을 만들어 채우지 않는다.

위 식을 원파일의 217개 lower 행에 직접 적용한 독립 계산 결과:

- Gamma 범위: 0.06534021644367108–1.5973483109163742.
- `c_s^2/c^2` 범위: 6.72999635734103e-6–0.0026519594679295807.
- raw epsilon/n에서의 변화: 0–0.004915606957411001 MeV.
- raw epsilon에 대한 최대 상대 변화: 5.245856306975331e-6.
- 마지막 lower mu: 938.9024559275274 MeV.

이 적분은 첫 에너지 경계값이 정한 유일한 표현이며 상대 SPG나 접합 성공을 보고
에너지 상수를 조절하지 않는다. source 전체에 대한 임의 에너지 gauge 조절과
구별된다. 압력을 1차 입력으로 고르는 것은 과학적 모델링 선택이며 표만으로
원저자의 완전한 microscopic EOS를 유일하게 복원했다는 뜻은 아니다.

## Lower–native bridge의 구성과 연속 구간 검사

공개 lower의 마지막 점을 고정하고 native SPG의 기존 point 중 높은 밀도의
endpoint를 선택할 수 있다. n,P,epsilon 양끝값을 그대로 쓰려면
`mu=(epsilon+P)/n`을 endpoint chemical potential로 사용한다.
양끝 mu 차이를 h>0, t=(mu-mu_L)/h, s=(P_R-P_L)/h라 두면

\[
P(t)=P_L+h n_Lt+a_2t^2+a_3t^3,
\quad a_2=h(3s-2n_L-n_R),
\quad a_3=h(n_L+n_R-2s).
\]

이때 `n=P_mu`, `epsilon=mu*n-P`로 정의하면 양끝 물리량과 first law가 보존된다.
`P_mumu`는 선형이므로 그 양끝 양성이 전 구간 양성과 동등하다:

\[
\frac{2n_L+n_R}{3}<s<\frac{n_L+2n_R}{3}.
\]

causality 함수 `C(mu)=mu*P_mumu-P_mu`는
`C'(mu)=mu*P_mumumu`이고 세 번째 미분이 상수다. mu>0에서 C는 단조이므로
양끝 C>=0 검사로 `c_s^2/c^2<=1`을 전 구간에서 확인할 수 있다.
매끄러운 endpoint의 sound speed가 양옆 source와 같아야 한다는 조건은 이 cubic의
계약에 포함되지 않는다. 필요한 경우 더 높은 매끄러움은 별도 sensitivity 표현이다.

위 압력 중심 lower의 마지막 점과 129점 native SPG 입력을 대조했을 때,
처음으로 convexity와 causality를 통과한 native endpoint는 다음과 같았다.
이는 생산 endpoint를 강제하는 선택 규칙이 아니라 구성 가능성에 대한 독립 확인이다.
최초 가능 점은 경계 조건에 가까울 수 있어 실제 baseline의 선택과 민감도를 별도로 남긴다.

| SPG | native 0-based row | n_R [fm^-3] | P_R [MeV fm^-3] | mu_R [MeV] | bridge max c_s^2/c^2 |
|---|---:|---:|---:|---:|---:|
| M1 | 36 | 0.00118949328287 | 0.00122188957159 | 940.155364524644 | 0.0338087583 |
| M2 | 34 | 0.00116166576203 | 0.00119915564904 | 940.138982897282 | 0.188926463 |
| M3 | 32 | 0.00101478087556 | 0.00106994645704 | 940.003841477676 | 0.0346900123 |
| M4 | 21 | 0.000661360769725 | 0.000795172757679 | 939.643201378794 | 0.0315305595 |
| M5 | 42 | 0.00154984847341 | 0.00152826647687 | 940.379415212785 | 0.0481001936 |

각 native 입력 경로는
`data/matched_eos_reference_20261117/campaign/pairs/SPG_Mi__RG_SLY4__quintic_P/work/crust_internal/`
이며, source metadata의 mn=938.9 MeV를 그대로 적용했다.

### M1의 유한 구간 선택 민감도

M1의 첫 허용 endpoint 바로 근처에서는 cubic의 낮은 쪽 compressibility가 작아
sound speed가 빠르게 변한다. 처음 가능한 점을 유일한 물리적 전이로 채택할 근거는
없다. 공개 outer 끝점을 고정하고 다음 세 endpoint를 비교했다.

| native row | n_R [fm^-3] | P_R [MeV fm^-3] | mu_R [MeV] | max c_s^2/c^2 | q=(s-n_L)/(n_R-n_L) |
|---|---:|---:|---:|---:|---:|
| 36 | 0.00118949328287 | 0.001221889571586 | 940.155364524644 | 0.0338087583 | 0.3352559553 |
| 37 | 0.00124457146555 | 0.001272203681634 | 940.196710739115 | 0.0161255555 | 0.3372622229 |
| 40 | 0.00142558585275 | 0.001443012745796 | 940.324777257681 | 0.00697434184 | 0.3417538543 |

convexity 경계는 q=1/3이다. 세 경우의 최대 sound speed는 낮은 endpoint에서 생긴다.
공통 밀도 구간 `[2.69021e-4, 0.00118949328287] fm^-3`의 101개 로그 균일 점에서
직접 cubic n(mu)를 역으로 풀어 비교했다. 이 비교의 표본 최댓값은 다음과 같다.

| 비교 | 압력 상대차의 절댓값 최댓값 | epsilon/n 차이의 절댓값 최댓값 [MeV] |
|---|---:|---:|
| row37 vs row36 | 0.00382179802 | 0.00184730289 |
| row40 vs row36 | 0.0128036004 | 0.00675535067 |

이 독립 검토의 제안은 row40을 **명시적인 유한 bridge variant**로 고정하고
row36/37을 선택 민감도로 남기는 것이다. 전이 밀도의 추정값이나 자유에너지
최소화 결과로 부르지 않는다. source의 에너지 상수와 핵 매개변수는 모든 variant에서
같다. 생산의 최종 선택은 해당 receipt에 명시한다.

## 표면 경계와 cutoff의 수렴

공개 첫 점 `n_a,P_a,epsilon_a` 아래에 명시적 Gamma_s=5/3 polytrope를
선택하면 다음과 같다.

\[
P=K n^{\Gamma_s},\qquad K=P_a/n_a^{\Gamma_s},\qquad
\epsilon=a n+\frac{P}{\Gamma_s-1},\qquad
a=\frac{\epsilon_a}{n_a}-\frac{P_a}{(\Gamma_s-1)n_a}.
\]

\[
\mu=a+\frac{\Gamma_s}{\Gamma_s-1}\frac{P}{n},\qquad
h(n)=\int_0^{P(n)}\frac{dP}{\epsilon+P}=\ln\frac{\mu(n)}{a}.
\]

Gamma_s>1, a>0이면 `n→0`에서 `P,epsilon→0`, `mu→a`, `h→0`이다.
Gamma_s=5/3은 저밀도 비상대론적 전자 압력의 지수를 참고한 표면 closure이며,
전리·격자·대기 미시물리를 새로 계산한 EOS로 해석하지 않는다. 다른 지수의
표면 closure는 선택 민감도로 구분한다.

첫 공개 점에서 `a=929.75301147835 MeV`, `h_a=5.171221145205719e-5`다.
cutoff에 남는 무차원 enthalpy는 다음과 같이 해석적으로 감소한다.

| n_cut/n_a | n_cut [fm^-3] | h(n_cut), Gamma_s=5/3 | h(n_cut), Gamma_s=4/3 |
|---:|---:|---:|---:|
| 1e-2 | 4.67962e-12 | 2.40032741424e-6 | 1.78265483066e-5 |
| 1e-4 | 4.67962e-14 | 1.11413456738e-7 | 3.84064026489e-6 |
| 1e-6 | 4.67962e-16 | 5.17135484178e-9 | 8.27442108485e-7 |
| 1e-8 | 4.67962e-18 | 2.40033029474e-10 | 1.78267056115e-7 |

Gamma_s=4/3 표면 표현은 같은 첫 공개 n,P,epsilon을 보존하며
`a=929.72416298185 MeV`, `h_a=8.274082196532061e-5`다. Gamma_s를 바꾸면
`h~n_cut^(Gamma_s-1)`이므로 cutoff 수렴 차수가 달라진다. 두 선택의 surface
enthalpy 차이와 각 선택 안의 cutoff 오차는 별도로 기록해야 한다.

이는 EOS 표면 표현과 적분 cutoff의 수렴이다. stellar M–R 및 mode의 실제 수렴은
원래 순서 D의 해당 소비자 계산에서 평가한다. 유한 양의 밀도 표와 함께 위 n=0
극한식을 제공하면 소비자가 임의 상수밀도 surface를 추가할 필요가 없다.

## 원자료 identity와 문헌

| 파일 | SHA256 |
|---|---|
| eos.nb | 46e5fb4504d3e53a6d61f210c1627234e34f8ed993fd54921299203756b0818c |
| eos.thermo | 20ef365b0fc9ff66b3b268c7f319d2d8dca0b782bd44a32612d3c4a82c8641e5 |
| eos.compo | 8ce704b9af8a3a7a74e6985d9edc5481a8dc7a8b03e50f7b4d39c3d292dfec27 |
| eos.pdf | 491c23709054bf8af3fa78b5f0563fbb5b1c3bc073dfea1747d3f6daa5e488d2 |

2026-11-17에 SciSpace/Parallel로 찾고 원 출처를 확인했다.

1. [BSk 공식 analytic-code 설명](https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html).
   독립적으로 fitting된 열역학 함수를 모두 동시에 사용하는 대신 항등식으로
   필요한 함수를 계산하는 consistency 처방을 제시한다. 여기의 구간별 압력 적분은
   그 일반 원칙에 따른 직접 구성이다. 본 문서의 수치와 endpoint는 자체 계산이다.
2. [BSk22 CompOSE 제출 문서](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK22/eos.pdf).
   outer mass input, inner ETFSI 및 모델 범위의 원 설명이다.
3. [SPG M2 공급자 페이지](https://compose.obspm.fr/eos/316)와 저장된 M1 PDF.
   BSk22 lower의 명시적 사용 범위와 model provenance를 확인했다.
4. [Pearson et al. 2018 논문 검색 레코드](https://scispace.com/papers/unified-equations-of-state-for-cold-non-accreting-neutron-dzcha15gut).
   통일된 BSk microscopic EOS와 이번 서로 다른 functional의 hybrid를 구별하는 근거다.
5. [Crust-core interface and bulk neutron star properties, Appendix A](https://arxiv.org/html/2102.10948v3).
   압력을 적분해 chemical potential을 구성하는 사례다. 이 논문의 core `n*Delta_mu`
   보정은 현재 source energy 보존 조건에 맞지 않으므로 이번 방법에 적용하지 않는다.
6. [CompOSE Reference Manual](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf).
   접합에서 압력·chemical potential의 공통 조건과 crust 선택 민감도를 설명한다.

확보한 결과는 원자료 진단과 조건부 barotrope 구성의 수학·수치 검산이다.
생산 산출물의 파일 왕복 읽기, origin map, 전체 retained 구간 검사 및 최종 선택은
별도 실행 기록과 연결한다.

## 추가 판정: 공개 outer nuclear-group fraction의 표준화

생산 native writer가 공개 M1의 composition을 거절한 원인을 공식 정의와 대조했다.
근거는 [CompOSE manual v3](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf)
§3.3, Eqs.(3.10)–(3.12), §4.2.3, Eqs.(4.19)–(4.22) 및
[출판본 DOI](https://doi.org/10.1140/epja/s10050-022-00847-y)다.

표준 particle fraction은 `Y_i=n_i/n_b`, baryon fraction은 `X_i=B_i*Y_i`다.
Nuclear group의 마지막 항도 baryon fraction이 아니라

\[
Y_I=\sum_{j\in M_I}Y_j
=\frac{\sum_{j\in M_I}A_jY_j}{A_I^{\rm av}}
\]

라는 합산 number fraction이다. 따라서 native parser의
`Yh=Aav*Y_I`는 올바른 표준 변환이다.

M1 공급자 PDF의 마지막 쪽은 composition을 neutron/proton/electron/muon 네 pair와
heavy nucleus 한 quadruple로 명시하고 phase 0을 outer, phase 7을 inner라고 정의한다.
직접 읽은 첫 217행은 모두 다음 조건을 만족했다.

- phase=0, 네 particle ID=10,11,0,1, nuclear group ID=1 하나.
- 정수 A,Z와 `0<Z<=A`; quadruple의 마지막 값은 정확히 1.
- free neutron, free proton, muon fraction이 모두 정확히 0.
- electron fraction은 Z/A를 소수점 다섯 자리로 출력한 값이다.
- `abs(Ye-Z/A)`의 최대는 `4.999999999977245e-6`이다.

이 조건에서는 baryon closure `A*Y_I=1`로부터 표준 fraction이 유일하게
`Y_I=1/A`로 결정된다. Raw last field=1을 표준 Y_I로 읽으면 baryon fraction이
A>1이 되므로 현재 export는 표준 number-fraction 정의와 맞지 않는다. Pure outer
조건 아래의 교정은 isotope abundance의 새로운 물리적 추정 없이 가능하다.

**판정: source identity와 위 조건이 확인된 217개 outer 행에 한정한 normalizer를
채택할 수 있다.** 내부 Yh를 1로 설정하고 표준 CoMPOSE writer로 `Y_I=1/A`를
출력한다. A,Z,Ye와 free fractions는 그대로 둔다. 교정 뒤 baryon residual은 0이며,
charge residual은 원본 Ye 출력 정밀도에서 온 최대 5e-6이다. 이 변환은 pressure,
epsilon 또는 chemical potential을 바꾸지 않는다.

Raw 값, 기존 parser에서 생긴 FAIL 상태, 변경 위치와 이유를 ledger에 보존한다.
해당 guard를 만족하지 않는 row에는 자동 변환을 적용하지 않는다. 특히 inner phase의
nuclear-group 의미는 이번 pure-outer 유도에서 결정되지 않으므로 그대로 미해결로 둔다.
이 교정과 pressure-primary 에너지 적분은 서로 다른 변경으로 기록한다.
