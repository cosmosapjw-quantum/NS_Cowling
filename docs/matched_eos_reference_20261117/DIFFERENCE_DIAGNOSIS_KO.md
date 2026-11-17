# GM1·SLy4 공개본 차이의 제한된 물리 진단

**실행 기준은 2026-11-17에 새로 내려받은 CompOSE 원본이다.** 저장소에 이전부터
있던 SLy4 표와 최신 원본은 동일하지 않다. 과거 표의 약 6% campaign 압력 차이와
진동하는 에너지 잔차를 최신 원본의 문제로 전용하지 않는다. 두 표가 달라진 경위가
제공자 revision인지 저장소 내부 변환인지는 이 진단에서 확정하지 않았다.

생성 EOS의 내부 일관성과 특정 제공자 표의 복제 정확도는 별도 관측량이다.
이 계산은 생산 Fortran이나 그 산출물을 바꾸지 않고, 공개본과의 차이에 대해
새로운 두 가지 수치 control을 수행한다. 정의된 다른 functional을 진단용으로
계산한 사실과 제공자가 그 functional을 썼다는 주장을 구분한다.

## 계산 계약

`scripts/diagnose_reference_differences_20261117.py`는 지정된 공개본의 실제 밀도점에서
기존 독립 Python evaluator를 직접 실행한다. 비교에 interpolation을 사용하지 않는다.
SLy4는 공개 `eos.compo`의 균일상 `Iphase=3`, GM1은 균일상이며 이번 campaign의
`0.08 ≤ n_B ≤ 0.5 fm^-3`에 있는 행을 선택한다. 실제 표의 밀도 index origin과
composition의 `(iT,iN,iYq)` identity를 사용한다. 없는 조성을 0으로 채우지 않는다.

`epsilon`은 정지에너지를 포함한다. `n_B`는 `fm^-3`, `epsilon`과 `P`는
`MeV fm^-3`, 화학퍼텐셜과 `epsilon/n_B`는 `MeV`다. 결과의 차이는
**model minus public reference**다. 파일 header로 복원한 물리량을 사용하며
임의의 정지에너지 이동을 적용하지 않는다. 입력 hash와 evaluator hash는
`data/matched_eos_reference_20261117/diagnosis/DIAGNOSIS.json`에 포함된다.

## SLy4: 공통 kinetic coefficient와 입자별 질량의 효과

현재 생산 모델은 공통 계수 `K=20.73553 MeV fm^2`와 별도 neutron/proton
정지질량을 사용한다. 진단 control은 interaction, 정지질량, lepton 처리와
beta equilibrium 조건을 유지하면서 kinetic term에만

\[
K_n=\frac{(\hbar c)^2}{2m_n c^2},\qquad
K_p=\frac{(\hbar c)^2}{2m_p c^2}
\]

를 사용한다. `hbar*c=197.3269804 MeV fm`, `m_n*c^2=939.5654 MeV`,
`m_p*c^2=938.2720 MeV`이므로 각각 `20.7212489912`, `20.7498130573 MeV fm^2`다.
이것은 앞선 공통 `K=20.734` 변화와 다른 control이다.

\[
\delta\epsilon=\delta K_n\tau_n+\delta K_p\tau_p,
\quad \tau_q=\frac35(3\pi^2)^{2/3}n_q^{5/3},
\]

\[
\delta P=\frac23\delta\epsilon,
\qquad
\delta\mu_q=\delta K_q(3\pi^2 n_q)^{2/3}.
\]

이는 고정된 조성에서 kinetic energy density가 밀도의 `5/3`차 동차함수라는
사실에서 나온다. 새 조성은 수정된 화학퍼텐셜로 beta equilibrium와 전하중성을
다시 풀어 구한다. 네 개의 독립 finite-difference control에서
`delta P=n_B^2 d(delta epsilon/n_B)/dn_B`를 검사한다.

공개 조성을 대입한 beta residual은 `mu_n(model)-mu_p(model)-mu_e(n_B*Ye_public)`로
정의한다. 공개본의 유한한 composition 출력 정밀도는 이 residual의 정확도를
제한한다. `Q5` 해석을 미리 가정하여 source의 물리적 비평형을 선언하지 않는다.
control이 조성 차이를 줄이는 것은 **numerically checked** 결과이며, 제공자가
입자별 kinetic mass를 사용했다는 귀속은 제공자의 정확한 구현을 확인하기 전까지
**unresolved**다. 압력·에너지 잔차가 함께 해소되는지도 별도로 검사한다.

### 실제 실행 결과

이번 campaign 구간에 있는 공개 원밀도점 20개(`0.08317638–0.4786301 fm^-3`)에서
다음 결과를 얻었다. 값은 보간된 native table의 비교와 구분되는, 같은 functional의
독립 Python evaluator를 원밀도점에서 직접 계산한 control이다.

| 관측량 | 현재 공통 K 모델 | 입자별 K 진단 control |
|---|---:|---:|
| 최대 `abs(delta P/P_public)` | 0.1827748% | 0.0112669% |
| 최대 `abs(delta epsilon/n_B)` | 0.0430150 MeV | 0.00288629 MeV |
| 최대 `abs(delta Yp)` | 8.94044e-5 | 5.75593e-6 |
| 공개 조성에서 beta residual RMS | 0.0632829 MeV | 0.00453769 MeV |

따라서 현재의 작은 차이 대부분은 **선언된 kinetic prescription의 변화로 수치적으로
재현되는 효과**다. 임의의 energy offset이나 데이터 보정으로 맞춘 결과가 아니다.
남은 최대 압력 차이 `0.01127%`와 `Yp` 차이는 정확한 provider parameter/convention을
확인하기 전까지 잔차로 보존한다. 이 control을 원저자의 입력 확인이나 생산 functional
변경으로 기록하지 않는다.

범위를 공개 균일상 전체 54행(`0.05248075–1.21 fm^-3`)으로 넓힌 보조 계산에서는
최대 상대 압력 차이가 `0.235022%`에서 `0.0126913%`로 감소했다. 이는 native matching
cohort의 core 상단 `0.5 fm^-3`를 확장한 산출물이 아니라 별도의 원밀도점 진단이다.
고정 조성 kinetic first-law control의 최대 상대 residual은 `3.40854e-11`이다.

## GM1: header mass와 인쇄 정밀도 설명의 한계

첫 control은 원래의 공통 kernel mass `939 MeV`를 공개 thermo header의 neutron
reference mass로 바꾸어 실제 균일상 해를 다시 구한다. Header mass는 에너지
인코딩의 reference이므로 이 control의 결과가 좋거나 나쁘다는 사실만으로
제공자의 kernel mass를 결정할 수 없다.

두 번째 control은 원전 Table II의 다섯 파라미터에 대해 마지막 인쇄 자리의
반 간격 안에서만 변화를 허용한다. `C_sigma,C_omega,C_rho,b,c`의 허용 반 간격은
각각 `0.005,0.0005,0.0005,5e-7,5e-7`이다. 매 여덟 번째 행과 마지막 행에서
선형화한 응답을 구하고 bounded least squares를 수행한 뒤, 얻어진 파라미터에서
실제 비선형 모델을 전체 밀도점에 다시 평가한다. 나머지 행은 held-out control이다.
목적함수의 scale은 `E/A`에 `1 MeV`, 압력에 `1 MeV fm^-3`로 명시한다.

이것은 생산 파라미터 재조정이 아니다. 이 계산이 충분히 맞지 않을 경우의 결론은
**이번 인쇄 정밀도 control로 설명되지 않았다**는 것이다. 전역 parameter box
최적화나 모든 가능한 provider 구현에 대한 불가능성 증명으로 확대하지 않는다.

최신 공개 원밀도점 132개에서 원래 모델의 최대 상대 압력 차이는 `0.486386%`,
최대 `epsilon/n_B` 차이는 `0.298580 MeV`다. Header mass를 kernel mass로 사용한
control은 각각 `0.634622%`, `0.507253 MeV`로 커진다. 따라서 header 값으로
생산 kernel mass를 교체할 근거가 없다.

18개 훈련 밀도와 114개 held-out 밀도의 인쇄 정밀도 control에서는 다섯 계수가
모두 허용 경계에 도달했다. 비선형 재계산의 최대 상대 압력 차이는 `0.178919%`,
최대 `epsilon/n_B` 차이는 `0.190246 MeV`로 남는다. Held-out 최대 에너지 차이도
`0.187327 MeV`다. 이 계산은 인쇄 정밀도가 일부 잔차에 영향을 준다는 사실을
보이지만 전체 불일치가 그 때문이라는 설명을 입증하지 않는다. GM1의 정확한
provider kernel과 파라미터에 관한 `UNRESOLVED_SOURCE_DIFFERENCE`는 유지한다.

## 근거와 종료 조건

원전·기존 독립 evaluator의 계약은
`docs/research_generator_20261116/core_literature.md`,
`reference_cores.md`, `skyrme_implementation.md`에 기록되어 있다.
GM1 원전은 [Glendenning–Moszkowski Table II](https://escholarship.org/uc/item/6sw424k1),
SLy4의 사용 파라미터는 [Douchin–Haensel Table 1](https://arxiv.org/abs/astro-ph/0006135),
공개 비교 대상은 [CompOSE GM1](https://compose.obspm.fr/eos/54)과
[cold SLy4](https://compose.obspm.fr/eos/134)다.

남는 에너지·압력 차이는 밀도별 값과 함께 `UNRESOLVED_SOURCE_DIFFERENCE`로
보존한다. 밀도에 따라 변하는 `epsilon/n_B` 잔차를 상수 offset으로 제거할 수 없고,
상수 offset은 압력에 영향을 주지 않는다. Source 내부의 interpolation 또는
다른 functional 선택은 정확한 제공자 입력 없이 확정하지 않는다.
이번 control 이후 동일한 원인 추정을 반복하지 않고 실제 matching 산출물과
해당 잔차를 함께 보고한다.

재실행은 새로운 출력 디렉터리를 지정해 다음과 같이 한다.

```sh
python scripts/diagnose_reference_differences_20261117.py \
  --source-root data/matched_eos_reference_20261117/references \
  --output build/new_reference_diagnosis
```
