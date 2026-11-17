# 공개 EOS의 최소 열역학 재구성

## 정의와 범위

이 계열의 입력은 공개 CoMPOSE 원본에서 문헌으로 확인된 core 또는 inner-crust
지지범위를 선택한 표다. 생성기의 핵 상호작용, 질량, 표면 매개변수 또는 입자 조성을
다시 계산하지 않는다. `public_original`과 `public_thermo_reconstructed`는 별도
계열이며 원본은 보존한다. native 생성 EOS를 공개본에 맞춘 결과도 이 계열에 넣지 않는다.

생산 계산은 `fortran_eos/app/reconstruct_public_thermo.f90`, 진단·변경 기록·파일 관리는
`scripts/reconstruct_public_thermo_20261117.py`가 담당한다. 기존 `reconstruct_master`는
조성, outer-crust 식과 에너지를 바꿀 수 있어 이 작업에 그대로 사용하지 않았다.

이 재구성의 명시적 선택은 **energy-primary**다. 공개 표의 밀도와 총 에너지 값을
기준으로 압력과 Euler 화학퍼텐셜의 모순을 수선한다. 불일치가 있는 표에서 어느 열이
원 저자의 의도를 더 정확히 표현하는지는 표만으로 유일하게 정해지지 않는다. 따라서
이 결과는 에너지 열을 보존한다는 선택을 붙인 재현 가능한 파생 EOS이지, 원 저자의
숨은 구현을 복원했다는 주장이 아니다. 특히 작은 차이를 저자 코드의 물리 오류로
귀속하지 않는다.

## 보존량과 식

단위는 \(n_B\)가 \(\mathrm{fm}^{-3}\), \(\epsilon,P\)가
\(\mathrm{MeV\,fm}^{-3}\), 화학퍼텐셜과 에너지/바리온이 MeV다. 표의 `mn`은
\(m_n^{\rm ref}c^2\)의 수치다. 에너지에는 정지 에너지가 포함된다.

\[
F_i=m_n^{\rm ref}c^2(1+Q_{7,i}),\qquad
x=\ln(n_B/n_0),\qquad \epsilon=n_BF.
\]

차가운, 검증된 beta-equilibrium 단일 barotrope에서는

\[
Q_1=\frac{dF}{dx},\quad P=n_B\frac{dF}{dx},\quad
\mu_E=\frac{d\epsilon}{dn_B}=F+\frac{dF}{dx},\quad
\frac{c_s^2}{c^2}=\frac{F'+F''}{F+F'}.
\]

각 연속 phase에서 source energy knots를 정확히 통과하는 C2 cubic을 사용한다.
양 끝 두 구간에서 세 번째 미분이 같다는 not-a-knot 조건을 고정하며, smoothing,
회귀, 매개변수 적합, 에너지 기준 이동은 하지 않는다. 미분 계산 중 상수를 빼는 것은
상쇄오차를 줄이는 수치 처리이며, 출력 `Q6`, `Q7`과 `mn`, `mp`는 입력과 정확히 같다.

지원되는 CoMPOSE 저장 규약에서는

\[
m_n^{\rm ref}c^2(1+Q_3)=
\begin{cases}
\mu_E-Y_qm_n^{\rm ref}c^2Q_4,& I_l=0,\\
\mu_E,& I_l=1,\ Q_5\simeq0.
\end{cases}
\]

여기서 `Q5`를 0으로 바꾸어 전제를 만들지 않는다. RG 공개본의 nonzero `Q5`와
정지 에너지 convention은 아직 분리하여 확인해야 하므로 해당 표는 이 단계에서
압력과 `Q3`까지 모두 그대로 둔다. provider가 beta-equilibrium이라고 기술했다는
사실만으로 숫자 열의 규약 문제가 해결되었다고 간주하지 않는다.

`n_B`, `Q2`, `Q4`, `Q5`, `Q6`, `Q7`, 추가 열, phase, 모든 조성 및 질량 기준은
보존한다. 바뀔 수 있는 것은 진단된 `Q1`과 이에 대응하는 `Q3`뿐이다. 코드에서
보존 여부를 배열 값 수준으로 검사하고 receipt에 기록한다. 원본 파일의 byte identity와
재직렬화된 표의 물리적 identity는 구분한다.

## 변경을 결정하는 진단

각 phase의 내부에서 독립적인 좌표 선택 두 가지를 비교한다. 주 표현은
\(F(\ln(n_B/n_0))\) cubic이며, 대조 표현은 \(F(n_B)\) cubic이다. 이들의 압력
미분 차이는 보간 선택에 대한 민감도다. 원본 토큰이 있으면 마지막 출력 자리의
반 단위를 사용해 인쇄 반올림 민감도도 계산한다. 고정 grid에서 cubic 미분은 에너지
샘플의 선형 함수이므로, 미분 행렬의 절댓값을 통해 에너지 반올림의 영향을 전달한다.
밀도 반올림은 일차 민감도로 포함한다. 이는 nuclear-model 불확실성이 아니다.

끝의 두 점씩은 변경 필요성 판정에서 제외한다. 내부 최소 3개 점에서

\[
|Q_{1,\mathrm{derived}}-Q_{1,\mathrm{source}}|
>5\,\delta Q_{1,\mathrm{coordinate}}
 +8\,\delta Q_{1,\mathrm{printing}}
 +10^{-7}|Q_{1,\mathrm{source}}|
\]

이면 해당 phase의 압력을 하나의 energy-primary 미분으로 재구성한다. 이 식의 계수는
고정된 진단 정책이다. 엄밀한 보간 오차 상계나 물리 모형의 신뢰구간을 뜻하지 않는다.
핵심 SPG inner-crust 불일치는 이전 독립 진단에서도 이 민감도보다 훨씬 크다.
압력 재구성이 필요 없고 `Q3`의 Euler 불일치만 출력 정밀도를 넘으면 `Q3`만 수선한다.
이미 일치하면 `UNCHANGED_THERMO`로 원래 숫자를 그대로 유지한다.

최소 3개의 내부 검사를 위해 7개 미만 점의 phase는 수선하지 않는다. phase 경계에는
어떤 보간도 넣지 않는다. SPG/GPPVA의 부착된 BSk22 outer crust는 `Iphase=0` 하나에
여러 핵종 전이가 들어 있으므로, 내부 접합 법칙이 별도로 표현되기 전까지 그 block을
통째로 smooth phase로 재구성하지 않는다. 이 캠페인의 inner-crust 입력은 공개본과
재구성본에서 동일한 phase-7 지지범위를 사용한다.

GPPVA의 추가 열은 source PDF가 정의한 pasta geometry 범주다. 캠페인 준비 단계는
`(source Iphase, geometry)`의 조합마다 working phase를 구분한 뒤 재구성을 호출한다.
원래 추가 열과 phase 표시는 `SOURCE_AUXILIARY.json`에 보존한다. 숫자로 코딩된
geometry를 연속량처럼 보간하거나, 그 열을 지운 뒤 서로 다른 geometry를 하나의
smooth phase로 합치지 않는다. 같은 준비 과정을 original/reconstructed 양쪽에 적용한다.

## 재구성과 Maxwell 채택은 별도

재구성된 표가 만들어졌다고 그 표가 전 구간에서 안정·인과적이라는 뜻은 아니다.
Maxwell 생산 kernel은 표를 다시 \(P(\mu_E)\) Hermite로 표현하며 각 knot에서
\(dP/d\mu_E=n_B\)를 사용한다. 이 표현과 재구성에 사용한 \(F(x)\) cubic은
격자 사이에서 같은 함수라고 가정하지 않는다. 생산 표현은 자체 first law와 각
source energy knot를 만족하며 그 자체의 전 구간 형상 검사를 통과해야 한다.

한 구간의 \(h=\mu_1-\mu_0>0\), \(s=(P_1-P_0)/h\)에 대해

\[
\left.\frac{dn_B}{d\mu_E}\right|_0
=\frac{6s-4n_0-2n_1}{h},\qquad
\left.\frac{dn_B}{d\mu_E}\right|_1
=\frac{-6s+2n_0+4n_1}{h}.
\]

양 끝의 양의 값을 검사하면 선형 \(dn_B/d\mu_E\)의 양성을 검사할 수 있다.
\(n_B\)와 \(\mu_E dn_B/d\mu_E-n_B\)는 이차식이므로 꼭짓점까지 검사한다.
후자가 음수가 아니어야 \(c_s^2/c^2\le1\)이다. Python은 이 검사로 불량 interval과
연속 유효 구간을 진단할 뿐이며, 실제 채택은 독립적으로 반복하는 Fortran Maxwell
검사가 결정한다. 불량 구간을 지우거나 연결하여 하나의 매끄러운 EOS로 만들지 않는다.

`READY`는 재구성 계산 및 변경 기록을 생성했다는 뜻이고 Maxwell 채택 상태가 아니다.
receipt의 `ready_for_matching`, `all_reconstructed_intervals_diagnostic_pass`,
`native_maxwell_source_admission`과 각 interval 기록을 함께 읽어야 한다.
소스 자체의 nonconvex/acausal 구간이 남으면 그대로 보존하고 보고한다. 이 경우에는
partner EOS나 root 위치를 보기 전에 정한 유효한 연속 source 구간만 유한 범위의
후보로 사용할 수 있다. 이를 whole-star EOS로 부르지 않는다.

## 직접 검증

Wolfram Language에서 first law, Euler 관계, sound-speed 식, Maxwell latent-energy
관계의 대수 residual 4개가 0이며, 위 Hermite endpoint 식도 독립적으로 얻었다.
코드와 실제 응답은 `validation/WOLFRAM_RECONSTRUCTION.json`에 보존했다.

native helper의 analytic controls는 정확한 cubic 에너지 미분, SciPy 대조, 압력만
오염된 표의 복구, 원 에너지 및 에너지 기준 보존, 정확한 pass-through, `Q3`만 수선,
phase 경계, nonzero `Q5` 보존, `Il=0` 전하 보정, 나쁜 형상 보존·보고, 따뜻한 표
거절, 부착된 outer-crust block 보존을 검사한다. 실제 실행 로그는
`validation/reconstruction_tests.log`에 남긴다. 최종 195조합 실행의 receipt가 실제
matched EOS 생성 여부와 가능한 유한 지지범위를 결정한다.

호출 예시는 다음과 같다. `prepared_public`에는 원본과 같은 물리량을 가진, 문헌에
따라 선택된 source domain을 넣는다.

```python
from reconstruct_public_thermo_20261117 import reconstruct_data
rebuilt, receipt = reconstruct_data(
    prepared_public, raw_source_dir=original_archive_directory, model=model_id
)
```

독립 helper에는 compiler가 필요하다. Python fallback은 없다. 이미 빌드한 파일을
쓰려면 `COMPOSE_PUBLIC_RECON_EXE`로 지정하며 source와 executable의 SHA-256을
receipt에 기록한다.
