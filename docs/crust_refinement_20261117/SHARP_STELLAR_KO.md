# G: sharp-interface TOV/Cowling reference

## 계산의 의미와 closure

입력은 `WholeCrustBarotrope`가 admission한 piece contract다. 수치 좌표는
\(h=\int_0^P dP/(\epsilon+P)=\log(\mu/\mu_{\rm vac})\)이며,
각 Maxwell interface에는 같은 압력/chemical potential과 서로 다른 \(n,\epsilon\)이 있다.
density gap에 임의의 homogeneous state를 만들지 않는다.
`SharpContractEnthalpyEOS`는 각 smooth interval의 원래 Hermite/outer/surface law를 직접 평가한다.

G의 stellar 입력은 E의 outer/surface contract를 그대로 재사용하고, source crust의 각 저장 knot를 낮은 밀도부터 검사하여 causal/stable Hermite lower bridge가 처음 허용되는 knot를 선택한다.
선택은 stellar 관측량 계산 전에 이루어진다. 다음 admissible knot를 lower-join sensitivity control로 보존한다.
이 bridge는 명시적 phenomenological lower closure다. native TF가 그 구간을 재생성했다는 뜻이 아니다.
실패한 lower join은 `FULL_CONTRACT_ADMISSION.json`에 모두 남긴다.

Bulk perturbation closure는 \(\Gamma_1=\Gamma_{\rm eq}\)인 barotropic perfect fluid다.
Frozen composition, elasticity, surface stress/energy, finite conversion time은 포함하지 않는다.
두 conversion limit은 모델 가정으로 따로 계산한다. EOS 표만으로 실제 conversion rate를 결정하지 않는다.

## TOV와 단위

Metric signature는 \((-,+,+,+)\),
\[
ds^2=-e^{2\Phi}c^2dt^2+e^{2\Lambda}dr^2+r^2d\Omega^2,
\qquad e^{-2\Lambda}=1-2m/r,\qquad m=GM/c^2 .
\]
코드의 \(r,m\)은 km, source \(P,\epsilon\)은 MeV fm\(^{-3}\)다.
`MEV_FM3_TO_GEOM_KM2`가 \(G/c^4\)와 SI–km 변환을 함께 적용한다.
\[
\frac{dr}{dh}=-\frac{r(r-2m)}{m+4\pi r^3P_g},\qquad
\frac{dm}{dh}=4\pi r^2\epsilon_g\frac{dr}{dh},
\]
\[
\Phi(h)=\tfrac12\log(1-2M_g/R)-h+h_s,\qquad
\Phi'=\frac{m+4\pi r^3P_g}{r(r-2m)} .
\]
각 interface에서 적분을 중단하고 같은 \(r,m\)으로 다음 one-sided branch를 시작한다.
Dense output도 개별 phase segment 안에서만 사용한다.
\(r,m,P,\Phi,\Lambda,\Phi'\)는 연속이고 \(\epsilon,n,c_s^2\)는 one-sided다.

중심 Taylor 초기값은
\[
m_0=\frac{4\pi}{3}\epsilon_{g,c}r_0^3,\qquad
h_c-h_0=2\pi(\epsilon_{g,c}/3+P_{g,c})r_0^2 .
\]
표면은 background에서 \(h_s=0\)인 해석적 vacuum endpoint까지 간다.
Mode integration은 singular endpoint를 피하도록 유한 \(h_{\rm mode}>0\)에서 \(\Delta P=0\)을 부과하고 cutoff sensitivity를 실제로 계산한다.
Finite scan의 mass crossing으로 고른 \(1.4M_\odot\) branch는 radial stability 판정과 구분한다.

## Cowling bulk equations와 sharp 접합 유도

\(\ell\ge2\), \(W=r^{\ell+1}w\), \(V=r^\ell v\),
\(\omega_g=2\pi f/(c\text{ in km s}^{-1})\)를 쓴다.
\[
\frac{dw}{dr}=\frac{\omega_g^2 e^{\Lambda-2\Phi}rv+\Phi'w}{c_s^2/c^2}
-\frac{(\ell+1)(\ell e^\Lambda v+w)}{r},
\]
\[
\frac{dv}{dr}=\left(2\Phi'-\frac{\ell}{r}\right)v-\frac{e^\Lambda w}{r}.
\]
Source의 smooth interval에서만 이를 적분하며 \(dh/dr=-\Phi'\)를 사용한다.
중심 regularity는 \((w,v)=(1,-1/\ell)\)이고 임의의 전체 normalization은 고유진동수에 영향을 주지 않는다.

\[
A=\omega_g^2e^{\Lambda-2\Phi}r,\qquad B=\Phi',
\qquad \xi^r=e^{-\Lambda}r^{\ell-1}w,
\]
\[
\Delta P=-(\epsilon+P)e^{-\Lambda}r^{\ell-1}(Av+Bw).
\]
여기서 \(A,B\)는 둘 다 km\(^{-1}\)이며 source/target side에서 동일하다.
표면층 없는 fluid interface에서는 \([\Delta P]=0\)이다.
원전의 nonradial junction conditions를 Cowling variables로 옮긴 결과는 다음과 같다.
[S1, Sec. II D; S2, Sec. III]

**Slow/material conversion:** fluid displacement도 연속이므로 \([w]=0\).
\(q=(\epsilon_s+P)/(\epsilon_t+P)\)라 하면
\[
\begin{pmatrix}w_t\\v_t\end{pmatrix}=
\begin{pmatrix}1&0\\(q-1)B/A&q\end{pmatrix}
\begin{pmatrix}w_s\\v_s\end{pmatrix}.
\]

**Rapid/equilibrium conversion:** interface displacement
\(\xi^r-\Delta P/P_0'\)가 연속이다.
\(P_0'=-(\epsilon+P)B\)를 대입하면
\[
\xi^r-\Delta P/P_0'
=-e^{-\Lambda}r^{\ell-1}(A/B)v
\]
이므로 \([v]=0\)이고
\[
\begin{pmatrix}w_t\\v_t\end{pmatrix}=
\begin{pmatrix}q&(q-1)A/B\\0&1\end{pmatrix}
\begin{pmatrix}w_s\\v_s\end{pmatrix}.
\]
두 행렬은 \(q\to1\)에서 identity이며,
역방향 transfer는 \(q\to1/q\)로 정확히 주어진다.
\(\Delta P\) 보존, inverse transfer, no-jump limit을 Wolfram으로 symbolic check하고 Python의 독립 물리 조건 테스트로 확인했다.
Source 논문의 hadron–quark conversion 동역학을 crust에 그대로 적용했다고 주장하지 않는다.
여기서의 사용은 표면층 없는 generic perfect-fluid interface의 두 이상화된 boundary condition이다.

고유진동수는 중심→matching point와 표면→matching point의 solution vectors가 선형 종속되는 normalized determinant의 영점이다.
Matching point는 interface와 겹치지 않는다.
Interface가 어느 integration 방향에 속하든 동일한 방향별 transfer를 적용한다.
Profile에는 interface 양쪽을 별도 행으로 보존한다. Discontinuous component의 sign change는 smooth radial node와 자동으로 동일시하지 않는다.
Finite frequency scan에서 검출한 root를 기록하며 complete spectrum이나 mode-family 식별을 자동 선언하지 않는다.

## 실제 검증 및 기록

- 최초 테스트는 module 부재로 수집 단계 exit 2였다. 구현 후 junction traction/kinematics, inverse/no-jump limit, density-gap rejection, smooth-background/mode limit, 별도의 pressure-coordinate TOV 비교가 **6 passed, exit 0**이었다.
- 독립 검토자는 pressure-coordinate TOV와 unscaled \(W,V\) Cowling equation 및 직접 2×2 junction solve를 구현하여 finite-jump slow/rapid를 별도로 비교했다. 결과는 `review/INDEPENDENT_SHARP_CONTROL.json`과 독립 검토 보고서를 따른다.
- `stellar/SYMBOLIC_CHECK.json`에는 정확한 Wolfram input/result와 선행 context helper의 internal error를 구분하여 보존했다.
- `stellar/PRODUCTION_RECEIPT.json`, `STELLAR_RESULTS.json`에는 실제 process exit, 각 full-contract hash, code hashes, scan, mode roots, 접합 residual, lower-join/mesh/numerical sensitivity가 기록된다.
- Stellar 수치 요약과 허용된 claim은 G의 `RESULTS_KO.md` 및 최종 검토 결과를 따른다. 이 문서의 방정식 검산만으로 실제 계산 성공을 대신하지 않는다.

## 원전과 구현

- **S1** H. Sotani, K. Tominaga, K. Maeda, *Density discontinuity of a neutron star and gravitational waves*, [arXiv:gr-qc/0108060](https://arxiv.org/abs/gr-qc/0108060), Sec. II D.
- **S2** L. Tonetto, G. Lugones, *Discontinuity gravity modes in hybrid stars: Assessing the role of rapid and slow phase conversions*, [arXiv:2003.01259](https://arxiv.org/abs/2003.01259), Sec. III.
- **Implementation:** `python_frontend/compose_eos/stellar_sharp.py`; old D/E `stellar_background.py`, `stellar_cowling.py`, and `whole_crust.py`의 smooth piece laws를 재사용한다.
- **Reproduction:** `scripts/prepare_sharp_stellar_20261117.py`, `scripts/run_sharp_stellar_20261117.py`; 입력은 `data/crust_refinement_20261117/stellar/INPUTS.json`.

### 중심 근접 interface의 경계 입력 수리

독립 검토에서 \(h_0\le h_t<h_c\)이면 중심 Taylor 영역 안에 interface가 숨는 입력을 확인했다.
이제 solver는 그 입력을 거절하고 더 작은 \(r_0\)를 요구한다.
최초 regression은 실제 실패했고, 수리 뒤 같은 상태에서 거절과 작은 \(r_0\)의 정상 interface 검출을 각각 확인한다.
이미 실행한 생산 사례에는 새 guard가 모두 false이므로 수치 경로가 같다.
원 생산 source snapshot과 guard만 제거했을 때 전체 module AST가 일치하는 검산을
CENTER_GUARD_REPAIR.json에 보존했다. 전체 생산을 새 코드로 다시 했다고 주장하지 않는다.
