# E: 실제 native SLy4 연장과 TOV 최대질량

## 결과와 판정

D의 입력 commit은 `300e5a0532dd137cb09738568fca0f35cfcc9282`다. B primary의 모든 기존 상태와 하부 연속 표현을 보존하면서, 같은 native SLy4 함수·파라미터·에너지 기준으로 \(n>1\;\mathrm{fm}^{-3}\) 영역을 실제 계산했다. 추가한 물리적 EOS 완성 모형이나 에너지 이동은 없다.

선택된 **B hybrid whole-star equilibrium barotrope**의 TOV 질량 최대점은 native 평형 barotropic 인과 경계보다 먼저 나온다. 검증된 수치 bracket은

\[
1.2060 < n_{c,\max} < 1.2061\;\mathrm{fm}^{-3}.
\]

대표 모델은 \(n_c=1.2060328050005111\;\mathrm{fm}^{-3}\), \(M=2.0503812906965106\,M_{\odot,\mathrm{nom}}\), \(R=10.006409438136822\;\mathrm{km}\), \(P_c=855.9246351636019\;\mathrm{MeV\,fm}^{-3}\)다. 추가 자릿수는 재현할 입력 모델을 식별한다. 과학적 요약에는 \(M_{\max}\simeq2.0503813\,M_{\odot,\mathrm{nom}}\), \(R\simeq10.0064\;\mathrm{km}\), 위의 중심밀도 bracket을 사용한다. 하부 crust와 matching 선택 자체의 모형 불확실성은 이 수치 오차에 포함하지 않는다.

최종 최대점 근거는 `data/stellar_extension_20261117/tov/MAXIMUM_REFINEMENT.json`이다. 초기 optimizer 결과와 진단 실패는 `TOV_RESULT_INITIAL.json`, `ENDPOINT_DERIVATIVE.json`, 원본 stdout에 보존했다. `TOV_RESULT.json`은 초기 후보와 최종 refinement 경로를 함께 담는다. 초기 optimizer 좌표만을 정밀 최대점으로 인용하지 않는다.

## EOS의 정의, 단위와 응답 범위

전 구간은 차가운 유효 beta-equilibrium barotrope다. \(n\)은 \(\mathrm{fm}^{-3}\), \(\epsilon,P\)는 \(\mathrm{MeV\,fm}^{-3}\), \(\mu_B=(\epsilon+P)/n\)는 MeV이며, \(\chi_{\rm eq}=c_{s,\rm eq}^{2}/c^2=dP/d\epsilon\)는 무차원이다. native source의 \(\hbar c=197.3269804\;\mathrm{MeV\,fm}\), kinetic coefficient \(20.73553\;\mathrm{MeV\,fm}^{2}\), rest masses와 SLy4 rounded coefficients를 변경하지 않았다. 정확한 사전은 `eos/SOURCE_BINDING.json`의 `native_definition.source`에 있다.

\[
d\epsilon=\mu_B\,dn,\qquad
\chi_{\rm eq}=\frac{n}{\mu_B}\frac{d\mu_B}{dn}.
\]

native analytic routine은 \(y=Y_p\) 및 공통 전자·뮤온 화학퍼텐셜로 lepton partition을 최소화한 \(\epsilon(n,y)\)의 Hessian을 사용한다. \(H_{yy}>0\)일 때

\[
\chi_{\rm eq}
=\frac{n}{\mu_B}\left(H_{nn}-\frac{H_{ny}^{2}}{H_{yy}}\right),\qquad
\chi_{Y_p}
=\frac{n H_{nn}}{\mu_B}.
\]

native 출력 이름 `cs2_frozen`은 두 번째 식이다. **이 값은 \(Y_p\)를 고정하면서 전자·뮤온 분배는 평형을 유지한다. 개별 모든 입자종을 동결한 완전한 frozen-species 응답은 아니다.** 이 구분과 호출 원자료는 `eos/RESPONSE_SCOPE.json`에 보존했다.

| 구분 | \(\chi=1\)의 native 중심밀도, \(\mathrm{fm}^{-3}\) | 최대질량 대표점의 \(\chi\) |
|---|---:|---:|
| equilibrium barotropic 응답 | 1.2650056715193052 | 0.9680741481811074 |
| 고정 \(Y_p\), lepton 평형 응답 | 1.2397422126050377 | 0.9799407903603997 |

따라서 최대점은 이 두 가지 검사된 응답의 인과 경계 안에 있다. 전체 microscopic causality나 full-crust frozen-species 응답은 이번 결과로 확정하지 않는다.

## 실제 native 연장과 경계 결정

`compose_eos_generate`는 D에서 사용한 빌드의 실제 실행파일을 재사용했다. 새로운 physics 코드는 작성하지 않았다. `native_endpoint_driver.f90`는 이미 존재하는 `solve_core_state`, `core_sound_speed_analytic`, `core_sound_speed_fd`만 호출하는 얇은 실행 진입점이다. 연결 object, 실행파일, D source와 입력의 SHA256을 기록했다.

1. \(n\in[1,1.3999]\;\mathrm{fm}^{-3}\)에서 81개 native 상태를 조사했다. \(P,\epsilon,\mu_B,\chi_{\rm eq}\), composition Hessian의 양수성을 확인했고 \(\chi_{\rm eq}\)가 증가했다.
2. 첫 표본 교차 bracket \([1.26493375,1.2699325]\;\mathrm{fm}^{-3}\) 안에서 analytic Hessian 응답의 \(\chi_{\rm eq}=1\)을 직접 native 호출로 정밀화했다.
3. 상대 간격 \(10^{-3},3\times10^{-4},10^{-4},2\times10^{-5},10^{-5}\)의 native 5점 유한차분으로 각각 경계를 다시 찾았다. analytic/FD 경계 최대 차이는 \(8.50\times10^{-12}\;\mathrm{fm}^{-3}\)였다. 가장 작은 간격에서 roundoff 영향이 보이므로 자릿수를 무한정 정확도 주장으로 해석하지 않는다.
4. 물리적 교차점의 상태는 \(\epsilon=1721.9997877909777\), \(P=976.9256097238558\;\mathrm{MeV\,fm}^{-3}\), \(\mu_B=2133.5282981564446\;\mathrm{MeV}\)다.
5. 공개한 연속 interpolation contract의 공통 상한은 실제 교차점보다 \(2\times10^{-5}\;\mathrm{fm}^{-3}\) 낮은 \(1.2649856715193053\;\mathrm{fm}^{-3}\)다. 이것은 보간의 인과 부등식까지 엄격하게 만족시키는 명시적 support 제한이다. \(\chi=1\) 이후는 native 진단 표본만 존재하며 TOV/Cowling 모형으로 적분하지 않았다.

81개 표본과 root refinement는 측정된 첫 교차의 근거이며, 모든 밀도에서의 전역 단조성에 대한 해석적 정리는 아니다. 현재 분기는 homogeneous native SLy4 functional branch다. 새 고밀도 상전이의 부재나 global ground state를 증명한 것으로 해석하지 않는다.

## 연속 표현과 보존성

기존 `extend_whole_crust`와 `WholeCrustBarotrope`를 변경 없이 사용했다. 새 upper 구간만 실제 native knot에서 cubic \(P(\mu_B)\) Hermite 표현으로 구성했다. 보간은 구간 내부에 국한되며 EOS 외삽은 없다. 각 구간의 \(n>0\), \(dn/d\mu_B>0\), \(\epsilon>0\), \(\mu_B dn/d\mu_B-n\geq0\)를 polynomial extrema로 연속 검사했다.

| contract | 새 native tail 구간 수 | 연장부 연속 최대 \(\chi_{\rm eq}\) | SHA256 |
|---|---:|---:|---|
| `eos/coarse/REPRESENTATION.json` | 64 | 0.9999918438563513 | `993e9e669f43e96e70bca55878c0473b1a51fe8a98c951e8fa572d38a08486b5` |
| `eos/medium/REPRESENTATION.json` | 128 | 0.9999902103496859 | `9ead5ff4548c23fba16c01600439175fff165a4f423ad0e13b9bd635d1950eef` |
| `eos/fine/REPRESENTATION.json` | 256 | 0.9999898002714080 | `f1bd5c3c41b45a07a7299e665594173473ac790af6f519bf53085db843bd936d` |

기존 표의 모든 row field와 composition mask, surface/outer/lower bridge, 기존 upper state 배열은 정확히 보존되었다. 기존 연속 \(n,P,\epsilon,\mu,h\) 값도 원래 knot와 추가 probe에서 bitwise 동일했다. 내부 구간의 \(\chi\)도 동일하다. 원래 끝점 \(n=1\)에서는 기존 마지막 구간의 왼쪽 미분과 새 첫 구간의 오른쪽 미분을 별도로 기록한다. 이는 보간 knot의 일방 미분 차이이며 상전이가 아니다. 원본 A/B 파일은 덮어쓰지 않았다.

## TOV 정의와 수치 검증

signature는 \((- ,+,+,+)\)다. D의 enthalpy DOP853 solver를 수정하지 않고 사용했다.

\[
m=\frac{GM}{c^2},\qquad
h=\int_0^P\frac{dP}{\epsilon+P}
=\log\frac{\mu_B}{\mu_{\rm vac}},\qquad
\frac{dm}{dr}=4\pi r^2\frac{G\epsilon}{c^4}.
\]

geometry는 km, 열역학 입력은 MeV fm\(^{-3}\)이며 \(c=299792458\;\mathrm{m\,s}^{-1}\), \(G=6.67430\times10^{-11}\;\mathrm{m^3\,kg^{-1}\,s^{-2}}\)를 유지한다. solar-mass 표시는 \(GM_{\odot,\mathrm{nom}}=1.3271244\times10^{20}\;\mathrm{m^3\,s^{-2}}\)의 기존 nominal convention이다. 중심 초기조건은 D의 \(O(r_0^2)\) enthalpy 및 \(O(r_0^3)\) mass 정칙 전개이며 surface는 \(h=0\)이다.

초기 21개 고밀도 별 sequence는 최대점 이후의 질량 하강까지 실제로 포함한다. 초기 optimizer는 \(n_c=1.2060130400651856\;\mathrm{fm}^{-3}\)를 선택했으나, 그 주변의 작은 one-sided 차분에서 약 \(10^{-6}\,M_\odot\)의 고립된 mass excursion이 확인됐다. 이를 물리적 구조로 해석하지 않았다. 원래 DOP853 조절값의 국소 미분 수렴 실패로 분류하고, 원자료를 보존한 채 \(\mathrm{rtol}=10^{-12}\), \(\mathrm{atol}=10^{-14}\), \(\Delta h_{\max}=5\times10^{-4},10^{-4}\)로 재계산했다.

64/128/256 tail grid 모두에서 \(n_c=1.2060\)의 중심차분 \(dM/dn_c\)는 \((3.72–3.85)\times10^{-5}\,M_\odot\,\mathrm{fm^3}\), \(n_c=1.2061\)에서는 \((-7.59–-7.47)\times10^{-5}\,M_\odot\,\mathrm{fm^3}\)였다. 차분 간격은 \(5\times10^{-4}\;\mathrm{fm}^{-3}\)다. 추가 \(10^{-3}\) 및 \(5\times10^{-4}\) 주변 질량의 국소 quadratic fit은 모두 음의 곡률 약 \(-1.14\,M_\odot\,\mathrm{fm^6}\)를 주었고, 대표 좌표는 \((1.2060328–1.2060339)\;\mathrm{fm}^{-3}\)에 모였다. 최종 보고 precision은 그보다 보수적인 실제 부호 bracket으로 고정했다.

대표점에서 source grid, ODE maximum step, 중심 반경 \(r_0=10,1,0.3\;\mathrm{m}\)를 바꾼 최대 관측 mass 차이는 \(2.07\times10^{-9}\,M_\odot\), radius 차이는 \(9.87\times10^{-7}\;\mathrm{m}\)였다. bracket 전체에 대한 국소 곡률 기반 mass 변화 추정은 \(2.59\times10^{-9}\,M_\odot\)다. 이는 검사한 오차 성분의 관측치이며 전체 model uncertainty나 엄밀한 전역 error bound는 아니다.

일반적인 평형 sequence의 turning point를 계산했으며 완전한 radial stability spectrum은 계산하지 않았다. Cowling의 두 고정 입력 모델은 `tov/MODE_SELECTION.json`에 있고, 두 번째 \(n_c=1.2060130400651856\)은 검증된 최대점 bracket 안의 선택된 별이다. mode 계산을 위해 이를 최종 대표 좌표로 몰래 교체하지 않았다.

## 실행·검증과 추적성

- C1: 정확한 D commit과 B primary를 고정하고 원본 바이트를 비교했다.
- C2: native source state → Euler chemical potential → \(P(\mu)\) contract → D enthalpy TOV 경로를 연결했다.
- C3: native source-grid, FD derivative, upper interpolation, ODE, 중심 전개 오차를 분리했다.
- C4: thread 수를 1로 고정하고 EOS → TOV → 선택 모델 handoff 순서로 실행했다.
- C5: 실행 가능한 runner와 실제 3개 연장 contract, endpoint 및 maximum 결과를 만들었다.
- C6: 기존 core-extension/TOV 대상 테스트 12개와 E 물리적 artifact regression 7개가 각각 실제 exit 0으로 통과했다.
- C7: 무에너지이동, 재계산된 \(n=1\) native anchor, cold Euler/first-law, compressibility, 응답 Hessian, 전체 보간 인과성을 확인했다.
- C8: 다섯 FD 간격, 세 source grid, 두 강화 ODE step과 세 중심 반경을 검사했고 초기 미분 실패를 보존·해결했다.

EOS 생성 실행은 exit 0, 25.36초, 초기 TOV pass는 exit 0, 20.26초, tightening 진단은 exit 0, 47.78초, 최종 refinement는 exit 0, 101.51초였다. 전 과정의 raw log와 구체적 native command/return code는 data 디렉터리에 보존했다. 원래 실행 runner의 SHA256과 일치하는 snapshot 및 refinement 실행 snapshot도 저장했다. 최종 source에는 두 가지 보고 산식/분기 정리가 추가됐으며, 결과의 metadata 정정은 `tov/POSTPROCESS_RECEIPT.json`에 기록했다.

재실행 예시는 다음과 같다. 기존 실행 디렉터리는 보존하며 새 output을 지정한다.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONPATH=python_frontend \
/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python \
scripts/extend_stellar_support_20261117.py \
--output data/stellar_extension_reproduce
```

EOS만 생성하려면 `--eos-only`, 이미 저장한 initial TOV 후보에 대한 정밀화만 실행하려면 `--refine-only`를 사용한다. 완료된 결과를 무조건 덮어쓰는 경로는 없다. `eos/native_endpoint_driver`는 재컴파일 가능한 로컬 실행물로 Git 대상에서 제외하고, 그 source와 실제 build receipt·hash를 보존한다.

부모 통합 담당자가 전체 테스트, 한 차례 독립 검토, 최종 commit/push 판정을 수행한다. 이 worker는 기존 D 코드나 원본 A/B 자료를 수정하지 않았고 Git commit/push를 수행하지 않았다.
