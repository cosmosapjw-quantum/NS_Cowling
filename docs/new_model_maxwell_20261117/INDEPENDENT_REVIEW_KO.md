# F 신규 RMF·sharp Maxwell 독립 검토

## 판정

**현 요청 범위에서 BLOCK_NOW 없음.** 고정된 다섯 입력의 균일 RMF 계산과 120개 모델·변형 조합의 sharp Maxwell 캠페인은 **선택한 유한 source/component 범위의 수치 결과로 수용**한다. 9개 primary 접합 출력과 신규 native core를 포함한 7개 coarse 출력의 실제 파일을 독립적으로 대조했다. 발견한 수렴 보고의 표현 오류 1건은 수치 계산을 변경하지 않고 닫았다.

이 판정은 original microscopic functional의 전역 최소, 연속 TF mesh 극한, 원전의 모든 핵물질 특성 재현, whole-star EOS나 별 관측량의 검증을 뜻하지 않는다. 특히 DDME-X 원전 supplemental J/L 불일치는 unresolved이며, 현재 provider 표본점 일치와 구분해야 한다.

검토 기준은 E 게시 commit `00bb5e1adcdbb8dfa0a0b2d5d7dd6f9037a4ddc3` 위의 F 작업이다. 하나의 fresh-context C9 검토로 새 driver, adapter, generator, matching code 및 위험한 물리·수치 가정을 검토했다. 기존 kernel·전체 역사 캠페인의 재감사는 수행하지 않았다. 부모가 실행한 전체 frontend/pipeline suite와 본 검토의 독립 수치 계산은 서로 다른 증거다.

## 1. 식과 구현의 대응

질량 변수는 rest energy $m_i c^2$를 MeV로 나타내고, 밀도는 fm⁻³, 에너지밀도와 압력은 MeV fm⁻³다. `ħc=197.3269804 MeV fm`을 유지한다. 물리계는 T=0의 npeμ, 중성미자 비포획, 국소 전하중성의 beta equilibrium이다.

검토한 원전 근거는 보존된 Long et al. `nucl-th/0311031`의 Lagrangian·self-energy·meson equations 및 Table I, Xia et al. `2208.12893`의 Table I–II, Taninah et al. `1910.13007`의 DD coupling equations와 supplemental Table II다. 원전 파일·추출문 및 바이트 identity는 `sources/SOURCE_LEDGER.json`에 결속되어 있다.

| 항목 | 직접 확인한 대응 | 판정 |
|---|---|---|
| scalar 부호 | 원전 $M^*=M+g_\sigma\sigma$, 코드 $\phi=-\sigma$, $M^*=M-g_\sigma\phi$ | 일관 |
| 비선형 scalar | $a_3=-g_2\hbar c$, $a_4=g_3$; quartic 부호를 추가 반전하지 않음 | 일관 |
| rho 정규화 | 원전 full vertex $\Gamma_\rho\tau_3$와 코드 $g_\rho\tau_3/2$, $g_\rho=2\Gamma_\rho$; neutron/proton isospin convention에 따른 field 부호도 함께 변환 | 일관 |
| omega | $c_3\ge0$, cross=0에서 양의 해에 대한 bracket; cubic의 도함수가 양수 | 이번 입력 범위에서 적절 |
| DD rearrangement | $\Sigma_R=g'_\omega\omega n+g'_\rho\rho(n_p-n_n)/2-g'_\sigma\phi n_s$ | chemical potentials에 공통 추가, P에 $n\Sigma_R$, ε에 독립 에너지 항을 추가하지 않음 |
| leptons | 전자·뮤온 rest energy, charge solve, muon threshold, beta residual | 실제 선택 상태에서 확인 |
| scalar branch | $0\le\phi<\min(M_n,M_p)/g_\sigma$, 양의 curvature를 갖는 검출 근 중 scalar+kinetic energy 최소 | 유한 scan의 검출 후보 범위 |
| beta endpoints | Yp=0/1의 KKT 부호와 interior solve의 방향 | 식은 일관; 생산 상태는 모두 interior |

On-shell energy와 pressure의 vector quartic 계수도 교차 확인했다. 에너지에는 $3c_3\omega^4/4$, 압력에는 $c_3\omega^4/4$가 들어간다. driver는 field source를 포함하는 동등한 식으로 ε를 계산하고, 독립 검산은 위의 on-shell 식을 사용했다. 두 표현이 동일한 결과를 주었다.

원천은 `fortran_eos/tools/rmf_candidate_driver.f90`, `python_frontend/compose_eos/rmf_candidate.py`, `scripts/generate_rmf_candidates_20261117.py`와 변경되지 않은 `compose_eos_tf_models.f90`의 호출 대상 함수다.

## 2. 독립 균일물질 계산

`independent_review/independent_check.py`는 `compose_eos`와 native RMF kernel을 import하지 않는다. 원전의 음의 sigma field와 full rho vertex를 사용하고, 192-node Gauss–Legendre 운동량 적분으로 각 Fermi gas의 ε, P, scalar density 및 mass derivative를 계산했다. DD coupling derivative는 complex step으로 얻었다. 선택된 beta 상태에서 4096-bin scalar scan과 별도의 Brent solve를 사용했다.

NL3, PK1, TM1, 두 DDME-X 입력 각각에서 fine-grid row 0, 160, 512, 736, 768을 선택했다. 총 25개 상태는 저밀도, $n=0.16\,\mathrm{fm}^{-3}$, 고밀도 및 $n=1\,\mathrm{fm}^{-3}$ 끝점을 포함한다. 각 상태의 composition 방향 국소 curvature를 검사했고, 각 입력의 $n=0.16\,\mathrm{fm}^{-3}$에서는 17개 Yp 표본으로 beta 함수의 부호와 단조성을 확인했다.

| 비교량 | 최대 scaled 잔차 |
|---|---:|
| ε | 8.0565×10⁻¹⁵ |
| P | 1.0640×10⁻¹³ |
| μn | 4.2296×10⁻¹⁵ |
| μp | 5.0461×10⁻¹⁵ |
| μe | 6.2003×10⁻¹⁶ |
| scalar field | 1.2396×10⁻¹⁴ |
| rearrangement | 1.2079×10⁻¹³ |
| scalar curvature | 4.2744×10⁻¹⁵ |

scaled 잔차는 해당 열의 지정 단위에서 `|reference-native|/max(1,|native|)`다. 독립적으로 계산한 beta residual의 최대 절댓값은 $1.4779\times10^{-12}\,\mathrm{MeV}$다. 25개 상태의 scalar root count와 국소 안정성은 생산 결과와 부합했다. 17개 composition 표본에서의 단조성은 연속 구간 전체의 전역 convexity 증명이 아니다.

10개 coarse/fine table에서 coarse의 모든 열은 fine의 짝수 index 열과 정확히 같았다. 이는 공통 밀도 계산의 일관성 증거다. source ledger의 87개 파일 바이트, generation에 결속된 namelist·density·실행파일·CompOSE hashes도 대조했다.

실제 실행 exit는 0이다. 최초 검산 시점에 production source hash와 현재 driver hash가 달라 검토 스크립트의 엄격한 byte assertion이 종료됐으며, 원인은 root-scan 범위를 명확히 한 주석 변경이었다. 정확한 원래 주석을 복원하면 production SHA가 재현됨을 확인했다. 이때의 로그와 주석 whitespace 재구성 오류 로그를 보존했고, 수치 식·대상 상태·허용오차는 바꾸지 않았다. 최종 증거는 `UNIFORM_SUBSET.json`, `UNIFORM_SUBSET_RUN.json` 및 log다.

## 3. Source identity와 남은 문헌 불일치

공개 PK1/TM1의 이름은 원 논문의 parameter vector와 actual thermodynamic curve를 함께 보아야 한다. 보존된 Xia Table I에는 canonical PK1/TM1 값의 이름 충돌이 있지만, 현재 raw 파일은 각각 같은 이름의 canonical model과 표본점에서 가깝다. 두 이름을 바꾸어 해석하면 pressure 차이가 약 13–15%인 반면 같은 이름의 P·ε·μ 비교는 약 10⁻⁶–5×10⁻⁶ 수준이다. `PROVIDER_COMPARISON.json`의 5×4 행렬과 raw header 복원 방식을 검토했다.

DDME-X는 사전에 고정된 두 n0 입력을 계속 구분한다.

| 입력 | current provider에 대한 max scaled P 차이 | 근거 범위 |
|---|---:|---|
| n0=0.152 fm⁻³ | 1.6575×10⁻³ | 명시된 rounded 입력의 forward evaluation |
| n0=0.151832 fm⁻³ | 3.9561×10⁻⁶ | 원래 고정된 대안 입력의 forward evaluation |

두 경우 physical nucleon mass 938.5 MeV 가정을 유지했다. parameter를 연속적으로 조정하거나 residual을 작게 만들기 위한 energy shift를 사용하지 않았다.

**unresolved:** ORIGINAL_TABLE_RHO0의 $J=32.2430\,\mathrm{MeV}$, $L=49.6520\,\mathrm{MeV}$는 original supplement의 32.865829, 47.806151 MeV와 다르다. 이 불일치의 원인은 입증되지 않았다. 따라서 현재 provider의 네 exact common density에서 열역학이 가깝다는 결과를 original supplemental nuclear-matter outputs 전체의 재현으로 확장할 수 없다.

자료에서 `mu_Q` 또는 `muQ`라고 부르는 비교량은 source Q3에서 복원한 baryon chemical potential이다. 전하 chemical potential과 혼동하지 않도록 본 검토에서는 μQ3와 μE=(ε+P)/n으로 구분했다.

## 4. Maxwell 캠페인과 실제 출력 대조

새 source parser는 phase 0을 terminal homogeneous core로, nonzero geometrical phase를 crust로 구분한다. TM1 Nadd=5와 다른 새 public source의 Nadd=7을 검증하며, Q1–Q7은 수선하지 않는다. 원문 auxiliary tokens·원문 thermo line·source phase를 sidecar에 모두 남긴다.

독립 readback은 34개 source/role/resolution record, 9개 public projection의 모든 selected Q1–Q7 및 모든 auxiliary line을 대조했다. native interval diagnostic에서 연속 valid interval을 직접 다시 묶어 준비된 component 목록과 비교했고, 각 component 안의 phase가 하나임을 확인했다. 이 과정에서 대조한 source·diagnostic·출력 파일 identity는 343개다.

| origin 범주 | 셀 | 실제 kernel 호출 | 수용 출력 |
|---|---:|---:|---:|
| 기존 native crust × 신규 native core | 65 | 305 | 6 |
| 신규 public crust × 기존 public GM1 | 5 | 28 | 1 |
| 신규 public crust × 신규 public core | 20 | 3696 | 0 |
| 명시적 reference-assisted hybrid | 30 | 168 | 2 |
| 합계 | 120 | 4197 | 9 |

성공한 여섯 native/native 조합은 SPG_M1–M5 및 GPPVA_DDME2 crust와 PK1_CANONICAL core다. 공개 TM1 crust는 public GM1과 한 건, native GM1과 별도 hybrid 한 건이 수용되었다. 공개 DDME-X crust와 native PK1의 hybrid가 나머지 한 건이다. 같은 GM1 이름이라도 public/native source origin을 합치지 않았다.

9개 primary 및 7개 coarse 출력에 대해 raw `RESULT.json`, `expected_internal`, `readback_internal`, 실제 CompOSE table을 함께 읽었다. 출력 hash, Q1–Q7, 모든 internal 열과 metadata가 일치했다. 두 interface row에서 P, μE, 양의 density jump와 `Δε=μt Δn`을 직접 재계산했다. **16개 실제 출력에서 이 세 상대 잔차의 최댓값은 8.9307×10⁻¹⁶**였다.

interface 밖의 retained rows를 source n·P·ε·Q3-derived μ와 대조했다. 서로 다른 mass header를 쓰는 입력에서도 물리량을 보존하는 재인코딩이었으며, physical energy shift는 0이다. 출력은 phase label만 제공하는 thermo-only EOS다. abundance를 이용한 C3/C4 검증이 불가능하다는 상태를 유지하므로, consistency PASS를 matched composition의 검증으로 읽어서는 안 된다.

독립 polynomial envelope 자료의 9개 결과도 대조했다. 같은 두 component의 전체 공통 μ 지원 범위에서 real root 하나, 왼쪽 `Pcore-Pcrust<0`, 오른쪽 `>0`, 전이 slope 차이 양수다. SciPy와 native root의 최대 차이는 $3.2173\times10^{-11}\,\mathrm{MeV}$다. 이는 선택한 finite Hermite 표현의 수치적 root enumeration이며, 다른 geometry나 제외 구간까지 포함한 global minimum 증명이 아니다.

주 증거는 `independent_review/MATCHING_READBACK.json`과 exit 0인 `MATCHING_READBACK_RUN.json`이다. 120셀 kernel 캠페인을 독립 검토에서 다시 실행하지 않았다.

## 5. 수렴 결과와 발견 사항

동일 crust component를 유지하고 신규 PK1 core의 385/769-row 두 격자를 바꾼 7개 성공에서 전이는 모두 유지된다.

| 양 | 최대 상대 차이 |
|---|---:|
| μt | 5.7713×10⁻¹² |
| Pt | 3.4730×10⁻⁹ |
| nlow | 2.3938×10⁻⁹ |
| nhigh | 2.1269×10⁻⁷ |
| Δn | 4.4848×10⁻⁷ |
| Δε | 4.4848×10⁻⁷ |

두 격자 비교는 core discretization sensitivity의 측정이다. asymptotic order, 엄밀한 error bound, public crust의 추가 정밀도 또는 기존 TF mesh 수렴을 제공하지 않는다.

### F-C9-001 — CLOSED, representation issue, 비차단

- **대상:** `convergence()`와 `CONVERGENCE.json`의 `relative_difference.maxwell_identity_residual`.
- **관측 결함:** fine residual이 0일 때 coarse 잔차를 1e-300으로 나누어 약 1e284의 의미 없는 상대차를 만들었다.
- **영향:** 실제 전이 위치·density jump·latent energy의 수렴 불량으로 오해할 수 있다. 전이 상태나 생성 table 자체는 변하지 않았다.
- **최소 수정:** 이 residual은 relative comparison에서 제외하고 fine/coarse 원값과 absolute difference로 기록했다.
- **종결 근거:** 7개 기존 결과의 원값이 보존됨을 직접 확인했다. production runner snapshot은 campaign SHA와 정확히 일치하며, AST 비교에서 변경 함수는 `convergence` 하나다. `prepare`, `prepare_record`, `execute`, `run_candidate`, `make_pairs`는 그대로다. 생산 계산을 재실행하지 않았다. 근거는 `matching_validation/CONVERGENCE_REPORTING_REPAIR.json`이다.

그 밖에 새로운 BLOCK_NOW defect는 발견하지 않았다. DDME-X NMP identity, source finite support, global branch claim은 아래의 미해결 범위로 남는다.

## 6. 실패 해석과 최종 claim 범위

4,188개 rejected component pair를 물리적 불가능으로 통합하지 않는다. 3,852건은 공통 μ 지원 부재, 328건은 공통 구간 내 pressure equality 부재, 2건은 nonpositive density jump, 2건은 multiple crossings다. raw UNKNOWN 네 건은 세 건의 retained-source-knot 부족과 한 건의 output C2 rejection으로 분류되었다.

마지막 C2 사례의 원본 public DDME-X crust에는 `max|μQ3-μE|=0.0041373185 MeV`가 남아 있다. 별도 source 진단에서 C2 상대 잔차 4.4023902×10⁻⁶가 strict 2×10⁻⁶ gate를 넘었다. source Q3·energy를 바꾸거나 gate를 넓혀 성공으로 만들지 않은 처리가 적절하다. 이 원인은 `REJECTION_DIAGNOSIS.json`에 raw 메시지와 함께 보존되었다.

부모의 전체 frontend/pipeline 실행은 1,415 passed, 9 skipped, 4 subtests passed, 실제 exit 0이며, 실행 중 408개 source byte identity가 안정적이었다. 이는 표현 수정 전 실행이다. 그 뒤 바뀐 부분은 위 수렴 aggregation이며, 해당 arithmetic와 AST 범위만 표적 검증했다. 기존 suite를 반복하지 않은 판단은 변경 범위와 일치한다.

**수용 가능한 결론:** 고정 parameter·고정 source를 이용한 신규 균일 RMF 계산, 유한 admissible component에서의 9개 sharp Maxwell 접합, 기록된 core 격자 민감도, physical energy를 보존하는 serialization과 provenance.

**계속 열려 있는 결론:** 원전 DDME-X NMP identity, 전역 microscopic minimum, 기존 crust calibration의 원전 동등성, 연속 mesh 극한, 제외된 source 구간의 물리성, whole-star/TOV/Cowling 응용. 실패 사례는 이들에 대한 강한 no-go를 증명하지 않는다.

전체 machine-readable 결과는 `data/new_model_maxwell_20261117/INDEPENDENT_REVIEW.json`에 있다. 이 검토가 직접 변경한 파일은 본 보고서와 독립 검토 JSON·스크립트·로그뿐이며, commit/push는 수행하지 않았다.
