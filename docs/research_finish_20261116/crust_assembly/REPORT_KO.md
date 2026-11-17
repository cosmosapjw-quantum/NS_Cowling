# Crust component construction: executed finite-domain results

상태는 **PARTIALLY_COMPLETE**다. 이 작업은 core–crust 결합이나 stellar 계산을 수행하지 않았다. 원래 13개 crust 대상의 변경하지 않은 inner functional을 독립 outer crust에 anchor한 BSk22 analytic representation과 같은 압력에서 비교했다. 네 GPPVA functional에는 실제 재최적화한 구형 TF 상태로 두 분지의 국소 Gibbs 교차가 존재한다. 다섯 SPG에는 별도로 정의한 species-rest-energy-restored 후보를 실제 beta 재최적화하여 교차와 부분 crust 테이블을 만들었다. 이 다섯 후보는 원본 SPG가 아니다.

## 실제 산출물과 범위

9개 native CompOSE 부분 crust 테이블, 총 **6,458행**을 write/readback했다. GPPVA DD2/TW/FSU2R/FSU2H는 각각 685/688/690/690행이며 구형 분지 n_B=.03 fm^-3까지만 포함한다. 복수 해상도 기존 stationary profile과 새 N384 interface를 사용했으므로 완전한 mesh-converged EOS가 아니다. SPG M1–M5 새 후보는 각 741행, n_B=.01 fm^-3까지만 포함한다. 이 cutoff는 core transition이 아니다. 표면 calibration, 고밀도 proton-bearing gas, pasta 및 homogeneous endpoint를 이 새 후보로 재계산하지 않았다.

| 모델 | 선언한 variant | n_lower (fm^-3) | n_inner (fm^-3) | P (MeV fm^-3) | Δμ_inner−lower (MeV) | 출력 cutoff (fm^-3) |
|---|---|---:|---:|---:|---:|---:|
| GPPVA_DD2 | unchanged RMF + analytic lower approximation | 0.001973484536 | 0.002170019775 | 0.002142467788 | -5.37e-10 | .03 |
| GPPVA_TW | unchanged RMF + analytic lower approximation | 0.00031735695 | 0.000572596747 | 0.0005587550856 | -1.12e-07 | .03 |
| GPPVA_FSU2R | unchanged RMF + analytic lower approximation | 0.0002842941979 | 0.0004672243068 | 0.0005188658083 | 1.41e-07 | .03 |
| GPPVA_FSU2H | unchanged RMF + analytic lower approximation | 0.0002728082878 | 0.0004353863803 | 0.0005014178466 | -2.33e-07 | .03 |
| SPG_M1 | species-rest-restored candidate | 0.002015448002 | 0.002163802468 | 0.002192197453 | 4.55e-13 | .01 |
| SPG_M2 | species-rest-restored candidate | 0.001903461566 | 0.001978024664 | 0.002060438248 | 5.68e-13 | .01 |
| SPG_M3 | species-rest-restored candidate | 0.001808943322 | 0.001963605896 | 0.001951639664 | -2.27e-13 | .01 |
| SPG_M4 | species-rest-restored candidate | 0.00142883172 | 0.001543344855 | 0.001537585238 | -1.14e-13 | .01 |
| SPG_M5 | species-rest-restored candidate | 0.001393628708 | 0.001512444822 | 0.00150122617 | 3.41e-13 | .01 |

`assembled_tf/MODEL/ASSEMBLY.json` 및 `restored_spg/MODEL/ASSEMBLY.json`에는 density jump 양측 row index와 `interpolation_across_density_gap=FORBIDDEN`이 저장된다. Native table만 임의의 연속보간기에 전달하면 이 계약이 보장되지 않는다. 모든 사용자는 companion manifest를 읽고 각 분지를 별도 취급해야 한다. 모든 row의 total composition이 계산되어 있지만, 공간 TF의 총 neutron/proton fraction을 free-particle fraction으로 바꾸지 않았다. TF native masks는 free/cluster 분해를 unavailable로 유지하고, `COMPOSITION.json`은 global row segment를 profile 및 total-composition companion에 연결한다. 전용 TF composition reader도 실제 실행했다.

## 물리 reference와 실제 변경

단위는 ε,P: MeV fm^-3, n: fm^-3, μ: MeV다. CompOSE header mn을 바꾸면서 Q3=μB/mn−1, Q4=μQ/mn, Q5=μL/mn, Q6=(F/A)/mn−1, Q7=ε/(n mn)−1로 다시 쓰는 것은 serialization normalization이다. 실제 ε,P,μ 및 조성은 변하지 않는다. 이 변환으로 상대 phase-energy 문제를 없앨 수 없다.

보존량 reference 변환 δε=a nB+b nQ에서는 δμn=a, δμp=a+b, δμe=−b이므로 beta 조건 μn−μp−μe가 불변이다. 중성계에서는 δε=a nB이며 모든 분지에 동일하게 적용해야 상대 μ 차이가 보존된다.

이번 SPG 후보는 그 변환과 다르다. 공통 kinetic mass M0=938.9 MeV를 유지하면서 ε′=ε+(Mn−M0)nn,total+(Mp−M0)np,total을 **조성 최소화 이전**에 추가했다. Mn=939.56542052, Mp=938.27208816 MeV를 사용한다. Electron functional은 그대로이므로 beta 방향 미분이 달라지고 실제 재평형이 필요하다. 이것은 명시적인 phenomenological energy-functional 변형이며 원래 covariant bulk action의 질량/결합 재보정이나 원저자 mass fit의 동등성을 주장하지 않는다. SurfaceParameters는 immutable한 기존 AME2408 candidate prescription으로 고정했고 그 JSON SHA를 보존했다. 이전 소수점으로 반올림한 질량 실험은 `restored_spg_rounded_mass_predecessor`에 보존하며 최종 표의 근거로 혼합하지 않는다.

## BSk22 표현의 지위

독립 outer-drip ε+P=nμ anchor에서 BSk22 published P(n)를 적분하여 d(ε/n)/dn=P/n²를 맞춘 기존 표현을 사용했다. 원래의 독립 P/E analytic fit residual은 별도 보존된다. 새 `bsk22_extended_representation.py`는 기존 `.002` API guard를 변경하지 않고, 별도 variant로 최대 .01까지 동일 published formula를 노출한다. Pearson et al. 및 erratum의 Eq.80/85는 inner-crust composition fit이고, Table8의 BSk22 Z=40 범위는 .034까지, corrected Table12의 proton drip은 .059이다. Zero free-proton convention은 저자 2023 analytic-code 수정에 따른다.

이 확장은 새 ETF/ETFSI 최소화가 아니다. 전체 table에서 analytic composition의 μn와 압력 적분의 μB 사이 residual은 약 −.0215..+.0255 MeV로 남는다. 따라서 열역학적 1차 법칙 및 두 분지 endpoint matching을 통과하더라도 microscopic thermochemical equivalence를 인증하지 않는다. 명시한 finite-error lower-piece approximation이다. SPG M1의 lower seam .002015 부근만 provider의 기존 .002 전환보다 약간 위에 있고, 나머지 구성은 그 아래에서 교차한다.

원전: [Pearson et al. + erratum](https://arxiv.org/pdf/1903.04981), [저자 analytic code 안내](https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html). 실제 수집 내용 및 SHA는 `data/research_finish_20261116/literature/retrievals_extension.json`, source contract는 `docs/research_finish_20261116/CALIBRATION_SOURCE_CONTRACT_KO.md`에 있다.

## 검증과 실패 분류

- **Derived / implementation-verified:** exact header rebase가 dimensional ε,P,μ,F를 보존하고 sparse masks를 유지한다. Interface helper는 finite positive tolerance, 양의 밀도 jump, positive P/μ, endpoint Euler identity와 P/μ matching을 검사한다. NaN/inf tolerance로 검사를 우회하는 결함을 독립 검토가 찾아 수정했고 refusal tests를 추가했다. Neutrality/beta 검증은 upstream physical solver의 계약이며 helper 자체가 sparse TF masks로 추론하지 않는다.
- **Numerically checked:** TF 네 interface의 |Δμ|<2.4e−7 MeV, P 상대 차이<2e−13, 실제 density perturbation의 |dε/dn−μ|<3.9e−6 MeV. N192→384의 E/A 변화 .00024–.00182 MeV를 보존하므로 mesh-converged physical transition으로 해석하지 않는다.
- **Numerically checked:** 다섯 SPG 후보의 3-seed 재최적화, insertion KKT, 65개 inner rows/모델, endpoint/중간/끝점 finite-difference first law. 독립 검토는 전체 325개 inner row의 Hessian을 추가 계산하여 양의 국소 곡률과 artificial-bound 비접촉을 확인했다. 최종 source fit/global minimum을 뜻하지 않는다.
- **Numerically checked:** `local_branch_ordering.json`은 실제 재최적화한 양측 rows에서 Δμ가 아래에서 양수, 위에서 음수임과 양의 dP/dn을 확인한다. dΔμ/dP=1/n_inner−1/n_lower<0과도 맞는다. 이 두 선택 분지 사이의 국소 순서만 확립하며 모든 nucleus/geometry/seed의 전역 lower convex envelope는 증명하지 않는다.
- **Preserved implementation/numerical failures:** 초기 runner의 keyword-only constructor 오호출은 `tf_interfaces_harness_red.log`; 엄격한 SLSQP ftol에서 status8을 낸 첫 시행은 `tf_interfaces_strict_optimizer_red/`와 log에 보존한다. 최종 runner는 실패 상태에서 재시작해 성공 여부와 동일 KKT/radial criteria를 재검사한다. 실패를 낮은 에너지 후보에서 조용히 제거하지 않는다.
- **Preserved cache defect and repair:** 독립 검토가 mutable configuration을 바꾼 뒤 stale cached state가 반환되고 반환 dict를 수정하면 cache가 오염되는 두 반례를 재현했다. 다섯 테스트 실패를 `restored_cache_tests_red.log`에 보존했다. Frozen configuration과 private cache/public deep copy로 수정했고, 채택 상태마다 artificial bound 및 free-coordinate Hessian도 검사한다. 원래 physical objective/optimizer는 그대로이며 `cache_repair_replay.json`은 전체 325개 archived inner row의 물리값을 다시 실행하여 bit identity를 확인한다. 원래 module bytes는 `PRE_CACHE_REPAIR_restored_spg_candidate.py`에 남겼다. 후속 심사에서 marginal t=1 경계의 inward saddle 반례를 찾아 두 테스트 실패를 `restored_critical_cone_red.log`에 보존했다. 최종 구현은 strict g_t<−1e−6인 경우에만 tangent subspace를 쓰며, marginal 경계에서는 one-sided inward t Hessian column까지 검사한다.
- **Current-code verification:** TF root 계산은 frozen SHA `11af1e007d997d9eea7465c4d42552440544e26f7070b0b3bd411f13990cb3f9`에서 수행했다. 이후 final electron SHA `8be3c6756f1b04a38cb59fdeeacb5c0aed5c98d8337c8ed14c64e7091501bdca`로 같은 profile을 재평가했고 |Δε|≤4.5e−16 MeV fm^-3, |Δμ|≤2.4e−12 MeV였다. 재최적화를 했다고 주장하지 않는다. 양쪽 exact SHA와 제약검사 결과는 `final_operator_profile_verification.json`에 저장했다.
- **Implementation-verified:** 신규 표적 테스트 **58 passed**. Native9 tables 및 TF profile reader 실제 실행. 독립 검토 상세는 `data/research_finish_20261116/independent_review/CRUST_TABLE_REVIEW.json`를 확인한다(부모 보고서의 최종 경로가 우선).

## 변경하지 않은 원래 13개에 대한 결과

`initial_original_gap_scan.json`은 BSk n≤.002, `extended_original_gap_scan.json`은 별도 published-fit extension n≤.01에서 수행했다. All13의 실제 기존 inner samples를 P에서 보간한 탐색이므로 무교차의 전역 증명은 아니다. SPG5, DDME2, FSU2, NL3wrL55, TM1e는 이 범위에서 μinner−μlower가 계속 음수이며 단순히 .002 이후로 BSk 범위를 늘려도 sampled crossing이 생기지 않았다. DD2/TW/FSU2R/FSU2H의 탐색 교차는 후속 실제 cell 재최적화로 검산했다.

이 결과는 원래 9모델을 "물리적으로 불가능"으로 분류하지 않는다. Provider lower-piece prescription과 absolute-energy functional의 관계, 원래 저자의 실제 matching/reference convention, candidate-model을 선택할지 여부를 분리해야 한다. Arbitrary relative energy offset으로 이 차이를 없애거나 .002에서 raw append하지 않았다. Source fidelity와 독립 thermodynamic implementation은 서로 다른 심사 항목이다.

## 재현 및 다음 최소 계산

저장소 루트, PYTHONPATH=python_frontend:scripts, OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1을 사용한다. Native export에는 COMPOSE_EOS_FORTRAN_ROOT=. 및 COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran을 지정한다.

```sh
python3 scripts/solve_crust_tf_interfaces.py
python3 scripts/build_restored_spg_crust.py
python3 scripts/assemble_tf_partial_crust.py
python3 scripts/verify_crust_interface_profiles.py
python3 scripts/check_crust_interface_ordering.py
python3 -m pytest -q python_frontend/tests/test_crust_assembly.py python_frontend/tests/test_bsk22_extended_representation.py python_frontend/tests/test_restored_spg_candidate.py
python3 scripts/seal_crust_assembly.py
```

새 candidate를 full crust로 확장하려면 exact-rest-energy functional에서 proton-bearing gas, surface alternatives, homogeneous endpoint까지 다시 최소화해야 한다. GPPVA는 mesh를 통일하고 geometry 전환을 실제 P/μ matching과 양측 profile로 계산해야 한다. Common analytic lower piece의 composition chemical residual은 full ETFSI 또는 명시적 일관성 모델로 닫아야 한다. 현재 source-original 원본 13개 전체의 독립 재생성 완료를 선언할 수 없다.

코드/입력 snapshot, native executable SHA 및 최종 산출물 SHA는 `PROVENANCE.json`과 `FINAL_SOURCE_INPUT_SNAPSHOT.json.gz`에 있다. 모델별 `SHA256SUMS.txt`는 native table, manifest 및 profile/readback receipts를 함께 묶는다.
