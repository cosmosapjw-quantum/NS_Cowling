# NS_Cowling: Fortran-native DVODE 성능·정확도 연구 루프 C1–C10

작업일 2026-11-18. 기준 소스: `990d07f4f0383f3777f65ed8d3454263b09e14f2`, 격리 브랜치: `research/fortran-dvode-perf-accuracy-20261118`.

## 정확한 주장 범위

**Implementation verified (scoped):** EOS `P,epsilon` 입력 평가, TOV, `l=2` tidal Love `y,k_2,\Lambda`, Maxwell sharp interface, 다중 중심압력의 batch 실행이 모두 standalone Fortran backend에서 진행된다. Python `native_tov_love.py`는 입력 파일·실행 환경·프로세스 호출·JSON parsing만 처리한다.

**Not migrated:** `python_frontend/compose_eos/stellar_background.py`, `stellar_cowling.py`, `stellar_sharp.py`의 기존 Python/DOP853 고유모드 `f,p_1` solver는 이번 구현에서 수정하거나 Fortran으로 대체하지 않았다. **NS_Cowling 전역이 Fortran-only라는 주장은 부적절하다.** 별의 nonradial `f,p_1` Fortran 이전은 `COWLING-NATIVE-01`이라는 별도 남은 코드 작업이다.

이 결과는 기존 EOS의 핵물질 functional·물성·원본 수치 데이터를 변경하지 않는다. `tov-love`는 EOS 입력에서 주어진 `P,epsilon` 표를 평가한다. 본 성능 자료의 4097행 PK1 표는 **저밀도 `n_B=0.015 fm^-3`에서 잘린 균일물질 core 모델**, 즉 NS surface까지 완비된 realistic crust+core 모델이 아니다. 여기서 얻은 `M,R,k2,\Lambda`는 **동일 입력·동일 물리식 알고리즘 검증 목적**이며 관측 예측값으로는 승인하지 않는다.

## 10개 코딩 단계: 실제 수행 내용과 acceptance 증거

| 단계 | 실행 내용 | Evidence·실제 상태 |
|---|---|---|
| **C1 기준 잠금** | 정확한 source HEAD, Fortran 13.3.0, 하드웨어 24 logical CPUs, PK1 native `eos.thermo` SHA와 4097/1025/257 동일 소스 표 잠금 | `C1_C4_PLAN_KO.md`, `BENCHMARK.json`. **완료** |
| **C2 수식·호출경로 감사** | `tov_rhs`의 잘못된 `(ε+p)c_s²`, `k2`의 1/2 계수(정답 8/5), RK4(5) stage times·누락된 k7·잘못된 trial acceptance, O(N) EOS bracket 탐색 발견. Python Cowling DOP853 별도 노출 확인 | Fortran 소스·baseline·독립 red log. **완료** |
| **C3 수치설계** | EOS 괄호 검색 O(log N), Adams MF10/BDF MF22의 적응형 IVP, 실제 수용/거절 RK45와 마지막 단계 기하학, surface dense-event interpolation, 각 Maxwell phase 분리 및 Love y-jump 보존 | Netlib 원전 `vode.f` 및 프로젝트 native Fortran. **완료** |
| **C4 실행/자원 DAG** | compiler source-signature build, 독립 DVODE process concurrency 최대6, 한 프로세스의 batch는 직렬, python은 frontend로만 사용. Netlib F77 COMMON solver state를 하나의 프로세스 내 thread 간 공유하지 않음 | `PARALLEL_STRESS.json` 및 build signature. **완료** |
| **C5 실제 코드 구현** | `compose_eos_tov.f90`의 Hinderer Love 물리식·DOPRI 오류 수정, binary bracket, native Netlib DVODE+LINPACK 원본 link, `tov-love batch-smooth` native endpoint, `NativeTovLove` Python I/O facade | 실제 Fortran 2 binary build 성공, source diff + vendor hash. **완료** |
| **C6 실행·계약검증** | 컴파일, native `tov-love` smooth/Maxwell, BDF/Adams/RK45, invalid input fail-closed, batch 4개 exact, tolerance controls | `test_native_tov_*_20261118.py`, Fortran `test_tov_love` native logs. **완료** |
| **C7 과학 검증** | Hinderer `k_2` 정규화 테스트는 구버전 실제 실패(0.3125×); 수정 후 일치. `Q(r)`의 독립 원전식·Fortran 구현 상대차 약 `1.12e-13`. Maxwell energy jump, (\Delta y=-4πr_t^3\Delta\epsilon/(m_t+4πr_t^3P_t))와 one-sided source admission 통과 | `TIDAL_RHS_PHYSICS.log`, `TOV_LITERATURE_RED.log`, Fortran native unit tests. **완료** |
| **C8 실제 정확도·성능** | 동일한 4097/1025/257 PK1 원자료 7회 측정, Adams/BDF/교정 RK/legacy 분리. 32개 native batch 결과 직렬값과 필드별 완전 일치. 6개 독립 Fortran 프로세스 병렬과 직렬 출력도 bit-level numeric identity | `BENCHMARK.json`, `BATCH_BENCHMARK.json`, `PARALLEL_STRESS.json`. **완료** |
| **C9 적대적 검토** | 물리·수치·구현·source·runtime 분리 감사. Missing tolerance, nonfinite, out-of-range branch, wrong density jump, wrong Love formula, parallel process isolation, false-publish claim 테스트. 불필요한 재귀 감사 없음 | `INDEPENDENT_AUDIT.json` 71 scoped checks PASS. **완료 (bounded)** |
| **C10 판정·인계** | 동일 commit에서 provenance·reproduction·known limitations와 scope를 보존, nonforce git push 및 remote SHA 일치로 마감 | `DELIVERY_RECEIPT.json`, 이 보고서, final full suite+push verification receipts. **Scoped TOV/Love release; Cowling migration separate** |

## C7 Love 방정식 수정의 물리적 의미

Metric convention ((-+++)), 내부 TOV 적분은 (G=c=1), 길이는 km, 압력·에너지밀도 입력은 MeV fm⁻³. 원래 구현의 (Q(r)) 항은 ((\epsilon+P)c_s^2)를 사용했는데, Hinderer 계열의 표준식은 (d\epsilon/dP=c_s^{-2})이므로

\[
Q(r)=4\pi e^\lambda\Bigl(5\epsilon+9P+(\epsilon+P)\frac{d\epsilon}{dP}\Bigr)
-\frac{6e^\lambda}{r^2}-(\nu')^2 .
\]

별표에서 사용하는 (k_2) 정규화는

\[
k_2=\frac{8}{5} C^5(1-2C)^2\,
\frac{2C(y_R-1)-y_R+2}{D(C,y_R)},\qquad
\Lambda=\frac{2k_2}{3C^5}.
\]

기존 Fortran 함수는 앞 계수를 (1/2)로 사용해 같은 (C,y_R) 기준 (k_2)에 정확히 (5/16=0.3125) 배율 오차가 있었다. Q 자체도 틀렸으므로 **기존 실제 (k_2,\Lambda)와 교정 값의 차이는 단순 3.2배 스케일링이 아니다.** Hinderer 정규화의 native red-to-green test와 독립 (Q(r)) 도함수 테스트를 수행했다. 이전 발표 수치를 자동으로 변경하지 않고 **re-evaluation required**로 표시한다.

## C8 성능: 직접 측정한 Fortran 실행

측정: 4097행 동일 PK1 native barotrope, 동일 `P_c=73.52922470749803 MeV fm⁻³`, 동일 CPU, 각 7회 별도 native executable invocation 중앙값.

| Backend | Median wall time [s] | 비교 의미 |
|---|---:|---|
| 구형 잘못된 Fortran RK | 0.823782 | 이전 시간 baseline, **물리식 결함 있으므로 정확도 공정 비교 아님** |
| 교정한 Fortran DOPRI RK45 | 0.046866 | 같은 교정 물리식의 independent numeric comparison |
| **Fortran DVODE Adams MF=10** | **0.016684** | 기본 solver; corrected RK보다 약 **2.81×** 빠름 |
| Fortran DVODE BDF MF=22 | 0.021336 | 강성 대안; Adams보다 느리지만 같은 물리량으로 수렴 |

구형 속도 대비 Adams는 약 **49.4×** 빠르다. 이 49.4×는 물리 정확도가 다르므로 마케팅성 정확도 동등 비교로 사용하지 않는다. 교정된 RK45 대비 비교가 이번 과제의 적합한 동일-물리식 speedup이다. 실제 순서는 4097행 바탕으로 관찰한 것이고 다른 EOS·컴파일러에서의 보편적 보장은 아니다.

동일 EOS 32개 별 TOV sequence는 개별 호출 중앙 `0.53293s`, native batch 중앙 `0.19263s`, **2.77× speedup**. 32개별 `R,M,k2,Lambda,y_surface,n_steps` 수치 결과가 정확히 동일함을 직접 대조했다. 별개의 DVODE 프로세스 6개를 병렬로 실행한 검증에서 직렬과 물리량이 수치적으로 정확히 같았고 속도는 `0.09888/0.02029≈4.87×`; 이 짧은 wall-clock 시간은 반복 측정 통계가 아니라 하나의 동시성 실험이다.

**EOS 해상도 수렴:** 257/1025/4097행은 동일 원 모형의 **다른 sampled barotropes**이다. 행 수에 따른 물리량 차이는 EOS 표 본래의 수치표현 오차를 포함하므로 `solver convergence`와 합산하지 않았다. DVODE 허용오차 2e-7/2e-9/2e-11의 독립 변화 및 Adams/BDF/RK 비교가 solver-error 검증이다.

## C9 중요 결함과 경계

- **과학적 오류 (교정됨):** Love Q의 (d\epsilon/dP) 누락, (k_2) 정규화, 잘못된 DOPRI stage-time 2·3·4·5와 7번째 embedded error 누락.
- **수치오류 (교정됨):** RK trial-step 무조건 수용, surface cutoff root overshoot, O(N) EOS pressure bracket 탐색; DVODE 옵션은 원전 의미대로 MF=10·MF=22를 적용.
- **Runtime/외부 코드:** 첫 DVODE link는 BLAS와 LINPACK 의존성이 빠져 실패했고, Netlib 원본 `vode.f` 및 `dgefa/dgesl/dgbfa/dgbsl`을 정확한 SHA로 고정해 해결했다.
- **Remaining scientific scope:** Python `stellar_cowling.py`의 (f,p_1) 계산은 native로 이식하지 않았다. Bulk equilibrium `\Gamma_1=\Gamma_{eq}`, slow/rapid phase conversion, elastic crust, frozen-species perturbations 및 full GR oscillations는 이 TOV/Love 성능 변경으로 검증되지 않았다.
- **Source domain:** benchmark PK1 table은 완전한 crust+core EOS가 아니라 homogeneous high-density source. 따라서 반경·질량의 절대값에 NS surface physics 보증을 부여하면 안 된다. Maxwell 테스트는 causal two-branch analytic fixture로 sharp jump implementation을 검증한 것이고, 7개 실제 matching의 global phase-stability를 새로 증명한 것이 아니다.
- **Speed claim:** 구형 49.4×는 wrong physics baseline, 교정된 RK와 2.81×가 동일-물리식 비교. DVODE와 RK가 완전히 동일한 제어 예산을 쓴 것은 아니며 실제 물리량 오차 허용 범위에서 비교했다.
- **Subagent:** host native subagent spawn tool unavailable; code executed in single parent research loop, six native scientific subprocesses are *not* called subagents.

## Reproduction

Native Fortran compile from repo worktree, with source signature guard:

```bash
fortran_eos/build.sh --ensure build/compose_eos_fortran
NS_COWLING_TOV_INTEGRATOR=dvode_adams NS_COWLING_TOV_RTOL=2e-9 \
  build/compose_eos_fortran/compose_eos_fortran tov-love smooth P_eps.csv P_c
NS_COWLING_TOV_INTEGRATOR=dvode_bdf \
  build/compose_eos_fortran/compose_eos_fortran tov-love smooth P_eps.csv P_c
NS_COWLING_TOV_INTEGRATOR=dvode_adams \
  build/compose_eos_fortran/compose_eos_fortran tov-love batch-smooth P_eps.csv centers.txt
NS_COWLING_TOV_INTEGRATOR=dvode_adams \
  build/compose_eos_fortran/compose_eos_fortran tov-love maxwell high.csv low.csv Pt delta_epsilon P_c
```

The accepted Python frontend is `python_frontend/compose_eos/native_tov_love.py`: subprocess/I/O without Python SciPy ODE. Fortran original external dependencies are documented under `fortran_eos/vendor/odepack/PROVENANCE.md`. Reproducible raw evidence is in `data/fortran_dvode_perf_20261118/{BENCHMARK.json,BATCH_BENCHMARK.json,PARALLEL_STRESS.json,INDEPENDENT_AUDIT.json,FULL_PYTEST.log}`.

Primary physical and numerical references:
- Hinderer, *Tidal Love Numbers of Neutron Stars*, https://arxiv.org/abs/0711.2420.
- Takátsy & Kovács, sharp density-jump tidal correction, https://arxiv.org/abs/2007.01139.
- Netlib DVODE Fortran original https://netlib.org/ode/vode.f ; exact SHA in vendor/PROVENANCE.md.
- LLNL ODEPACK method overview https://computing.llnl.gov/projects/odepack.

**Next science/code ID:** `COWLING-NATIVE-01`: native Fortran (\ell=2) Cowling background + ODE eigenvalue f/p1 shooting, including normalization, node labels, slow/rapid Maxwell interface and mode residuals. Must be implementation- and independently benchmark-verified before replacing the existing Python reference. Not claimed or shipped in this branch.
