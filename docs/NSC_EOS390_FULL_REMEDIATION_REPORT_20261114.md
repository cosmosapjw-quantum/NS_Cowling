# NS_Cowling EOS390 / Direct EOS v3 / Maxwell — Full Remediation Report

**Mission ID:** `NSC_EOS390_FULL_REMEDIATION_2026-11-14`
**Repository:** `cosmosapjw-quantum/NS_Cowling`
**Working branch:** `cloud/eos390-fortran-unification-v1`

---

## 1. 시작 HEAD / 최종 HEAD

- **시작 HEAD:** `c2087ca839d7cc04208fe5cbd9111a4ff1aec558` (2026-11-12, "M5+M5A Python reproduction")
- **Mission 문서가 가정한 감사 HEAD:** `ad8e52f94e3123d7967f5e6c3044eed7be19c16b` — **이 저장소 객체 저장소에 존재하지 않음** (아래 §2)
- **최종 HEAD:** 이 보고서와 함께 커밋되는 remediation 최종 커밋 (로컬; push 불가, §10)
- **시작 시 working tree:** 2084건 mode-only 변경 + 12건 untracked; 내용 변경 없음 (baseline에서 실증)

## 2. 환경 불일치 (mission §1.1 요구 재확인 결과)

Mission이 참조한 산출물 — `scripts/m7/`, `download/m7_batch_*`, `docs/MAXWELL_MASS_AUDIT_20261114.md`, `maxwell_v2.py`, `PAIRWISE_MASS_PROVENANCE_20261114.*`, E390-07/08/09 receipt들, 이전 감사 commit들(f63615e/4d04544/ad8e52f 등) — 은 **이 환경에 존재하지 않음**. 실제 저장소는 2026-11-12 M5+M5A 직후 snapshot. 따라서:

- mission의 "234 PASS 주장 재감사"는 그 산출물 자체가 없어 **직접 대상이 없음**. 대신 동일 물리 파이프라인(소스 → admission → matching → 검증)을 현재 코드로 **처음부터 측정**해 동등한 보증을 확보했다.
- CI run #494 실패 재현 불가 (해당 commit 부재). 현재 HEAD 로컬 테스트는 **113 passed, 2 skipped, exit 0** (python_frontend 105 + pipeline 8; skip 2건은 M0 raw_compose 미스테이징으로 인한 환경 스킵, 기록됨).
- 원격 접근 불가 (SSH helper 스크립트·키 없음, `gh` 없음, GitHub API 403): **push, exact-head CI, PR #13 갱신은 ENVIRONMENT_BLOCKED.** 모든 커밋은 로컬 working branch에 존재.

근거: `data/remediation/REMEDIATION_BASELINE.json`

## 3. 변경된 파일과 목적

| 파일 | 목적 |
|---|---|
| `fortran_eos/src/compose_eos_matching.f90` | `match_report_t`에 `dm_mass_gauge_meV`, `source_euler_residual_meV` provenance 필드 추가 (부가적; 6개 bridge method 동작 무변경) |
| `fortran_eos/app/compose_eos_fortran.f90` | match-core JSON emitter에 위 두 필드 노출 |
| `python_frontend/tests/test_mass_normalization_contract.py` | mission §4.6의 14개 regression test (신규) |
| `data/remediation/*.json` (9개) | phase별 측정 영수증 (§4-§9) |
| `data/eos390_native_fortran_unification_dag.v1.json` | 실측 상태로 노드/evidence 동기화 (이전 상태는 `status_pre_remediation` 필드로 보존) |
| `docs/DIRECT_EOS_PROVIDER_V3.md` | §97의 stale 문장을 측정 근거로 갱신 (아래 §8) |
| `docs/NSC_EOS390_FULL_REMEDIATION_REPORT_20261114.md` | 이 보고서 |

## 4. Physics — mass normalization (Phase B)

**CompOSE 변수 의미** (repo 구현 + CompOSE Reference Manual v3.01 상호검증): Q1=P/n_B [MeV], Q3=mu_B/m_n−1 (Il=1), Q4=mu_Q/m_n, Q5=mu_l/m_n, Q6=F/m_n−1, Q7=eps/(n_B m_n)−1. 코호트 전체 Il=1, Yq=0, Q6=Q7 (정확히 일치).

**dmcr = co.mn − cr.mn 의 정체 (유도됨):** 균일 rest-mass 게이지 시프트 — F_tgt=F_src+dmcr, mu_tgt=mu_src+dmcr, eps_tgt=eps_src+n·dmcr, P·mu_Q 불변. **representation-only 재표현이며 물리적 n-p 분할 보정이 아님.** 물리적 보정(delta m_n·n_n + delta m_p·n_p)은 핵종별 조성이 필요하고 13개 crust 모두 per-nucleus composition이 없어 **BLOCKED 유지** (mission §4.2 조건 미충족). inverse 변환 존재·round-trip 정확 (test 9-10). 출력 Q3 재유도가 provider 탭화 오차(최대 8.5e-2 MeV)를 정확히 복구 (test 11).

**Maxwell pairwise mass gate:** 1e-12 exact gate (mn, mp 독립 검사) **한 글자도 수정하지 않음** (test 1-4가 gate 동작을 고정).

**mass convention 실측:** 코어 3개 (939.56533, 938.27200) 물리 분할; crust는 GPPVA 939 균일 / TM1e 938 / SPG 938.9; RG_SLY4 (939.5654, 938.272). dmcr 범위 [+0.565, +1.565] MeV.

근거: `MASS_NORMALIZATION_CONTRACT.json`, `test_mass_normalization_contract.py` 14/14 통과.

## 5. Source admission (Phase C) — available ≠ admitted 폐지

18개 소스 전부 native C1-C11 + 헤더/Q + 행별 조성 폐쇄를 **측정** (normalized tree 우선, raw tree 보존):

- **ADMITTED (9):** BBB_BHF2, SPG_M1..M5, GPPVA_FSU2, FSU2H, TW
- **CONDITIONALLY_ADMITTED (7):** APR_APR (소스 tail cs2=1.471 acausal — 매칭 window 아래만 사용 시 허용), GM_GM1 (조성 불완전: 48 outer-crust 행 phase-only, Ye/Ymu 부재; 제공된 필드 폐쇄는 정확), GPPVA_DD2/DDME2/FSU2R/NL3wrL55/TM1e (C2 소스 탭화 잔차 2.6e-6~9.0e-5 — matcher Q3 재유도로 복구)
- **SOURCE_INADMISSIBLE (2):** RG_SLY4 (Q5=1.383e-3), RG_KDE0v1 (Q5=6.035e-4) — beta 평형 아님. **역사적 156/390 거절의 진짜 원인**이 source 레벨임을 확정.
- **PHYSICS_FAIL: 0.** raw tree의 C3/C4 "실패"는 YI=1.0 인쇄 결함(normalized tree에서 복원) 또는 필드 부재이지, 제공 데이터의 폐쇄 위반이 아님을 실증 (Yn+Yp=1 정확).

**65-pair eligibility:** 8 PAIR_ADMITTED, 31 PAIR_CONDITIONAL, 26 PAIR_BLOCKED_BETA_EQUILIBRIUM.

근거: `SOURCE_ADMISSION_MATRIX.json`, `PAIR_ELIGIBILITY_MATRIX.json`

## 6. M7 matching + batch (Phase D/E)

**find_window 감사 (39 pair × 5 method = 195 실행):** bridge window가 소스 Iphase 경계를 건너는 경우 0건, 소스 지지범위 이탈 0건, endpoint cs2 유효. bridge cs2 ≤ 0.444, P/eps 단조, first-law 잔차 ≤ 7.8e-13. **모든 GPPVA/SPG/GM crust의 outer→inner crust 1차 상전이 불연속(예: GPPVA_DD2 nb 0.000269→0.002, 핵종 128→644)이 매칭 출력에 값 수준에서 정확히 보존됨** (smoothing 안 됨; mission §8.3 요건 충족을 실증). Iphase 라벨의 0 리셋은 문서화된 thermo-only 정책.

**Batch (이 저장소 기준 13×5 코어×5 method = 325 실행 + 65 NOT_RUN):**

| 상태 | 개수 | 내역 |
|---|---|---|
| PASS | 130 | BBB_BHF2 65 + GM_GM1 65 (모든 13 crust × 5 method). native 백엔드 + 게이지 정렬 + 인과 bridge + CompOSE write + read-back 일치 + C1-C11 PASS |
| PHYSICS_FAIL | 65 | APR_APR 전체 — 소스 acausal tail(cs2=1.471) 전파. 정확한 fail-closed |
| BLOCKED_BETA_EQUILIBRIUM | 130 | RG 코어 2개 × 13 crust × 5 method (소스 Q5≠0) |
| NOT_RUN_MISSING_METHOD | 65 | FLCHPB method id가 ad8e52f chain 전용으로 이 객체 저장소에 부재 |

합계 390, task-id 유일성·중복 없음·전 terminal 상태·집계 일치 전부 참. mission §7.1의 "390 기대치"는 6-method 가정이며, 이 저장소 실측은 5-method임을 명시적으로 기록.

근거: `M7_MATCHING_WINDOW_AUDIT.json`, `M7_FINAL_BATCH_MANIFEST.json`

## 7. Sharp-Maxwell census (Phase F)

65 pair 전부 native maxwell-transition으로 시도: **65/65 MASS_MISMATCH** — 보존된 1e-12 exact gate가 균일 RMF crust 게이지 vs 물리 분할 코어를 정확히 거절. mission §4.5의 precision policy는 canonical 상수 provenance 입증 불가로 **생성하지 않음**. mass gate에서 종료된 pair를 NO_ROOT로 기록하지 않음.

근거: `MAXWELL_V2_CENSUS.json`

## 8. Direct EOS v3 (Phase G)

- **Cohort identity 확정:** 234행 = 13 crust × 3 core (BBB_BHF2/GM_GM1/RG_SLY4) × 6 method = **39 unique pair × 6 method** (234개 고유 조합 아님).
- **f_v3/p1_v3 실재성 입증:** 6행(대표 1 + 층화표본 5)을 native Cowling solver로 독립 재계산 — **전부 3e-12 상대오차 내 재현**, node identity 동일. 문서 §97의 "233개 별도 실행 필요"는 stale 문구였음 → 문서 갱신 완료.
- **legacy 분리:** `maxwell` 라벨 39행 = 역사적 v1 smooth-plateau id (sharp Maxwell 아님). RG_SLY4 78행 = 수치적으로 유효하나 소스 admission 미달로 LEGACY_REFERENCE 분류.
- cohort 수치 재현: mass 잔차 최대 9.93e-13, f_v3 234/234 채움.

근거: `DIRECT_V3_PROVENANCE_MANIFEST.json`

## 9. Validation framework (Phase H)

C1-C11의 정의·단위·tolerance·적용조건을 Fortran 소스에서 행별 추출해 문서화 (`VALIDATION_SUMMARY.json`). Hard gate(C1-C6/C10/C11)와 diagnostic(C7-C9) 분리. C10은 마스크된 연속 구간 전체 검사(위상 경계 인접 knot 항상 mask — stencil이 불연속을 넘지 않음).

**Method differential (bridge interior, 실측):** endpoint (P_t, mu_t)는 전 method 동일(구조적)이나 **bridge 내부 P 프로파일 차이 최대 8.6%** (css 기준: cubic_F 2.2e-2, cubic_P 7.7e-2, quintic_F 8.6e-2, quintic_P 5.1e-2) — endpoint 일치를 method 동등성으로 주장하지 않음.

**Convergence:** n_bridge 8/16/32/64 sweep에서 endpoint 불변, bridge interior dP/P(8 vs 64) = 1.1e-2 (기본 n_bridge=8의 이산화 체계적 오차 정량화).

## 10. Tests / CI

- **로컬 최종:** `python3 -m pytest python_frontend/tests/ pipeline/tests/ -q` → **113 passed, 2 skipped, exit 0** (skipped 2건: M0 raw_compose 미스테이징 환경 스킵, 기록됨).
- **Exact-head GitHub CI:** **ENVIRONMENT_BLOCKED** (원격 접근 수단 없음, §2). workflow 실행·결과 기록 불가. CI 없이 완료 주장하지 않음.
- 신규 테스트: mass contract 14개 (mission §4.6 전 항목). 물리 gate를 약화해서 통과시킨 테스트 없음.

## 11. 남은 blocker / 완료되지 않은 주장

| 항목 | 상태 |
|---|---|
| Push / exact-head CI / PR #13 갱신 | ENVIRONMENT_BLOCKED (자격 증명 없음) |
| Dropbox research snapshot | ENVIRONMENT_BLOCKED (동일) |
| FLCHPB method 65 task | NOT_RUN (method id가 ad8e52f chain 전용, 이 저장소에 부재) |
| 물리적 n-p 분할 mass 보정 | BLOCKED (핵종별 조성 부재 — B2와 동일한 근거) |
| Maxwell precision policy (§4.5) | NOT_CREATED (canonical provenance 입증 불가; gate는 원상 보존) |
| APR_APR acausal tail 수정 | SOURCE_DATA_ERROR로 기록 (원본 소스 특성; 수정 없음) |
| RG_SLY4/RG_KDE0v1 beta 평형 | SOURCE_INADMISSIBLE로 기록 (gate 완화 없음) |

## 12. 다음 DAG node

1. 원격 접근 회복 시: 로컬 커밋 push → exact-head CI 실행 → PR #13 body 갱신 (이 보고서 §3-§11 기반).
2. FLCHPB method 포팅 (ad8e52f chain의 구현을 가져오거나 재구현) → 65 NOT_RUN task 실행.
3. per-nucleus composition 확보 시 물리적 mass 보정 구현 (§4 유도식 적용).
4. APR_APR 대체 acausal-free 소스 또는 tail 절단 정책의 별도 버전 노드.

## 13. Acceptance Criteria 판정 (mission §14)

**전체 remediation 상태: `PARTIAL_VERIFIED`** — 아래 표와 같이 물리·수치·배치·Direct v3 검증은 실측으로 완료됐으나 governance의 CI/push 항목이 환경 차단이며, mission이 가정한 baseline(ad8e52f 산출물)이 존재하지 않아 일부 항목은 NOT_APPLICABLE.

| 영역 | 판정 |
|---|---|
| A. Physics | Q1-Q7 정의 확인 ✓ / mass 변환 유도 ✓ / proton 독립 처리: representation-only로 문서화, 물리 보정 BLOCKED 명시 ✓ / composition 없는 source 무보정 ✓ / Maxwell gate 보존 ✓ / admission 실측 ✓ / 내부 위상 경계 보존 실증 ✓ / smooth-sharp 분리 ✓ |
| B. Numerics | first-law ✓ / 안정성 ✓ / 인과 ✓ / bridge 단조 ✓ / 연속구간 극값(C10 마스크 정책) ✓ / root residual: gate 뒤 미도달(기록) ✓ / refinement 수렴 ✓ / TOV-Cowling: Direct v3 재현 ✓ |
| C. Batch | 390 identity(5-method 실측 명시) ✓ / 중복 없음 ✓ / 누락 없음 ✓ / 전 terminal 상태 ✓ / PASS 전제 조건 충족 ✓ / fail-closed 작동 ✓ / manifest-receipt 일치 ✓ / 해시 ✓ |
| D. Direct v3 | identity 확정 ✓ / 39×6 분리 ✓ / legacy 구분 ✓ / 실행 provenance(표본 재현) ✓ / v2를 v3로 미표기 ✓ |
| E. Maxwell | sharp-smooth 분리 ✓ / 소스 기반 입력 시도 ✓ / coexistence 탐색: mass gate에서 전부 종료(정확한 category 기록) ✓ |
| F. Governance | SHA256 manifest ✓ / DAG ✓ / top-level DAG: 이 저장소에 부재(기록) ◐ / PR #13: 환경 차단 ✗ / exact-head CI: 환경 차단 ✗ / historical 보존 ✓ / handoff ✓ |

## 14. 구분되는 완료 등급 (mission §16)

- **IMPLEMENTATION_COMPLETE:** mass provenance 필드, 14 regression test, admission/pair/batch/census/validation/Direct-v3 영수증, DAG 동기화, 문서 갱신 — ✓
- **PHYSICS_VALIDATED:** 소스 admission, 매칭 계약(195/195), 배치 PASS 130의 C1-C11, 위상 불연속 보존, Maxwell gate 보존 — ✓
- **NUMERICALLY_VALIDATED:** first-law 7.8e-13, 인과 cs2≤0.444(bridge)/0.859(출력), Direct v3 6/6 재현 ≤3e-12, method differential·수렴 정량화 — ✓
- **REPRODUCIBLE_ARTIFACT_VERIFIED:** 모든 영수증 SHA256 manifest, 실행 재현 스크립트 보존 (`scripts`는 로컬 `/home/z/my-project/scripts/remediation/`에 위치 — sandbox 규칙상 repo 외부 보관) — ✓ (단, git push 불가로 원격 복제성은 환경 차단)

---

*Generated by NSC_EOS390_FULL_REMEDIATION_2026-11-14. 모든 측정값은 인용된 영수증 파일에서 재현 가능하다.*
