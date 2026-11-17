# A. Whole-crust 연결 실행 기록

기준 commit: `8c917baf3f5975a6109c3b6e59177079c285c4e0`.
권위: 현재 사용자의 A 실행 요청, `docs/matched_eos_reference_20261117/PLAN_KO.md`의 A–D 순서.

## 이번 범위

공개 BSk22 outer crust와 기존 코드 생성 SPG_M1–SLy4 quintic_P 결과를 연결하는
별도 `reference-assisted hybrid` baseline을 실제로 생성한다. M1 CSS 및 surface
gamma 변경은 구성 민감도 대조다. 기존 세 계열 195개 sharp-Maxwell census는 보존한다.
최소 완료 조건은 표면 처방, outer/inner/core 출처, 경계 조건, 연속 구간의 열역학 검증,
CoMPOSE readback을 갖춘 하나 이상의 사용 가능한 baseline이다.

## 실제 진단으로 정한 방법

1. 공개 outer 원본의 에너지·압력·밀도를 모두 고정한 cubic P(mu_E) 표현은
   216개 interval 중 40개에서 양의 압축률에 실패했다. 원자료와 실패를 보존한다.
2. A에 한정해 별도 pressure-primary outer representation을 사용한다. n/P knot는
   유지하고 첫 원본 energy/baryon에 고정한 제1법칙 적분으로 energy를 재구성한다.
   이는 앞선 energy-primary reconstructed family와 별도 variant다.
3. 원본의 17개 isotope switch는 명시적인 tabulation-gap regularization으로 표시한다.
   그 구간 내부의 composition은 미정이며 sharp nuclear transition으로 인증하지 않는다.
4. outer의 마지막 실제 knot와 native inner의 가능한 knot 사이에 causal convex
   P(mu) bridge를 구성한다. 기존 source energy를 상대 이동하지 않는다.
5. 첫 outer knot 아래는 gamma=5/3의 analytic surface polytrope를 사용하며 gamma=4/3
   및 positive-density export cutoff를 바꿔 surface 처방의 영향을 계산한다.

## 실제 생산에서 확인한 조성 규약 보정

공개 원본의 첫 217개 phase-0 행은 free neutron/proton/muon이 없고 하나의 정수 A,Z
핵종만 있다. 그러나 `eos.compo`의 nuclear group 마지막 값은 모두 1이다. CoMPOSE
매뉴얼의 number-fraction 규약으로 읽으면 Yh=A가 되어 바리온 합이 맞지 않는다.
원본 hash, 단일 group inventory, field presence, Ye≈Z/A를 먼저 확인한 뒤 이 범위만
Yh=1, 표준 nuclear number fraction Yi=1/A로 명시적으로 변환한다. 원본 Ye의 인쇄
잔차 최대는 약 5e-6이며 원본 Ye 자체는 바꾸지 않는다. 이 adapter는 열역학량을
바꾸지 않으며, energy 재구성과 별도 행별 ledger에 기록한다. Inner-crust group의
의미를 이 결과로 추정하지 않는다. 최초 실패한 생산 경로도 별도로 보존한다.

## 구현과 독립 검증

새 surface, outer reconstruction, lower bridge 행은 standalone Fortran kernel로 계산한다.
Python의 독립적인 해석적 표현은 연속 평가, 구간 극값, first law 및 native 대조에 사용한다.
기존 matching/generator의 유효 실행 기록은 identity를 확인해 재사용한다.
공개 원본, 수정량, 미사용 source prefix, 각 row origin, 방법별 실패를 보존한다.

실행 순서: 입력 진단 → 구성식/구현 → native 생산 → 연속/직렬화 검증 → 산출물/게시.
다음 원래 순서는 **B. 필요한 core 지원 확장과 source-variant 정책**이다.
B는 앞선 세 core의 0.02–1.0 fm^-3 생성 결과를 재사용하고, A baseline에 실제 필요한
core support와 source variant를 구분해 평가한다. A의 최초 baseline 상단은 기존0.5다.
