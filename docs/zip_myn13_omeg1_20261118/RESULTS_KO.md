# ZIP-only MYN13–OMEG1 Maxwell 재평가 · 2026-11-18

**Source scope:** 사용자가 별도 제공한 원본 ZIP의 `EoS.d`, `RHF-TF.data`, `Unified-EoS.d`, `Crust-core-trans.d`만 물리 데이터로 사용했다. 파일의 원본 bytes는 변경하지 않았다. **이 브랜치는 공개 가능한 계산 요약과 SHA-256만 보존한다. 원본 사용자 ZIP·첨부 PDF와 수치 표 자체는 GitHub에 게시하지 않았다.** 전체 재현 코드, 4쪽 PDF 감사보고서, M–R/Cowling 그래프와 raw 계산은 대화의 다운로드 묶음에서 제공한다.

## 확인된 물리적 결론

- `Crust-core-trans.d`의 구분 밀도는 n_B=0.079 fm^-3. `Unified-EoS.d`는 마지막 crust n_B=0.075 fm^-3 (P=0.398130 MeV fm^-3, mu_B=951.215067 MeV) 다음에 첫 core n_B=0.079 fm^-3 (P=0.6781606351 MeV fm^-3, mu_B=955.886470 MeV)를 붙인 **hard splice**다. 채택된 두 support는 P와 mu_B에서 모두 겹치지 않는다.
- 첨부 Maxwell 검산 PDF는 MYN13 **비균질 crust의 더 좁은 경계** n_B=0.06757 fm^-3를 사용했으며, 해당 행은 P=0.305670 MeV fm^-3, mu_B=950.668492 MeV다. 원래 unified splice의 마지막 crust 행 n_B=0.075와 혼동하면 안 된다.
- 원본 `EoS.d`에는 n_B<0.079인 **실제로 계산된 균일물질 가지**가 있다. 이를 메타안정성 검증이 필요한 uniform branch로 취급하면 공통 chemical potential 구간은 약 939.917389–951.215067 MeV가 된다. 이 방법은 미지의 외삽을 도입하지 않지만 채택된 core support를 물리적으로 확장한 것은 아니다.

## 첨부 PDF 결과의 독립 재계산

첨부 PDF의 계약을 그대로 구현: log epsilon(log P)를 로그선형 / cubic Hermite / quintic Hermite로 보간한 뒤 d ln(mu_B)=dP/(epsilon+P)를 각 source n_B=0.06757 (MYN13), n_B=0.079 (OMEG1)에서 정규화하고 n=(epsilon+P)/mu로 재구성. **채택 범위 0교점 / MYN13 hetero crust + OMEG1 전체 1교점 / 전체 원자료 4교점**이 세 방법 모두 재현됐다.

| 전체 원자료 교점 | 로그선형 P_t | cubic P_t | quintic P_t (독립) | 물리적 거절 이유 |
|---|---:|---:|---:|---|
| #1 | 0.001104263 | 0.001094616 | 0.001094646 | uniform core 후보 n≈0.001825, adopted core 시작 n=0.079 아래 |
| #2 | 0.394911474 | 0.394509504 | 0.394492818 | 음의 density jump, hetero crust 범위 밖 |
| #3 | 10.1941057 | 10.1980094 | 10.1981626 | MYN13 측이 비균질 crust의 범위 밖 |
| #4 | 69.407893 | 69.4055091 | 69.4051829 | 음의 density jump와 crust 범위 밖 |

압력 단위 MeV fm^-3. 로그선형과 cubic 수치는 첨부 PDF 표시 정밀도와 일치. **Quintic의 작은 차이는 PDF에서 2차 도함수 추정법을 완전히 지정하지 않은 method convention에 의존**하므로 모든 유효숫자까지 검증했다고 주장하지 않는다.

각 방법의 첫 번째 교점에서 c_s^2/c^2는 양측 모두 양수이며 1 미만이다. 그러나 국소 안정성은 클러스터 불안정성/전역 에너지 최소를 증명하지 않는다. **물리적으로 채택할 strict Maxwell crust–core 접합은 0개**다. 첫-node 정규화로 바꾸어도 교점 수 4는 유지되나 첫 교점 P_t는 변한다.

별도 진단으로 원자료 각 knot의 mu=(epsilon+P)/n을 전부 고정해 P(mu)를 직접 Hermite 보간하면 (PDF의 적분 재구성과 다른 thermodynamic closure) root 수가 달라진다. 이 현상은 PDF 자체 오류 증명이 아니라 source precision/closure/normalization 문제다.

## ZIP-only M–R 및 Cowling 산출

물리적으로 승인된 Maxwell-star 곡선이 없으므로 **Maxwell-star M–R 및 f/p1 그래프는 결과로 승인하지 않았다**. 기존 ZIP-only 계산은 1.4 M_sun에서 다음과 같다. CSS는 현상론적 sharp junction, smooth Hermite는 no-jump 연결, unified는 원본 hard-splice 대조군이라는 명시적 차이가 있다.

| EOS | R (km) | f (Hz) | p1 (Hz) |
|---|---:|---:|---:|
| CSS+cubic, slow | 12.790524 | 2156.020 | 5798.866 |
| CSS+quintic, slow | 12.790525 | 2156.024 | 5798.866 |
| smooth+cubic | 12.790523 | 2156.022 | 5797.822 |
| smooth+quintic | 12.790546 | 2155.978 | 5796.333 |
| unified hard splice + cubic | 12.771275 | 2157.295 | 5828.750 |
| unified hard splice + quintic | 12.771215 | 2157.305 | 5828.795 |

기존 수치 근은 l=2, perfect-fluid equilibrium Cowling으로 계산했으며 정확한 microscopic phase/elasticity/frozen composition을 해결하지 않는다. 표면 이하 외삽도 조건부 closure이다.

## 검증 상태 및 남은 일

대화 실행에서 원본 ZIP byte-identity, unified splice, PDF root counts/물리적 거절 등을 **58개 검사 모두 PASS**했다. 기존 TOV/Cowling raw evidence는 별도 240개 점검을 통과했으며 두 점검의 범위는 다르다. **이 GitHub summary 자체에는 원본 소스 표를 포함하지 않으므로 repo 단독 재실행은 불가능**하며, 사용자가 제공한 동일 SHA-256 원본 ZIP이 필요하다.

원 functional로 n~0.0018 fm^-3 OMEG1 homogeneous branch의 cluster/spinodal/pasta 전역 경쟁을 확인하는 것이 다음 최소 물리작업이다. 원자료를 넘어선 새 모델을 보태지 않았다.
