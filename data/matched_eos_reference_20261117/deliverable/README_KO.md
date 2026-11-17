# Core–crust matched EOS 20개

2026-11-17 campaign에서 확보한 SPG M1–M5 × GM1·SLy4의 10개 모델 쌍이다.
각 쌍에 quintic_P와 CSS 출력이 하나씩 있어 총 20개다. 19개는 이번 campaign의 새 matching 결과이며,
SPG_M1–GM_GM1 CSS 한 개는 기존 산출물의 동일성을 확인하고 현재 native readback·consistency·독립 admission으로 재검증한 결과다.

## 파일

- `eos/<case>/compose/`에 각 EOS의 CompOSE 원본 5개 파일이 있다.
- `eos/<case>/PAIR.json`과 `evidence/`에 해당 실행·검증 기록을 보존했다.
- `CATALOG.json`과 `CATALOG.csv`에는 20개 EOS와 10개 Maxwell 종료 결과를 모두 기록했다.
- `rejected/<case>/PAIR.json`은 EOS가 생성되지 않은 Maxwell 실행의 원 기록이다.
- `SHA256SUMS.txt`에는 ZIP 안의 각 파일 SHA-256이 있다.

`eos.nb`는 baryon density n_B [fm^-3]다. 이 출력의 CompOSE 표현에서는 P=n_B Q1,
epsilon=n_B m_n(1+Q7), mu_B=m_n(1+Q3)이며 m_n은 각 eos.thermo header의 MeV 단위 reference energy다.
epsilon은 rest energy를 포함한다. 음속 항 chi는 c_s^2/c^2=dP/d epsilon이다.

## 범위와 검증

모든 표는 생성된 inner crust부터 core의 n_B=0.5 fm^-3까지를 포함한다.
저밀도 시작점은 모델별로 다르므로 CATALOG의 n_min_fm3를 확인해야 한다.
표면 P=0까지의 EOS와 outer crust 연결은 이 묶음에 포함되어 있지 않다.
matching은 thermo-only이며 eos.compo는 phase 정보를 담는다. 누락된 abundance를 물리적 0으로 해석하면 안 된다.

quintic_P 10개는 native smooth construction, 전체 수치·메타데이터 readback과 serialized consistency를 통과했다.
CSS 10개에는 이 검증에 더해 serialized rows에서 계면 조건을 다시 계산한 독립 CSS admission 기록이 있다.
두 검증 범위를 CATALOG의 verification_type에 구분했다. Smooth의 diagnostic_only는 10개 모두 false다.

CSS 중 SPG_M1__GM_GM1__css, SPG_M3__GM_GM1__css, SPG_M5__GM_GM1__css의 Seidov small-core stability flag는 false다.
이 결과는 삭제하지 않았다. 열역학적 matching admission과 작은 core 출현 시의 안정성 조건은 별도 판단이다.
다른 CSS의 local flag가 true여도 별 전체의 안정성을 증명하지 않는다. 이 묶음으로 TOV 또는 Cowling 관측량을 계산하거나 검증하지 않았다.

## Maxwell 결과

30개 요청 중 Maxwell 10개에는 출력 EOS가 없다. 8개는 현재 공통 mu_B 지원구간에서 압력 교차를 찾지 못했고,
2개(SPG_M2–SLy4, SPG_M5–SLy4)는 두 입력의 mu_B 지원구간이 겹치지 않았다.
이는 이번에 제공된 유한 입력구간의 결과이며, 모든 밀도에서 해당 모델 조합의 공존이 불가능하다는 주장이 아니다.

## 물리 모델의 provenance

SPG crust는 `public_composition_inferred_20261116` candidate surface calibration을 사용한 생성물이다.
원저자의 surface calibration 또는 CompOSE 공개표와의 동일성을 주장하지 않는다.
각 source request·receipt·artifact hash는 CATALOG에 보존했다. 공개표와의 정량 차이 설명은 같은 저장소의
`data/matched_eos_reference_20261117/comparison/` 및 `docs/matched_eos_reference_20261117/` 기록과 함께 읽어야 한다.

포장 과정에서는 모든 EOS 파일을 campaign output inventory와 SHA-256으로 대조한 뒤 ZIP에서 다시 읽어 검증한다.
새로운 생성·matching·별 구조 계산은 수행하지 않는다.
