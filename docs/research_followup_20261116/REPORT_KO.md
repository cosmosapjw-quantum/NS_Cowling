# 독립 EOS 후속: SPG surface와 DD2 geometry 표본

기준 HEAD `ad1ff42fa596a2db4c78b90d764d876de722ee68`. 전체 목표는 **PARTIALLY_COMPLETE**다. 기존 계산을 재사용했으며 새 core–crust 결합이나 stellar/mode 계산은 없다.

## SPG 표면: 추가 원전 구분과 고정 셀 감도

[Scurto et al.](https://arxiv.org/html/2402.15548v1)의 M1–M5별 원저자 surface fit 계수 및 실행 처방은 확인되지 않았다. [Zenodo 11084025](https://zenodo.org/records/11084025)의 공개 목록에는 EOS, 모델 매개변수, NMP 및 README가 있으나 별도 surface 계수나 fit 코드가 나열되어 있지 않다. 원저자 fit cohort와 2408개 score cohort의 관계, 질량 변환·유한핵 밀도조건·항·가중치·최적화 및 각 모델 계수는 여전히 `SOURCE_CALIBRATION_INPUT_MISSING`이다. 이것은 기존 후보 모델의 수치 실패 판정이 아니다.

Scurto가 관련 문헌으로 인용한 [Dinh Thi et al.의 유한핵 질량식](https://arxiv.org/html/2109.13638v2)과 [후속 논문](https://arxiv.org/html/2111.04374)은 surface 외 curvature 항도 포함한다. 그 논문은 비대칭 핵물질의 평형밀도를 사용한다고 설명하면서 에너지밀도에 대한 미분이 0인 조건을 인쇄한다. 이는 에너지/입자의 정상점 또는 압력 0 조건과 구분해야 한다. 표기와 실제 구현 중 무엇이 원저자의 의도인지는 미확정이다. 다른 모형을 위한 표의 계수는 SPG M1–M5의 원저자 계수로 대체하지 않았다.

기존 SPG 후보 셀 2,565개에서 **dense-phase** 양성자 분율을 추출하여, 기존에 실행한 네 가지 9/8핵종·면적 처방의 surface 장력과 미분을 2408핵종 후보 fit과 비교했다. 실제 셀의 y 범위는 모델별 대략 0.0428–0.3257이다. 최대 절대 상대 차이는 sigma **32.9746%**, sigma' **32.8736%**이며 둘 다 M4의 geometric-area nine 대안에서 발생했다. 비교는 기존 후보의 동일한 셀 y에 함수를 평가한 고정 셀 민감도이고 EOS를 재최적화하지 않았다. 모두 p=3이므로 이 수치는 저-y 전체 불확실도 대역이나 원저자 fit 동등성이 아니다. 입력 및 실행 코드 SHA와 각 대안의 최대점은 `SPG_FIXED_CELL_SURFACE.json`에 있다.

## DD2 geometry 직접 계산

기존 검토의 N192 midpoint/quarter 자료를 기반으로, 기존 0.060 droplet, 0.065 rod, 0.070 slab의 세 profile을 초기값으로 0.0625와 0.0675 fm^-3를 직접 재최적화했다. 기존 및 신규 계산의 공유 solver SHA가 일치한다. N192를 고정했으므로 이 작업은 density 표본 정련이지 N192/N384 공간격자 수렴 검사가 아니다. 기존의 공간격자 결과는 이전 receipt를 참조한다.

첫 실행에서는 0.0625의 rod 초기값이 screened-electron Newton line search에서 예외를 냈음에도 이전 `select`가 예외에 에너지가 없다는 이유로 droplet을 선택해 native 출력까지 만들었다. 이것은 **불완전 geometry 비교의 잘못된 선택**으로 분류하며, 첫 출력은 과학적 결과로 사용하지 않는다. 당시 스크립트 원본은 receipt의 SHA와 바이트 단위로 일치하도록 복원해 함께 보존했다. 반례 테스트를 추가하고, 각 시드에서 적어도 한 수렴 결과가 없으면 `HOLD_UNRESOLVED_SEED`를 반환하도록 새 runner를 고쳤다. 실패 warm seed 후 cold seed로 한 차례 재시도했지만 같은 오류가 반복되어 0.0625는 HOLD다. 이 점을 보간하거나 0으로 채우지 않았다.

0.0675에서는 세 시도가 수렴했고 rod가 가장 낮은 관측 후보였다. slab와 droplet은 각각 2.55750e-5, 3.12100e-4 MeV/baryon 높았다. 이는 finite seed 비교다. 기존 0.065 rod 및 0.070 slab와 함께, rod–slab 표본 bracket은 **[0.0675,0.070] fm^-3**로 좁아졌다. 작은 에너지 차이, 유한 시드 및 미정련 공간격자를 고려해 정확한 전이 root나 전역 최솟값을 주장하지 않는다. 0.0625 HOLD 때문에 droplet–rod bracket은 추가로 좁히지 않는다.

0.065, 0.0675, 0.070 세 행만 별도로 native write/readback했다. 0.0625는 이 산출물 밖에 명시적으로 HOLD로 남는다. 0.0675–0.070에는 geometry 변화가 있으므로 세 readback 행은 위상 경계를 가로지르는 연속 보간의 승인이 아니다. 원자료, 원래 failed selection, 재시도 및 세 행의 조성 mask/profile을 모두 보존한다.

## 검증과 후속

재현 명령, 코드/seed/출력 SHA와 상태는 `data/research_followup_20261116/RECEIPT.json`과 각 작업 디렉터리 receipt에 있다. 테스트에서 예외 가지가 성공 가지에 숨겨지는 반례는 구현 전 실패했고, 수정 후 unresolved seed와 낮은 미수렴 에너지 가지 모두 HOLD임을 검사했다. 새 geometry 샘플은 원본 EOS 수치 일치를 승인 조건으로 사용하지 않았다.

미완료: 0.0625 rod의 screened-electron solver 원인 분류와 안전한 수리, DD2 이외 geometry와 밀도 정련, 밀도 보간과 공간격자 오차의 독립 평가, SPG 원저자 fit과 저-y 모델 불확실성, 공통 lower crust와 각 inner functional의 질량·에너지·조성 기준 정합. BSk22 원전 analytic fit, 경계 고정 first-law 표현 및 full ETFSI는 별개의 모델이다. 실패점의 보간이나 미정의 조성의 0 대체는 허용하지 않는다.

최종 local 회귀는 frontend/pipeline **775 passed, 9 skipped**, generator **142 passed**이며 generator에 기존 Matplotlib 경고가 있다. 독립 검토는 세 native 행의 phase `[22,22,23]`, 미정의 자유 입자·A/Z mask, profile 전하 잔차 <=1.25e-17 및 에너지 재계산 차이 <=1.43e-14를 확인했다. SPG 고정 셀 평가의 40개 극값도 별도 p=3 식과 최대 6.7e-16 차이였다. 이는 **PASS_SCOPED**이며 density 전 구간의 수렴·원저자 surface calibration을 승인하지 않는다.
