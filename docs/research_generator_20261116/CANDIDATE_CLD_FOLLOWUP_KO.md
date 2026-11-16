# SPG 후보 CLD 유한 표본 실행 — 2026-11-16

기준 HEAD `674aa04bd3dc90b8d3437ffa8401950ccc43861c`. 판정 **PARTIALLY_COMPLETE**, 독립 검토 **PASS_SCOPED**(아래 유한 표본의 후보 셀에 한함). 사용자 목표인 5 core·13 crust의 독립 재생성은 완료되지 않았다.

## 구현 및 물리 범위

`candidate_gdfm_cld.py`는 저자 GDFM 후보와 기존에 실행된 AME2408 vacuum-droplet 후보 fit을 연결한다. 순수 중성자 기체, 균일 전자, 구형 sharp interface를 가정한다. 원전 식과 같은 형태의 bulk, Coulomb, surface, 전자 전체 에너지와 고정 총밀도에서 네 좌표 `(f,R,n_i,y_p)`의 미분을 사용한다. GDFM nucleon rest energy는 공통 938.9 MeV다. 원전 surface fit과의 동등성, 외피, 전역 최소, 실제 원전 crust 재현은 주장하지 않는다. 저밀도 양성자 기체와 뮤온을 억제한 가정을 각 셀의 삽입 안정성 여유도로 점검했다. 원전은 [Scurto et al., Eqs. 20–25](https://arxiv.org/html/2402.15548v1)이고, 기존 proxy는 변경하지 않았다.

최초 미분 구현에서는 `y_p`를 독립 좌표로 고정하면서 `n_i` 미분에 잘못된 surface 연쇄항을 넣었다. SPG M1/M3/M5의 독립 차분에서 세 건 실패했고, 그 항을 제거한 뒤 검사가 통과했다. 첫 missing-module 실패, 미분 반례와 수정 후 출력은 receipt 디렉터리의 로그에 보존했다. 기체 양성자 삽입 미분은 별도 유한차분으로 검사했다.

## 실제 실행 및 한계

각 SPG 모델에서 원본 phase-7 격자의 약 0.005, 0.015, 0.03, 0.045, 0.055 fm^-3에 가장 가까운 다섯 점을 선택하고 각각 세 초기값을 실행했다. 총 25점 중 24점이 설정된 탐색 범위 내부의 정지 후보이고 M3 최저 에너지 시도 한 점은 optimizer 실패로 **HOLD**다. 모든 초기값의 상태와 원본 비교가 `CANDIDATE_CLD.json`에 있다. 탐색 경계, 혼합 단위 좌표 미분의 수치 기준, input 및 실제 실행 코드 해시는 `RECEIPT.json`에 있다. 이 기준은 새 진단용이며 기존 tolerance를 변경하지 않는다. HOLD에서는 기체 화학퍼텐셜을 이용한 압력을 기재하지 않았다.

24점의 압력 원본 대비 최대 상대 차이 **12.8151%**, 전자 분율 최대 상대 차이 **12.2717%**다. 25점 모두 비교한 npemu 균일상보다 후보 셀 에너지가 낮았다. 이는 해당 점과 finite starts에서의 결과다. 클러스터 최소 양성자 분율은 0.0951로, 원자핵 fit 학습의 최소 0.3077보다 낮다. 저 Yp 외삽 정확도는 미검증이다. 조성은 클러스터와 자유 중성자·전자 분율로 기록되지만, CompOSE 조성 테이블 출력 및 native readback은 아직 없다.

독립 재검토는 25개 셀을 재계산하고 조성 보존 잔차 <=1.11e-16, 채택 24개 셀의 국소 Hessian 최소 고유값 0.0522741 MeV/A, 양성자 삽입 여유도 >=2.54837 MeV, 뮤온 삽입 여유도 >=31.7305 MeV를 기록했다. 세 모델의 추가 밀도 재최적화는 envelope 압력과 최대 8.70e-10 MeV fm^-3 차이를 보였다. 이 국소 표본은 밀도 격자 수렴, 전 구간 인과성, 전역 최소 또는 물리적 crust 완성의 증명이 아니다. 상세 수치는 `INDEPENDENT_REVIEW.json`에 있다.

실제 명령과 SHA는 `data/research_generator_20261116/candidate_cld_followup/RECEIPT.json`을 참조한다. 전체 post-change frontend/pipeline **641 passed, 9 skipped**, generator **142 passed**, 기존 Matplotlib 경고 1건. 실패와 skip은 보존한다.

## 남은 범위

5점/model은 독립 crust 표 전체가 아니다. 밀도 격자 정련, 모든 구간의 안정성·인과성, 독립 조성 테이블 출력·readback, 원본 조성과의 정량 검증, 원전 fit/저 Yp 외삽 검증과 outer crust가 남았다. GPPVA의 meson/Coulomb 공간장과 geometry 최적화는 이 작업에서 구현되지 않았으며 SPG fit에 의존시키지 않는다. APR Gibbs 이후 provider 차이, BBB 저자 압력 미분 처방과 Skyrme 저밀도 차이도 미해결이다. 새 core–crust 결합이나 stellar 계산은 실행하지 않았다.
