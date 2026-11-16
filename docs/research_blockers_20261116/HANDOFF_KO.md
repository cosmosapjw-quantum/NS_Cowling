# Local 인계

다음 프롬프트를 사용한다. 실행 시 최종 답변에 제공된 게시 SHA를 PUBLISHED_HEAD에 넣는다. 이 문서 자체의 git commit을 현재 origin 상태와 혼동하지 않는다.

~~~text
PUBLISHED_HEAD=<최종 게시 SHA>
repo=cosmosapjw-quantum/NS_Cowling
branch=remediation/2026-11-15-m5a-maxwell-eos-contracts
PR=15 (draft/open; merge 승인 없음)

원래 목표인 core5 + crust13의 재생성과 독립 물리 검증을 이어라.
먼저 git status, worktree list, origin branch SHA와 PR/CI를 확인하고
PUBLISHED_HEAD 대비 새 변경만 읽어라. dirty checkout을 reset하거나 덮어쓰지 말라.
중앙 인계는 docs/research_blockers_20261116/{REPORT_KO.md,REPRODUCE.md,INTEGRATION_REVIEW.md}
및 data/research_blockers_20261116/{BLOCKER_MATRIX.json,GOAL_DAG.json,MANIFEST.json}이다.

Fortran은 production physics/최적화/phase선택/생성/readback을 맡고
Python은 frontend/transport/provenance를 맡는다.
원전 충실도, 독립 물리, 수치 지원구간, CompOSE 비교를 분리한다.
CompOSE 정확 일치나 full BHF/ETFSI 재계산을 모든 생성의 필수조건으로 추가하지 말라.
반대로 후보 구현 성공을 원본 calibration 동일성이나 전구간 물리 지원으로 승격하지 말라.

재사용:
- 18/18 native 모델 실행 및 이전807core/985SPG/590TFgrid 증거
- 새 TF50구간190개 N32/N64 +4reuse, 97끝점485continuation, N1289점
- core 해석/차분 인과성 경계와 새API지원구간 검사
- BBB 원전Table1/3 대조와 source-only9행 경계 재구성
- BSk 저자 권고 열역학적 적분 및 code–paper C18 불일치
- SPG760행 조건부계수,15새forward평형,fullsigma′ 대 (1−f)sigma′ 잔차 패턴
실패·provenance race·미확인 결과도 원형 보존한다.

다음 실제 물리 작업:
1. FSU2R/TM1e slab–tube 공통접선 후보를 상세TF보고서에서 읽고
   동일 local-neutrality 가정 아래 mesh·branch 정련을 수행하라.
   common tangent는 global two-charge Gibbs/전역최소 증명이 아니다.
2. SPG 원본 조성의 표면미분 잔차 가설을 판별하라.
   추정sigma0,b는 조건부 연구입력이며 production 기본값으로 자동교체하지 말라.
   공개 원전/코드/데이터 재사용 후 정말 없는 입력만 특정하라.
3. 미표본 phase/segment 구간을 자동 보간하지 말고
   명시한 representation과 support에 맞춰 core-crust 조합으로 연결하라.
4. 실제 변경 범위에 필요한 검사 후 commit/push하고 원격CI를 확인하라.
   소스변경 없는 대규모 캠페인을 다시 실행하지 말라.
~~~

새 공개 출력은 sound_speed.csv와 attachment_diagnostics.json이다. 등록된 진단 파일이 사라지는 경우 checksum 검사는 거부한다. source_prescription_valid는 구간 평형·보간 승인과 다르다.

자료 확보 상태: 현재 새 원문 PDF 요청은 필요 없다. 원저자 SPG 계수/export 규약 및 BBB smooth-fit 계수는 정확 동등성을 더 주장할 때 필요한 별도 입력이다. 해당 입력 부재가 선언된 독립 EOS 생성 자체를 막지는 않는다.
