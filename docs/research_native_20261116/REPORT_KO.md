# Fortran native EOS 생성 연구 루프

기준 HEAD는 `4de49430560930be84d9a4147d71c3c9841da2db`이다. 이번 작업은 기존 계산 결과를 Fortran 실행으로 재표기하지 않고, core 5종·SPG 5종·spatial TF 8종의 실제 생성 경로를 새로 구현하고 실행했다. Python은 입력·실행·파일 검증·표시를 담당하며, functional·평형·field·gradient·최소화·후보 선택·EOS 생성은 Fortran에서 수행한다.

**18개 모델 이름에 대한 native 생성 경로와 원저자의 완전한 EOS 재현은 서로 다른 완료 조건이다.** 공개 CLI의 유한 검증에서는 18/18개가 성공했지만, 원 SPG 표면 보정의 동등성, BSk ETFSI 자체의 재계산, BBB microscopic BHF kernel, 연속적인 mesh 수렴·전역 최소를 완료했다고 선언하지 않는다.

## 이론 연구 결과

원전 충실도, 독립 물리 검증, 수치 지원구간, CompOSE 비교를 별도 상태로 유지했다. 13개 crust 공급자 PDF는 모두 명목 밀도 0.002 fm^-3 아래에 BSk22를 사용한다고 명시한다. 따라서 새 common-tangent 해가 존재해야만 그 처방을 재생성할 수 있다는 기준은 원전 요구가 아니다. 반면 원 자료에는 공통 outer 최종점 0.000269021 fm^-3와 각 inner 시작점 사이의 저장 행이 없고, 그 양끝 압력도 같지 않다. 이를 Maxwell 전이로 해석하거나 미해결 구간을 연속 물리 지원으로 표시하지 않았다.

SPG에서는 원전의 같은 밀도 cluster/uniform 에너지 crossing을 native로 구현했다. 별도의 common-tangent 처방과 혼동하지 않는다. 표면 계수는 명시적으로 선택한 후보이며, 2,408개 핵종 목록의 동일성은 fit 가중치·곡률 처방·원저자 계수의 동일성을 입증하지 않는다. Scurto 2025 학위논문 확인 후에도 원저자 sigma0,b 동등성은 미해결이다.

Native outer drip와 BSk source-fit 시작점은 약 2.90e-8 fm^-3 겹친다. outer→lower source-fit 차이는 ΔP=-3.34991421748e-6 MeV fm^-3, ΔμB=-0.0508156514111 MeV이다. Source-fit 조성 μ와 열역학적 μ 차이의 이번 표본 최대값은 0.0869130 MeV이다. 이것은 이전 outer-anchored 처방의 약 0.0255 MeV와 다른 계산이다. 각 처방과 부호를 보존했고, 상대 에너지 이동을 넣어 접합을 강제하지 않았다.

## 구현과 직접 계산

| 범위 | 실행·검증 | 의미와 한계 |
|---|---|---|
| Core 5종 | 807행 native 생성 및 독립 구현 대조 | APR two-charge Gibbs; BBB published-energy reconstruction |
| SPG 5×5 | native 985행; 기본 후보는 모델별 129점, 대안은 17점 | 원저자 surface fit 동등성 미해결 |
| TF 8종 | N16, 3밀도×6 family=144개 계열 정지 후보 | 유한 후보; 전역 최소·연속 mesh 한계 미증명 |
| Native lower/outer | 8,481개 핵종 목록; outer 99행·17 전이; BSk analytic | ETFSI 재계산 아님; 접합 잔차 보존 |
| 공개 CLI 전체 | 18/18개, native 테이블 44개, 3,048행 readback | main 916행+각 crust의 별도 lower/outer |

N16 검증은 저해상도 실행 가능성 검증이다. 별도로 보존된 FSU2R N384 profile을 native field/functional로 재평가했고, 독립 유한차분 gradient, charge/baryon closure, field 잔차와 비교했다. 이전 156개 고해상도 캠페인의 producer와 이번 native 계산을 혼합하지 않았다.

Native SLSQP는 SciPy v1.15.3의 허가된 Fortran 원본을 exact commit과 SHA로 고정했다. Callback 실패·비유한값·반복 한도·infeasible equality를 구분하고, optimizer 성공 자체를 물리 정지 조건으로 사용하지 않는다. TF family 실패와 모든 재시도 profile을 보존하며, 미해결 family 또는 더 낮은 미수렴 후보가 있으면 HOLD한다.

통합 독립 검토에서는 TF의 모형별 rest-mass header, uniform의 빈 cluster carrier, 실제 실패 원인 보존을 수리했다. 성공 receipt 문자열만으로 잘못된 파일이 승인되던 결함도 재현 후 수리했다. 현재는 별도의 실제 Fortran readback, 요청 밀도·mesh·APR 추가점·composition mask·출력 identity를 검사한다. C1–C11은 Fortran 진단으로 기록하며 일괄 PASS 조건으로 쓰지 않는다.

최종 전체 회귀는 frontend/pipeline **1,074 passed, 9 skipped**, generator **142 passed**, native optimizer **10개 PASS**다. 기존 325개 코드·빌드·CI 파일은 실행 전후 동일했다. 독립 cohort runner가 추가되어 최초 source-set 검사가 실패한 기록과 그 정확한 차이는 REGRESSION_RAW.json 및 검토한 REGRESSION.json에 보존했다.

추가로 TF 8종의 원 배포 inner 밀도 590점 전체를 N16에서 계산했다. 3,641개 시도 중 87개 내부 field 평가 실패를 보존했으며, 최종 590점은 모두 후보 선택·native readback이 가능했다(HOLD 0). 일부 geometry는 N192/N384 결과와 다르므로 mesh 수렴 주장으로 사용하지 않는다. 별도 lower/outer를 포함한 이 캠페인은 공개 readback 테이블 24개, 1,902행이다.

## 판단

구현 검증은 원저자 입력 회복을 대신하지 않으며, 같은 표면·functional 구현의 재생성 성공도 모든 보간 구간의 물리적 유효성을 뜻하지 않는다. 공개 생성기는 모델별 선언된 representation을 재생성할 수 있다. 원전과 동등한 완전 EOS라는 더 강한 최종 목표는 아래 독립 항목이 남아 있어 부분 완료다.

1. SPG 원저자 surface fit 계수·가중치·curvature 설정의 확인 또는 동등성 검증.
2. BSk analytic lower와 원 ETFSI 계산의 구별 및 source attachment의 미확인 상세 처방.
3. 필요할 경우 BBB microscopic BHF 재계산을 별도 목표로 정의; 현재 published-energy reconstruction과 구별.
4. 실제 사용하는 TF 밀도 구간에서 mesh·seed·radius 정련 및 phase 경계별 지원구간 확정.
5. 접합 부근의 물리적 연결 처방 선택. 문헌 prescribed boundary를 그대로 생성하는 것과 thermodynamic seamless EOS는 별도 결과다.

재현 명령은 REPRODUCE.md, 모델별 상태와 정확한 실행 결과는 data/research_native_20261116의 receipt 및 MODEL_STATUS.json을 기준으로 한다. 최종 회귀와 게시 SHA는 EXECUTION.json 및 Git commit에 연결된다.
