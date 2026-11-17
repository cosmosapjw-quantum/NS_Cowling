# EOS generator blocker 연구 — 2026-11-16

기준은 원격 HEAD a49aaf9f82ddfda484e0f965aefdfe73aa7918f2다. 여섯 작업을 병렬로 수행해 원전 확보, 변분식 유도, native Fortran 계산과 공개 생성 경로 통합을 마쳤다. 이전 18/18 모델 생성 실행 증거는 재사용했다. 이 수치는 원저자 입력의 완전 회수나 모든 밀도에서 물리적 지원을 뜻하지 않는다.

물리 계산은 Fortran, Python은 입력·실행·readback·provenance를 담당한다. 외부 CAS와 Python 수치 적분은 독립 연구 검산에만 사용했다. 원전 충실도, 독립 물리 검증, 수치 지원구간, CompOSE 비교를 따로 기록한다. CompOSE 수치 동일성은 통과 조건이 아니다.

## 결과와 원인 분류

| 대상 | 실제 새 결과 | 원인과 남은 범위 |
|---|---|---|
| TF 압력 감소 | 나머지 50구간 끝점 N32/N64 190개 신규 + 4개 재사용. 97개 N64 끝점에서 485회 branch continuation | N64 잔여 3구간 중 하나는 N128에서 해소. 두 slab–tube 구간은 finite branch 비볼록성으로 분리하여 공통접선 계산 |
| APR/KDE0v1/SLy4 | 해석 Hessian, 압력 차분, 화학퍼텐셜 차분으로 인과성 경계 확인 | 변경하지 않은 functional의 고밀도 한계 |
| BBB | 원 TeX Table1·3 대조, 35개 full-source 및 21개 source-only native 상태 | 공개 에너지 재구성 가능. 저자 smooth fit 동일성은 별도 미해결. .7<n≤.8 보간이 .9 외삽 노드에 의존 |
| BSk lower | source-energy anchor 압력 적분, 41점 독립 검산, code–paper 비교 | analytic fit의 독립 근사 차이. 열역학 재구성은 저자 권고 경로이나 조성 stationarity까지 복구하지 않음 |
| SPG 직접 입력 | arXiv source 38파일과 공개 선행 코드 조사 | 원저자 계수 payload 미회수. APS 접근 차단은 supplement 부재 증거가 아님 |
| SPG 역복원 | 760개 조성에서 조건부 표면 계수 복원, 15개 새 native 평형 | 표면식은 정밀하게 복원되나 원 조성 평형 잔차가 남음. 변분적으로 일관된 대안 계산 |

### TF: 유한 branch에서 실제 공존 후보 구성

N64에서 남은 세 압력 감소 중 TM1e rod–slab 구간은 N128 세 점 모두 slab로 선택되며 사라졌다. 두 slab–tube 구간에서는 각각의 평형 profile을 native로 다시 최소화해 P와 μ_B가 같은 두 밀도를 찾았다.

| 모델 | n_slab (fm⁻³) | n_tube (fm⁻³) | P_t (MeV fm⁻³) |
|---|---:|---:|---:|
| FSU2R | 0.07275856095 | 0.07431969379 | 0.2342521711 |
| TM1e | 0.08122389021 | 0.08237448523 | 0.3275351091 |

차분 간격을 절반으로 한 최종 ΔP는 각각 5.53e-8, −8.63e-8 MeV fm⁻³로 지정한 2e-7 목표 안에 있다. TM1e의 초기 검산 실패는 보존하고 더 작은 step·엄격한 optimizer 설정으로 정련했다. 내부 세 표본에서는 공통 chord가 두 branch보다 낮았다. 세 점 research CSV를 Fortran에서 생성·재읽기했으며 μ_Q나 조성을 임의로 채워 CompOSE 표로 만들지 않았다.

이는 국소 전하중성·β평형인 두 canonical branch의 Maxwell-type 공존 후보다. 전역 전하를 허용한 두 전하 Gibbs 해, 전역 최소, mesh 극한을 증명하지 않는다. raw production 표와 기존 source 전이 처방도 보존한다.

### Core: 모델 한계와 미분 한계

냉각·중성·β평형 경로에서 c_eq²/c²=dP/dε=n(dμ_B/dn)/μ_B를 쓴다. ε는 정지에너지 포함, n은 fm⁻³, P·ε는 MeV fm⁻³, μ는 MeV다. 조성 최소화의 Hessian Schur complement를 유도했고 APR 혼합상은 두 전하 Gibbs susceptibility로 별도 확인했다.

| 모델 | c_eq=c의 native 수치 부호 bracket (fm⁻³) |
|---|---|
| APR | [0.8689531515166163, 0.8689531516283750] |
| KDE0v1 | [1.4109658931940792, 1.4109658933803437] |
| SLy4 | [1.2650056713819504, 1.2650056715682150] |

여섯 차분 간격에서 bracket이 같고 해석식과 c² 비율 차이는 최대 8.19e-12다. 이는 엄밀한 interval-arithmetic 오차막대가 아니다. 1e-8 fm⁻³ 운영 buffer는 관측된 수치 오차만 겨냥하며 계수·모델 오차를 포함하지 않는다. BBB와 GM1은 조사한 표본에서 인과적이었다. 인과성은 전체 유효성의 충분조건이 아니다.

생성기는 sound_speed.csv를 내보낸다. muon onset·phase 경계·BBB PCHIP knot를 가로지르거나 모델 지원 밖으로 나가는 차분은 상태 메시지와 빈 수치 필드로 남긴다. 실패한 계산의 초기값 0을 물리적 음속으로 내보내지 않는다. 원 EOS 값은 보존한다.

### BBB: 재현 가능한 표현과 원저자 fit

공개 SNM/PNM 에너지와 parabolic isospin interpolation의 재구성은 정의된 EOS다. 이 구현에 모든 microscopic BHF kernel 재계산이 필요하다는 조건은 제거한다. 그러나 반올림 표와 현재 PCHIP로 저자의 미공개 smooth-fit 계수나 Table3 압력이 복구되었다고 주장하지 않는다.

.8 fm⁻³까지의 source-only 구현은 첫 9행으로 endpoint slope를 다시 구성한다. .8에서 압력은 247.887661→250.152005 MeV fm⁻³(+0.91346%), .75에서 −0.23597% 바뀐다. 기존 표현은 보존한다. 반올림 corner와 cubic 보간 비교도 실행했지만 엄밀 오차 상한은 아니다. 큰 Table3 압력 차이를 단순 반올림으로 설명할 근거는 얻지 못했다.

### Lower crust: 독립 fit의 열역학적 재구성

Ioffe 저자 페이지는 독립 fit의 모든 출력을 동시에 쓰기보다 열역학 항등식으로 n 또는 ρ를 계산하도록 권고한다. 같은 source-energy anchor를 고정해 다음 식을 구현했다.

e_I(n)=e_source(n0)+∫[n0,n] P_source(s)/s² ds, ε_I=n e_I, μ_I=e_I+P_source/n.

branch를 맞추려고 에너지 상수를 fit하지 않았다. 41점 독립 적분과 e 차이는 최대 1.51e-11 MeV, 첫 법칙 차분 잔차는 1.21e-8 MeV다. source-fit e 변화 최대 0.06416 MeV와 조성 μ 잔차 0.05729 MeV는 적분 오차와 구분한다. 압력의 pointwise envelope를 가정할 때만 성립하는 조건부 오차 전파식을 제공하며 논문의 경험적 정확도를 그 상한으로 대체하지 않았다.

Pearson 원 계산에도 outer HFB와 inner ETFSI의 근사 차이가 있다. 전체 ETFSI 재계산은 seam 일치를 자동 보장하지 않는다. 인쇄 C18의 n,+259n²와 배포 코드의 n_free,−259n_free²도 다르며 조사 구간 μ 차이는 최대 0.09154 MeV였다. 기존 native는 배포 코드를 따른다.

공개 생성기는 attachment_diagnostics.json에 실제 outer/lower와 lower/inner의 gap·overlap 및 P·μ 차이를 추가한다. source_prescription_valid는 명목 경계와 끝점의 구조적 관계만 뜻한다. 공존 평형이나 구간 보간의 승인이 아니다. source-fit 출력과 연구용 적분 표현을 혼합하지 않았다.

### SPG: 조건부 복원과 독립 평형 검증

구형·양성자 gas 없음·geometric A 가정에서 y=Z/A, f=(Ye/y)n_B/n_I, R=(3A/4πn_I)^(1/3)를 사용한다. virial과 기계적 평형을 결합해 표면 계수 없이 n_I와 σ(y)를 복원한 뒤 두 계수 표면식에 맞췄다. 원 테이블 gas 칸 0은 실제 gas 부재로 해석하지 않았다. gas-excess A, proton-drip 행 제외, 입력 반올림, 합성 forward→inverse 검사를 따로 보존했다.

복원 계수는 원저자 값으로 확정하지 않았다. held-out σ 오차는 작으나 화학·β평형 잔차는 반올림 효과보다 크다. β 잔차가 f·3σ′/(R n_I)를 따르는 패턴은 source table이 표면 미분항에 1−f를 적용한 처방에 가깝다는 가설을 지지한다. 원 논문·thesis의 에너지 변분식과 현재 구현은 full derivative를 요구한다. 저자 코드 없이 특정 코드 버그의 증명으로 부르지 않는다.

복원 계수를 현재 변분식에 넣은 5모델×3밀도 모두 stationary였다. 원 조성과 차이는 Ye 최대2.46%, A1.87%, Z2.68%다. 이는 조건부 입력으로 만든 물리적으로 일관된 대안이며 원본 동일성 선언이 아니다. 생산 기본 표면 처방을 이 계수로 바꾸지 않았다.

## 검증과 게시

통합 Fortran build, 공개 CLI core5와 DD2+lower 생성·native readback을 실행했다. 최종 frontend/pipeline 1078 passed, 9 skipped, generator 142 passed였다. 독립 검토에서 seal된 optional 진단 파일 삭제를 놓치는 문제를 수리했다. 수정 후 receipt·architecture·dispatch 65 passed에 두 삭제 회귀가 포함된다. 직접 호출하는 sound-speed API의 지원범위 검사도 수리했고 전용 20개 경계 검사를 통과했다. 최종 source identity는 실행 전후 같았다. 진단 추가 전후 45개 native 테이블 파일은 byte-identical이며 기존 물리 dependency 24개도 기준 HEAD와 같았다. 원격 CI는 게시 인계에서 따로 확인한다.

native TF 두 실행은 child 성공 뒤 동시 source 변경의 상위 provenance 검사에서 실패했다. 원 실패 receipt를 보존하고 실제 TF dependency·binary identity와 native readback으로 재사용 근거를 별도 작성했다. 이를 물리 solver 실패나 원 실행 PASS로 바꿔 쓰지 않았다.

## 재개 지점

원본 동등성의 남은 입력은 SPG 계수·bulk 연결·export convention, BBB smooth-fit/미분 규약이다. 공개 원문은 확보했으며 지금 사용자에게 새 논문 PDF를 요청할 필요는 없다. 이 입력 없이도 선언된 EOS 표현과 독립 물리 검증은 진행 가능하다.

다음 물리 작업은 TF 공존 후보의 mesh·branch 정련과 미표본 교차 구간 지원 결정, SPG 표면 미분 처방 가설의 추가 판별이다. 동일성 미해결을 모든 core/crust 생성 불가능으로 확대하지 않는다. 유한 실행 성공도 모든 밀도의 물리적 지원으로 확대하지 않는다.

세부 근거는 [core](core/CAUSAL_DOMAIN.md), [BBB](bbb/RESOLUTION.md), [lower](lower/REPORT_KO.md), [SPG 공개 입력](spg_sources/SOURCE_RECOVERY_KO.md), [SPG 역복원](spg_inverse/REPORT.md) 및 [TF](TF_REPORT_KO.md) 보고서와 [통합 검토](INTEGRATION_REVIEW.md)에 있다. 기계 판정과 원자료는 data/research_blockers_20261116에 보존한다.
