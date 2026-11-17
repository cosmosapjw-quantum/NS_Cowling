# C 독립 검토: 범위 내 수용

검토 시점은 2026-11-17T10:16:03.700595+00:00, 기준 commit은 `894cbf1c8315f07ee926c29b83766ffdab2ac7b7`다. 선택된 C의 새 runner, native continuation adapter, 실제 profile, source seals와 390-cell 합집합을 한 번의 fresh-context pass로 검토했다. **현재 범위의 blocking defect는 발견하지 않았다.** 결과는 branch-specific finite-mesh TF와 source-bound matching 근거로 수용한다. 전역 geometry 최소나 연속 TF EOS, fully frozen composition closure의 승인은 아니다.

## 직접 확인한 사항

- Native TF 핵심4개와 legacy matching3개 파일은 기준 commit과 exact-byte 동일하다. TF 최종 identity9개와 retained 초기 state6개의 seal도 실제 파일과 일치한다. 새 driver는 보존 nn/np를 기존 native API로 전달하며 density profile과 log-radius를 재최적화한다. 공개 Ye는 초기 평균에만 적용되고 endpoint N512의 최종 Ye=0.07005257970은 공개 초기0.06988과 다르다.
- 64개 실제 profile에서 직접 계산한 가중 baryon residual 최대는2.1654e-15 fm^-3, 전하 residual은5.4210e-20 fm^-3, 종 밀도의 최솟값은0이다. 45 stationary와19 mode8 원본이 구분되어 남아 있다. Retained Python functional의64개 잔차 근거를 읽었으며 epsilon 최대차7.33e-15 MeV fm^-3, chemical potential 최대차1.69e-8 MeV다. 이는 별도 구현의 동일 방정식 검증이다.
- 실제 profile의 proton-excess radial centroid와 경계 형태는 보고된 droplet/rod label에 부합한다. 같은 N256에서 rod−droplet=0.1414414910 MeV/baryon이고 관측 mesh 민감도 합(rod128→256, droplet256→512)의73.1배다. 이 두 branch의 에너지 순서는 관측 오차척도보다 크지만, 미확보 geometry branch에 대한 최소성은 여전히 unresolved다.
- Droplet/rod 각각 반경 시작값±10%의4회 재최적화는 기존 branch로 돌아왔고 최대 에너지 차이는1.18e-13 MeV fm^-3다. Endpoint의 재최적화3점 density derivative를 실제 저장값에서 다시 계산하면 P_h=0.0019557976816415375 MeV fm^-3로 기록과 정확히 일치한다. Euler 압력과 상대차는−1.5607e-7이다.
- 신규 matching80개, 전체390 unique ID 및 source/output seal506건을 직접 확인했다. 합집합은165 kernel-evidence IDs,105 table IDs,25 independent method-admitted IDs다. 신규 smooth70은 construction-only, CSS10은 independent admission이다. 미완료225는 TF non-Maxwell200과 SPG→APR25로 분리된다. Maxwell427 component invocations는65 unique cell의 근거이며 full-source TF certification으로 바뀌지 않는다.
- B domain audit의 exact hash를 재사용했고 새 표의 실제 core 입력 범위0.08–0.5 fm^-3와 B의0.02–1.0 범위를 혼동하지 않는다. BBB에는 자체 sealed single-phase source/native verification을 적용한다. 물리적 epsilon+Cn 이동이나 public field/radius/cluster 정보의 생성은 없다.
- 공개복원행 endpoint 비교는 diagnostic으로만 읽었다. 공개행 Euler residual0.08244845 MeV/baryon이 보존되어 전체 native−public energy 차이를 모델 효과로 해석하지 않는다. 해당 비교의 입력 두 seal을 직접 확인했다.

## 판정과 남은 의무

현재 수치 결과와 TF_RESULTS_KO/MATCHING_RESULTS_KO의 제한된 선언은 일치한다. Mode8는 numerical optimizer failure, mesh-label 충돌은 수정된 implementation identity 오류, WolframContext 실패는 tool/environment 오류로 구분되어 있다. 성공한 warm restart가 원 실패를 덮어쓰지 않는다. 현재 추가 수리를 요구할 blocking finding은0개다.

Smooth C8 최대3.94%는 공개된 finite-difference diagnostic이며 analytic first-law 실패 또는 독립 smooth-method admission으로 오인하지 않는다. 동일 pass의 종료 전에 추가 bridge80/160 최종 receipt와 갱신된 본문을 읽었다. 같은 endpoint/analytic bridge에서 C8은0.03935953→0.01166746→0.00334835로 감소하고 공통220개 source/endpoint 행은 exact equality다. C7의 최대는 unchanged source 첫 행에 남는다. 이는 sampling diagnostic의 해석을 지지하며 추가2회 실행은 coverage0개·독립 admission0개다. Reviewer는 이 추가 native 계산을 중복 실행하지 않았다. 후속 최상위 CURRENT_OUTCOME/보고서 조립은 root가 맡으며 그 이후 파일은 이 시점의 byte snapshot과 구별한다.

실제 검토 대상 SHA-256, 입력·출력 identity, 직접 재계산 수치와 보존한 한계는 `data/tf_pair_expansion_20261117/INDEPENDENT_REVIEW.json`에 있다. 기존2개 TF/15개 matching-focused 회귀 로그를 읽고 재사용했으며 전체 historical suite나 production을 중복 실행하지 않았다. 이 reviewer는 구현/source/result를 수정하지 않았고 두 검토 파일만 작성했다. C의 commit/push/PR은 수행하지 않았다.

`skill_context_used`: PhysMath Research Loop1.3.0 main과 POLICY_INTEGRATION을 실제 읽고 selected-stage-only, one fresh pass, typed identity, 실제 수치·실패 보존, 기존 근거 재사용을 적용했다. 하위 reviewer는 만들지 않았다.

최종 종료 통지 반영: 2026-11-17T10:16:44.061440+00:00에 같은 pass를 닫았다. TF runner의 환경변수 output-root override는 출력 위치만 바꾸며 기본 경로·물리를 보존함을 읽었다. Parent 최종 full suite는 처음 C native EOS executable 미존재로 실패한 환경/build 기록을 보존하고 CI build 뒤 실행 중이었다. 이 검토는 해당 full suite PASS를 주장하지 않으며 최종 결과는 parent validation receipt가 담당한다.
