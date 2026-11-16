# GPPVA 고정 공간격자 밀도 세분화 독립 검토

`N=192`개의 공간 유한체적 셀을 고정하고, 최종 동일 solver campaign 안에서만 밀도 보간을 비교했다. 모델마다 12개의 coarse 밀도와 실제 재최적화한 중간 밀도 11개를 사용했다. 저밀도 끝점 `n_B=0.002 fm^-3`의 production 선택은 N384이므로, 보존된 동일 밀도의 수렴 N192 시도를 사용했다. 다른 solver 버전이나 공간 격자의 결과를 밀도 수렴으로 섞지 않았다.

측정량은 압력 P, 전체 에너지밀도 ε, 전자분율 Ye다. 선형 보간, PCHIP(n,value), exp[PCHIP(log n, log value)]를 비교했다. 로그 보간은 양의 값에서만 정의한다. 위상 전이가 걸린 구간은 별도로 표시하며 그 구간의 보간 오차를 매끄러운 단일 branch의 오차로 해석하지 않는다.

## 실제 추가 계산

같은 geometry를 유지하는 구간 중 로그 PCHIP 압력 오차가 가장 큰 구간을 모델마다 골라 1/4 및 3/4 밀도에서 새 최적화를 실행했다. 각 점에는 왼쪽·중간·오른쪽의 서로 다른 세 수렴 profile을 초기값으로 사용했다. 모든 시도와 실패·polish 기록을 보존하며 기준을 완화하지 않았다. 16개의 새 점 모두 수렴 후보로 받아들여졌다. 이 새 점들을 공통 정답 표본으로 두고 coarse 12점 보간과 refined 23점 보간을 비교하므로, 다음 표는 동일한 검증점에서 측정한 개선이다.

| 모델 | 추가 밀도 구간 [fm^-3] | coarse 최대 P 차이 | refined 최대 P 차이 |
|---|---|---:|---:|
| GPPVA_DD2 | 0.002–0.005 | 0.686602% | 0.196245% |
| GPPVA_DDME2 | 0.01–0.02 | 0.36979% | 0.0498856% |
| GPPVA_FSU2 | 0.01–0.02 | 3.49725% | 0.394765% |
| GPPVA_FSU2H | 0.002–0.005 | 0.586584% | 0.246429% |
| GPPVA_FSU2R | 0.002–0.005 | 0.823093% | 0.226503% |
| GPPVA_NL3wrL55 | 0.01–0.02 | 0.442499% | 0.0949965% |
| GPPVA_TM1e | 0.002–0.005 | 0.792424% | 0.199229% |
| GPPVA_TW | 0.002–0.005 | 0.749048% | 0.243395% |

새 quarter 표본에서 refined 로그 보간의 최대 상대 ε 차이는 `2.19389e-06`, 최대 절대 Ye 차이는 `0.000126076`이다. 압력 보간 차이는 여덟 모델 모두 감소했다. 이것은 이 유한 표본에서의 개선이며 보편적인 오차 상계나 점근적 수렴 차수의 증명은 아니다.

FSU2의 최초 `0.01–0.02 fm^-3` 중간점 압력 차이는 로그 PCHIP에서도 4.692%였으나, 추가 quarter 표본에서는 coarse 3.497%에서 refined 0.3948%로 감소했다. 세 profile 초기값이 각각 같은 에너지에 약 10^-10 MeV/baryon 수준으로 도달하고 반경이 32.0021/30.1944 fm으로 연속적으로 변했다. 이 검사에서는 다른 더 낮은 droplet branch를 찾지 못했다. 이 근거는 보간 격자 해상도 부족과 일치하며 functional 오류를 입증하지 않는다.

## 위상과 의미론

기존 88개 실제 중간 밀도 검사 중 geometry가 바뀌는 구간은 24개다. 대체 geometry의 에너지 차이가 절대값 10^-4 MeV/baryon보다 작은 경우는 3개이며 이는 설명용 flag이지 새 통과 threshold가 아니다. 원시 에너지 차이와 chemical KKT residual을 AUDIT.json에 보존했다. 이 자료로 보간함수를 위상 경계의 물리적 해로 승인하지 않는다.

CompOSE와의 일치 여부는 이 검토의 조건이 아니다. 표면/저밀도 source 선택, 원전 모델 동일성 및 공간격자 수렴은 별도 receipt의 주장이다. 새 profile은 독립 검토 자료에 저장했으며 production 테이블이나 solver를 수정하지 않았다.

## 재현

```sh
python3 data/research_completion_20261116/tf_density_review/audit_density.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 data/research_completion_20261116/tf_density_review/run_quarters.py
```

quarter_probes 출력 디렉터리는 새 경로여야 한다. 이미 보존된 데이터를 덮어쓰지 않는다. 재실행 전 기존 `quarter_probes/`를 다른 이름으로 보존하면 같은 명령이 새 디렉터리에 결과를 생성한다. 기존 실행 결과와 SHA manifest는 함께 보관한다. `AUDIT.json`은 입력 SHA, `quarter_probes/SUMMARY.json`은 불변 코드 SHA 및 실행 중 불변 여부, 개별 gzip은 모든 profile/시도, `SHA256SUMS.txt`는 보존물 identity를 기록한다.
