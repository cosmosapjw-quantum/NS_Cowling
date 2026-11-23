# 재구성 EOS와 32개 코어의 접합 시험

기존 source-gate census의 672개 행과 새 결과의 대응 키는 `census_cell_id`다. 각 행은 코어 *자료원*과 방법의 조합이다. BPS 후보질량 재구성과 두 가지 고정 N,Z 프로파일을 원래 문헌 식별자에 각각 대응시켰다. 미시적 원본에 대한 동일성은 미확인이다.

| 입력 및 추가 조건 | 조합 | 유한 표 통과 | 공개 EOS |
|---|---:|---:|---:|
| BPS, MTW1964 + MS1966 재구성 | 224 | 132 | 127 |
| NV 인쇄 계수, 고정 N,Z | 224 | 0 | 0 |
| NV D03×100, 고정 N,Z | 224 | 0 | 0 |
| 위 두 프로파일의 화학응답만 제거하고 원래 고정 Yq 축 유지 | 448 | 0 | 0 |
| 위 두 프로파일을 Yq=0 단일축 표로 각각 축약 | 448 | 276 | 264 |
| D03×100 유한 기저 궤도 해, 11개 밀도, 계산된 화학응답 | 224 | 0 | 0 |
| 같은 궤도 해의 단일축 축약 | 224 | 138 | 132 |
| 균일물질 β평형 해 | 224 | 138 | 132 |
| BPS + 유한 합성 보간 + NV1973 저자 fit | 224 | 138 | 132 |
| 합계 | 2464 | 822 | 787 |

**시험의 의미.** 각 조합에 대해 동결된 native matcher에 실제로 입력했으며, 허용된 phase-view 쌍을 합하면 6,016회다. 과거 양성·음성 혼합 control 7개와 번들 재현 시험은 본표에서 제외했다. 단일축 축약은 n, P, ε와 상 레이블을 보존하면서 독립적인 화학좌표를 제거한다. 입력은 영온 스칼라 열역학이므로 Q3=(ε+P)/(n mₙ)−1, Q4=Q5=0으로 놓고, 완전한 축약에서는 Yq=0 단일축 관례를 사용했다. 원래 고정 조성의 β평형과 이 조건은 별개다. 두 고정 조성 프로파일의 원래 Yq=32/982 관례는 코어와 달라 native 게이트에서 먼저 반려됐다. 궤도 해는 원래 Q5와 β 메타데이터로 반려됐다. 유한 후보군 기준으로 BPS의 통과 수는 first-law Hermite 17, cubic/quintic F와 P 각각 23, CSS 23, Maxwell 0이다. 전체 12상 BPS가 남은 경우는 30개다.

**통과 판정.** 커널은 사용 가능한 source view에서 처음으로 통과한 후보를 반환한다. 조립된 표의 native 검산, CompOSE 왕복 검산, phase transport 및 방법별 interface 게이트를 거친 후 실제로 남은 지각 행에 별도의 Fortran 유한 구간 필요조건 검사를 적용했다. 음성 결과는 해당 입력과 수치 탐색에 한정된다. 전 별 구조, 연속적인 모든 구간의 충분조건, 모든 조성에 대한 전역 최소는 별도 문제다. 복합 지각의 통과 138개 중 115개는 BPS 623행, 합성 구간 내부 63행, 실제 NV fit 220행을 모두 남기며 23개 CSS 결과는 외부 BPS 구간에서 접합된다.

**파일.**

- `ALL_CELLS.csv`에 2,464개 조합의 방법, 범위 및 결과를, `MODEL_METHOD_COUNTS.csv`에 각 입력과 방법별 수치를 담았다. `CORE_KEYS.csv`는 32개 입력 식별자의 공개 대응표다. 공개 CompOSE 자료원 식별자가 같지 않은 경우 같은 모델명도 별개 소스다.
- `PUBLIC_MATCHED_ROWS.zip`은 787개 공개 후보를 물리단위 CSV로 제공한다. 경로는 `<reconstruction>/<census_cell_id>.csv`이며 열은 n [fm⁻³], ε 및 P [MeV fm⁻³], Iphase, 구간명이다. ε에는 정지질량이 포함된다.
- 각 자료원 하위 폴더의 `PUBLIC_MATCHED_EOS.zip`에는 공개 통과 후보의 전체 CompOSE 파일과 `MATCH.json`이 있다. `PUBLIC_RECEIPTS.jsonl.gz`에는 해당 217개 공개 셀의 성공·실패 native stage, 방법, source view, 수치 및 오류 원인이 있다. `SOURCE_AUDIT.json`에는 재구성 지각 자료의 전체 구간 필요조건 검사가 있다. 비공개 코어에는 현황과 수치 없는 범위 표시만 공개했다.
- `BPS_MTW1964_MS1966_NV1973_Hermite/COMPOSITE_RETENTION.csv`는 합성지각 후보가 실제 NV 원자료를 갖는지 구분한다. `LEGACY_CONTROL.json`은 이전 NV1973 저자 fit × D1MS 코어 × 일곱 방법의 native 보고서 일치 결과다.
- `SUMMARY.json`의 각 입력 항목에는 밀도 범위, 행 수, 원자료 해시, matcher·auditor 해시 및 출력 해시가 있다. 부가 실험의 결과를 기존 672개 독립 조합에 추가되는 새로운 코어·방법 발견 수로 해석해서는 안 된다.

## 재실행

저장소 루트에서 Python, NumPy, GNU Fortran을 사용한다. 원 세션은 동결된 `compose_eos_fortran` SHA-256 `d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588`로 실행했다. `MATCHER`와 `RUNTIME`은 각각 이 실행 파일과 저장소 밖 작업 디렉터리다. 공식 32개 소스 재현에는 기존 OMEG1의 선택행과 메타데이터가 실제로 있는 core catalog를 `CORES`로 지정한다. 공개 31개 자료원만 재실행할 때는 저장소에 포함된 `data/eos_matching_reaudit_20261123/SOURCE_CATALOG_EXPANDED.json`과 `--public-only`를 사용한다.

```sh
python scripts/run_remaining_reconstruction_matching.py \
  --model BPS_MTW1964_MS1966 \
  --compose data/eos_reconstructed_matching_20261124/generated_sources/BPS_MTW1964_MS1966/generated \
  --cores "$CORES" --matcher "$MATCHER" --work "$RUNTIME"
gfortran -O2 -std=f2008 -ffree-line-length-none \
  fortran_eos/tests/source_prefix_audit.f90 -o "$RUNTIME/source_prefix_audit"
python scripts/audit_remaining_reconstruction_matching.py \
  --work "$RUNTIME" --model BPS_MTW1964_MS1966 --auditor "$RUNTIME/source_prefix_audit" \
  --out "$RUNTIME/public-projection"
```

고정 N,Z 프로파일의 엄격한 시험은 원래 CompOSE n, P, ε, 화학열 및 Yq 축을 사용한다. 추가 시험의 옵션은 화학열 축약에 `--scalar-projection`, Yq까지 단일좌표 관례로 축약할 때 `--scalar-projection --reduce-charge-axis`다. 궤도 해의 원래 축 좌표값은 이미 0이다. 공개 CSV를 다시 생성하는 명령은 다음과 같다.

```sh
python scripts/export_remaining_matching_rows.py \
  --release data/eos_reconstructed_matching_20261124/remaining_matches \
  --out "$RUNTIME/PUBLIC_MATCHED_ROWS.zip"
```
