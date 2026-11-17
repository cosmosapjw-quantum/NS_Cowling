# 세 계열 sharp Maxwell 산출물의 전달 계약

패키징은 실행이 끝난 `CAMPAIGN.json`과 그 파일이 hash로 지목한 preparation
manifest를 읽는다. 새로운 EOS 계산이나 source 보정은 수행하지 않는다.
`native_generated`, `public_original`, `public_reconstructed` 각각에 독립 ZIP을
만들고, 세 계열의 crust 13 × core 5 조합을 합친 195개 cell을 JSON과 CSV에 기록한다.

## 포함 조건과 정확한 개수

`matched_eos/`에는 candidate의 상태가
`GENERATED_VERIFIED_FINITE_DOMAIN`이고, admission receipt가 `ADMITTED`, 실제
serialized output의 native consistency가 `PASS`인 경우만 포함한다. 모든 EOS
member는 생성 결과의 SHA256과 대조한 후 원본 bytes 그대로 ZIP에 넣는다.
실패 후 디렉터리에 남아 있는 compose 파일은 포함하지 않는다.

같은 crust/core 모델 쌍의 독립 source component들에서 여러 candidate가 성공할
수 있다. 따라서 `pairs_with_output`은 모델 쌍의 개수이고 `generated_tables`는
실제 전달된 candidate 표의 개수다. 이를 유일한 전역 평형 EOS의 개수로 해석하지
않는다. 결과가 없는 계열도 65개 cell의 실패·미지원 상태와 source evidence를
담은 독립 ZIP을 제공한다. source 데이터 파일의 수를 생성된 matched EOS 수에
더하지 않는다.

## 유한 지원구간과 source tail 보존

검증된 중앙 산출물은 해당 source의 실제로 허용된 연결 component를 사용한
유한 구간 EOS다. 공개 crust에 원래 붙어 있던 저밀도 자료, 허용되지 않은 보간
구간, core의 고밀도 원본 자료까지 전부 검증됐다는 뜻은 아니다.

각 ZIP의 `sources/`에는 원본 전체 EOS, 준비된 thermodynamic-only 입력,
source-interval 진단, 재구성 변경 기록, composition 및 auxiliary projection의
receipt를 별도로 보존한다. 공개 원본은 그 EOS bytes와 provider 자료를 유지한다.
reconstructed 계열의 공개 원본은 그 재구성의 부모 자료이며, 서로 다른 계열의
EOS를 접합한 것으로 분류하지 않는다.

각 성공 candidate에는 추가로 `diagnostics/source_retaining_thermo.csv`를 제공한다.
이 표는 source crust의 앞부분, 검증된 matched EOS, source core의 뒷부분을 밀도
순서로 기록한다. 행마다 출처와 원래 row index가 명시된다. 재구성 계열에서는
preparation이 명시한 source row 범위에만 준비된 값을 사용하고, 그 밖의 공개
저밀도·고밀도 원본은 그대로 둔다.

이 CSV의 상태는 `DIAGNOSTIC_SOURCE_RETAINING_NOT_ADMITTED`다. prefix/suffix
연결 경계에서 보간이나 외삽을 수행하지 않으며, 그 경계의 평형 조건이나 전체
인과성·열역학적 admission을 주장하지 않는다. 이 검토용 표를 accepted EOS의
추가 개수로 세지 않는다. `SOURCE_RETENTION.json`은 중앙 EOS가 제외한 crust
prefix/core suffix의 실제 row 수와 밀도 범위를 기록한다.

## 단위와 기준 에너지

CSV의 밀도는 fm⁻³, 압력과 에너지 밀도는 MeV fm⁻³, chemical potential은 MeV다.
CoMPOSE header에 기록된 neutron reference rest energy를
\(E_{n,\mathrm{ref}}=m_{n,\mathrm{ref}}c^2\)라고 쓰면, source 물리량은

\[
P=n_B Q_1,\qquad
\epsilon=n_B E_{n,\mathrm{ref}}(1+Q_7),\qquad
\mu_{B,\mathrm{stored}}=E_{n,\mathrm{ref}}(1+Q_3)
\]

로 읽는다. 서로 다른 source의 header를 강제로 같게 바꾸지 않는다.
\((\epsilon+P)/n_B\)는 별도 진단 열로 제공하므로 source의 stored chemical
potential과의 불일치도 숨겨지지 않는다. 패키징은 물리적인
\(\epsilon\mapsto\epsilon+C n_B\) 변환을 적용하지 않는다.

## 재현과 검증

```bash
python scripts/package_maxwell_three_families_20261117.py \
  --campaign data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json \
  --output data/maxwell_three_families_20261117/deliverable_v2
```

output 디렉터리가 이미 존재하면 덮어쓰지 않고 실패한다. 각 ZIP에는
`ARCHIVE_MANIFEST.json`이 들어 있으며 자신의 manifest를 제외한 모든 member의
SHA256과 byte 수를 기록한다. ZIP 생성 후 모든 member를 다시 열어 CRC, 중복
이름 여부, byte 수와 SHA256을 실제로 확인한다.

네 targeted regression은 195-cell census와 ZIP 내부 EOS byte identity, 실패한
candidate compose 제외, reconstructed selected-row overlay의 범위, 변경된
verified EOS의 거부, cross-family source reference의 거부, 명시적 density index에
따른 행 재정렬과 0-based index 및 Fortran D/d exponent 처리를 확인한다. 중복·누락
index는 거부한다. 이들은
패키징 계약 테스트이며 합성 fixture를 물리적 admission 계산으로 주장하지 않는다.
실제 EOS의 물리적 검증 근거는 각 campaign candidate의 native receipt다.
