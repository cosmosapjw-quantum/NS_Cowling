# Whole-crust A 결과 사용법

## 결과의 정의

이번 A의 주 baseline은 공개 BSk22 outer의 **pressure-primary 재구성**과 코드 생성
SPG_M1–SLy4의 기존 `quintic_P` 연결을 사용하는 `reference-assisted hybrid`다.
공개 원본 (n,P\) knot와 첫 energy/baryon을 고정한 뒤 제1법칙을 적분하며,
나머지 outer energy 및 Euler chemical potential의 수정량을 기록했다.
앞선 세 계열 sharp-Maxwell 결과에 이 hybrid를 합산하지 않는다.

전체 생산은 **주 baseline 1개와 구성 민감도 표 7개**, 합계 8개의 표로 정의된다.
모든 표가 같은 SPG_M1 및 SLy4 계열을 사용하므로 8개 독립적인 핵 모형으로 세지 않는다.

| 변화 | 의미 |
|---|---|
| `quintic_P` baseline | 주 smooth core–crust 연결; lower endpoint native row 40; surface \(\gamma=5/3\) |
| `lower_endpoint_036` | lower bridge의 오른쪽 endpoint를 native row 36으로 변경 |
| `lower_endpoint_037` | lower bridge의 오른쪽 endpoint를 native row 37로 변경 |
| `surface_gamma_4over3` | lower endpoint row 40을 유지하고 surface 지수만 \(4/3\)으로 변경 |
| `css`의 같은 네 가지 설정 | 기존 CSS core–crust construction을 사용하는 민감도 대조 |

정확한 결과 상태, 행 수, 밀도·압력 범위, 각 파일의 상대 경로는 bundle의 `CATALOG.json`
및 각 모델의 `RESULT.json`을 따른다. `whole_crust_primary.zip`은 primary 하나를 포함하고,
`whole_crust_all_variants.zip`은 여덟 표와 검증·진단 기록을 포함한다.

## 바로 평가하기

ZIP을 새 폴더에 풀고 해당 폴더에서 Python 3, NumPy, SciPy를 사용한다.
제공된 연속 consumer의 실행에는 Fortran 컴파일러나 native executable이 필요하지 않다.

```bash
python evaluate.py --list
python evaluate.py --pressure 0 0.0001 1
python evaluate.py --pressure 0 0.0001 1 --output states.json
```

`--output` 파일이 이미 있으면 덮어쓰지 않고 실패한다. `--model`을 생략하면 primary가
선택된다. 다른 variant는 `--list`에 표시된 정확한 ID로 선택한다.

```bash
python evaluate.py \
  --model BSk22pp__SPG_M1__RG_SLY4__css__baseline \
  --pressure 0.0001 0.1 1
```

반환하는 물리량은 다음과 같다.

| 키 | 의미 | 단위 |
|---|---|---|
| `n` | baryon number density | fm^-3 |
| `p` | 압력 | MeV fm^-3 |
| `epsilon` | rest energy를 포함한 energy density | MeV fm^-3 |
| `mu` | cold barotropic chemical potential | MeV |
| `cs2_over_c2` | \(dP/d\epsilon=c_s^2/c^2\) | 무차원 |
| `enthalpy` | \(h=\int_0^P dP/(\epsilon+P)\) | 무차원 |

압력은 `--pressure`에 **MeV fm^-3**로 넣는다. 출력 범위 밖의 양수 압력, 음수,
NaN/Infinity는 거절한다. 현재 A의 실제 upper core 밀도는 **0.5 fm^-3**이며 고밀도
외삽은 정의하지 않았다.

## 표면과 경계

각 모델의 `REPRESENTATION.json`은 surface, outer pressure law, lower bridge,
upper interpolation의 입력을 갖는다. `WholeCrustBarotrope.from_contract`가 이를 복원한다.
이 consumer는 **(P=0\)에서 (n=P=\epsilon=0\)**인 analytic surface를 제공한다.
`compose/`의 CoMPOSE 표에는 양의 밀도만 기록되어 있으므로, 첫 행을 별의 정확한 표면으로
간주하지 않는다. Export cutoff의 영향은 누락되는 dimensionless enthalpy로 기록했다.

Surface는

\[
P=K n^\gamma,\quad \epsilon=a n+\frac{P}{\gamma-1},\quad
h=\log\!\frac{\mu}{a}
\]

를 사용한다. \(\gamma=5/3\) baseline과 \(4/3\) 대조는 명시적인 surface 모델 선택이다.

공개 outer 원본의 **17개 동위원소 tabulation gap**은 명시적인 effective regularization이다.
그 내부 조성이나 microscopic sharp coexistence를 계산한 것으로 해석하지 않는다.
CoMPOSE에 보존된 isotope 조성은 원본 knot의 참고 조성이고, 재구성된 energy functional에서
새로 최소화한 조성은 아니다. `OUTER_COMPOSITION_NORMALIZATION.json`에는 원본의 heavy-group
표기를 number-fraction 형식으로 해석한 규칙과 보존한 원본을 기록했다.

구체적으로 공개 217개 pure-isotope 행의 raw heavy-group 값은 `Yi=1`이다. 이를 표준
number-fraction으로 그대로 읽으면 frontend의 heavy-nucleus baryon fraction이 `Yh=A`가
된다. 이번 작업은 원본의 순수 동위원소 의미에 맞춰 **`Yh=1`, 출력 CoMPOSE `Yi=1/A`**로
명시적으로 정규화한다. 이 semantic correction은 pressure-primary의 열역학적 energy
재구성과 별도 ledger에 기록하며, 다운로드 원본은 그대로 보존한다.

CSS의 sharp interface 압력에서는 저밀도·고밀도 상태를 구분한다. 해당 transition pressure는
각 모델의 `RESULT.json`과 `REPRESENTATION.json`에서 확인한다.

```bash
python evaluate.py --model MODEL_ID --pressure TRANSITION_PRESSURE --transition-side low
python evaluate.py --model MODEL_ID --pressure TRANSITION_PRESSURE --transition-side high
```

`MODEL_ID`와 `TRANSITION_PRESSURE`를 실제 값으로 바꿔 사용한다. 두 상태 사이의 density gap을
일반 spline으로 연결하지 않는다. Smooth와 CSS 모두 전체 CoMPOSE 표를 한 개의 임의 spline으로
다시 보간하면 제공된 구간별 함수·미분과 달라질 수 있으므로, 기록된 consumer를 사용한다.

## 파일 구조와 추적성

| 경로 | 내용 |
|---|---|
| `models/<ID>/compose/` | 실제 native writer로 출력하고 다시 읽어 검증한 CoMPOSE 표 |
| `models/<ID>/REPRESENTATION.json` | 정확한 구간별 연속 평가 계약 |
| `models/<ID>/ROW_PROVENANCE.json` | 각 출력 행의 surface/public-derived/bridge/native origin |
| `models/<ID>/RESULT.json` | 실제 결과, 경계·인과성·독립 대조 및 readback 결과 |
| `models/<ID>/UNUSED_NATIVE_PREFIX.csv` | 이번 baseline에서 사용하지 않은 native 저밀도 prefix의 진단용 view |
| `OUTER_RECONSTRUCTION.json` | 공개 원본 대비 energy/chemical-potential 수정량 |
| `sources/` | 원본 공개 자료와 이전 native matched source·receipt |
| `code/python_frontend/compose_eos/` | 실제 생산에 사용된 consumer의 Python 소스 |
| `code/fortran_eos/app/whole_crust_kernel.f90` | 새로운 thermodynamic piece를 계산한 native source |
| `figures/` | 연속 consumer로 계산한 그림과 plotting data |
| `MANIFEST.json`, `SHA256SUMS.txt` | 파일 크기·해시 |

`SHA256SUMS.txt`와 `MANIFEST.json`은 두 자체 기술 파일을 제외한 모든 ZIP member를 포함한다.
ZIP 전체의 해시와 실제 압축 해제 후 CLI 검증 결과는 외부 `PACKAGING_RECEIPT.json`에 있다.
전체 bundle에는 과거 실패·원자료 진단도 들어 있지만, 채택된 EOS 목록은 `CATALOG.json`으로
구분된다. Native 실행 명령의 절대 경로는 당시 실행의 provenance로 남겨져 있으며,
휴대 가능한 평가 경로는 모두 bundle 내부의 상대 경로다.

## 다음 순서

다음 원래 순서는 **B. 필요한 core 지원 확장과 source-variant 정책**이다.
이미 생성한 GM1·KDE0v1·SLy4의 0.02–1.0 fm^-3 자료를 재사용 가능한 입력으로 검토하고,
현재 A baseline의 core를 실제 별 모형에 필요한 central-density 범위로 연결해야 한다.
이번 bundle은 EOS의 조성 의존 부력, frozen-composition closure, TOV 또는 Cowling 결과를
포함하지 않는다. 그 관측량과 경계조건의 검증은 이후의 해당 단계에서 수행한다.
