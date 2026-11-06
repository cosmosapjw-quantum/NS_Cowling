# 제공된 방법으로 별을 계산하기

TOV를 구현해 보았고 기본 Cowling 이론을 공부했다면 수치 적분기를 다시 작성할 필요는 없다. 먼저 완성된 물리 예제를 호출하여 두 모드와 수렴표를 얻고, 본문에서 EOS·경계조건·오차의 의미를 확인한다. 자신의 물리 함수는 별도의 명시적 연결로 사용한다.

실행에 필요한 것은 GNU Fortran, fpm, Python과 `requirements.txt`의 Python 패키지다. 확인한 조합은 GNU Fortran 14.2.0, fpm 0.13.0, Python 3.13.5의 Linux 환경이다. `fpm`이 PATH에 없으면 `NS_FPM`에 실제 실행파일의 절대경로를 지정한다. `fpm.toml`만 복사하는 것은 fpm 설치가 아니다. 문서를 읽거나 계산하는 데 LaTeX 설치는 필요하지 않다.

```bash
python -m venv ../ns-env
source ../ns-env/bin/activate
python -m pip install -r requirements.txt
```

## 1. 제공 EOS 방법을 확인한다

```bash
python tools/run_assignment.py A20 --plan
python tools/run_assignment.py A20 --run --output ../A20.json
```

기본 workspace는 `calculation`이다. 소스부터 빌드하여 한 구간의 열역학과 압력 역함수 등을 검사한다. 이 성공이 별의 주파수 수렴을 뜻하지는 않는다.

## 2. 두 별을 계산한다

패키지 루트에서 Python으로 실행한다.

```python
from tools.stellar_methods import solve_pair
pair = solve_pair("../first_pair", level="fine")
print(pair["delta_f_hz"], pair["delta_p1_hz"])
print(pair["physical_claim_status"])
```

저밀도 변형을 준 별과 기준 별을 같은 질량에 맞추고 두 모드의 고유함수를 비교한다. 차이는 변형 별에서 기준 별을 뺀 값이다. 단일 pair는 `pending_additional_checks`다. 결과의 부호를 미리 정답으로 강제하지 않는다.

## 3. 전체 수렴과 민감도를 확인한다

```bash
python tools/stellar_methods.py --output ../my_study --max-new-runs 2
```

두 단일별 계산을 수행한 뒤 멈춘다. 이후 같은 명령을 반복하거나 다음 명령으로 남은 계산을 끝낸다.

```bash
python tools/stellar_methods.py --output ../my_study
```

`study.json`은 판정과 남은 사례, `convergence.csv`는 각 별의 절대 주파수와 차이의 수열이다. 네 정밀도, 다른 적분 단계법, 중심·표면·matching 위치, 질량 잔차와 무변형 대조가 모두 필요하다. 기본 preset은 21개 단일별 실행이다. 완료 여부와 해의 물리적 신뢰성을 구별하며, `numerically_resolved_in_fixed_model`은 이 고정 모델에서의 수치적 분해를 뜻한다.

중단된 시도는 성공으로 재사용하지 않는다. 같은 입력을 다시 시도할 때에는 `--retry-failed`를 지정하면 원래 시도를 `failed_attempts/`에 보존한다. 코드를 바꾸거나 다른 물리 입력을 선택하면 새 출력 폴더를 사용한다. 완료한 계산은 입력과 소스, 출력 내용이 맞는지 확인한 뒤 재사용하므로 과거 숫자로 조용히 대체하지 않는다.

## 4. 자기 물리 함수를 연결한다

직접 구현을 선택하면 `learner/src`의 TOV·Cowling 우변·중심상태·표면잔차를 완성한다. 같은 공개 함수 인터페이스를 지켜야 한다. 단위나 lapse를 맞추지 않은 임의의 TOV 코드를 그대로 연결할 수 있다는 뜻은 아니다.

```bash
python tools/make_calculation.py --physics learner/src \
  --output ../my_physical_calculation
```

네 파일을 모두 명시적으로 선택하고 바이트 그대로 새 계산공간으로 옮긴다. 누락·중복·기존 목적지는 거절한다. 미완성 함수를 완성 예제로 대체하지 않는다. 새 위치에서는 다시 같은 `stellar_methods.py` 명령을 사용한다.

제공된 물리 예제를 그대로 새 위치에 옮기려면 `--physics` 없이 `make_calculation.py --output ../my_calculation`을 사용한다. `solutions/`는 포함하지 않는다. 어떤 물리 함수를 선택했는지는 `CALCULATION.json`에 기록된다.

## 읽을 순서

제3장의 EOS·TOV에서 단위·에너지·lapse를 맞춘다. 제4장의 해석적 양방향 예제와 제5장의 Cowling 유도·중심·표면조건을 읽는다. 제6장에서 모드 이름과 수렴을 구분한다. 제7·9장의 한 함수 보간과 구간 적분을 거쳐 제10장에서 전체 비교를 해석한다. A20과 A21은 제9·10장의 본 과제이며 나중에 별도 수치해석 구현을 요구하는 추가 과제가 아니다.

`documents/Lecture_Notes.pdf`, `Course_Plan.pdf`, `Assignments.pdf`가 전체 강의와 16회 과정, 28개 과제의 큐시트다. 개별 과제는 `documents/cues/`, 빠른 함수 안내는 `documents/Precision_Study.pdf`, 빌드 구성은 `documents/FPM_Note.pdf`에 있다.

## 해 전체의 검사와 중심 시작

전체 `stellar_methods.py` 실행은 마지막에 `variational_checks.json`을 만든다. 이는 주파수를 다시 대입하지 않은 scalar 적분과 유한 표면항을 검사한다. 선택적인 `Model(r0_m=100.0, center_order=2)`는 중심 다음 항을 사용한다. 기본 1 m 결과와 섞지 말고 별도의 출력 경로에서 비교한다. 물리적 의미와 사용조건은 본문 제3·5·6장에 있다.


## 설명을 읽고 계산으로 확인하기

과제 큐시트의 첫 질문에 짧게 답한 뒤 본문과 제공 예제를 읽는다. 처음 예측의 정답 여부는 채점하지 않는다. 계산 뒤에는 어느 식이나 결과가 판단을 바꾸었는지 적고, 다음 만남에서 다른 입력으로 같은 원리를 설명한다. 필요한 함수는 본문과 큐시트가 가리키는 제공 방법을 사용한다.

지도와 피드백의 예는 `documents/Mentor_Guide.pdf`, 과제별 평가 기준은 `docs/ASSESSMENT.md`에 있다. 문헌의 결과와 이 수업에서 선택한 적용은 지도 안내에서 구분한다.
