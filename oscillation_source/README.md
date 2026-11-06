# 중성자별 구조와 비방사 진동

상태방정식에서 평형 별을 구하고 상대론적 Cowling 근사에서 유체 고유진동을 계산한다. 차갑고 비회전인 구대칭 완전유체, 평형 barotropic 압축응답과 명시한 양의 표면압력을 사용한다. 기본 비교는 ell=2의 f와 첫 압력모드 p1이다.

`START_HERE.md`에서 시작한다. 완성된 계산 경로는 `calculation/`, 제공 수치 방법은 `methods/`, 공통 상수·자료형은 `support/`에 있다. 적응 적분·구간별 EOS 평가·양방향 matching을 직접 구현하지 않고 사용할 수 있다. 물리식은 `calculation/src/physics/`에 명시적으로 드러나 있다.

```python
from tools.stellar_methods import solve_pair, run_study
pair = solve_pair("../first_pair")
study = run_study("../full_study")
```

한 pair와 전체 수렴 판정은 다르다. `run_study`는 21개 단일별 계산과 모드·경계·질량·단계법 검사를 연결한다. 같은 명령을 반복하면 같은 소스·입력의 완료 사례만 재사용한다. 아직 필요한 검사가 남으면 `pending`이고, 성공한 전체 비교에서만 `numerically_resolved_in_fixed_model`을 반환한다. 엄밀한 오차 상한이나 관측 정확도의 보장은 아니다.

직접 구현 과제는 `learner/`, 풀이 자료는 전체 패키지의 `solutions/`에 별도로 있다. 기본 계산은 이 경로를 자동 호출하지 않는다. 자기 물리 함수를 연결할 때에는 `make_calculation.py --physics learner/src`로 명시하고, 미구현 함수가 실패하는 것을 풀이로 숨기지 않는다.

강의 노트의 제3·7·9장은 EOS 보간을 열역학적 함수의 정의로 다룬다. 제4·5·6장은 단방향 예제에서 양방향 경계값 문제로 연결한다. 제10장은 이 방법으로 실제 작은 주파수 차이를 분해한다. 전체 강의·계획·과제 PDF는 `documents/`, 원문은 `notes/`와 `curriculum/`, LaTeX 소스는 `documents/latex/`에 있다. `Precision_Study.pdf`는 본문을 대신하는 추가 이론이 아니라 함수 사용 안내다.
