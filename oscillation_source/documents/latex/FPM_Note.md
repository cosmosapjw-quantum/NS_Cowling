---
title: "fpm.toml의 역할"
subtitle: ""
date: ""
lang: ko-KR
---

`fpm.toml`은 물리식을 더하는 파일이 아니라 소스·실행 대상·의존성을 정하는 빌드 설정이다. 기본 계산은 다음 세 묶음을 사용한다.

- `support/`: 수치 자료형, 상수, 입출력과 작은 공통 연산.
- `methods/`: 적응 적분·근 찾기·EOS 구간 표현과 역함수. `support`를 사용한다.
- `calculation/`: TOV·Cowling 물리식과 제공된 구간 적분·양방향 연결. 앞의 두 묶음을 사용한다.

따라서 기본 계산에 필요한 manifest는 세 개다. `calculation`에서 빌드하면 의존 라이브러리도 함께 컴파일된다. 학생이 별도의 적분기 파일을 작성하거나 세 디렉터리를 각각 수동 링크할 필요는 없다.

```bash
cd calculation
fpm build --profile release
fpm test cell_checks --profile release
```

전체 자료에는 직접 구현용 `learner/`와 풀이 `solutions/`의 manifest도 있다. 이 둘은 같은 Fortran 모듈 이름을 쓰므로 동시에 하나의 library로 링크하지 않는다. 기본 `calculation`은 이 두 경로를 의존하지 않는다. 자기 물리 함수를 연결할 때는 지정한 파일을 새 계산공간에 명시적으로 옮긴다.

```bash
python tools/make_calculation.py --physics learner/src \
  --output ../my_physical_calculation
```

Makefile로 같은 소스와 의존성을 지정하는 것도 가능하다. fpm은 빌드를 정리할 뿐 물리식이나 수렴을 보증하지 않는다. 이 구성은 fpm manifest의 library·executable·local dependency 기능을 사용한다. 상세 형식은 `https://fpm.fortran-lang.org/spec/manifest.html`을 참고한다.
