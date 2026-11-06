# 별과 두 모드를 계산하기

완성된 물리 예제를 제공된 수치 방법으로 적분한다. `src/physics/`의 식과
`src/models/`의 저밀도 변형은 모두 공개되어 있다. 이 경로의 성공은 `learner/`의
직접 구현 완료를 의미하지 않는다. 자신의 물리식을 사용할 때는 해당 함수의
입출력·단위를 유지하고 수정된 소스에서 다시 계산한다. 다른 풀이로 자동 전환하지 않는다.

`use ns_stellar_methods`로 구간 적분과 양방향 matching을 호출한다.
전체 수렴 실험은 상위 디렉터리에서 `python tools/stellar_methods.py --output ../my-study`로 실행한다.
