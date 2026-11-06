# 기본 수치 커널의 수식표 계약

## Scalar callback
`subroutine f(x,value,status)`: 입력 x의 단위와 출력 F의 단위는 문제 정의에서 고정한다. 실패 status 또는 비유한 F는 유효한 bracket 점이 아니다. 상위 호출은 mutable global model 대신 수명이 유지되는 contained callback을 사용할 수 있다.

`bisect(f,a,b,root_options,result)`: 순서가 있는 구간 a<=b, 유한·비음수 xtol/rtol/ftol, max_iter>0. 같은 끝점은 그 점에서 계산된 F=0인 경우만 성공한다. exact 잔차, 양의 ftol의 잔차 기준, 양의 좌표 기준 중 사용한 종료 사유를 기록한다. 좌표기준은 반구간 폭<=xtol+rtol max(|lo|,|hi|). 반복 초과·callback 실패·정체는 result.valid=false이며 root를 성공값으로 내보내지 않는다. 중첩 호출은 별도 지역 상태를 쓴다.

## ODE callback과 한 단계
`subroutine rhs(x,y,dydx,status)`: [dydx_i]=[y_i]/[x]. 배열 크기는 같아야 한다. 모든 성분을 채우고 status를 명시한다.

`euler_step`, `rk4_step`, `dopri_trial`는 intent(in)의 x,y를 바꾸지 않는다. 성공할 때만 step_result.y를 채운다. 제공 DOPRI 표는 7개의 RHS 평가를 쓰고 FSAL cache는 적용하지 않는다. error는 h sum_i (b_i-bhat_i) k_i로 계산한다.

## Adaptive 기본값과 실패
오차식과 단계 제어는 본문 §2.3을 따른다. atol_i>0, rtol>=0. 초기 단계·최소 단계·최대 단계는 양의 크기이며 방향은 x1-x0에서 정한다. h_min<=h_max, max_steps는 성공·거절을 합친 시도 예산이다. 기본값은 기록용 설정이며 모든 문제의 정확도를 보장하지 않는다.

계산 완료 상태는 result.valid=true와 status=ok의 조합이다. x/y profile은 initial+accepted, trace는 오차가 평가된 trial만 포함한다. callback 도중 실패한 시도는 attempts에 세되 error trace를 꾸며 넣지 않는다. 끝점 또는 지정 정지점까지의 마지막 성공 단계만 h_min보다 작을 수 있다. 그 단계가 거절되면 더 작은 시도를 강행하지 않는다. callback의 domain failure는 임의의 smoothing·재시도로 감추지 않는다.

`event_schedule(:)`은 진행 방향에 따라 엄격히 정렬된 유한 내부 좌표이다. 중복, 구간 밖과 양 끝점의 중복 선언은 거절한다. x0=x1일 때 유효한 설정이면 callback 없이 초기 상태를 반환한다. 상태 종속 표면 탐색, dense output, 계면 전달, stiffness 전환은 이 API의 현재 기능이 아니다.

## 계약에 대한 추가
수치 상태는 max_iterations, stagnation, step_too_small, max_steps와 nonfinite_result를 구별한다. 각 상태에서는 마지막으로 유효했던 입력과 구간을 기록한다.
