# 해석적 shooting의 수식표: 제4장 / 

## 정의와 단위

$x\in[0,L]$, $L>0$이며 $x,L,y$의 단위는 m이다. $A=y'(0)\ne0$는 무차원이다. $k\ge0$의 단위는 m$^{-1}$이며 이 문제의 고유값은 $k^2$이다. 검색 결과의 $k$를 Hz 또는 중성자별의 $f,p_1$이라고 부르지 않는다.

\[
y''+k^2y=0,\qquad y(0)=0,\quad y(L)=0.
\]

사격법(shooting)의 초기값 문제는 $y'(0)=A$로 정규화한다. 독립적인 해석해는

\[
y(x;k)= A\frac{\sin(kx)}{k}\quad(k>0),\qquad y(x;0)=Ax.
\]

따라서 $k_n=n\pi/L$, $n=1,2,\ldots$이며 $k=0$은 비자명한 고유해가 아니다. 이 해석식은 pytest 및 해석용 도구의 oracle에만 사용하고 수치 solver에는 넣지 않는다.

## 실제 계산과 잔차

$u=y/A$, $v=y'/A$를 적분한다.

\[
u'=v,\qquad v'=-k^2u,\qquad (u(0),v(0))=(0,1).
\]

각 함수의 입력과 결과는 같은 SI 길이 단위를 사용한다. RK4 또는 adaptive DP5(4)를 호출하고 $x=L$ 도달과 유한 상태를 확인한다. 상태를 해석적 sin으로 대체하지 않는다. 출력에서는 $y=A u$, $y'=A v$를 복원하고 표현 범위를 검사한다.

\[
\mathcal R(k)=\frac{y(L;k)}{|A|L}=\operatorname{sgn}(A)\frac{u(L;k)}L.
\]

잔차는 무차원이다. 양의 진폭 배율에는 불변이며 음의 배율에는 전체 부호가 뒤집힌다. 이 방향은 한 스캔 동안 고정한다. 진폭의 절댓값을 바꾸어 오차 제어가 달라지지 않도록 실제 IVP는 unit-amplitude 변수에서 푼다. $\mathcal R(0)=\operatorname{sgn}(A)$이다.

$y(L)/y'(L)$는 이 문제의 대체 잔차로 사용하지 않는다. $kL=\pi/2$에서 분모가 0이 되어 root가 아닌 pole이 부호를 바꿀 수 있다.

## Scan 계약

`linear_scan_grid(lo,hi,n,grid,status)`는 유한하고 엄격히 증가하는 검색 좌표를 만든다. `build_brackets(trials,zero_tol,result)`는 실제 수치 solver와 분리된 순수한 기록 해석 단계이다. `trial_record_t`의 `parameter,residual,orientation,valid,status`를 모두 사용한다.

유효 상태이면서 성공했고 유한한 잔차를 가진 **인접** trial만 sign bracket을 만든다. 부호 비교에 함수값의 곱을 사용하지 않는다. 실패/NaN/Inf 행을 삭제한 뒤 이웃을 다시 정하지 않는다. 이 행은 barrier로 보존한다.

작은 잔차는 root 후보이다. 이미 sign-change cell의 끝점에 들어간 near-zero sample은 따로 중복 등록하지 않는다. 정확한 zero sample의 양옆이 성공했고 반대 부호면 두 cell에 걸친 bracket을 만든다. 그 조건이 없으면 `CAND_SAMPLE`로 보존하고 자동으로 고유값을 선언하지 않는다. 이 처방은 검색 경계의 root, 접하는 root, 같은 cell 안의 여러 root를 모두 찾는 완전성 증명이 아니다.

서로 이웃한 두 sign bracket은 일반적으로 합치지 않는다. 공유하는 sample 양옆에 실제로 서로 다른 root가 있을 수 있다. Generic scan의 성공은 기록 처리의 성공이며 모든 trial의 성공이나 spectrum completeness와는 다른 값이다.

## Refinement 계약

`refine_string(L,candidate,options,result)`는 근 찾기 전에 실제 IVP 잔차를 다시 평가한다. 외부 기록의 부호를 곧바로 믿지 않는다. `root.ftol=0`으로 좌표 bracket 종료를 사용하며, 최종 root에서 IVP를 다시 적분해 별도의 dimensionless `residual_tol`도 확인한다. `xtol`은 m$^{-1}$, `rtol`은 무차원이다.

`sample_candidate`는 bracket 없이 자동 승격하지 않는다. Root 폭, 경계 잔차, IVP 설정에 대한 안정성은 다른 검사이다. 이 함수의 유효한 반환은 지정한 수치 문제의 결과이지 엄밀한 참 고유값의 enclosure가 아니다.

`scan_string`은 모든 trial과 candidate별 refinement 상태를 보존한다. 실패가 있는 spectrum을 전체 `ok`로 보고하지 않으며 성공한 부분은 typed record의 진단으로 남긴다. CLI에서는 전체 실패 후 `k1_per_m` 같은 정상 결과 열을 내보내지 않는다. 반환 root의 순서는 요청 구간 내부 순서이고, 구간이 $k_1$보다 위에서 시작하면 첫 결과는 기본 모드가 아니다.

## API 경계

- `support/src/ns_scan_types.f90`: 자료형만 제공. 스캔이나 물리 solver는 없다.
- `solutions/src/modes/ns_scan.f90`, `ns_shooting_toy.f90`: 완성 해설.
- `learner/src/modes/`: 동일 signature의 직접 구현 슬롯. `not_implemented`를 반환한다.
- `ns_shooting_cases.f90`: 공통 입출력과 합성 결함 사례. 항상 선택된 workspace의 커널을 호출한다.

각주파수를 탐색 변수로 사용할 때에는 입력 단위를 별도로 명시한다.
