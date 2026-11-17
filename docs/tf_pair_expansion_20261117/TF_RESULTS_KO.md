# C: DD2 native TF의 제한된 저밀도 continuation

기준 commit은 `894cbf1c8315f07ee926c29b83766ffdab2ac7b7`다. 기존 native N16, n_B=0.003 fm^-3 해의 실제 nn/np profile과 cell radius에서 출발하여 0.0025, 0.0020000001 fm^-3로 continuation했다. 두 새 밀도에서 N512 droplet 정지해를 확보했다. 마지막 mesh 차이는 각각 에너지 0.000115891, 0.000148122 MeV/baryon, 압력 0.01810%, 0.001941%다. 이는 유한 branch와 mesh의 계산 증거이며 연속 EOS, 전체 geometry 전역 최소, whole-crust 지원 승인으로 확장하지 않는다.

## 정의와 보존한 물리

생성은 기존 `compose_eos_spatial_tf.f90`와 `compose_eos_tf_models.f90`의 native Fortran으로 수행했다. Python은 보존 profile의 비음수 보간·가중평균 복원, 실행, 독립 검산에만 사용했다. 기존 native physics와 optimizer source의 SHA-256은 시작/종료에 동일하다. 새 `tf_continuation_20261117.f90`는 nn/np 입력을 기존 API에 전달하고 실제 상태를 출력하는 adapter다.

T=0, local density-dependent DD2 RMF, screened electrons, 셀 평균 baryon 보존과 전하 중성, beta equilibrium을 유지했다. 미지수는 nn_i,np_i,log R이며 meson/electron fields는 매 functional 평가에서 푼다. Yp=Ye는 최소화 결과다. 공개 DD2 row217의 실제 n_B=0.0020000001 fm^-3, Ye=0.06988, Ymu=0을 endpoint의 초기 평균조성에만 사용했다. 최종 N512 Ye=0.07005257970은 이 초기값에 고정되지 않았다. 없는 public meson fields, radius, nuclear A/Z/abundance를 채우지 않았다. Phase7은 inner-crust 표시이며 geometry label로 사용하지 않았다.

에너지밀도는 rest energy를 포함한다. n은 fm^-3, R은 fm, ε와 P는 MeV fm^-3, chemical potential은 MeV다. ħc=197.3269804 MeV fm, DD2 nucleon rest energies mn c²=939.56536 MeV, mp c²=938.27203 MeV를 보존한다. Coupling과 rearrangement는 local nn+np에서 평가한다. 경계는 중심 regularity와 outer Neumann, Coulomb mean-zero gauge다. Muon threshold crossing은 별도 저장·거절하며 이번 모든 상태는 crossing=false다.

공간 평균을 ⟨·⟩, dimensionless positive radial operator를 L로 쓰면 native 계산은

ε=⟨εFn+εFp+εFe+Uσ/(ħc)^3+gωωn+gρρ(np−nn)/2−Uωρ/(ħc)^3⟩+Gσ−Gω−Gρ+EC,

Ga=⟨aLa⟩/[2ħc R²], EC=½⟨(np−ne)VC⟩,

∂ε/∂log R=2(−Gσ+Gω+Gρ+EC)

를 사용한다. 여기서 fermion energies는 effective nucleon mass와 electron rest mass를 포함한다. 목적함수의 `(ε−mn c²⟨n⟩)/nbar`는 조건수 개선용이며 반환 ε에 임의 Cn 보정을 하지 않는다. KKT는 양수 밀도에서 μn=μp+μe=μB, 영 밀도에서 해당 부등식이다. 압력은 P=n_B μB−ε이며 parent의 `SYMBOLIC_CHECK.json`과 아래 재최적화 density derivative가 별도로 연결한다.

## 실제 지원과 mesh 계산

출발점은 `data/research_native_20261116/tf/family_campaign/DD2_0.0030_d{1,2,3}_{droplet,hole,uniform}.out` 중 실제 존재하는 여섯 initialization family다. 각 파일 identity와 native provenance는 `INPUT.json`에 있다. 기존 native receipt의 밀도 노드는 .003/.03/.075, mesh는 N16이었다. 이번은 이 유한 native 증거의 아래쪽에 필요한 두 노드를 추가한다. 과거 Python 생산 해나 public table row를 새 native 해로 집계하지 않았다.

각 새 밀도의 시작 family를 비교했다. d1-droplet, d2-hole, d3-uniform 시작은 near-uniform으로 수렴했고 d2-droplet은 rod, d3 시작은 droplet을 만들었다. 초기 family명과 반환 geometry는 별도로 보존했다. coarse N16 에너지가 훨씬 낮아 보이는 현상은 mesh를 늘리면 사라지는 underresolution이었다. 서로 다른 mesh의 에너지를 phase 경쟁으로 비교하지 않았다.

Endpoint n_B=0.0020000001 fm^-3의 droplet:

| N | ε (MeV fm^-3) | P (MeV fm^-3) | R (fm) | Ye | max chemical KKT (MeV) |
|---:|---:|---:|---:|---:|---:|
|64|1.88065861941301|0.00197098930149|42.75740466|0.06959364652|8.42e-6|
|128|1.88066508195265|0.00195517707387|42.53582878|0.07026287826|1.30e-5|
|256|1.88066624226361|0.00195579800230|42.57803009|0.07008168910|1.59e-4|
|512|1.88066653850851|0.00195576004780|42.58456400|0.07005257970|6.34e-5|

N128→256과 N256→512의 Δ(ε/n_B)는 각각 0.00058015545, 0.00014812244 MeV/baryon이다. 이 세 mesh에서 추정한 관측 에너지 수렴차수는 1.96965다. 이 유한 자료의 비율은 약 2차 거동과 부합하지만 엄밀한 continuum error bound는 아니다. 압력 변화는 부호가 바뀌므로 동일한 차수를 강제하지 않고 마지막 실제 차이를 제시한다.

n_B=0.0025 fm^-3에서는 endpoint의 실제 N256 droplet을 density continuation한 뒤 N512로 재최적화했다. N512 ε=2.35138626364943, P=0.00251873421999, R=41.25776996 fm, Ye=0.06164334093이다. N256→512 Δ(ε/n_B)=0.00011589080 MeV/baryon, 상대 ΔP=−0.00018101755다. 최대 KKT=3.55e-4 MeV, radius residual=3.75e-8 MeV로 기존 acceptance 안이다. N512가 자동으로 더 작은 KKT를 보장하지 않는 점도 수치에 드러난다.

같은 N256 endpoint에서 rod−droplet 에너지 차이는 +0.1414414910 MeV/baryon이다. 계산된 두 branch 중 droplet이 낮다. N256 rod는 ε=1.88094912525978, P=0.00173405295862, R=65.41752284 fm이다. Endpoint droplet/rod의 시작 R을 각각 0.9배,1.1배로 바꾸어 다시 최소화한 네 결과는 모두 기존 branch로 돌아왔고, |Δε|≤1.18e-13 MeV fm^-3였다. 이는 관측된 basin의 cell-size 안정성 검사이며 모든 radius/geometry의 전역 최적성을 증명하지 않는다. 모든 geometry의 별도 fine branch가 확보된 것은 아니므로 전역 geometry 선택은 UNRESOLVED로 남는다.

## 실패와 수리의 구분

실제 최소화 64회 중 45회가 stationary로 수용되고 19회는 optimizer mode8로 남는다. 원 mode8 기록은 삭제하거나 성공으로 덮어쓰지 않았다. 해당 반환 profile의 baryon 평균을 정확하게 복원하고 SLSQP state를 새로 초기화한 warm minimization을 별도 실행했다. ftol=1e-14, maxiter=1800, KKT<0.02 MeV, radius residual<2e-4 MeV와 native field 조건을 유지했다. 각 새로운 실패 노드의 반환-profile 재시작을 최대 한 번 수행하는 제한된 episode였다. 일부 coarse family는 그 자리에서 끝냈으며 다른 실제 fine branch에서 계속했다.

Mode8은 이 계산에서 field-PDE 실패와 구별된다. 작은 KKT, cell residual, 물리 constraint에도 optimizer success=false인 원 상태를 보존했다. 내부 SLSQP 원인 자체를 이번 계산에서 재계측하지 않았으므로 과거의 h3 원인 진단을 새 실행의 직접 관측처럼 쓰지 않는다. 새 warm solve의 성공이 실제 수리 근거다.

실행기에는 N32 다음 N64 출력명을 N32로 다시 만들던 결함이 한 번 있었다. 기존 디렉터리 덮어쓰기 방지가 실행을 중단했고, 결과는 그대로 남았다. `test_mesh_label_changes_after_first_refinement` RED→GREEN으로 수정했고 원 refine.log, repaired log를 구분했다. native kernel 결함이나 물리 실패가 아니다. 최종 실행기 회귀 2개는 returncode0, 2 passed다.

## 독립 검산과 열역학

64개 모든 반환 profile을 독립적으로 유지된 Python functional로 평가했다. 최대 차이는 ε 7.33e-15 MeV fm^-3, local chemical potential 1.69e-8 MeV, radius derivative 7.80e-14 MeV fm^-3다. 별도 electron Poisson residual 최대 1.30e-12, charge residual 5.42e-20 fm^-3, baryon residual 2.17e-15 fm^-3다. Python에서 다시 계산한 scalar/vector residual 최대는 각각 1.89e-13,1.68e-14이며 KKT도 다시 계산하여 보존했다. 이 검산은 동일 방정식의 별도 구현 비교이며 모델 자체의 원저자 fidelity 증명은 아니다. 원 실패 profile도 검산에 포함되므로 PDE 확인과 optimizer admission을 혼동하지 않는다.

N256 endpoint droplet에서 h=1e-6 fm^-3의 n0+h,n0+2h,n0+4h를 실제 native 재최적화했다. 모두 같은 droplet branch의 stationary state다. e(n)=ε(n)/n에 대해 전진 3점 미분 P_h=n0²[−3e0+4e(n0+h)−e(n0+2h)]/(2h)를 사용했다.

- Euler P=0.001955797986885388 MeV fm^-3.
- P_h=0.001955797681641538, P_2h=0.001955796529652695 MeV fm^-3.
- |P_h−P_Euler|=3.05e-10 MeV fm^-3, 상대차 1.56e-7.
- |P_h−P_2h|=1.15e-9 MeV fm^-3.

Frozen-profile derivative가 아니다. Native finite KKT와 radius residual의 기여, 최적화 roundoff, 미분 truncation을 단일 정확한 원인으로 귀속하지 않는다. 관측 차이는 마지막 endpoint mesh pressure 차이보다 작다. 이 국소 미분 검산이 전체 밀도구간의 continuous EOS admission을 대신하지 않는다.

## 재현·상태

`tf_continuation/`의 INPUT, FINAL_SOURCE_IDENTITY, 각 STATE/input/stdout/stderr, COARSE/REFINE/REPAIR/FINE/ENDPOINT, STATES, SUMMARY, INDEPENDENT_RESIDUALS, THERMODYNAMIC_CHECK가 실제 근거다. 전체 native 시간은 272.17초, 가장 긴 호출은76.69초였다. Thread 수는 OMP/OPENBLAS/MKL 모두1이다. 소스와 모든 증거 파일의 최종 SHA-256은 MANIFEST.json에 있다. 원자료 INPUT의 script hash는 repair 전 실행기이고 최종 수정본 identity는 별도로 보존했다.

빌드는 `bash fortran_eos/tests/blocker_tf_build.sh build/tf_continuation_20261117` 후 같은 native objects에 새 driver를 연결한다. `REPRODUCE.sh`에 실제 compile/run 명령이 있다. 예: `bash data/tf_pair_expansion_20261117/tf_continuation/REPRODUCE.sh /tmp/tf_C_replay_unique`. 출력 디렉터리가 이미 있으면 중단하며 원 증거를 덮어쓰지 않는다. 선택 범위 안의 producer, validator, 수리 스크립트들은 모두 보존했다. Publication/commit/push는 이 lane에서 수행하지 않았다.

skill_context_used: PhysMath Research Loop1.3.0 `SKILL.md`, `references/POLICY_INTEGRATION.md`, `references/CODING_LOOP.md`를 실제 읽었다. 적용한 지시는 selected project stage C, definition→computation→validation 추적, precompute checkpoint, 실패와 productive repair 구분, typed identity, 한 primary verification pass다. 별도 reviewer tree는 만들지 않았다.

두 추가 밀도에서 branch-specific native 계산과 mesh/cell-size 검산은 완료했다. 전역 geometry 최소, 이산 노드 사이 밀도 보간오차, 전체 TF EOS 구성·승인은 이 결과로 해결되지 않는다. 이 결과가 자동으로 TF matching source seal을 만들거나 support 문자열을 바꾸지 않는다.
