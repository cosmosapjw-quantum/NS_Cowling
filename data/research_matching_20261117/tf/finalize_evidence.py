"""Finalize the finite TF report from completed native records; no solver here."""
import hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return json.loads((BASE/name).read_text())
def rows(path,name):return [line.split()for line in path.read_text().splitlines()if line.startswith(name+' ')]
summary=read('SUMMARY.json')
assert (BASE/'reproduction/exit_code.txt').read_text().strip()=='3'
assert (BASE/'continuation_repaired/exit_code.txt').read_text().strip()=='0'
prior=ROOT/'data/research_coexistence_20261116/tf/output_n256/FSU2R.stdout'
current=BASE/'reproduction/FSU2R.stdout'
a,b=rows(prior,'STATE'),rows(current,'STATE')
assert a==b and len(b)==30
assert all(x[3]=='STATIONARY'for x in b[:-1]) and b[-1][3]=='REQUIRES_REFINEMENT'
exact=rows(BASE/'original_mode_probe/stdout.txt','STATE')[0]
assert exact[3:]==b[-1][5:]
diag=rows(BASE/'original_mode_probe/stdout.txt','SLSQP_DIAGNOSTICS')[0]
restart=rows(BASE/'restart_probe/stdout.txt','SLSQP_DIAGNOSTICS')[0]
assert diag[1]=='8'and restart[1]=='0'
assert float(diag[12])>float(diag[13])and float(restart[12])<float(restart[13])
identity=read('SOURCE_IDENTITY.json')
for record in identity['files']:
 assert record['byte_identical']and sha(ROOT/record['path'])==record['current_sha256']
for name,state in [('REPRODUCTION.json','PRESERVED_EXACT_FAILURE_REPRODUCED'),('REPAIR.json','FINITE_CHORD_COMPLETED')]:
 d=read(name);d.update(state=state,finished_unix=time.time())
 if name=='REPRODUCTION.json':
  d.update(exit_code=3,states=30,stationary_states=29,root_and_halfstep_stationary_states=28,
   prior_state_rows_byte_identical=True,first_failed_optimizer_mode=8)
 else:d.update(chord_exit_code=0,chord_stationary_states=10,minimum_sampled_interior_gap_MeV_fm3=summary['models']['FSU2R']['minimum_sampled_interior_gap_MeV_fm3'])
 (BASE/name).write_text(json.dumps(d,indent=2)+'\n')
evidence={'scope':'finite two-neutral-branch TF coexistence, no global or continuum certification',
 'baseline_head':'478c7b5de1514ac512f91156f36551478f36f688','state':'COMPLETE_FINITE_RESEARCH_REPAIR',
 'production_sources_modified':False,'tolerances_relaxed':False,
 'original_failure':{'mode':8,'iterations':int(diag[2]),'hessian_resets':int(diag[3]),'evaluations':812,
 'normalized_constraint_error':float(diag[5]),'max_bound_violation':float(diag[6]),
 'h3':float(diag[12]),'termination_threshold':float(diag[13]),'ftol':float(diag[15]),
 'state_values_identical_to_prior_failure':True},
 'fresh_restart':{'mode':0,'h3':float(restart[12]),'termination_threshold':float(restart[13]),'iterations':int(restart[2]),'evaluations':286,
 'note':'Fresh optimizer initialization from returned alternate failed state; actual mode0, no acceptance override'},
 'retained_prior_state_rows_identical':30,'root_h_half_stationary_states':28,
 'completed_chord_native_states':10,'completed_minimization_attempts_this_lane':45,'stationary_attempts_this_lane':41,'mode8_attempts_this_lane':4,'aborted_ieee_diagnostic_probes':1,'native_field_readback_profiles':10,'reproduction_build_exit':int((BASE/'reproduction_build.exit').read_text()),
 'external_indicator':'IEEE_INVALID traced to deliberate vendor unused-multiplier NaN at slsqp_optmz.f:787',
 'source_identity':identity,'created_unix':time.time()}
evidence['files']={str(p.relative_to(ROOT)):sha(p)for p in BASE.rglob('*')if p.is_file()and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.o','.mod','.so') and p.read_bytes()[:4]!=b'\x7fELF' and p.name not in ('EXECUTION.json','finalization.log','finalization.exit')}
(BASE/'EXECUTION.json').write_text(json.dumps(evidence,indent=2,allow_nan=False)+'\n')
lines=[]
for model,o in summary['models'].items():
 lines.append(f"| {model} | {o['n_low_fm3']:.12g} | {o['n_high_fm3']:.12g} | {o['mu_t_MeV']:.12g} | {o['P_t_MeV_fm3']:.12g} | {o['minimum_sampled_interior_gap_MeV_fm3']:.8g} |")
fs=summary['models']['FSU2R'];tm=summary['models']['TM1e']
readback=fs['retained_profile_field_readback']
report=f"""# FSU2R N256 finite coexistence repair — 2026-11-17

기준 HEAD는 478c7b5de1514ac512f91156f36551478f36f688이다. 기존 FSU2R N256 실패의 실제 원인을 재현했고, 동일한 물리 functional·경계·optimizer 허용오차로 재최적화하여 slab/tube 공존 후보의 유한 chord 검증을 마쳤다. TM1e N256 결과는 변경되지 않은 물리 의존성을 확인해 재사용했다. 전체 EOS 또는 연속 mesh 극한을 승인한 결과가 아니다.

## 실제 수치 결과

각 phase는 셀 평균 전하 중성을 만족하는 T=0 TF 정지해이다. 전자장은 screened Poisson 방정식으로 풀며, 점별 n_p=n_e 조건을 가정하지 않는다. 두 branch의 공존 조건은 P_low=P_high, mu_B,low=mu_B,high이며 혼합 chord는 epsilon=mu_t n_B-P_t이다. 거시적 계면 비용/부피는 무시한다. global two-charge Gibbs 해와 동일하다고 주장하지 않는다.

| 모형 | n_low (fm^-3) | n_high (fm^-3) | mu_t (MeV) | P_t (MeV fm^-3) | 최소 표본 내부 gap (MeV fm^-3) |
|---|---:|---:|---:|---:|---:|
{chr(10).join(lines)}

FSU2R root 잔차는 Delta mu={fs['root_delta_mu_MeV']:.10g} MeV, Delta P={fs['root_delta_P_MeV_fm3']:.10g} MeV fm^-3이다. h/2 검산은 각각 {fs['halfstep_delta_mu_MeV']:.10g}, {fs['halfstep_delta_P_MeV_fm3']:.10g}이며 기존 2e-6 MeV, 2e-7 MeV fm^-3 기준을 유지했다. 원 실행의 30개 STATE 행은 재현 실행과 byte-identical이며, 처음 28개는 root/h/2 정지 계산이고 마지막 두 개는 최초 chord의 slab 성공과 tube 실패다.

복구 계산은 다섯 밀도에서 두 branch를 실제로 다시 최소화했다. 10개 모두 optimizer mode0, active radius/density bound 없음, STATIONARY 상태다. 세 내부 밀도점의 두 branch 에너지는 chord보다 높았다. 이는 조사한 두 유한 branch/mesh/초기값에 대한 결과이며, 모든 geometry의 전역 최소나 optimizer 에너지 오차의 엄밀한 상한을 증명하지 않는다. N128→N256의 관측 변화는 SUMMARY.json에 보존했다.

## 실패 원인과 구성적 수리

원 실패 profile을 만드는 입력을 별도의 진단용 Fortran 실행으로 그대로 반복했다. 원 실패와 물리값이 일치했고, SLSQP mode8(positive directional derivative), 223 iterations/812 evaluations, 6 Hessian resets였다. bounds 위반은 0이지만 내부 정규화 equality 잔차 h3={float(diag[12]):.14g}가 변경하지 않은 종료 문턱 10*ftol={float(diag[13]):.14g}를 넘었다. 이 때문에 작은 field/KKT/constraint 물리 잔차에도 optimizer success가 false였다. 이는 장 방정식 실패나 모델의 물리적 불안정 증거가 아니다.

다른 보존 초기값에서도 같은 mode8 메커니즘을 재현했다(h3=1.0502709812953981e-13). 그 반환 profile에서 SLSQP를 새로 초기화한 실제 재최적화는 mode0, 41 iterations/286 evaluations, h3={float(restart[12]):.14g}로 종료했다. 실패 상태의 플래그를 덮어쓰지 않았고 ftol=1e-14, maxiter=1800 및 기존 물리 판정은 유지했다. 그 실제 정지 profile에서 forward continuation을 수행해 위 10개 chord 표본을 얻었다. 원 실패, 대체 초기값 실패, 진단 반복, 성공 재시작, 복구 continuation 로그를 각각 보존했다. 이번 lane의 완료된 최소화 호출은 45회(정지 41회, mode8 4회)이며 같은 실패의 진단 반복을 포함한 실행 횟수다. 별도 IEEE trap 진단 1회는 중간 중단되었고, profile field 재평가 10회는 최소화 횟수에 넣지 않았다.

production solve_tf_family_native에는 반환 profile을 사용한 별도의 warm refinement와 모든 시도 보존이 이미 구현되어 있다. 따라서 production optimizer/backend를 바꿀 이유는 없었다. 즉시 중단하던 연구용 공존 driver의 continuation 절차만 보완했다.

## IEEE flag와 profile readback

IEEE_INVALID는 성공 상태에서도 나타났다. -ffpe-trap=invalid 진단은 vendor/slsqp/slsqp_optmz.f:787의 사용하지 않는 bound-multiplier NaN 생성(0/0)을 직접 가리켰다. 이 flag를 TF 물리 실패와 혼동하지 않았고 vendor를 수정하지 않았다.

보존된 10개 nn/np/radius profile을 native field solver로 다시 평가했다. scalar/vector/electron 잔차 최댓값은 각각 {readback['scalar_pde_residual_max']:.8g}, {readback['vector_pde_residual_max']:.8g}, {readback['electron_pde_residual_max']:.8g}로 기존 기준을 만족했다. epsilon 재평가 최대 차이는 {readback['maximum_absolute_re_evaluation_difference']['epsilon_MeV_fm3']:.8g} MeV fm^-3이다. 이는 동일 native 방정식의 저장 profile readback이며 독립 물리 구현이라고 부르지 않는다. mu_B 비교는 solve_spatial_tf와 같은 active-density KKT multiplier 관례로 수행했다. evaluator의 raw neutron-volume 평균은 다른 유한-KKT 관측량이며, 그 차이를 별도로 기록했다. root의 작은 수치 잔차를 물리적 mu 오차 상한으로 해석하지 않는다.

두 후보 각각 endpoint/midpoint의 3행 혼합 chord를 Fortran으로 작성하고 Fortran으로 다시 읽었다. candidate_chords/ 아래 CSV는 연구용 열역학 표다. 혼합상의 mu_Q 관례와 핵 조성이 제공되지 않았으므로 abundance와 presence mask를 만들지 않았다. 생산 EOS를 교체하지 않았다.

## 재현과 근거

- SOURCE_IDENTITY.json: 기존 8개 producer/physics 파일 hash 동일.
- REPRODUCTION.json, reproduction/: 원 실패까지 실제 다시 실행한 전체 기록과 각 profile.
- mode_probe/, original_mode_probe/, restart_probe/: 내부 종료조건 및 성공 재시작.
- continuation_repaired/: 실제 10개 정지해와 모든 profile.
- SUMMARY.json, profile_readback/, candidate_chords/: root/h/2, 유한 chord, N128→N256 비교와 native I/O.
- EXECUTION.json: 종료 상태, 실행 기록과 결과 파일 checksum.

소스 checkout에서 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1을 설정한다. bash fortran_eos/tests/tf_finish_build.sh build/tf_finish로 재현 실행파일들을 만든다. 각 실행의 입력과 명령은 REPRODUCTION.json/CONTINUATION.json 및 probe RECEIPT.json에 있다. 내부 SLSQP 진단은 별도 build 디렉터리에 bash fortran_eos/tests/tf_finish_diagnostic_build.sh를 실행한다. 원자료를 덮어쓰지 않는 새 작업 디렉터리에서 실행한다.

선택된 연구 경로는 실패 재현→내부 원인 확인→실제 warm 재최적화→유한 chord 계산→native 저장/readback이다. PhysMath Research Loop1.3 main과 부모가 제공한 변경 없는 정책/단계 문맥을 사용했으며 추가 검토 트리를 만들지 않았다. 생산 계산은 Fortran, Python은 실행·자료 검증·표현만 수행했다.
"""
dest=ROOT/'docs/research_matching_20261117/tf/REPORT_KO.md'
dest.write_text(report)
print(json.dumps({'state':evidence['state'],'report':str(dest),'models':summary['models'],'report_sha256':sha(dest)},indent=2))
