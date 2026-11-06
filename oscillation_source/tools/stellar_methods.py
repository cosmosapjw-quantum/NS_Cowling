#!/usr/bin/env python3
"""Call the supplied Fortran methods; keep the physical problem explicit.

solve_pair is an exploratory calculation. Only run_study with all required
checks can return numerically_resolved_in_fixed_model. All paths are local.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import csv
from dataclasses import dataclass
import math
from pathlib import Path
import sys
import time
import uuid

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from tools import precision_study as kernel
from tools.atomic_io import write_json,read_json,digest

LEVELS={'coarse':(1e-7,.2),'medium':(1e-9,.1),'fine':(1e-11,.05),'audit':(1e-13,.025)}
REQUIRED_GROUPS=('main','method','boundary','mass','null')

@dataclass(frozen=True)
class Model:
    """SI controls of the supplied continuous barotropic model family."""
    target_mass_kg: float=kernel.MASS
    p_cut_pa: float=1e20
    r0_m: float=1.0
    split_index: int=69
    match_index: int=72
    table_path: Path|None=None
    center_order: int=0

    def parameters(self)->dict:
        if not all(math.isfinite(v) and v>0 for v in (self.target_mass_kg,self.p_cut_pa,self.r0_m)):
            raise ValueError('Mass, cutoff and central radius must be finite and positive')
        if type(self.split_index) is not int or type(self.match_index) is not int or min(self.split_index,self.match_index)<2:
            raise ValueError('EOS indices must be integers >= 2')
        if type(self.center_order) is not int or self.center_order not in (0,2):
            raise ValueError('center_order must be 0 or 2')
        r={'target_mass_kg':self.target_mass_kg,'p_cut_pa':self.p_cut_pa,'r0_m':self.r0_m,
           'split_index':self.split_index,'match_index':self.match_index}
        if self.center_order: r['center_order']=self.center_order
        if self.table_path is not None: r['table_path']=str(Path(self.table_path).resolve())
        return r

def solve_star(output:Path|str,*,alpha:float=0.0,level:str='fine',method:str='dopri',
               model:Model=Model(),workspace:str='calculation',timeout:float=180)->dict:
    """Solve one declared model. This alone does not establish convergence."""
    if level not in LEVELS: raise ValueError('Choose coarse, medium, fine or audit')
    if method not in {'dopri','rk4'}: raise ValueError('Choose dopri or rk4')
    if not math.isfinite(alpha): raise ValueError('alpha must be finite')
    rtol,cap=LEVELS[level]
    parameters={**model.parameters(),'alpha':float(alpha),'method':method,
                'rtol':rtol,'mode_rtol':rtol,'step_cap':cap if method=='dopri' else .005}
    return kernel.run_one(Path(output),parameters,workspace,timeout)

def solve_pair(output:Path|str,*,alpha:float=.1,level:str='fine',model:Model=Model(),
               workspace:str='calculation',timeout:float=180)->dict:
    """Compute a same-mass pair and its mode identity, not an error bound."""
    if level not in LEVELS or not math.isfinite(alpha): raise ValueError('Invalid level or alpha')
    model.parameters();out=Path(output).resolve()
    solve_star(out/'reference',alpha=0.,level=level,model=model,workspace=workspace,timeout=timeout)
    solve_star(out/'variant',alpha=alpha,level=level,model=model,workspace=workspace,timeout=timeout)
    r=kernel.load_pair(out/'reference',out/'variant')
    r['physical_claim_status']='pending_additional_checks'
    write_json(out/'pair.json',r)
    return r

def absolute_convergence(rows:list[dict])->dict:
    """Reject common-mode drift that cancels in a paired difference.

The last two changes are compared with a root/roundoff scale, not with the
physical signal. This diagnostic is empirical, not a convergence theorem.
"""
    if len(rows)<3: return {'status':'insufficient_levels'}
    result={'status':'absolute_convergent','modes':{}}
    for side in ('reference','variant'):
        for prev,nxt in zip(rows,rows[1:]):
            if any(not (0<nxt[side][k]<prev[side][k]) for k in ('rtol','mode_rtol','step_cap')):
                return {'status':'invalid_refinement_schedule'}
        for mode in ('f','p1'):
            vals=[r[side][mode+'_hz'] for r in rows]
            if not all(math.isfinite(v) and v>0 for v in vals): return {'status':'invalid_numerical_data'}
            jumps=[abs(b-a) for a,b in zip(vals,vals[1:])]
            floor=64*sys.float_info.epsilon*max(vals)+max(r[side][mode+'_bracket_hz'] for r in rows)
            growing=jumps[-1]>max(jumps[-2],floor)
            result['modes'][side+'_'+mode]={'frequencies_hz':vals,'changes_hz':jumps,
                'algorithmic_scale_hz':floor,'convergent':not growing}
            if growing:result['status']='absolute_convergence_unresolved'
    return result

@contextmanager
def _study_lock(out:Path):
    # flock is released by the OS if the interpreter is interrupted.
    try: import fcntl
    except ImportError as e: raise kernel.StudyError('This resumable runner requires a POSIX system') from e
    with (out/'.study.lock').open('a+') as stream:
        try: fcntl.flock(stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as e: raise kernel.StudyError('Another process is using this study directory') from e
        try:yield
        finally:fcntl.flock(stream,fcntl.LOCK_UN)

def _recipe(workspace:str)->dict:
    return {'workspace':workspace,'source_sha256':kernel.source_digest(workspace),
            'method_api_sha256':digest(__file__),'runner_sha256':digest(kernel.__file__),
            'eos_sha256':digest(ROOT/'data/normalized/compose_69.csv'),
            'integral_checker_sha256':digest(ROOT/'tools/variational_check.py'),
            'integral_eos_sha256':digest(ROOT/'tools/native_eos_reference.py'),
            'groups':{k:kernel.experiment_groups()[k] for k in REQUIRED_GROUPS}}

def _finite_tree(value):
    kernel._require_finite(value)

def _report(out:Path,recipe:dict)->dict:
    pending=[];done=[]
    for group in REQUIRED_GROUPS:
        for name,params in recipe['groups'][group]:
            p=out/name
            if not (p/'execution.json').is_file() or read_json(p/'execution.json').get('status')!='complete':
                pending.append(name);continue
            rec=read_json(p/'execution.json')
            if rec['workspace']!=recipe['workspace']:raise kernel.StudyError('Different physics workspace in study')
            for key,value in params.items():
                if rec['parameters'].get(key)!=value: raise kernel.StudyError('Experiment controls changed: '+name+' '+key)
            kernel._verified_result(p);done.append(name)
    rows=[]
    for i in range(4):
        names=[f'L{i}_a{a}' for a in (0.0,0.1)]
        if all(n in done for n in names):rows.append(kernel.load_pair(*(out/n for n in names)))
    r=kernel.compare_sequence(rows);r['paired_status']=r['status'];r['levels']=rows
    r['absolute_convergence']=absolute_convergence(rows)
    # A partial late experiment is preserved, not passed to a completed-pair reader.
    if len(rows)==4 and not pending:r.update(kernel.supplementary_checks(out,rows))
    else:r['physical_claim_status']='pending_additional_checks'
    if r['absolute_convergence']['status'] not in {'absolute_convergent','insufficient_levels'}:
        r['physical_claim_status']=r['absolute_convergence']['status']
    if pending:r['physical_claim_status']='pending_additional_checks'
    if not pending and r['physical_claim_status']=='numerically_resolved_in_fixed_model':
        from tools.variational_check import audit_study
        integral=audit_study(out,r['numerical_budget'])
        r['integral_check']={k:v for k,v in integral.items() if k!='levels'}
        if integral['status']!='integral_consistent':r['physical_claim_status']=integral['status']
    r['status']='pending' if pending else r['physical_claim_status']
    r['completed_runs']=done;r['pending_runs']=pending
    r['limitations']='Empirical fixed-model numerical evidence; not an interval enclosure, EOS uncertainty or observational accuracy.'
    _finite_tree(r);write_json(out/'study.json',r)
    write_json(out/'paired_report.json',r)
    if rows:
        with (out/'convergence.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['level','rtol','step_cap','reference_f_hz','variant_f_hz','delta_f_hz','reference_p1_hz','variant_p1_hz','delta_p1_hz'])
            for i,q in enumerate(rows):w.writerow([i,q['reference']['rtol'],q['reference']['step_cap'],q['reference']['f_hz'],q['variant']['f_hz'],q['delta_f_hz'],q['reference']['p1_hz'],q['variant']['p1_hz'],q['delta_p1_hz']])
    return r

def run_study(output:Path|str,*,group:str='all',workspace:str='calculation',
              max_new_runs:int|None=None,retry_failed:bool=False,timeout:float=180)->dict:
    """Run/resume the fixed comparison, recording each completed star separately.

Completed runs are content-verified before reuse. A failed attempt is never
converted to a success. retry_failed=True preserves it under failed_attempts
before attempting the same input again. Changing sources or controls requires
a new output directory. max_new_runs bounds one session, not the tolerances.
"""
    if group not in (*REQUIRED_GROUPS,'all'):raise ValueError('Unknown study group')
    if max_new_runs is not None and (type(max_new_runs) is not int or max_new_runs<0):raise ValueError('max_new_runs must be nonnegative')
    out=Path(output).resolve()
    if out==ROOT or out.is_relative_to(ROOT):raise ValueError('Output must be outside the source tree')
    out.mkdir(parents=True,exist_ok=True)
    with _study_lock(out):
        recipe=_recipe(workspace);plan=out/'study_plan.json'
        # JSON normalization prevents tuple/list differences across a resumed call.
        import json
        recipe=json.loads(json.dumps(recipe))
        if plan.exists():
            if read_json(plan)!=recipe:raise kernel.StudyError('Study source or controls changed; choose a new output directory')
        else:
            if any(p.name!='.study.lock' for p in out.iterdir()):raise kernel.StudyError('Nonempty output has no study plan')
            write_json(plan,recipe)
        count=0
        try:
            for g in REQUIRED_GROUPS:
                if group!='all' and g!=group:continue
                for name,params in recipe['groups'][g]:
                    dest=out/name;recfile=dest/'execution.json'
                    complete=recfile.is_file() and read_json(recfile).get('status')=='complete'
                    if not complete and max_new_runs is not None and count>=max_new_runs:continue
                    if dest.exists() and not complete:
                        if not retry_failed:raise kernel.StudyError('Incomplete attempt retained: '+name+'; use retry_failed to try again')
                        failed=out/'failed_attempts';failed.mkdir(exist_ok=True)
                        dest.rename(failed/(name+'_'+uuid.uuid4().hex))
                    kernel.run_one(dest,params,workspace,timeout)
                    if not complete: count+=1
                    _report(out,recipe)
                    print(name, 'verified reuse' if complete else 'computed',flush=True)
            result=_report(out,recipe);result['new_runs_this_call']=count
            write_json(out/'study.json',result)
            return result
        except (Exception,KeyboardInterrupt) as e:
            write_json(out/'last_interruption.json',{'error':str(e),'type':type(e).__name__,'completed_new_runs':count,'time_unix':time.time()})
            raise

def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--group',choices=(*REQUIRED_GROUPS,'all'),default='all')
    ap.add_argument('--workspace',choices=('calculation','learner','solutions'),default='calculation')
    ap.add_argument('--max-new-runs',type=int)
    ap.add_argument('--retry-failed',action='store_true')
    a=ap.parse_args(argv)
    try:
        result=run_study(a.output,group=a.group,workspace=a.workspace,max_new_runs=a.max_new_runs,retry_failed=a.retry_failed)
        print(result['status'], 'remaining:',len(result['pending_runs']),flush=True)
        return 0 if result['status'] in {'pending','numerically_resolved_in_fixed_model'} else 3
    except (OSError,ValueError,kernel.StudyError) as e:print(str(e),file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
