#!/usr/bin/env python3
"""Run native-cell Fortran calculations and examine matched mode differences.

Python handles files, verification and mode-shape diagnostics, not stellar ODEs.
All runs select one workspace explicitly. Interrupted directories are preserved.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from tools.atomic_io import write_json, read_json
MASS = 2.644606257493485e30
class StudyError(RuntimeError): pass

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source_digest(workspace: str) -> str:
    if workspace not in {'learner','solutions','calculation'}: raise StudyError('Unknown workspace')
    h=hashlib.sha256()
    for base in [ROOT/workspace,ROOT/'support'] + ([ROOT/'methods'] if workspace=='calculation' else []):
        if not base.is_dir(): raise StudyError(f'Workspace is absent: {base.name}')
        files=[base/'fpm.toml',*base.glob('src/**/*.f90'),*base.glob('app/*.f90')]
        for p in sorted(files):
            h.update(str(p.relative_to(ROOT)).encode());h.update(p.read_bytes())
    return h.hexdigest()

def _fpm() -> str:
    value=os.environ.get('NS_FPM') or shutil.which('fpm')
    if not value or not Path(value).is_file(): raise StudyError('fpm was not found; set NS_FPM to its executable')
    return str(Path(value).resolve())

def run_one(output: Path, parameters: dict, workspace: str='solutions', timeout: float=120) -> dict:
    output=Path(output).resolve()
    if output==ROOT or output.is_relative_to(ROOT):
        raise StudyError('Choose an output directory outside the source tree')
    defaults={'task':'spectrum','shooting':'matched','nmode':2,'split_index':69,'match_index':72,
              'table_path':str(ROOT/'data/normalized/compose_69.csv'),'alpha':0.1,
              'target_mass_kg':MASS,'rtol':1e-11,'mode_rtol':1e-11,'step_cap':.05}
    defaults.update(parameters)
    table=Path(defaults['table_path']).resolve()
    if not table.is_file(): raise StudyError(f'EOS table is absent: {table}')
    defaults['table_path']=str(table)
    for k,v in defaults.items():
        if not k.replace('_','').isalnum(): raise StudyError('Invalid parameter name')
        if isinstance(v,(float,int)) and not math.isfinite(v): raise StudyError('Nonfinite parameter')
    def fmt(v): return "'"+v.replace("'","''")+"'" if isinstance(v,str) else repr(v)
    nml='&precision\n'+''.join(f' {k}={fmt(v)}\n' for k,v in sorted(defaults.items()))+'/\n'
    fc=os.environ.get('FPM_FC') or os.environ.get('FC') or 'gfortran'
    try: compiler=subprocess.check_output([fc,'--version'],text=True,timeout=10).splitlines()[0]
    except (OSError,subprocess.SubprocessError) as e: raise StudyError('Compiler identity unavailable') from e
    identity={'fpm_sha256':sha(Path(_fpm())),'compiler':compiler,
              'build_environment':{k:os.environ.get(k,'') for k in ['FPM_FC','FC','FPM_FFLAGS','FFLAGS','FPM_CFLAGS','CFLAGS']},
              'source_sha256':source_digest(workspace),'table_sha256':sha(table),
              'runner_sha256':sha(Path(__file__)),
              'input_sha256':hashlib.sha256(nml.encode()).hexdigest(),'workspace':workspace}
    receipt=output/'execution.json'
    if output.exists():
        if not receipt.is_file(): raise StudyError('Existing output has no completed execution record')
        old=json.loads(receipt.read_text())
        if old.get('status')!='complete' or any(old.get(k)!=v for k,v in identity.items()):
            raise StudyError('Existing result is partial, failed or has different inputs; choose a new directory')
        for name,digest in old['outputs'].items():
            p=output/name
            if not p.is_file() or sha(p)!=digest: raise StudyError('Existing result content does not match its record')
        return _verified_result(output)
    output.mkdir(parents=True)
    (output/'input.nml').write_text(nml)
    cmd=[_fpm(),'run','cell_solve','--profile','release','--',str(output/'input.nml'),str(output/'result')]
    rec={**identity,'status':'running','command':cmd,'parameters':defaults}
    write_json(receipt,rec)
    start=time.monotonic()
    try:
        p=subprocess.run(cmd,cwd=ROOT/workspace,capture_output=True,text=True,timeout=timeout)
        (output/'stdout.log').write_text(p.stdout);(output/'stderr.log').write_text(p.stderr)
        rec.update(returncode=p.returncode,seconds=time.monotonic()-start)
        if p.returncode!=0: raise StudyError(f'Fortran execution failed ({p.returncode}); see {output}')
        result=json.loads((output/'result.json').read_text())
        _require_finite(result)
        if result.get('status')!='ok': raise StudyError('Fortran did not close the requested calculation')
        rec['outputs']={p.name:sha(p) for p in sorted(output.glob('result*')) if p.is_file()}
        rec['status']='complete';write_json(receipt,rec)
        return result
    except (subprocess.TimeoutExpired,OSError,ValueError,StudyError) as exc:
        rec['status']='partial' if isinstance(exc,subprocess.TimeoutExpired) else 'failed'
        rec['error']=str(exc);rec['seconds']=time.monotonic()-start
        if isinstance(exc,subprocess.TimeoutExpired):
            for name,value in [('stdout.log',exc.stdout),('stderr.log',exc.stderr)]:
                (output/name).write_text(value.decode(errors='replace') if isinstance(value,bytes) else (value or ''))
        write_json(receipt,rec)
        raise StudyError(str(exc)) from exc

def _require_finite(value):
    if isinstance(value,float) and not math.isfinite(value):
        raise StudyError('Nonfinite numerical data')
    if isinstance(value,dict):
        for v in value.values(): _require_finite(v)
    elif isinstance(value,(list,tuple)):
        for v in value: _require_finite(v)

def _verified_result(path:Path)->dict:
    rec=json.loads((path/'execution.json').read_text())
    if rec.get('status')!='complete': raise StudyError('An incomplete run cannot enter a comparison')
    if rec['source_sha256']!=source_digest(rec['workspace']): raise StudyError('Source changed since this calculation')
    if rec.get('runner_sha256')!=sha(Path(__file__)): raise StudyError('Execution interface changed')
    if sha(path/'input.nml')!=rec['input_sha256']: raise StudyError('Input hash mismatch')
    if sha(Path(rec['parameters']['table_path']))!=rec['table_sha256']: raise StudyError('EOS input changed')
    for name,digest in rec['outputs'].items():
        if Path(name).name!=name or not (path/name).is_file() or sha(path/name)!=digest:
            raise StudyError('Output hash mismatch: '+name)
    r=json.loads((path/'result.json').read_text());_require_finite(r)
    if r.get('status')!='ok': raise StudyError('Non-success result')
    # Reports must correspond to the actual requested problem, not only a filename.
    for key in ['alpha','target_mass_kg','rtol','mode_rtol','step_cap','p_cut_pa','r0_m','match_index','split_index','method','shooting','task','center_order']:
        if key in rec['parameters'] and key in r and r[key]!=rec['parameters'][key]:
            raise StudyError('Result disagrees with requested '+key)
    return r

def _integral(y,x): return float(np.sum((y[1:]+y[:-1])*.5*np.diff(x)))

def mode_profile(path:Path,index:int)->tuple[np.ndarray,np.ndarray,dict]:
    r=_verified_result(path);R=r['radius_m']
    rec=json.loads((path/'execution.json').read_text())
    names=[f'result_mode{index}_inner.csv',f'result_mode{index}_outer.csv']
    for name in names:
        if sha(path/name)!=rec['outputs'][name]: raise StudyError('Mode profile hash mismatch')
    a,b=[np.genfromtxt(path/n,delimiter=',',names=True) for n in names]
    def wave(z):
        return np.column_stack((z['w_over_r3']*(z['r_m']/R)**3,z['v_over_r2']*(z['r_m']/R)**2))
    wa,wb=wave(a),wave(b)
    den=float(np.dot(wb[-1],wb[-1]))
    if not den>0: raise StudyError('Zero outer matching direction')
    gain=float(np.dot(wa[-1],wb[-1])/den)
    radius=np.r_[a['r_m'],b['r_m'][-2::-1]]
    wv=np.concatenate((wa,gain*wb[-2::-1]))
    if not np.all(np.isfinite(wv)) or not np.all(np.isfinite(radius)): raise StudyError('Nonfinite mode profile')
    dif=np.diff(radius)
    if np.any(dif<0): raise StudyError('Mode profile radius is not ordered')
    duplicate=np.where(dif==0)[0]
    if len(duplicate):
        if np.max(np.abs(wv[duplicate+1]-wv[duplicate]))>1e-10*max(1.,float(np.max(np.abs(wv)))):
            raise StudyError('A discontinuity in W,V cannot be deleted as a duplicate')
    keep=np.r_[True,dif>0];x=radius[keep]/R;wv=wv[keep]
    norm=math.sqrt(_integral(wv[:,0]**2+6*wv[:,1]**2,x))
    if norm<=0: raise StudyError('Zero mode normalization')
    wv/=norm;counts=[]
    for xmax in [.99,.999,.9999]:
        for threshold in [1e-8,1e-7,1e-6]:
            pair=[]
            for k in [0,1]:
                yy=wv[(x>=.001)&(x<=xmax),k]
                yy=yy[np.abs(yy)>threshold*np.max(np.abs(yy))]
                pair.append(int(np.sum(yy[:-1]*yy[1:]<0)))
            counts.append({'xmax':xmax,'threshold':threshold,'nodes':pair})
    return x,wv,{'node_tests':counts,'outer_scale':gain,
        'matching_vector_relative_error':float(np.linalg.norm(wa[-1]-gain*wb[-1])/np.linalg.norm(wa[-1]))}

def load_pair(reference:Path,variant:Path)->dict:
    a,b=_verified_result(reference),_verified_result(variant)
    ea=json.loads((reference/'execution.json').read_text());eb=json.loads((variant/'execution.json').read_text())
    for key in ['source_sha256','table_sha256','workspace','fpm_sha256','compiler','build_environment']:
        if ea[key]!=eb[key]: raise StudyError(f'Pair differs in {key}')
    for key in ['target_mass_kg','p_cut_pa','r0_m','shooting','method','match_index','rtol','mode_rtol','step_cap']:
        if a[key]!=b[key]: raise StudyError(f'Pair differs in control {key}')
    if a.get('center_order',0)!=b.get('center_order',0): raise StudyError('Pair differs in control center_order')
    if not a.get('core_function_unchanged') or not b.get('core_function_unchanged'): raise StudyError('Core function changed')
    if a['target_mass_kg']<=0: raise StudyError('A same-mass pair requires a positive common mass target')
    for r in [a,b]:
        if abs(r['mass_kg']/r['target_mass_kg']-1)>1e-12: raise StudyError('Mass constraint did not close')
        # A gross consistency check, not a global error bound. Refinement and
        # matching-location sensitivity still determine the numerical allowance.
        limit=max(10*max(r['rtol'],r['mode_rtol']),256*np.finfo(float).eps)
        for mode in ['f','p1']:
            mismatch=max(abs(r[mode+'_match_radius_mismatch_m'])/r['radius_m'],
                         abs(r[mode+'_match_mass_mismatch_kg'])/r['mass_kg'])
            if mismatch>limit: raise StudyError('Inner/outer background does not agree at matching point')
        if max(abs(r['f_residual']),abs(r['p1_residual']))>1e-8: raise StudyError('Matching residual did not close')
        if max(abs(r['f_surface_residual']),abs(r['p1_surface_residual']))>1e-14: raise StudyError('Surface condition failed')
    info={};identity=True
    for j,label in [(1,'f'),(2,'p1')]:
        xa,u,ca=mode_profile(reference,j);xb,v,cb=mode_profile(variant,j)
        expected=[0,0] if j==1 else [1,1]
        ok=all(c['nodes']==expected for c in ca['node_tests']+cb['node_tests'])
        x=np.linspace(.001,.999,10001)
        uu=np.column_stack([np.interp(x,xa,u[:,k]) for k in [0,1]])
        vv=np.column_stack([np.interp(x,xb,v[:,k]) for k in [0,1]])
        dot=lambda p,q:_integral(p[:,0]*q[:,0]+6*p[:,1]*q[:,1],x)
        overlap=dot(uu,vv)/math.sqrt(dot(uu,uu)*dot(vv,vv))
        identity &= ok and abs(overlap)>.99
        info[label]={'reference':ca,'variant':cb,'overlap':overlap}
    return {'delta_f_hz':b['f_hz']-a['f_hz'],'delta_p1_hz':b['p1_hz']-a['p1_hz'],
            'delta_radius_m':b['radius_m']-a['radius_m'],'identity_consistent':bool(identity),
            'core_function_unchanged':True,'reference':a,'variant':b,'mode_identity':info}

def compare_sequence(rows:list[dict])->dict:
    if len(rows)<3: return {'status':'insufficient_levels'}
    try: _require_finite(rows)
    except StudyError: return {'status':'invalid_numerical_data'}
    if not all(r['identity_consistent'] for r in rows): return {'status':'identity_unresolved'}
    if not all(r.get('core_function_unchanged',False) for r in rows): return {'status':'controls_failed'}
    if all('reference' in r and 'variant' in r for r in rows):
        controls=['alpha','target_mass_kg','p_cut_pa','r0_m','method','shooting','match_index','center_order']
        for r in rows[1:]:
            if any(r[side].get(k)!=rows[0][side].get(k) for side in ['reference','variant'] for k in controls):
                return {'status':'controls_failed'}
    ans={'status':'paired_convergent','modes':{}}
    for label in ['f','p1']:
        values=np.array([r[f'delta_{label}_hz'] for r in rows]);jumps=np.abs(np.diff(values))
        scale=max([abs(r.get(side,{}).get(label+'_hz',0.)) for r in rows for side in ['reference','variant']]+[1.])
        root=max([sum(r.get(side,{}).get(label+'_bracket_hz',0.) for side in ['reference','variant']) for r in rows])
        roundoff=64*np.finfo(float).eps*scale+root
        # Growth above the independently computed root/roundoff floor is not convergence.
        growing=jumps[-1]>jumps[-2] and jumps[-1]>roundoff
        if growing: ans['status']='convergence_unresolved'
        floor=float(max(max(jumps),roundoff))
        verdict='resolved' if abs(values[-1])>10*floor else 'below_resolution'
        if growing: verdict='convergence_unresolved'
        ans['modes'][label]={'delta_hz':float(values[-1]),'refinement_scale_hz':floor,
                             'root_roundoff_scale_hz':float(roundoff),'verdict':verdict}
    return ans

def experiment_groups()->dict:
    groups={}
    groups['main']=[(f'L{i}_a{a}',{'alpha':a,'rtol':rt,'mode_rtol':rt,'step_cap':cap})
        for i,(rt,cap) in enumerate([(1e-7,.2),(1e-9,.1),(1e-11,.05),(1e-13,.025)]) for a in [0.,.1]]
    groups['method']=[(f'rk4_a{a}',{'alpha':a,'method':'rk4','rtol':1e-12,'mode_rtol':1e-12,'step_cap':.005}) for a in [0.,.1]]
    groups['boundary']=[(f'{name}_a{a}',{'alpha':a,'rtol':1e-12,'mode_rtol':1e-12,'step_cap':.025,**change})
        for name,change in [('center',{'r0_m':.5}),('cutoff',{'p_cut_pa':1e22}),('match',{'match_index':71})] for a in [0.,.1]]
    groups['null']=[('null_a0',{'alpha':0.0,'rtol':1e-13,'mode_rtol':1e-13,'step_cap':.025})]
    groups['response']=[(name,{'alpha':a,'rtol':1e-12,'mode_rtol':1e-12,'step_cap':.025})
        for name,a in [('half',.05),('opposite',-.1)]]
    groups['mass']=[(f'mass_{sign}_a{a}',{'alpha':a,'target_mass_kg':MASS*(1+sign*1e-6),
        'rtol':1e-12,'mode_rtol':1e-12,'step_cap':.025}) for sign in [-1,1] for a in [0.,.1]]
    return groups

def supplementary_checks(out:Path,rows:list[dict])->dict:
    answer={'additional_groups':{},'physical_claim_status':'pending_additional_checks'}
    if len(rows)<4: return answer
    ref=rows[-1];checks={}
    null_path=out/'null_a0'
    null_complete=(null_path/'execution.json').is_file()
    if null_complete:
        q=load_pair(out/'L3_a0.0',null_path)
        if not q['identity_consistent'] or max(abs(q['delta_f_hz']),abs(q['delta_p1_hz']))>1e-9:
            raise StudyError('Repeated alpha=0 null did not close')
        answer['null_check']={'complete':True,'delta_f_hz':q['delta_f_hz'],'delta_p1_hz':q['delta_p1_hz']}
    else: answer['null_check']={'complete':False}

    for name,keys in [('method',['rk4']),('boundary',['center','cutoff','match']),('mass',['mass_-1','mass_1'])]:
        pairs={}
        for key in keys:
            paths=[out/f'{key}_a{a}' for a in [0.,.1]]
            if all((p/'execution.json').is_file() for p in paths):
                q=load_pair(*paths)
                if not q['identity_consistent']: raise StudyError('Supplementary mode identity failed')
                # Re-evaluate all records with the current source identity, not merely status strings.
                for p in paths:
                    e=json.loads((p/'execution.json').read_text())
                    if e['source_sha256']!=source_digest(e['workspace']): raise StudyError('Stale supplementary source')
                pairs[key]=q
        checks[name]=pairs;answer['additional_groups'][name]={'complete':len(pairs)==len(keys),'pairs':pairs}
    if not null_complete or not all(answer['additional_groups'][k]['complete'] for k in ['method','boundary','mass']): return answer
    seq=compare_sequence(rows)
    if seq['status']!='paired_convergent':
        answer['physical_claim_status']=seq['status'];return answer
    budgets={}
    for mode in ['f','p1']:
        k=mode+'_hz';dk='delta_'+k
        method=checks['method']['rk4']
        components={'paired_refinement':seq['modes'][mode]['refinement_scale_hz'],
             'step_method_absolute_disagreement':sum(abs(method[side][k]-ref[side][k]) for side in ['reference','variant'])}
        for name,q in checks['boundary'].items(): components[name]=abs(q[dk]-ref[dk])
        mass_floor=0.;derivatives={}
        for side in ['reference','variant']:
            lo=checks['mass']['mass_-1'][side];hi=checks['mass']['mass_1'][side]
            slope=(hi[k]-lo[k])/(hi['target_mass_kg']-lo['target_mass_kg']);derivatives[side]=slope
            mass_floor+=abs(slope)*abs(ref[side]['mass_kg']-ref[side]['target_mass_kg'])
        components['mass_residual']=mass_floor
        allowance=10*sum(components.values())
        budgets[mode]={'delta_hz':ref[dk],'components_hz':components,'df_dM_hz_per_kg':derivatives,
                      'empirical_allowance_hz':allowance,
                      'verdict':'resolved' if abs(ref[dk])>allowance else 'below_resolution'}
    answer['numerical_budget']=budgets
    answer['physical_claim_status']='numerically_resolved_in_fixed_model' if all(q['verdict']=='resolved' for q in budgets.values()) else 'below_resolution'
    return answer

def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--group',choices=[*experiment_groups(),'all'],default='main')
    ap.add_argument('--workspace',choices=['solutions','learner','calculation'],default='solutions')
    args=ap.parse_args(argv);out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    groups=experiment_groups()
    selected=groups.values() if args.group=='all' else [groups[args.group]]
    try:
        for group in selected:
            for name,params in group:
                r=run_one(out/name,params,args.workspace)
                print(name,r['f_hz'],r['p1_hz'],flush=True)
        rows=[]
        for i in range(4):
            paths=[out/f'L{i}_a{a}' for a in [0.,.1]]
            if all((p/'execution.json').is_file() for p in paths):rows.append(load_pair(*paths))
        report=compare_sequence(rows);report['levels']=rows
        report['limitations']='Empirical fixed-model numerical evidence, not an interval proof, EOS error bar or observation.'
        report.update(supplementary_checks(out,rows))
        write_json(out/'paired_report.json',report)
        print('paired report:',report['status'],flush=True)
        return 0
    except (StudyError,OSError,ValueError) as exc:
        print(str(exc),file=sys.stderr);return 2
if __name__=='__main__': raise SystemExit(main())
