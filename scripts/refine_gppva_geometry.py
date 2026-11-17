#!/usr/bin/env python3
"""Recompute spatial TF geometry competitors on density and radial meshes.

The search compares the six declared Wigner--Seitz seed families. It records
their final physical geometry, including collapse into another branch. A family
is covered by a stationary recomputation, not by the success of every arbitrary
initial guess. Missing families and potentially lower unresolved states remain
visible. No CompOSE target thermodynamics enters this calculation.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from compose_eos.spatial_tf_equilibrium import solve_screened_cell
from compose_eos.tf_selection import select_geometry_candidates,feasibility
from compose_eos.tf_search import weak_cluster_seed
from build_gppva_tf import MODELS,dump,export_consecutive,physical_geometry
from refine_gppva_tf import good
from tf_polish import profile_initial,polish_profile

FAMILIES={'droplet':(3,'droplet'),'rod':(2,'droplet'),'slab':(1,'droplet'),
          'tube':(2,'hole'),'bubble':(3,'hole'),'uniform':(3,'uniform')}
SOURCE=ROOT/'data/research_completion_20261116/gppva_tf/final'
DEPENDENCIES=['scripts/refine_gppva_geometry.py','scripts/build_gppva_tf.py',
    'scripts/tf_polish.py','scripts/refine_gppva_tf.py',
    'python_frontend/compose_eos/tf_selection.py','python_frontend/compose_eos/spatial_tf.py',
    'python_frontend/compose_eos/spatial_tf_equilibrium.py',
    'python_frontend/compose_eos/gppva_rmf.py','python_frontend/compose_eos/cold_leptons.py',
    'python_frontend/compose_eos/tf_composition.py','python_frontend/compose_eos/tf_search.py',
    'python_frontend/compose_eos/tf_morphology.py']

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def family_select(attempts,expected=tuple(FAMILIES),energy_tolerance=1e-6):
    """Select a finite-family candidate, preserving missing/lower failures."""
    decision=select_geometry_candidates(attempts,[FAMILIES[g] for g in expected],
                                        energy_tolerance=energy_tolerance)
    reverse={v:k for k,v in FAMILIES.items()}
    missing=[reverse[tuple(f)] for f in decision['unresolved_families']]
    return decision['selected'],decision['status'],missing

def branch_minima(attempts):
    out={}
    for r in attempts:
        if not good(r) or feasibility(r,r['n_b'])!='FEASIBLE':continue
        g=physical_geometry(r)
        if g not in out or r['epsilon']<out[g]['epsilon']:out[g]=r
    return out

def geometry_windows(rows,*,recorded=False):
    """Return observed changes in the original finite-density profile sequence."""
    ordered=sorted((r for r in rows if r['n_b']>=.002),key=lambda r:r['n_b'])
    geometry=lambda row:row['geometry'] if recorded and 'geometry' in row else physical_geometry(row)
    return [{'left':a['n_b'],'right':b['n_b'],
             'geometries':[geometry(a),geometry(b)]}
            for a,b in zip(ordered,ordered[1:]) if geometry(a)!=geometry(b)]

def choose_profile_seeds(candidates,family,continuations=()):
    """Retain family continuation even when it relaxed to another morphology."""
    if family=='uniform':return []
    same=[r for r in candidates if physical_geometry(r)==family]
    selected=(same+[r for r in candidates if physical_geometry(r)!=family])[:2]
    for old in continuations:
        if old is not None and not any(old['n_b']==r['n_b'] and old['points']==r['points'] and old['radius_fm']==r['radius_fm'] for r in selected):selected.append(old)
    return selected

def all_profiles(value):
    if isinstance(value,dict):
        if 'profile' in value and 'model_id' in value:yield value
        else:
            for item in value.values():yield from all_profiles(item)
    elif isinstance(value,list):
        for item in value:yield from all_profiles(item)

def scale_initial(old,n,points):
    initial=profile_initial(old,points)
    for key in ('nn','np','ne'):initial[key]*=n/old['n_b']
    return initial

def solve_attempt(model,n,points,family,old=None,label='COLD'):
    dim,seed=FAMILIES[family]
    metadata={'requested_family':family,'initial_kind':label,
              'initial_density':old['n_b'] if old else None}
    try:
        row=solve_screened_cell(model,n,points=points,dimension=dim,seed=seed,
            radius=30.,radius_bounds=(5.,180.),maxiter=900,ftol=1e-13,
            initial=scale_initial(old,n,points) if old else None)
        row.update(metadata,geometry=physical_geometry(row))
        attempts=[row]
        if not good(row):
            try:
                fixed=polish_profile(row,points)
                for retry in fixed.get('restart_attempts',[]):
                    retry.update(metadata,geometry=physical_geometry(retry),initial_kind=label+'_POLISH_RESTART')
                    attempts.append(retry)
                fixed.update(metadata,geometry=physical_geometry(fixed),initial_kind=label+'_POLISHED')
                attempts.append(fixed)
            except Exception as exc:
                attempts.append(dict(metadata,model_id=model,n_b=n,points=points,
                    dimension=dim,seed=seed,execution='FAILED',exception=type(exc).__name__,message=str(exc)))
        return attempts
    except Exception as exc:
        return [dict(metadata,model_id=model,n_b=n,points=points,dimension=dim,seed=seed,
                     execution='FAILED',exception=type(exc).__name__,message=str(exc))]

def transition_evidence(windows,records,meshes):
    lookup={(r['n_b'],r['points']):r for r in records};out=[]
    for window in windows:
        g0,g1=window['geometries'];samples=[]
        densities=sorted({r['n_b'] for r in records if window['left']<=r['n_b']<=window['right']})
        for n in densities:
            values={}
            for mesh in meshes:
                row=lookup.get((n,mesh));branches=branch_minima(row['attempts']) if row else {}
                if g0 in branches and g1 in branches:
                    values[str(mesh)]={'delta_E_MeV_per_baryon':branches[g1]['energy_per_baryon']-branches[g0]['energy_per_baryon'],
                        'geometries':[g0,g1],
                        'energies_per_baryon':[branches[g0]['energy_per_baryon'],branches[g1]['energy_per_baryon']],
                        'pressures':[branches[g0]['pressure'],branches[g1]['pressure']],
                        'mu_B':[branches[g0]['mu_B'],branches[g1]['mu_B']]}
            fine=values.get(str(meshes[-1]));coarse=values.get(str(meshes[-2])) if len(meshes)>1 else None
            shift=abs(fine['delta_E_MeV_per_baryon']-coarse['delta_E_MeV_per_baryon']) if fine and coarse else None
            roundoff=(64*np.finfo(float).eps*max(1.,*(abs(x) for x in fine['energies_per_baryon']))) if fine else None
            samples.append({'n_b':n,'mesh_values':values,'adjacent_mesh_deltaE_change_MeV_per_baryon':shift,
                            'fine_sign_exceeds_mesh_change':bool(fine and shift is not None and abs(fine['delta_E_MeV_per_baryon'])>shift),
                            'roundoff_scale_MeV_per_baryon':roundoff,
                            'roundoff_near_degeneracy':bool(fine and abs(fine['delta_E_MeV_per_baryon'])<=roundoff)})
        brackets=[]
        for a,b in zip(samples,samples[1:]):
            af=a['mesh_values'].get(str(meshes[-1]));bf=b['mesh_values'].get(str(meshes[-1]))
            if af and bf and af['delta_E_MeV_per_baryon']*bf['delta_E_MeV_per_baryon']<0:
                brackets.append({'n_left':a['n_b'],'n_right':b['n_b'],
                    'both_signs_exceed_mesh_change':a['fine_sign_exceeds_mesh_change'] and b['fine_sign_exceeds_mesh_change'],
                    'both_endpoints_exceed_roundoff_scale':not a['roundoff_near_degeneracy'] and not b['roundoff_near_degeneracy']})
        out.append({**window,'samples':samples,'fine_mesh_sign_change_brackets':brackets,
            'interpretation':'Finite radial-grid energy ordering; a sign bracket is not automatically a continuum transition.'})
    return out

def run_model(model,output,meshes,subdivisions,only_window=None,warm_start=None):
    start=time.monotonic();out=Path(output)/model;out.mkdir()
    source=SOURCE/model/'selected_states.json.gz'
    selected_old=json.loads(gzip.decompress(source.read_bytes()))
    windows=geometry_windows(selected_old,recorded=True)
    if only_window:windows=[w for w in windows if w['left']>=only_window[0] and w['right']<=only_window[1]]
    if not windows:raise ValueError('no observed geometry windows')
    candidates=list(selected_old)
    attempt_source=SOURCE/model/'attempts.json.gz'
    candidates.extend(all_profiles(json.loads(gzip.decompress(attempt_source.read_bytes()))))
    warm_paths=sorted((Path(warm_start).resolve()/model).glob('attempts_N*.json.gz')) if warm_start else []
    for path in warm_paths:candidates.extend(all_profiles(json.loads(gzip.decompress(path.read_bytes()))))
    templates={}
    for g in FAMILIES:
        templates[g]=[r for r in candidates if good(r) and physical_geometry(r)==g]
    ns=sorted({float(n) for w in windows for n in np.linspace(w['left'],w['right'],subdivisions+1)})
    records=[];previous_mesh={};current_guess={}
    for points in meshes:
        rows=[]
        for n in ns:
            attempts=[]
            for family in FAMILIES:
                warm=[]
                key=(n,family)
                if key in previous_mesh:warm.append(previous_mesh[key])
                if family in current_guess:warm.append(current_guess[family])
                warm.extend(sorted(templates[family],key=lambda r:(abs(r['n_b']-n),abs(r['points']-points)))[:2])
                unique=[]
                for old in warm:
                    if not any(old['n_b']==q['n_b'] and old['points']==q['points'] and old['radius_fm']==q['radius_fm'] for q in unique):unique.append(old)
                # Two neighboring independent profile continuations for a clustered
                # branch. The uniform seed is evaluated directly from its equations.
                seeds=choose_profile_seeds(unique,family,[previous_mesh.get(key),current_guess.get(family)])
                family_attempts=[]
                for k,old in enumerate(seeds):
                    family_attempts.extend(solve_attempt(model,n,points,family,old,'PROFILE_'+str(k)))
                if not seeds or not any(good(r) for r in family_attempts):
                    family_attempts.extend(solve_attempt(model,n,points,family,None,'COLD_RECOVERY' if seeds else 'COLD'))
                if family!='uniform' and not any(good(r) for r in family_attempts):
                    uniform_trials=solve_attempt(model,n,points,'uniform',None,'UNIFORM_FOR_WEAK_RECOVERY')
                    attempts.extend(uniform_trials)
                    uniform_rows=[r for r in uniform_trials if good(r)]
                    if uniform_rows:
                        uniform=min(uniform_rows,key=lambda r:r['epsilon']);sign=-1 if FAMILIES[family][1]=='hole' else 1
                        for amplitude in [.03,.3]:
                            for radius in [15.,30.,60.]:
                                initial=weak_cluster_seed(uniform,FAMILIES[family][0],sign*amplitude,radius)
                                family_attempts.extend(solve_attempt(model,n,points,family,initial,f'WEAK_RECOVERY_{sign*amplitude}_{radius}'))
                                if any(good(r) for r in family_attempts):break
                            if any(good(r) for r in family_attempts):break
                attempts.extend(family_attempts)
                accepted=[r for r in family_attempts if good(r)]
                if accepted:
                    chosen=min(accepted,key=lambda r:r['epsilon'])
                    previous_mesh[key]=chosen;current_guess[family]=chosen
            best,status,missing=family_select(attempts)
            record={'n_b':n,'points':points,'status':status,'missing_families':missing,
                    'selected_geometry':physical_geometry(best) if best else None,'attempts':attempts}
            records.append(record)
            if best:rows.append(best)
            dump(out/'progress.json',{'points':points,'completed':len([r for r in records if r['points']==points]),
                                      'total':len(ns),'last_density':n,'last_status':status})
            dump(out/f'attempts_N{points:03}.json.gz',[r for r in records if r['points']==points])
        meshout=out/f'N{points:03}';meshout.mkdir()
        dump(meshout/'selected_states.json.gz',rows)
        dump(meshout/'RECEIPT.json',{'export':export_consecutive(model,rows,ns,meshout),
            'status':'FINITE_FAMILY_SEARCH','held_densities':[r['n_b'] for r in records if r['points']==points and r['selected_geometry'] is None]})
    transitions=transition_evidence(windows,records,meshes)
    final=[r for r in records if r['points']==meshes[-1]]
    summary={'model_id':model,'windows':windows,'meshes':meshes,'subdivisions':subdivisions,
        'requested_density_nodes':len(ns),'fine_selected_nodes':sum(r['selected_geometry'] is not None for r in final),
        'fine_held_nodes':[{'n_b':r['n_b'],'status':r['status'],'missing_families':r['missing_families']} for r in final if r['selected_geometry'] is None],
        'transitions':transitions,'seed_input_sha256':{str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):digest(p) for p in [source,attempt_source,*warm_paths]},
        'seconds':time.monotonic()-start,'global_minimum':'Not established; finite families and starts only.',
        'source_table_matching':'No target EOS values used or exact-match gate.'}
    dump(out/'SUMMARY.json',summary)
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    p.add_argument('--models',nargs='+',choices=MODELS,default=list(MODELS))
    p.add_argument('--meshes',nargs='+',type=int,default=[192,384])
    p.add_argument('--subdivisions',type=int,default=2);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--window',nargs=2,type=float)
    p.add_argument('--warm-start',type=Path,help='Prior campaign profiles are used only as initial guesses; every final cell is recomputed.')
    a=p.parse_args()
    if len(a.meshes)<2 or any(np.diff(a.meshes)<=0) or min(a.meshes)<8 or a.subdivisions<2:p.error('two increasing meshes and at least two density subdivisions required')
    a.output.mkdir(parents=True,exist_ok=False)
    hashes={name:digest(ROOT/name) for name in DEPENDENCIES}
    dump(a.output/'INPUT.json',{'head':__import__('subprocess').check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'models':a.models,'meshes':a.meshes,'subdivisions':a.subdivisions,'window':a.window,'warm_start':str(a.warm_start) if a.warm_start else None,'source_sha256':hashes})
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        jobs={pool.submit(run_model,m,a.output,a.meshes,a.subdivisions,a.window,a.warm_start):m for m in a.models}
        for job in as_completed(jobs):
            try:r=job.result();results.append(r);print(json.dumps({k:r[k] for k in ('model_id','requested_density_nodes','fine_selected_nodes','seconds')}),flush=True)
            except Exception as exc:
                r={'model_id':jobs[job],'status':'EXECUTION_FAILED','exception':type(exc).__name__,'message':str(exc)};results.append(r);print(json.dumps(r),flush=True)
    unchanged=all(digest(ROOT/name)==value for name,value in hashes.items())
    dump(a.output/'CAMPAIGN.json',{'models':results,'source_sha256':hashes,'source_unchanged_during_execution':unchanged})
    if not unchanged:raise RuntimeError('scientific source changed during campaign')
    if any(r.get('status')=='EXECUTION_FAILED' or r.get('fine_held_nodes') for r in results):
        raise SystemExit(2)
if __name__=='__main__':main()
