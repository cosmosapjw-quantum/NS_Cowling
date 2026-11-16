#!/usr/bin/env python3
"""Complete missing fine-mesh density nodes in parallel with frozen TF physics."""
from __future__ import annotations
import argparse,gzip,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from refine_gppva_geometry import (FAMILIES,DEPENDENCIES,SOURCE,digest,good,solve_attempt,
    family_select,geometry_windows,transition_evidence)
from compose_eos.tf_search import weak_cluster_seed
from compose_eos.tf_selection import select_geometry_candidates
from build_gppva_tf import MODELS,dump,export_consecutive,physical_geometry
from export_gppva_geometry import additional_attempts

def read(path):return json.loads(gzip.decompress(Path(path).read_bytes()))

def atomic_dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.stem+'.tmp'+path.suffix)
    dump(temporary,value);temporary.replace(path)

def extend_supplement(attempts,extras):
    seen={(r.get('supplemental_attempt_source'),r.get('supplemental_attempt_index')) for r in attempts if 'supplemental_attempt_source' in r}
    for row in extras:
        key=(row.get('supplemental_attempt_source'),row.get('supplemental_attempt_index'))
        if key not in seen:attempts.append(row);seen.add(key)

def solve_node(model,n,points,plans,uniform,extras,checkpoint_dir,source_hashes):
    start=time.monotonic();attempts=list(extras);completed=[]
    checkpoint_path=Path(checkpoint_dir)/model/f'n{n:.12f}_N{points}.json.gz'
    def checkpoint(rows,current_family):
        atomic_dump(checkpoint_path,{'model_id':model,'n_b':n,'points':points,'status':'RUNNING',
            'completed_families':list(completed),'current_family':current_family,'attempts':rows,
            'source_sha256':source_hashes,'seconds':time.monotonic()-start})
    checkpoint(attempts,None)
    for family in FAMILIES:
        existing=[r for r in attempts if r.get('requested_family')==family]
        if any(good(r) for r in existing):
            completed.append(family);checkpoint(attempts,None);continue
        trials=[]
        for index,old in enumerate(plans[family]):
            trials.extend(solve_attempt(model,n,points,family,old,f'PARALLEL_COARSE_PROFILE_{index}'))
            checkpoint(attempts+trials,family)
        if family=='uniform' or not any(good(r) for r in trials):
            trials.extend(solve_attempt(model,n,points,family,None,'PARALLEL_COLD'))
            checkpoint(attempts+trials,family)
        if family!='uniform' and not any(good(r) for r in trials):
            sign=-1 if FAMILIES[family][1]=='hole' else 1
            for amplitude in [.03,.3]:
                for radius in [15.,30.,60.]:
                    seed=weak_cluster_seed(uniform,FAMILIES[family][0],sign*amplitude,radius)
                    trials.extend(solve_attempt(model,n,points,family,seed,f'PARALLEL_WEAK_{sign*amplitude}_{radius}'))
                    checkpoint(attempts+trials,family)
                    if any(good(r) for r in trials):break
                if any(good(r) for r in trials):break
        attempts.extend(trials);completed.append(family);checkpoint(attempts,None)
    best,status,missing=family_select(attempts)
    result={'model_id':model,'n_b':n,'points':points,'status':status,'missing_families':missing,
        'selected_geometry':physical_geometry(best) if best else None,'attempts':attempts,
        'execution_provenance':'INDEPENDENT_DENSITY_JOB_FROM_EXPLICIT_COARSE_PROFILES',
        'seconds':time.monotonic()-start,'source_sha256':source_hashes}
    atomic_dump(checkpoint_path,result);return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('baseline',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--additional',nargs='*',type=Path,default=[]);p.add_argument('--workers',type=int,default=12)
    p.add_argument('--equivalence-receipt',type=Path,help='Reviewed numerical-equivalence evidence when the field solver source changes; bulk functional must remain identical.')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);extras,input_hash=additional_attempts(a.additional)
    hashes={n:digest(ROOT/n) for n in DEPENDENCIES+['scripts/complete_gppva_geometry.py','scripts/export_gppva_geometry.py','python_frontend/compose_eos/tf_search.py']}
    previous_sources=json.loads((a.baseline/'INPUT.json').read_text())['source_sha256']
    bulk='python_frontend/compose_eos/gppva_rmf.py'
    if previous_sources.get(bulk)!=hashes[bulk]:raise ValueError('bulk functional/input identity changed; retained nodes must be recomputed')
    solver_names=['python_frontend/compose_eos/spatial_tf.py','python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/cold_leptons.py']
    changed_solvers={name:{'before':previous_sources.get(name),'after':hashes[name]} for name in solver_names if previous_sources.get(name)!=hashes[name]}
    if changed_solvers and not a.equivalence_receipt:raise ValueError('changed field solver requires explicit reviewed numerical-equivalence receipt')
    equivalence={'changed_solver_sources':changed_solvers,'receipt':str(a.equivalence_receipt) if a.equivalence_receipt else None,
        'receipt_sha256':digest(a.equivalence_receipt) if a.equivalence_receipt else None,
        'scope':'Caller-supplied reviewed evidence; byte comparison alone does not prove numerical equivalence.'}
    dump(a.output/'DEPENDENCY_COMPARISON.json',equivalence)
    dump(a.output/'INPUT.json',{'baseline':str(a.baseline),'additional':list(map(str,a.additional)),
        'workers':a.workers,'source_sha256':hashes,'scope':'Reuse completed exact-solver nodes and independently recompute every missing fine-density node.'})
    snapshot={name:{'sha256':sha,'content':(ROOT/name).read_text()} for name,sha in hashes.items()}
    dump(a.output/'SOURCE_SNAPSHOT.json.gz',snapshot)
    bymodel={};jobs=[];windows={}
    for model in MODELS:
        out=a.output/model;out.mkdir();coarse_path=a.baseline/model/'attempts_N192.json.gz'
        coarse=read(coarse_path);input_hash[str(coarse_path)]=digest(coarse_path)
        for r in coarse:
            r['execution_provenance']={'baseline_file':str(coarse_path),'sha256':input_hash[str(coarse_path)]}
            r['status_before_supplement']=r['status']
            extend_supplement(r['attempts'],extras.get((model,round(r['n_b'],12),192),[]))
            best,r['status'],r['missing_families']=family_select(r['attempts']);r['selected_geometry']=physical_geometry(best) if best else None
        fine_path=a.baseline/model/'attempts_N384.json.gz';fine=read(fine_path) if fine_path.exists() else []
        if fine_path.exists():input_hash[str(fine_path)]=digest(fine_path)
        for r in fine:
            r['execution_provenance']={'baseline_file':str(fine_path),'sha256':input_hash[str(fine_path)]}
            r['status_before_supplement']=r['status']
            extend_supplement(r['attempts'],extras.get((model,round(r['n_b'],12),384),[]))
            best,r['status'],r['missing_families']=family_select(r['attempts']);r['selected_geometry']=physical_geometry(best) if best else None
        bymodel[model]=coarse+fine
        old_path=SOURCE/model/'selected_states.json.gz';windows[model]=geometry_windows(read(old_path),recorded=True);input_hash[str(old_path)]=digest(old_path)
        for record in coarse:
            n=record['n_b']
            if any(abs(r['n_b']-n)<1e-12 for r in fine):continue
            plans={}
            for family in FAMILIES:
                candidates=[r for r in record['attempts'] if r.get('requested_family')==family and 'profile' in r]
                candidates.sort(key=lambda r:(not good(r),r['energy_per_baryon']))
                chosen=candidates[:1]
                neighbors=[r for q in coarse if q is not record for r in q['attempts'] if r.get('requested_family')==family and good(r)]
                neighbors.sort(key=lambda r:(abs(r['n_b']-n),physical_geometry(r)!=family,r['energy_per_baryon']))
                if neighbors:chosen.append(neighbors[0])
                elif len(candidates)>1:chosen.append(candidates[1])
                plans[family]=chosen if family!='uniform' else []
            uniform=next(r for r in record['attempts'] if r.get('requested_family')=='uniform' and 'profile' in r)
            jobs.append((model,n,384,plans,uniform,extras.get((model,round(n,12),384),[]),str(a.output/'nodes'),hashes))
        for mesh in [192,384]:dump(out/f'attempts_N{mesh}.json.gz',sorted([r for r in bymodel[model] if r['points']==mesh],key=lambda r:r['n_b']))
    dump(a.output/'PROVENANCE.json',{'input_sha256':input_hash,'reused_nodes':sum(len(r) for r in bymodel.values()),'redispatched_nodes':len(jobs)})
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures={pool.submit(solve_node,*job):job for job in jobs}
        for future in as_completed(futures):
            job=futures[future];model=job[0]
            try:r=future.result()
            except Exception as exc:r={'n_b':job[1],'points':job[2],'status':'EXECUTION_FAILED','selected_geometry':None,'missing_families':list(FAMILIES),'attempts':[],'message':str(exc)}
            bymodel[model].append(r)
            atomic_dump(a.output/model/'attempts_N384.json.gz',sorted([x for x in bymodel[model] if x['points']==384],key=lambda x:x['n_b']))
            print(json.dumps({'model_id':model,'n_b':r['n_b'],'status':r['status'],'seconds':r.get('seconds')}),flush=True)
    summaries=[]
    for model,records in bymodel.items():
        coarse=sorted([r for r in records if r['points']==192],key=lambda r:r['n_b']);ns=[r['n_b'] for r in coarse]
        for mesh in [192,384]:
            stage=sorted([r for r in records if r['points']==mesh],key=lambda r:r['n_b']);rows=[]
            for record in stage:
                if not record['attempts']:continue
                row=select_geometry_candidates(record['attempts'],FAMILIES.values(),n_b=record['n_b'])['selected']
                if row:rows.append(row)
            out=a.output/model/f'N{mesh}';out.mkdir();dump(out/'selected_states.json.gz',rows)
            dump(out/'RECEIPT.json',{'export':export_consecutive(model,rows,ns,out)})
        fine=[r for r in records if r['points']==384]
        summary={'model_id':model,'meshes':[192,384],'windows':windows[model],'requested_density_nodes':len(ns),
            'fine_selected_nodes':sum(r['selected_geometry'] is not None for r in fine),
            'fine_held_nodes':[{'n_b':r['n_b'],'status':r['status']} for r in fine if r['selected_geometry'] is None],
            'transitions':transition_evidence(windows[model],records,[192,384])}
        dump(a.output/model/'SUMMARY.json',summary);summaries.append(summary)
    unchanged=all(digest(ROOT/n)==h for n,h in hashes.items())
    dump(a.output/'CAMPAIGN.json',{'models':summaries,'source_sha256':hashes,'source_unchanged_during_execution':unchanged,
        'baseline':'Completed baseline nodes retain original source identity; supplementary attempts retain source file SHA.'})
    if not unchanged:raise RuntimeError('source changed during parallel completion')
    if any(r['fine_held_nodes'] for r in summaries):raise SystemExit(2)
if __name__=='__main__':main()
