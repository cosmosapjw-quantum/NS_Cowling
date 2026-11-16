#!/usr/bin/env python3
"""Execute explicit upper-density TF extension with unchanged production solver."""
from __future__ import annotations
import argparse, gzip, json, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'python_frontend'),str(ROOT/'scripts')]
from complete_gppva_geometry import solve_node, atomic_dump
from refine_gppva_geometry import FAMILIES, DEPENDENCIES, good, digest, family_select, transition_evidence
from compose_eos.tf_morphology import geometry_label
from build_gppva_tf import dump, export_consecutive

def read(path):
    return json.loads(gzip.decompress(Path(path).read_bytes()))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('primary',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--model',default='GPPVA_NL3wrL55')
    p.add_argument('--densities',nargs='+',type=float,default=[.0875,.09])
    p.add_argument('--anchor',type=float,default=.085)
    p.add_argument('--meshes',nargs='+',type=int,default=[192,384])
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    if any(n<=a.anchor for n in a.densities) or len(set(a.densities))!=len(a.densities):
        p.error('distinct extension densities must exceed anchor')
    a.output.mkdir(parents=True,exist_ok=False)
    out=a.output/a.model;out.mkdir()
    dependencies=list(dict.fromkeys(DEPENDENCIES+['scripts/complete_gppva_geometry.py',
        'scripts/export_gppva_geometry.py',str(Path(__file__).resolve().relative_to(ROOT))]))
    hashes={n:digest(ROOT/n) for n in dependencies}
    snapshot={n:{'sha256':h,'content':(ROOT/n).read_text()} for n,h in hashes.items()}
    dump(a.output/'SOURCE_SNAPSHOT.json.gz',snapshot)
    records=[];jobs=[];input_hash={}
    for mesh in a.meshes:
        source=a.primary/a.model/f'attempts_N{mesh:03}.json.gz'
        rows=read(source);input_hash[str(source)]=digest(source)
        anchor=next(r for r in rows if abs(r['n_b']-a.anchor)<1e-12)
        # Primary anchor is reused by identity; extension states are all new solves.
        anchor=dict(anchor,extension_anchor_provenance={'file':str(source),
            'record_index':rows.index(anchor),'sha256':input_hash[str(source)]})
        best,status,missing=family_select(anchor['attempts'])
        if best is None:raise ValueError('anchor is unresolved')
        anchor.update(status=status,missing_families=missing,selected_geometry=geometry_label(best))
        records.append(anchor)
        for n in sorted(a.densities):
            plans={}
            for family in FAMILIES:
                eligible=[(r,x) for r in rows for x in r['attempts']
                    if x.get('requested_family')==family and good(x)]
                eligible.sort(key=lambda pair:(abs(pair[0]['n_b']-n),
                    geometry_label(pair[1])!=family,pair[1]['energy_per_baryon']))
                seeds=[]
                for r,x in eligible:
                    if not seeds or abs(x['n_b']-seeds[0]['n_b'])>1e-12:
                        seeds.append(x)
                    if len(seeds)==2:break
                plans[family]=seeds if family!='uniform' else []
            uniform=next(x for x in anchor['attempts']
                if x.get('requested_family')=='uniform' and good(x))
            jobs.append((a.model,n,mesh,plans,uniform,[],str(a.output/'nodes'),hashes))
    dump(a.output/'INPUT.json',{'primary':str(a.primary),'model_id':a.model,
        'densities':sorted(a.densities),'anchor':a.anchor,'meshes':a.meshes,
        'workers':a.workers,'new_nodes':len(jobs),'reused_anchor_nodes':len(a.meshes),
        'source_sha256':hashes,'input_sha256':input_hash,
        'scope':'All six initial families; source profiles are initial values, not fitted targets.'})
    for mesh in a.meshes:
        atomic_dump(out/f'attempts_N{mesh:03}.json.gz',[r for r in records if r['points']==mesh])
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures={pool.submit(solve_node,*job):job for job in jobs}
        for future in as_completed(futures):
            r=future.result();records.append(r)
            mesh=r['points']
            atomic_dump(out/f'attempts_N{mesh:03}.json.gz',
                sorted([q for q in records if q['points']==mesh],key=lambda q:q['n_b']))
            print(json.dumps({k:r[k] for k in ['model_id','n_b','points','status','seconds']}),flush=True)
    ns=sorted([a.anchor,*a.densities]);windows=[]
    for mesh in a.meshes:
        stage=sorted([r for r in records if r['points']==mesh],key=lambda r:r['n_b']);selected=[]
        for r in stage:
            row,status,missing=family_select(r['attempts'])
            if row:selected.append(row)
        mout=out/f'N{mesh:03}';mout.mkdir()
        dump(mout/'selected_states.json.gz',selected)
        dump(mout/'RECEIPT.json',{'export':export_consecutive(a.model,selected,ns,mout)})
        for left,right in zip(stage,stage[1:]):
            lg,rg=left['selected_geometry'],right['selected_geometry']
            if lg and rg and lg!=rg:
                w={'left':left['n_b'],'right':right['n_b'],'geometries':[lg,rg]}
                if w not in windows:windows.append(w)
    summary={'model_id':a.model,'windows':windows,'meshes':a.meshes,'subdivisions':None,
        'requested_density_nodes':len(ns),'new_mesh_density_nodes':len(jobs),
        'reused_anchor_nodes':len(a.meshes),'fine_selected_nodes':sum(
            r['selected_geometry'] is not None for r in records if r['points']==a.meshes[-1]),
        'fine_held_nodes':[{'n_b':r['n_b'],'status':r['status']} for r in records
            if r['points']==a.meshes[-1] and r['selected_geometry'] is None],
        'transitions':transition_evidence(windows,records,a.meshes),
        'seconds':time.monotonic()-started,'scope':'Four new extension nodes plus two unchanged primary anchors.'}
    dump(out/'SUMMARY.json',summary)
    unchanged=all(digest(ROOT/n)==h for n,h in hashes.items())
    inputs_unchanged=all(digest(Path(n))==h for n,h in input_hash.items())
    dump(a.output/'CAMPAIGN.json',{'models':[summary],'source_sha256':hashes,
        'source_unchanged_during_execution':unchanged,'inputs_unchanged':inputs_unchanged,
        'primary_new_node_count':len(jobs)})
    if not unchanged or not inputs_unchanged:raise RuntimeError('source or input changed')
    if any(r['selected_geometry'] is None for r in records):raise SystemExit(2)
if __name__=='__main__':main()
