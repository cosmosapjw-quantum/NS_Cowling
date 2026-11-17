#!/usr/bin/env python3
"""Execute native TF driver finite starting-radius probes; no Python physics."""
import argparse,hashlib,json,os,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DRIVER=ROOT/'build/native_tf_test/native_tf_driver'
POINTS={'FSU2R':(.048,.048999999),'FSU2H':(.050999999,.052000001),'TM1e':(.052000001,.052999999)}
FAMILIES=((2,'droplet'),(2,'hole'),(3,'hole'))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def task(out,model,n,dim,seed,radius,timeout):
    name=f'{model}_n{n:.9f}_d{dim}_{seed}_r{radius:02d}'
    cmd=[str(DRIVER),'solve','GPPVA_'+model,'64',str(dim),format(n,'.17g'),str(radius),seed,'1800','1e-13']
    meta={'id':name,'command':cmd,'model':model,'density_fm3':n,'mesh':64,'dimension':dim,'seed':seed,'initial_radius_fm':radius,'status':'RUNNING'}
    dest=out/(name+'.json');dest.write_text(json.dumps(meta,indent=2)+'\n')
    started=time.monotonic()
    try:
        r=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=timeout,check=False)
        (out/(name+'.stdout')).write_bytes(r.stdout);(out/(name+'.stderr')).write_bytes(r.stderr)
        meta.update(status='COMPLETED' if r.returncode==0 else 'PROCESS_FAILURE',returncode=r.returncode,
                    stdout_sha256=sha(out/(name+'.stdout')),stderr_sha256=sha(out/(name+'.stderr')))
    except subprocess.TimeoutExpired as e:
        (out/(name+'.stdout')).write_bytes(e.stdout or b'');(out/(name+'.stderr')).write_bytes(e.stderr or b'')
        meta.update(status='TIMEOUT',stdout_sha256=sha(out/(name+'.stdout')),stderr_sha256=sha(out/(name+'.stderr')))
    meta['seconds']=time.monotonic()-started;tmp=dest.with_suffix('.tmp');tmp.write_text(json.dumps(meta,indent=2)+'\n');tmp.replace(dest)
    return {'id':name,'status':meta['status'],'seconds':meta['seconds']}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('--workers',type=int,default=3);p.add_argument('--timeout',type=int,default=300);a=p.parse_args()
    if not 1<=a.workers<=4: p.error('workers must be between 1 and 4')
    a.output.mkdir(parents=True,exist_ok=False)
    inputs=[DRIVER,ROOT/'fortran_eos/tests/native_tf_driver.f90',ROOT/'fortran_eos/src/compose_eos_spatial_tf.f90',ROOT/'fortran_eos/src/compose_eos_tf_models.f90',ROOT/'fortran_eos/src/compose_eos_optimize_native.f90',ROOT/'fortran_eos/vendor/slsqp/slsqp_optmz.f',Path(__file__)]
    identities={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    summary={'status':'RUNNING','producer':'native Fortran driver','source_sha256':identities,'mesh':64,'ftol':1e-13,'maxiter':1800,'scope':'three explicit seed families and three initial radii; not global optimization','results':[]}
    (a.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        jobs=[pool.submit(task,a.output,m,n,d,s,r,a.timeout) for m,ns in POINTS.items() for n in ns for d,s in FAMILIES for r in (30,45,70)]
        for f in as_completed(jobs):
            result=f.result();summary['results'].append(result);print(json.dumps(result),flush=True)
            (a.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    summary['source_unchanged']=all(sha(ROOT/p)==h for p,h in identities.items())
    summary['status']='COMPLETED' if all(r['status']=='COMPLETED' for r in summary['results']) and summary['source_unchanged'] else 'COMPLETED_WITH_FAILURES'
    (a.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
