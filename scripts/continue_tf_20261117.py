#!/usr/bin/env python3
"""Bounded native DD2 continuation; Python only builds seeds and validates results."""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('NS_TF_CONTINUATION_OUTPUT',str(ROOT/'data/tf_pair_expansion_20261117/tf_continuation')))
EXE=ROOT/'build/tf_continuation_20261117/tf_continuation_20261117'
KEYS='epsilon pressure mu_B mu_e scalar_residual vector_residual electron_residual d_epsilon_d_log_radius kkt radius_stationarity Y_p radius iterations evaluations optimizer_mode baryon_error charge_error optimizer_success radius_bound density_bound muon_threshold_crossed'.split()
FAMILIES=[(3,'droplet'),(2,'droplet'),(1,'droplet'),(2,'hole'),(3,'hole'),(3,'uniform')]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def weights(n,d):return np.diff(np.linspace(0,1,n+1)**d)
def refine_label(label,n):return f"N{n:03d}_"+label.split("_",1)[1]
def initial(old,n,target,yp=None,rscale=1):
    p=np.asarray(old['profile']);x=(np.arange(n)+.5)/n;xo=(np.arange(len(p))+.5)/len(p)
    nn=np.interp(x,xo,p[:,0]);pr=np.interp(x,xo,p[:,1]);w=weights(n,old['dimension'])
    if yp is None:
        fac=target/(w@(nn+pr));nn*=fac;pr*=fac
    else:
        nn*=target*(1-yp)/(w@nn);pr*=target*yp/(w@pr)
    return nn,pr,old['radius']*rscale

def run(old,n,target,label,yp=None,rscale=1):
    nn,pr,r=initial(old,n,target,yp,rscale);d=old['dimension'];seed=old['seed']
    dest=OUT/label;dest.mkdir(exist_ok=False)
    inp=f'DD2 {n} {d} {seed} {target:.17g} {r:.17g}\n'+'\n'.join(' '.join(format(z,'.17g') for z in a) for a in [nn,pr])+'\n'
    (dest/'input.txt').write_text(inp)
    rec={'id':label,'n_b':target,'points':n,'dimension':d,'seed':seed,'parent_id':old['id'],'initial_mean_Yp':float(weights(n,d)@pr/target),'initial_radius':r,'public_Yp_initial_only':yp,'status':'RUNNING','command':[str(EXE)],'stdin_sha256':sha(dest/'input.txt')}
    dump(dest/'STATE.json',rec);t=time.monotonic()
    try:
        p=subprocess.run([str(EXE)],input=inp,text=True,capture_output=True,timeout=1200,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
        (dest/'stdout.txt').write_text(p.stdout);(dest/'stderr.txt').write_text(p.stderr);rec.update(returncode=p.returncode,seconds=time.monotonic()-t)
        lines=p.stdout.splitlines()
        if p.returncode==0 and lines[0]=='ok':
            st,geom=lines[1].split();vals=list(map(float,lines[2].split()));prof=np.loadtxt(lines[3:]).tolist()
            rec.update(dict(zip(KEYS,vals)));rec.update(status=st,geometry=geom,profile=prof)
            rec['accepted_stationary']=st in ('STATIONARY','UNIFORM_STATIONARY') and rec['optimizer_success']==1 and rec['muon_threshold_crossed']==0 and rec['baryon_error']<1e-10 and rec['charge_error']<1e-10
        else:rec.update(status='PROCESS_OR_FIELD_FAILURE',accepted_stationary=False)
    except subprocess.TimeoutExpired as e:
        (dest/'stdout.txt').write_bytes(e.stdout or b'');(dest/'stderr.txt').write_bytes(e.stderr or b'');rec.update(status='TIMEOUT',accepted_stationary=False,seconds=time.monotonic()-t)
    dump(dest/'STATE.json',rec);print(json.dumps({k:rec.get(k) for k in ['id','status','geometry','epsilon','pressure','kkt','radius','seconds']}),flush=True)
    return rec

def starts():
    ans=[]
    for d,s in FAMILIES:
        p=ROOT/f'data/research_native_20261116/tf/family_campaign/DD2_0.0030_d{d}_{s}.out';lines=p.read_text().splitlines();v=list(map(float,lines[1].split()))
        assert lines[0]=='ok' and v[15]==1
        ans.append({'id':str(p.relative_to(ROOT)),'source_sha256':sha(p),'n_b':.003,'points':16,'dimension':d,'seed':s,'radius':v[11],'epsilon':v[0],'profile':np.loadtxt(lines[2:]).tolist()})
    return ans

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['coarse','refine','validate'],required=True);a=ap.parse_args()
    if a.phase=='coarse':
        OUT.mkdir(parents=True,exist_ok=True)
        old=starts();files=[Path(__file__),EXE,ROOT/'fortran_eos/tests/tf_continuation_20261117.f90',*[ROOT/f'fortran_eos/src/{s}.f90' for s in ['compose_eos_spatial_tf','compose_eos_tf_models','compose_eos_optimize_native']],ROOT/'fortran_eos/vendor/slsqp/slsqp_optmz.f',ROOT/'data/matched_eos_reference_20261117/reference_support/SEED_CATALOG.json']
        dump(OUT/'INPUT.json',{'base_commit':'894cbf1c8315f07ee926c29b83766ffdab2ac7b7','source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'starting_states':[{k:v for k,v in s.items() if k!='profile'} for s in old],'units':{'density':'fm^-3','radius':'fm','energy':'MeV fm^-3','hbar_c':'197.3269804 MeV fm'},'closure':'T=0 beta-equilibrium npe; screened electrons, global neutrality; Yp optimized; no muons allowed beyond threshold'})
        rows=[]
        for s in old:
            f=f'd{s["dimension"]}_{s["seed"]}'
            for i,t in enumerate([.0025,.0020000001]):
                q=run(s,16,t,f'N016_{f}_n{i}',yp=.06988 if i==1 else None);rows.append(q)
                if not q.get('accepted_stationary'):break
                s=q
        dump(OUT/'COARSE.json',rows)
    elif a.phase=='refine':
        rows=[]
        coarse=json.loads((OUT/'COARSE.json').read_text())
        # Refine all six converged families at both new density nodes, preserving geometry diagnostics.
        for old in coarse:
            if not old.get('accepted_stationary'):continue
            for n in (32,64):
                label=refine_label(old['id'],n)
                saved=OUT/label/'STATE.json'
                q=json.loads(saved.read_text()) if saved.exists() else run(old,n,old['n_b'],label)
                rows.append(q)
                if not q.get('accepted_stationary'):break
                old=q
        dump(OUT/'REFINE.json',rows)
    else:validate()

def validate():
    from compose_eos.spatial_tf import CellFunctional,RadialGrid,ELECTRON_MASS,fermions,E2
    allrows=[json.loads(p.read_text()) for p in sorted(OUT.glob('*/STATE.json'))]
    checks=[]
    for s in allrows:
        if 'profile' not in s:continue
        p=np.array(s['profile']);g=RadialGrid(s['points'],s['dimension']);f=CellFunctional('DD2',g,s['n_b'])
        ref=f.evaluate(p[:,0],p[:,1],p[:,2],s['radius'])
        threshold=max(1e-10,s['n_b']*1e-7)
        masks=[p[:,0]>threshold,p[:,1]>threshold]
        chem=[ref['mu_n'],ref['mu_p']+s['mu_e']]
        coeff=np.concatenate([g.weights[m] for m in masks])
        mub=float(coeff@np.concatenate([c[m] for c,m in zip(chem,masks)])/sum(coeff))
        kk=[np.max(abs(c[m]-mub),initial=0) for c,m in zip(chem,masks)]+[np.max(np.maximum(mub-c[~m],0),initial=0) for c,m in zip(chem,masks)]
        checks.append({'id':s['id'],'independent_KKT_MeV':float(max(kk)),'independent_mu_B_difference':abs(mub-s['mu_B']),'scalar_residual':float(ref['scalar_pde_residual']),'vector_residual':float(ref['vector_pde_residual']),'native_field_oracle_difference':float(np.max(abs(np.column_stack([ref['phi'],ref['omega'],ref['rho'],ref['coulomb']])-p[:,6:10]))),'energy_abs_difference':abs(ref['epsilon']-s['epsilon']),'chemical_max_difference':float(np.max(abs(np.column_stack([ref['mu_n'],ref['mu_p'],ref['mu_e']])-p[:,3:6]))),'radius_derivative_abs_difference':abs(ref['d_epsilon_d_log_radius']-s['d_epsilon_d_log_radius']),'baryon_residual':abs(g.mean(p[:,0]+p[:,1])-s['n_b']),'neutrality_residual':abs(g.mean(p[:,1]-p[:,2])),'electron_poisson_residual':float(np.max(abs(g.minus_laplacian@p[:,9]/s['radius']**2-4*np.pi*E2*(p[:,1]-p[:,2])))),'electron_equilibrium_residual':float(np.max(abs(fermions(p[:,2],ELECTRON_MASS)['ef']-p[:,9]-s['mu_e'])))})
    dump(OUT/'INDEPENDENT_RESIDUALS.json',{'oracle':'Independent retained Python functional; native Fortran remains producer','rows':checks,'pass':all(c['energy_abs_difference']<2e-9 and c['chemical_max_difference']<2e-7 and c['radius_derivative_abs_difference']<1e-9 and c['electron_poisson_residual']<2e-10 and c['neutrality_residual']<2e-11 for c in checks)})
    print(json.dumps({'independent_rows':len(checks),'max_energy_abs_difference':max(c['energy_abs_difference'] for c in checks)}),flush=True)
if __name__=='__main__':main()
