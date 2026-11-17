#!/usr/bin/env python3
"""Stage B2: bounded native homogeneous-core support audit, no extrapolation.

The continuous certificate applies to the declared cubic P(mu) interpolant.
Actual functional evidence consists of reused native knots, native interior
states/derivatives, and an independent Skyrme Python oracle; it is not a global
ground-state proof or a validated-interval proof for the underlying functional.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, platform, subprocess, sys, time
from pathlib import Path
from collections import Counter
import numpy as np
from numpy.polynomial import Polynomial as Poly
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.skyrme_uniform import beta_equilibrium, parameter_provenance
MODELS=('GM_GM1','RG_SLY4','RG_KDE0v1')
BASE=ROOT/'data/maxwell_three_families_20261117/native_core_extended_v2'
PINS={
'GM_GM1':'fefdbfdfd6a53041038533ce5812f0ceeb2503bdd4c5486fb85c6a24a56a5f3b',
'RG_SLY4':'b3cfc9e54324880f9a298559344b720adb60a00e6e5a8eb30693a34721a66d6d',
'RG_KDE0v1':'b44daaca894370bda0263a336b5c65814f79723a66882fa4211ed5257891e122'}
REP=ROOT/'data/whole_crust_20261117/production_v2/models/BSk22pp__SPG_M1__RG_SLY4__quintic_P__baseline/REPRESENTATION.json'
FORTRAN=r"""
program core_support_audit
 use compose_eos_types,only:wp
 use compose_eos_core_native,only:core_state_t,solve_core_state
 use compose_eos_core_domain,only:core_sound_speed_analytic
 implicit none
 type(core_state_t)::s,t(5)
 character(len=32)::model
 character(len=256)::status,ast
 real(wp)::n,h,de,dp,dm,cs,csm,fl,ca,fr,hy,be,w(5)
 integer::ios,k,j,shift,attempt
 logical::same
 write(*,'(a)') 'model,n_b,epsilon,pressure,mu_n,mu_p,mu_e,Yp,Ye,Ymu,cs2_fd,cs2_mu,firstlaw_relative,fd_h,cs2_analytic,fixed_proton_fraction_lepton_equilibrated_curvature_ratio,hyy,beta_gradient,phase,fd_status,analytic_status'
 do
  read(*,*,iostat=ios) model,n
  if(ios/=0)exit
  call solve_core_state(trim(model),n,s,status)
  if(status/='ok')then
   write(0,*) model,n,trim(status)
   stop 21
  endif
  h=n*1.e-4_wp
  same=.false.
  do attempt=1,24
   shift=-2;w=[1._wp,-8._wp,0._wp,8._wp,-1._wp]
   if(n-2*h<.02_wp)then
    shift=0;w=[-25._wp,48._wp,-36._wp,16._wp,-3._wp]
   else if(n+2*h>1._wp)then
    shift=-4;w=[3._wp,-16._wp,36._wp,-48._wp,25._wp]
   endif
   do k=1,5
    j=shift+k-1
    call solve_core_state(trim(model),n+j*h,t(k),status)
    if(status/='ok')then
     write(0,*) model,n+j*h,trim(status)
     stop 22
    endif
   enddo
   same=all(t%phase==s%phase).and.all((t%ymu>0).eqv.(s%ymu>0))
   if(same)exit
   h=h/2
  enddo
  if(.not.same)then
   write(0,*) 'could not keep derivative stencil within a composition branch',model,n
   stop 23
  endif
  de=sum(w*t%epsilon)/(12*h);dp=sum(w*t%pressure)/(12*h);dm=sum(w*t%mu_n)/(12*h)
  cs=dp/de;csm=n*dm/s%mu_n;fl=abs(de/s%mu_n-1)
  call core_sound_speed_analytic(trim(model),n,ca,fr,hy,be,ast)
  write(*,'(a,",",17(es24.16e3,","),i0,",ok,",a)') trim(model),n,s%epsilon,s%pressure,s%mu_n,s%mu_p,s%mu_e,s%yp,s%ye,s%ymu,cs,csm,fl,h,ca,fr,hy,be,s%phase,trim(ast)
 enddo
end program
"""
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v): Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def rows(p,strings=()):
    out=[]
    for r in csv.DictReader(Path(p).open()):
        out.append({k:(v if k in strings or v=='' else float(v)) for k,v in r.items()})
    return out
def csvwrite(p,records):
    if not records:return
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
def real_roots(p):
    c=np.trim_zeros(p.coef,'b')
    if len(c)<2:return []
    return [float(z.real) for z in Poly(c).roots() if abs(z.imag)<1e-9 and -1e-11<z.real<1+1e-11]
def extremum(p):
    points=[0.,1.]+real_roots(p.deriv())
    vals=[float(p(x)) for x in points]
    return min(vals),max(vals)
def ratio_extremum(a,b):
    points=[0.,1.]+real_roots(a.deriv()*b-a*b.deriv())
    vals=[float(a(x)/b(x)) for x in points]
    return min(vals),max(vals)
def interval_certificate(rr):
    result=[]
    for i,(l,r) in enumerate(zip(rr[:-1],rr[1:])):
        m0,m1=l['mu_n'],r['mu_n'];dm=m1-m0
        if dm<=0:raise ValueError('nonincreasing native chemical potential')
        delta=r['pressure']-l['pressure']
        p=Poly([l['pressure'],dm*l['n_b'],3*delta-dm*(2*l['n_b']+r['n_b']),-2*delta+dm*(l['n_b']+r['n_b'])])
        n=p.deriv()/dm;dn=n.deriv()/dm;mu=Poly([m0,dm])
        nrange=extremum(n);drange=extremum(dn);margin=extremum(mu*dn-n)
        cs=ratio_extremum(n,mu*dn);compress=ratio_extremum(n,dn)
        en=mu*n-p
        eerr=max(abs(en(0)-l['epsilon']),abs(en(1)-r['epsilon']))
        erel=eerr/max(1.,abs(l['epsilon']),abs(r['epsilon']))
        result.append(dict(interval=i,n_left=l['n_b'],n_right=r['n_b'],
            mu_left=m0,mu_right=m1,n_min=nrange[0],dn_dmu_min=drange[0],
            dn_dmu_max=drange[1],cs2_min=cs[0],cs2_max=cs[1],
            dP_dn_min=compress[0],dP_dn_max=compress[1],
            causal_margin_mu_dn_dmu_minus_n_min=margin[0],
            endpoint_energy_error=eerr,endpoint_energy_relative_error=erel,
            status='PASS' if nrange[0]>0 and drange[0]>0 and cs[0]>0 and cs[1]<1 and erel<1e-11 else 'FAIL'))
    return result
def analytic_onset(model):
    fun=lambda n:beta_equilibrium(model,float(n))['mu_e_mev']-105.6583755
    return float(brentq(fun,.02,1.,xtol=2e-14,rtol=2e-14))
def scalar_max(rr,key):return max(abs(float(r[key])) for r in rr)
def run(cmd,log,**kw):
    t=time.monotonic()
    p=subprocess.run(cmd,capture_output=True,text=True,**kw)
    Path(log).write_text(p.stdout+p.stderr)
    rec=dict(argv=[str(x) for x in cmd],exit_code=p.returncode,seconds=time.monotonic()-t,log=str(Path(log).relative_to(ROOT)),log_sha256=sha(log))
    if p.returncode:raise RuntimeError(json.dumps(rec))
    return p,rec
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='data/core_support_20261117/domain_audit')
    args=parser.parse_args();out=(ROOT/args.output).resolve();out.mkdir(parents=True,exist_ok=True)
    report={'schema':'ns-cowling.core-domain-audit.v1','stage':'B2','status':'IN_PROGRESS',
        'repository_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'assumptions':{'temperature_MeV':0,'neutrinos':'untrapped','geometry':'uniform spin-saturated homogeneous nuclear matter',
            'units':{'n':'fm^-3','epsilon_P':'MeV fm^-3','mu':'MeV','sound_speed_ratio':'c_s^2/c^2; c retained in definition'},
            'energy_reference':'original source rest energies; no shift or fitted calibration',
            'selection':'neutral beta-equilibrium roots searched on 0<=Yp<=1/2; finite root search is not a global minimum proof',
            'continuous_representation':'cubic Hermite P(mu), P and dP/dmu=n fixed at both native knots; epsilon=mu*n-P',
            'derivative':'d epsilon/dn=mu; dP/dn=n/(dn/dmu); c_s^2/c^2=n/(mu dn/dmu)',
            'proof_scope':'all stationary candidates of each binary64 polynomial interval are evaluated; no outward-rounded interval arithmetic'},
        'source_sha256':{},'models':{},'execution':[],
        'limits':['No TOV, maximum-mass, radial-stability or Cowling claim.',
                  'No global inhomogeneous ground-state or crust-transition-density determination.',
                  'No spin, gradient, finite-wavelength or pairing stability claim.',
                  'No source-provider equivalence claim; unresolved GM1 calibration/provider residual is preserved.',
                  'No support beyond actual n=1.0 fm^-3 certified; metadata upper bounds are computational capability only.',
                  'Positive beta-path compressibility is not full microscopic stability.',
                  'Native interior derivative evidence and polynomial certificate are distinct from a rigorous continuum theorem for the functional.',
                  'GM1 is intentionally npe; excluded muons are not added silently.',
                  'The native frozen diagnostic fixes proton fraction while leptons equilibrate at combined susceptibility; it is not a fully frozen-species buoyancy closure.']}
    sources=[ROOT/'fortran_eos/src'/x for x in ('compose_eos_types.f90','compose_eos_core_native.f90','compose_eos_core_domain.f90')]
    sources += [ROOT/'python_frontend/compose_eos/skyrme_uniform.py',ROOT/'python_frontend/compose_eos/cold_leptons.py',Path(__file__).resolve(),REP,ROOT/'python_frontend/tests/test_core_support_domain_audit.py']
    report['source_sha256']={str(p.relative_to(ROOT)):sha(p) for p in sources}
    native={};target=[];onsets={}
    for model in MODELS:
        p=BASE/model;receipt=json.loads((p/'RECEIPT.json').read_text())
        assert sha(p/'RECEIPT.json')==PINS[model], 'canonical receipt mismatch'
        rr=rows(p/'generated/states.csv');native[model]=rr
        assert len(rr)==385 and rr[0]['n_b']==.02 and rr[-1]['n_b']==1.
        assert np.all(np.diff([r['n_b'] for r in rr])>0)
        expected=receipt['provenance']['evidence']['output_sha256']
        file_checks={name:sha(p/'generated'/name)==value for name,value in expected.items()}
        assert all(file_checks.values()), 'historical generated file identity mismatch'
        provenance_sources=receipt['provenance']['source_sha256']
        source_checks={}
        for f in sources[:3]:
            matched=[h for old,h in provenance_sources.items() if old.endswith('/'+str(f.relative_to(ROOT)))]
            source_checks[str(f.relative_to(ROOT))]={'current':sha(f),'historical':matched,'match':matched==[sha(f)]}
        assert all(v['match'] for v in source_checks.values()),'native implementation drift'
        delta=max(abs(r[k]-s[k]) for r,s in zip(rr,receipt['rows']) for k in r)
        assert delta==0,'receipt rows and states CSV differ'
        intervals=interval_certificate(rr);csvwrite(out/(model+'_INTERVALS.csv'),intervals)
        ss=rows(p/'generated/sound_speed.csv',strings=('fd_status','analytic_status'))
        usable=[r for r in ss if r['fd_status']=='ok']
        original={key:[min(r[key] for r in usable),max(r[key] for r in usable)] for key in ('cs2_over_c2_fd','cs2_over_c2_mu','firstlaw_relative')}
        mm=dict(species=receipt['species'],source_representation=receipt['representation'],
            receipt_sha256=sha(p/'RECEIPT.json'),input_sha256={str(f.relative_to(ROOT)):sha(f) for f in p.glob('generated/*') if f.is_file()},
            existing_generated_sha256_checks=file_checks,historical_native_source_identity=source_checks,
            native_rows=385,actual_generated_domain_fm3=[.02,1.],
            declared_computational_interval_fm3=receipt['computational_interval_fm3'],
            sampled_beta_residual_max_MeV=max(abs(r['mu_n']-r['mu_p']-r['mu_e']) for r in rr),
            charge_residual_max_fm3=max(abs(r['n_b']*(r['Yp']-r['Ye']-r['Ymu'])) for r in rr),
            baryon_residual_max=max(abs(r['Yn']+r['Yp']-1) for r in rr),
            euler_residual_max_MeV_fm3=max(abs(r['epsilon']+r['pressure']-r['n_b']*r['mu_n']) for r in rr),
            original_fd_ranges=original,original_fd_status_counts=dict(Counter(r['fd_status'] for r in ss)),
            original_analytic_status_counts=dict(Counter(r['analytic_status'] for r in ss)),
            original_first_law_step_refinement=receipt['first_law_independent_differences'],
            interval_certificate={'status':'PASS' if all(r['status']=='PASS' for r in intervals) else 'FAIL','intervals':len(intervals),
                'density_domain_fm3':[.02,1.],
                'cs2_range':[min(r['cs2_min'] for r in intervals),max(r['cs2_max'] for r in intervals)],
                'dn_dmu_range':[min(r['dn_dmu_min'] for r in intervals),max(r['dn_dmu_max'] for r in intervals)],
                'dP_dn_range':[min(r['dP_dn_min'] for r in intervals),max(r['dP_dn_max'] for r in intervals)],
                'minimum_causal_margin_fm3':min(r['causal_margin_mu_dn_dmu_minus_n_min'] for r in intervals),
                'endpoint_energy_error_max_MeV_fm3':max(r['endpoint_energy_error'] for r in intervals),
                'endpoint_energy_relative_error_max':max(r['endpoint_energy_relative_error'] for r in intervals),
                'endpoint_energy_relative_tolerance':1e-11,
                'scope':'continuous defined interpolant; homogeneous functional between audit samples is not rigorously enclosed'})
        if model!='GM_GM1':
            onsets[model]=analytic_onset(model);mm['parameter_provenance']=parameter_provenance(model)
            mm['muon_onset_density_fm3']=onsets[model]
            mm['muon_onset_class']='continuous free-lepton threshold, not a sharp nuclear first-order interface'
            cross=[r for r in intervals if r['n_left']<=onsets[model]<=r['n_right']]
            mm['muon_onset_native_interval']=cross[0]['interval']
        else:
            mm['muon_policy']='npe only even when mu_e exceeds muon rest energy; separate physical species model'
        report['models'][model]=mm
        for n in [.02]+[float(np.sqrt(a['n_b']*b['n_b'])) for a,b in zip(rr[:-1],rr[1:])]+[1.]:
            target.append((model,n))
        print(model,'native rows reused, continuous polynomial audit',mm['interval_certificate']['status'],flush=True)
    # Supplemental native states resolve interval behavior and both endpoint derivatives.
    build=out/'build';build.mkdir(exist_ok=True)
    src=out/'native_domain_probe.f90';src.write_text(FORTRAN)
    exe=build/'native_domain_probe'
    cmd=['gfortran','-std=f2008','-O2','-Wall','-Wextra','-fcheck=all','-ffree-line-length-none','-J',str(build),'-I',str(build)]
    cmd += [str(p) for p in sources[:3]]+[str(src),'-o',str(exe)]
    _,rec=run(cmd,out/'NATIVE_BUILD.log',cwd=build);report['execution'].append(rec)
    request=''.join(f'{m} {n:.17g}\n' for m,n in target)
    (out/'NATIVE_INTERIOR_INPUT.txt').write_text(request)
    proc,rec=run([str(exe)],out/'NATIVE_INTERIOR_RUN.log',input=request,cwd=ROOT)
    report['execution'].append(rec)
    (out/'NATIVE_INTERIOR.csv').write_text(proc.stdout)
    probe=rows(out/'NATIVE_INTERIOR.csv',strings=('model','fd_status','analytic_status'))
    assert len(probe)==len(target)==1158
    report['native_probe_identity']={'source_sha256':sha(src),'executable_sha256':sha(exe),'input_sha256':sha(out/'NATIVE_INTERIOR_INPUT.txt'),
        'purpose':'384 geometric interval interiors plus both endpoints/model; generated only to diagnose actual-domain continuum and endpoint derivatives',
        'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0]}
    aliases={'epsilon':'energy_density_mev_fm3','pressure':'pressure_mev_fm3','mu_n':'mu_n_mev','mu_p':'mu_p_mev','mu_e':'mu_e_mev','Yp':'Y_p','Ye':'Y_e','Ymu':'Y_mu'}
    oracle=[]
    for model in MODELS:
        mm=report['models'][model];pr=[r for r in probe if r['model']==model]
        mm['native_interval_probe']={'samples':len(pr),'phase_codes':dict(Counter(str(int(r['phase'])) for r in pr)),
            'cs2_fd_range':[min(r['cs2_fd'] for r in pr),max(r['cs2_fd'] for r in pr)],
            'minimum_dP_dn_MeV':min(r['cs2_fd']*r['mu_n'] for r in pr),
            'max_firstlaw_relative':max(r['firstlaw_relative'] for r in pr),
            'max_cs2_fd_minus_mu':max(abs(r['cs2_fd']-r['cs2_mu']) for r in pr),
            'beta_residual_max_MeV':max(abs(r['mu_n']-r['mu_p']-r['mu_e']) for r in pr),
            'charge_residual_max_fm3':max(abs(r['n_b']*(r['Yp']-r['Ye']-r['Ymu'])) for r in pr),
            'all_derivative_samples_within_actual_domain':True,
            'stencil':'5-point, relative h=1e-4, one-sided at actual endpoints; shrink to avoid muon/phase threshold crossing'}
        assert all(0<r['cs2_fd']<1 and r['firstlaw_relative']<1e-7 for r in pr)
        if model!='GM_GM1':
            mm['native_interval_probe'].update(analytic_cs2_range=[min(r['cs2_analytic'] for r in pr),max(r['cs2_analytic'] for r in pr)],
                composition_hessian_min_MeV_fm3=min(r['hyy'] for r in pr),
                max_cs2_fd_minus_analytic=max(abs(r['cs2_fd']-r['cs2_analytic']) for r in pr),
                beta_gradient_max_MeV_fm3=max(abs(r['beta_gradient']) for r in pr))
            assert all(r['hyy']>0 and abs(r['cs2_fd']-r['cs2_analytic'])<1e-7 for r in pr)
            diffs={k:[] for k in aliases};candidate=[]
            for origin,seq in [('reused_native_knot',native[model]),('native_interval_probe',pr)]:
                for r in seq:
                    py=beta_equilibrium(model,r['n_b']);candidate.append(py['candidate_count'])
                    er={k:abs(py[v]-r[k]) for k,v in aliases.items()}
                    for k,v in er.items():diffs[k].append(v)
                    oracle.append(dict(model=model,origin=origin,n_b=r['n_b'],candidate_count=py['candidate_count'],**{k+'_abs_error':v for k,v in er.items()}))
            mm['independent_python_oracle']={'rows':len(candidate),'max_abs_error':{k:max(v) for k,v in diffs.items()},
                'candidate_counts':dict(Counter(candidate)),
                'scope':'independent implementation, shared published functional and conventions; no independent provider/calibration validation'}
            assert max(diffs['epsilon'])<1e-8 and max(diffs['pressure'])<1e-8 and max(diffs['Yp'])<1e-10
        print(model,'native midpoint/endpoint derivatives complete',flush=True)
    csvwrite(out/'SKYRME_ORACLE.csv',oracle)
    upper=json.loads(REP.read_text())['upper'];n=upper['n'][-1];py=beta_equilibrium('RG_SLY4',n)
    report['A_upper_anchor_identity']={'density_fm3':n,'A_P_MeV_fm3':upper['p'][-1],'native_oracle_P_MeV_fm3':py['pressure_mev_fm3'],
        'A_epsilon_MeV_fm3':upper['epsilon'][-1],'native_oracle_epsilon_MeV_fm3':py['energy_density_mev_fm3'],
        'pressure_difference_MeV_fm3':py['pressure_mev_fm3']-upper['p'][-1],
        'energy_difference_MeV_fm3':py['energy_density_mev_fm3']-upper['epsilon'][-1],
        'interpretation':'difference ledger only; extension join policy remains B3 responsibility; no energy shift or retuning'}
    report['GM1_source_identity']={'mass_MeV':939.,'species':'npe','couplings':{'Cs2_fm2':11.79,'Cv2_fm2':7.149,'Crho2_fm2':4.411,'b':.002947,'c':-.001070},
        'provider_calibration_residual':'UNRESOLVED_PRESERVED; this audit does not fit or repair source identity'}
    report['domain_decision']={'SLy4':'[0.02,1.0] native homogeneous branch numerically checked; B3 may use validated upper portion to 1.0 subject to join identity',
        'KDE0v1':'[0.02,1.0] same finite scope; 1.8 metadata is not actual support',
        'GM1':'[0.02,1.0] npe branch finite scope; no beyond-1.0 capability, analytic Hessian unavailable',
        'upper_endpoint_class':'actual input/generated-support endpoint; not a causal crossing, maximum mass, or microscopic validity boundary',
        'subsaturation':'uniform state may be metastable against cluster formation; use only as explicitly scoped homogeneous branch, not as replacement for A crust'}
    report['environment']={'python':sys.version,'numpy':np.__version__,'platform':platform.platform()}
    report['status']='COMPLETE_BOUNDED_B2'
    generated_names=[m+'_INTERVALS.csv' for m in MODELS]+['native_domain_probe.f90','NATIVE_BUILD.log','NATIVE_INTERIOR_INPUT.txt','NATIVE_INTERIOR_RUN.log','NATIVE_INTERIOR.csv','SKYRME_ORACLE.csv']
    report['artifacts_sha256']={name:sha(out/name) for name in generated_names}
    dump(out/'REPORT.json',report)
    lines=['# B2 native core support audit','',
        'Status: COMPLETE_BOUNDED_B2. All three actual native inputs contain 385 states on [0.02, 1.0] fm^-3. The higher Skyrme metadata limits are not generated support.','',
        '## Numerical evidence','',
        '| Model | species | continuous cubic P(mu) cs2/c2 | native probe cs2/c2 | max relative first-law residual |',
        '|---|---|---|---|---|']
    for m,v in report['models'].items():
        lines.append(f"| {m} | {v['species']} | {v['interval_certificate']['cs2_range']} | {v['native_interval_probe']['cs2_fd_range']} | {v['native_interval_probe']['max_firstlaw_relative']:.6e} |")
    lines += ['', 'Every one of the 384 intervals/model was evaluated at endpoints and every real stationary candidate of the relevant polynomial or rational derivative. This certifies the declared interpolant numerically; it does not prove the full underlying functional is stable at every real density. Native geometric midpoints and one-sided endpoint derivative stencils provide additional functional evidence. No densities outside [0.02,1.0] were sampled by the native probe.',
        '', '## Thermodynamic construction', '',
        'With x=(mu-mu0)/Delta_mu, P(x)=P0+Delta_mu*n0*x+[3*Delta_P-Delta_mu*(2*n0+n1)]*x^2+[-2*Delta_P+Delta_mu*(n0+n1)]*x^3. Define n=dP/dmu and epsilon=mu*n-P. Then d epsilon/dn=mu, dP/dn=n/(dn/dmu), and cs2/c2=n/(mu*dn/dmu). Native source states and energies are retained at knots to floating point roundoff.',
        '', '## Interpretation and exclusions', '']
    lines += ['- '+v for v in report['limits']]
    lines += ['', 'The upper endpoint 1.0 fm^-3 is an actual-input boundary, not a measured causal crossing. Subsaturation homogeneous states can be metastable relative to nonuniform matter; a beta-path mechanical check does not determine the crust transition. Skyrme muon onset is a continuous lepton threshold, not a first-order nuclear transition. GM1 retains its intentional npe species convention.', '',
        'Full exact hashes, residuals, onset locations, native probe identities and A-anchor differences are in REPORT.json. Reproduce with:',
        '', '    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python scripts/audit_core_support_20261117.py', '']
    (out/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({'status':report['status'],'output':str(out),'A_anchor':report['A_upper_anchor_identity'],'models':{m:{'continuous':v['interval_certificate']['cs2_range'],'native':v['native_interval_probe']['cs2_fd_range']} for m,v in report['models'].items()}},indent=2))
if __name__=='__main__':
    main()
