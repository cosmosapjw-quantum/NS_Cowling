from __future__ import annotations
import sys, json, hashlib, csv, math
from pathlib import Path
from collections import Counter
import numpy as np
sys.path.insert(0,'/tmp/eos_revised_adversarial_20261106/DELIVERY_v3_3/compose_eos_pkg')
from compose_eos.io import read_compose
from compose_eos.causality import causality_profile
BASE=Path('/tmp/eos_revised_adversarial_20261106')
MATCH=BASE/'matched_eos_v3_2'
CLASS=json.loads((BASE/'release_v3_3/MATCHING_FAILURE_CAUSE_CLASSIFICATION_v3_3.json').read_text())
OUT=Path('/tmp/ns_oscillation_followup_20261106')
DAT=OUT/'normalized_tables'
DAT.mkdir(parents=True,exist_ok=True)
HEADER='source_row,i_t,i_nb,i_yq,n_b_fm3,n_b_m3,pressure_mev_fm3,pressure_pa,energy_mev_fm3,energy_j_m3,mu_b_mev,mu_b_j,q1_p_over_nb_mev,q2_entropy_per_baryon,q3_mu_b_over_mn_minus_1,q4_mu_q_over_mn,q5_mu_l_over_mn,q6_f_over_nb_mn_minus_1,q7_e_over_nb_mn_minus_1'
MEVFM3_TO_JM3=1.602176634e32
MEV_TO_J=1.602176634e-13
rows=[]
count=Counter()
for rec in CLASS['results']:
    d=rec['directory']
    meta_path=MATCH/d/'MATCHING_METADATA.json'
    if not meta_path.exists():
        count['MISSING_METADATA']+=1;continue
    meta=json.loads(meta_path.read_text())
    if rec.get('status')!='PASS' or meta.get('status')!='PASS':
        count['NOT_JOIN_PASS']+=1;continue
    if meta.get('diagnostic_only',False):
        count['DIAGNOSTIC_ONLY']+=1;continue
    p=MATCH/d
    try:
        eos=read_compose(str(p))
        n=np.asarray(eos.nb,float)
        P=np.asarray(eos.p,float)
        eps=np.asarray(eos.eps,float)
        mu=np.asarray(eos.mub,float)
        if not (eos.n>=3 and np.all(np.isfinite(n+P+eps+mu))):
            raise ValueError('nonfinite/short table')
        # Independent cold barotropic first-law identity from CompOSE Q3/Q7.
        ident=np.abs(eps+P-mu*n)/np.maximum(np.abs(eps+P),1e-300)
        identity=float(np.max(ident))
        if identity>2e-12: raise ValueError(f'first-law identity residual={identity:.3e}')
        if np.any(np.diff(n)<=0) or np.any(np.diff(P)<=0) or np.any(np.diff(eps)<=0):
            raise ValueError('density/pressure/energy not strictly increasing')
        # Direct sampled slope is a prefilter only; native P(mu) cell audit follows.
        cs2=np.diff(P)/np.diff(eps)
        if np.any(~np.isfinite(cs2)) or np.any(cs2<=0):
            raise ValueError('invalid sampled cs2')
        source_hashes={f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in ['eos.nb','eos.thermo','eos.t','eos.yq']}
        out=DAT/(d+'.csv')
        with out.open('w',newline='') as fh:
            w=csv.writer(fh);w.writerow(HEADER.split(','))
            for i in range(eos.n):
                vals=[
                  i+1,eos.iT,eos.inb_start+i,eos.iYq,n[i],n[i]*1e45,
                  P[i],P[i]*MEVFM3_TO_JM3,eps[i],eps[i]*MEVFM3_TO_JM3,
                  mu[i],mu[i]*MEV_TO_J,
                  eos.Q1[i],eos.Q2[i],eos.Q3[i],eos.Q4[i],eos.Q5[i],eos.Q6[i],eos.Q7[i]
                ]
                w.writerow([f'{v:.17g}' if isinstance(v,(float,np.floating)) else v for v in vals])
        rows.append({
          'directory':d,'core':meta.get('source_core'),'crust':meta.get('source_crust'),
          'method':meta.get('method'),'rows':int(eos.n),'table_path':str(out),
          'target_Msun':1.4,'mn_MeV':eos.mn,'mp_MeV':eos.mp,'Il':eos.Il,
          'beta_equilibrated':bool(eos.beta_equilibrated),
          'first_law_rel_max':identity,'cs2_sample_min':float(np.min(cs2)),
          'cs2_sample_max':float(np.max(cs2)),
          'n_max_fm3':float(n[-1]),'P_max_MeV_fm3':float(P[-1]),'eps_max_MeV_fm3':float(eps[-1]),
          'mu_min_MeV':float(mu[0]),'mu_max_MeV':float(mu[-1]),
          'matching_metadata_sha256':hashlib.sha256(meta_path.read_bytes()).hexdigest(),
          'eos_source_hashes':source_hashes,
          'composition_available':bool(eos.composition_available),
        })
        count['PREPARED']+=1
    except Exception as e:
        count['COMPOSE_INTAKE_REJECTED']+=1
        rows.append({'directory':d,'core':meta.get('source_core'),'crust':meta.get('source_crust'),
                     'method':meta.get('method'),'preparation_status':'REJECTED','reason':f'{type(e).__name__}: {e}'})
(OUT/'prepared_eos_manifest.json').write_text(json.dumps({
 'schema':'neutron-star-oscillation-prepared-eos-v1',
 'source_matrix':'/tmp/eos_revised_adversarial_20261106/matching_matrix_v32.json',
 'failure_classifier':'/tmp/eos_revised_adversarial_20261106/release_v3_3/MATCHING_FAILURE_CAUSE_CLASSIFICATION_v3_3.json',
 'n_classified':CLASS['n_results'],'counts':dict(count),
 'conversion':'CompOSE Q1..Q7; P=nb*Q1, eps=nb*mn*(1+Q7), muB=mn*(1+Q3); SI factor 1 MeV fm^-3=1.602176634e32 J m^-3, n conversion 1e45, mu conversion 1.602176634e-13',
 'identity':'eps+P=muB*nb checked rowwise to max relative residual 2e-12',
 'entries':rows
},indent=2))
print('COUNTS',dict(count))
print('PREPARED',sum('table_path' in x for x in rows))
print('REJECTED',sum(x.get('preparation_status')=='REJECTED' for x in rows))
for x in rows:
 if x.get('preparation_status')=='REJECTED': print('REJECT',x['directory'],x['reason'])
