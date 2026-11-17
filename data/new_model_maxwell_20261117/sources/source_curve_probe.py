#!/usr/bin/env python3
"""Read-only source curve reconnaissance; no EOS generation or Maxwell construction.
All interpolation is exploratory. Never bridges phase boundaries or invalid edges.
"""
from pathlib import Path
import json, hashlib
import numpy as np
ROOT=Path('/home/cosmosapjw/ns_cowling_wt/stellar_extension_20261117')
OUT=Path(__file__).resolve().parent

def load(path):
    n=np.loadtxt(path/'eos.nb',skiprows=2)
    rows=np.loadtxt(path/'eos.thermo',skiprows=1)
    mn,mp=map(float,(path/'eos.thermo').read_text().splitlines()[0].split()[:2])
    p=n*rows[:,3];e=n*mn*(1+rows[:,9]);mu=(e+p)/n
    phases=np.array([int(x.split()[3]) for x in (path/'eos.compo').read_text().splitlines() if x.strip()]) if (path/'eos.compo').exists() else np.zeros(len(n),int)
    return dict(n=n,p=p,e=e,mu=mu,phases=phases,mn=mn,mp=mp,muq=mn*(1+rows[:,5]),q5=rows[:,7],rows=rows)

def segment(d,a,b,identity,role):
    return {'id':identity,'role':role,'start':int(a),'stop':int(b),'n':d['n'][a:b],'p':d['p'][a:b],'mu':d['mu'][a:b]}

def own_segments(d,identity,role):
    sel=(d['phases']==0) if role=='core' else ((d['phases']!=0)&(d['n']>=0.002))
    good=sel[:-1]&sel[1:]&(np.diff(d['phases'])==0)&(np.diff(d['p'])>0)&(np.diff(d['mu'])>0)&(np.diff(d['e'])>0)&(np.diff(d['p'])<=np.diff(d['e']))
    blocks=[];a=None
    for i,v in enumerate(np.r_[good,False]):
        if v and a is None:a=i
        if not v and a is not None:
            if i+1-a>=3:blocks.append(segment(d,a,i+1,identity,role))
            a=None
    return blocks

def probe(l,h):
    a=max(l['mu'][0],h['mu'][0]);b=min(l['mu'][-1],h['mu'][-1])
    if b<=a:return None
    # Union knots resolves every sign change of the two piecewise linear source curves.
    grid=np.unique(np.r_[a,b,l['mu'][(l['mu']>a)&(l['mu']<b)],h['mu'][(h['mu']>a)&(h['mu']<b)]])
    dp=np.interp(grid,h['mu'],h['p'])-np.interp(grid,l['mu'],l['p'])
    bracket=[]
    for j in np.flatnonzero(dp[:-1]*dp[1:]<0):
        mid=.5*(grid[j]+grid[j+1]);nl=float(np.interp(mid,l['mu'],l['n']));nh=float(np.interp(mid,h['mu'],h['n']))
        bracket.append({'mu_bracket_MeV':grid[j:j+2].tolist(),'deltaP_bracket_MeV_fm3':dp[j:j+2].tolist(),'orientation':'CORE_OVERTAKES' if dp[j]<0 else 'REVERSE','n_crust_mid_fm3':nl,'n_core_mid_fm3':nh,'delta_n_mid_fm3':nh-nl})
    return dict(crust=l['id'],core=h['id'],crust_rows=[l['start'],l['stop']],core_rows=[h['start'],h['stop']],mu_overlap_MeV=[float(a),float(b)],deltaP_range_MeV_fm3=[float(dp.min()),float(dp.max())],min_abs_deltaP=float(np.min(np.abs(dp))),brackets=bracket)

manifest=ROOT/'data/maxwell_three_families_20261117/prepared_v2/MANIFEST.json'
old=json.loads(manifest.read_text());oldsegments=[]
for record in old['source_records']:
    if record['family'] not in ['native_generated','public_original'] or not record.get('components'):continue
    d=load(ROOT/record['compose_path'])
    for c in record['components']:
        a,b=c['start_zero_based'],c['stop_exclusive']
        oldsegments.append(segment(d,a,b,record['id'],record['role']))
newsegments=[];inventory=[]
for path in sorted(OUT.glob('XMLSLZ_*')):
    if not (path/'eos.thermo').exists():continue
    d=load(path);segs=[]
    for role in ['crust','core']:segs+=own_segments(d,'public_original:'+path.name,role)
    newsegments+=segs
    inventory.append({'model':path.name,'mass_header_MeV':[d['mn'],d['mp']],'q5_max_abs':float(np.max(np.abs(d['q5']))),'muQ_vs_muE_max_abs_MeV':float(np.max(np.abs(d['muq']-d['mu']))),'components':[{'role':s['role'],'rows':[s['start'],s['stop']],'n_fm3':[float(s['n'][0]),float(s['n'][-1])],'mu_MeV':[float(s['mu'][0]),float(s['mu'][-1])]} for s in segs]})
results=[]
for n in newsegments:
    for o in oldsegments:
        if n['role']==o['role']:continue
        l,h=(n,o) if n['role']=='crust' else (o,n)
        r=probe(l,h)
        if r:r['comparison']='RAW_PUBLIC' if o['id'].startswith('public_original') else 'CROSS_REPRESENTATION_EXPLORATORY';results.append(r)
report={'status':'EXPLORATORY_SOURCE_CURVE_SCREEN_ONLY','scientific_admission':'NOT_RUN','generation':'NOT_RUN','maxwell_construction':'NOT_RUN','energy_changes':'NONE','mu_convention':'mu_E=(epsilon+P)/n; exact source rest-energy normalization retained','new_segment_rule':'partner-independent phase blocks with positive source secants dP/dn,dmu/dn,de/dn and 0<dP/de<=1; not continuous derivative certification','old_segment_rule':'reuse prepared_v2 admitted components','manifest_path':str(manifest),'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'inventory':inventory,'comparisons':results}
(OUT/'SOURCE_CURVE_SCREEN.json').write_text(json.dumps(report,indent=2)+'\n')
for i in inventory:print(json.dumps(i))
for r in results:
    if r['brackets']:print(json.dumps(r))
print('summary',json.dumps({'comparisons':len(results),'sign_bracket_cases':sum(bool(r['brackets']) for r in results),'core_overtakes':sum(sum(b['orientation']=='CORE_OVERTAKES' for b in r['brackets']) for r in results),'file':str(OUT/'SOURCE_CURVE_SCREEN.json')}))
