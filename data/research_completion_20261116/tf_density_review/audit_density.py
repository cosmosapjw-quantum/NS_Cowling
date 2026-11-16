"""Independent density interpolation audit at fixed N=192; no production writes."""
from pathlib import Path
import argparse,gzip,hashlib,json
import numpy as np
from scipy.interpolate import PchipInterpolator
COARSE=np.array([.002,.005,.01,.02,.03,.04,.05,.06,.07,.08,.09,.10])
FIELDS=('pressure','epsilon','Y_e')
METHODS=('linear','pchip','log_pchip')
def read(path):
    return json.loads(gzip.decompress(path.read_bytes())) if path.suffix=='.gz' else json.loads(path.read_text())
def good(row):
    return row.get('numerical_status') in ('STATIONARY','UNIFORM_STATIONARY') and not row.get('muon_threshold_crossed',True)
def geometry(row):
    if row['baryon_density_contrast']<1e-4:return 'uniform'
    p=row['profile'];dense=p['nn'][0]+p['np'][0]>p['nn'][-1]+p['np'][-1]
    return {1:'slab',2:'rod' if dense else 'tube',3:'droplet' if dense else 'bubble'}[row['dimension']]
def all_attempts(directory):
    result=[]
    for name in ('attempts.json.gz','repair_attempts.json.gz'):
        path=directory/name
        if path.exists():
            for group in read(path):result.extend(group['attempts'])
    return result

def fixed_mesh_states(directory):
    selected=read(directory/'selected_states.json.gz');attempts=all_attempts(directory)
    available=[a for a in attempts if a.get('points')==192 and good(a)]
    states=[]
    for n in np.sort(np.r_[COARSE,(COARSE[:-1]+COARSE[1:])/2]):
        rows=[s for s in selected if abs(s['n_b']-n)<1e-13 and s['points']==192]
        if not rows:rows=[s for s in available if abs(s['n_b']-n)<1e-13]
        if not rows:raise ValueError(f'{directory.name}: no accepted N192 state at {n}')
        chosen=min(rows,key=lambda s:s['epsilon']);assert good(chosen)
        states.append(chosen)
    return states,attempts

def audit(directory):
    states,attempts=fixed_mesh_states(directory)
    coarse=states[::2];mids=states[1::2]
    nodes=np.array([s['n_b'] for s in coarse]);xn=np.array([s['n_b'] for s in mids]);predictions={}
    for field in FIELDS:
        y=np.array([s[field] for s in coarse])
        if np.any(y<=0):raise ValueError('log interpolation requires positive values')
        predictions[field]={'linear':np.interp(xn,nodes,y),'pchip':PchipInterpolator(nodes,y)(xn),'log_pchip':np.exp(PchipInterpolator(np.log(nodes),np.log(y))(np.log(xn)))}
    rows=[]
    for j,s in enumerate(mids):
        left,right=coarse[j:j+2];g=[geometry(v) for v in (left,s,right)];competitors={}
        for a in attempts:
            if a.get('points')==192 and good(a) and abs(a['n_b']-s['n_b'])<1e-13:
                kind=geometry(a);competitors[kind]=min(competitors.get(kind,float('inf')),a['energy_per_baryon'])
        gaps={kind:e-s['energy_per_baryon'] for kind,e in competitors.items() if kind!=geometry(s)}
        record={'n_b':s['n_b'],'interval':[left['n_b'],right['n_b']],'geometry_triple':g,
                'phase_change':len(set(g))>1,'alternative_geometry_energy_gaps_MeV_per_baryon':gaps,
                'near_degenerate_competitor':any(abs(v)<1e-4 for v in gaps.values()),
                'max_chemical_KKT_residual_MeV':s['max_chemical_KKT_residual_MeV'],'errors':{}}
        for field in FIELDS:
            record['errors'][field]={'actual':s[field]}
            for method in METHODS:
                delta=float(predictions[field][method][j]-s[field])
                record['errors'][field][method]={'signed_absolute':delta,'absolute':abs(delta),'relative':abs(delta)/max(abs(s[field]),1e-30)}
        rows.append(record)
    summary={}
    for group,filtered in [('all',rows),('same_geometry',[r for r in rows if not r['phase_change']]),('phase_change',[r for r in rows if r['phase_change']])]:
        summary[group]={'count':len(filtered)}
        if filtered:
            for field in FIELDS:
                summary[group][field]={method:{'max_relative':max(r['errors'][field][method]['relative'] for r in filtered),'max_absolute':max(r['errors'][field][method]['absolute'] for r in filtered)} for method in METHODS}
    return {'model':directory.name,'source_density_nodes':23,'independently_optimized_midpoints':11,'mesh_points':192,'summary':summary,'midpoints':rows}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--campaign',type=Path,default=Path('data/research_completion_20261116/gppva_tf/final'));parser.add_argument('--output',type=Path,default=Path(__file__).with_name('AUDIT.json'));a=parser.parse_args()
    inputs=[a.campaign/'SUMMARY.json',a.campaign/'INPUT.json']+[p for d in sorted(a.campaign.glob('GPPVA_*')) for p in d.iterdir() if p.name in ('selected_states.json.gz','attempts.json.gz','repair_attempts.json.gz')]
    result={'status':'DENSITY_INTERPOLATION_AUDIT_COMPLETED','fixed_mesh_points':192,'scope':'within final unchanged solver campaign, source interval0.002..0.1fm^-3; no external EOS agreement gate',
      'density_convergence_claim':'One true midpoint refinement; measured interpolation differences, not asymptotic convergence proof.',
      'degeneracy_policy':'Alternative accepted geometry within1e-4MeV/baryon is flagged descriptively; not an acceptance threshold.',
      'source_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
      'audit_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'models':[audit(d) for d in sorted(a.campaign.glob('GPPVA_*'))]}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for m in result['models']:
        s=m['summary'];print(m['model'],'P_pchip_all',s['all']['pressure']['pchip']['max_relative'],'P_pchip_same_geometry',s['same_geometry']['pressure']['pchip']['max_relative'],'phase_changes',s['phase_change']['count'],flush=True)
if __name__=='__main__':main()
