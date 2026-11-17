"""Unmodified source-inner Gibbs diagnostics against the anchored BSk22 branch."""
from pathlib import Path
import json,gzip,hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator
from compose_eos.bsk22_thermodynamic_reconstruction import BoundaryAnchoredBSk22
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
root=ROOT/'data/research_completion_20261116'
d=json.loads((root/'outer_crust/generated/RECEIPT.json').read_text())['drip']
b=BoundaryAnchoredBSk22(outer_density=d['n_b'],outer_pressure=d['pressure'],outer_mu_b=d['g'],outer_epsilon=d['epsilon'])
print('BSK',b.state(b.n0)['pressure'],b.state(.002)['pressure'],b.state(.002)['mu_b_euler'])
results={}
for name in [f'SPG_M{i}' for i in range(1,6)]+['GPPVA_DD2','GPPVA_DDME2','GPPVA_TW','GPPVA_FSU2','GPPVA_FSU2R','GPPVA_FSU2H','GPPVA_NL3wrL55','GPPVA_TM1e']:
    if name.startswith('SPG'):
        path=root/f'spg/{name}/all_attempts.json.gz'
        rows=json.loads(gzip.decompress(path.read_bytes()))['fine']
        n=np.array([r['n_B'] for r in rows]);p=np.array([r['best']['pressure_MeV_fm3'] for r in rows]);mu=np.array([r['best']['mu_B_MeV'] for r in rows])
    else:
        path=root/f'gppva_tf/final/{name}/selected_states.json.gz'
        rows=json.loads(gzip.decompress(path.read_bytes()))
        n=np.array([r['n_b'] for r in rows]);p=np.array([r['pressure'] for r in rows]);mu=np.array([r['mu_B'] for r in rows])
    order=np.argsort(p);p,n,mu=p[order],n[order],mu[order]
    assert np.all(np.diff(p)>0)
    ng=np.geomspace(b.n0,.002,65);bs=[b.state(x) for x in ng];pg=np.array([s['pressure'] for s in bs])
    good=(pg>=p[0])&(pg<=p[-1]);pg=pg[good];ng=ng[good];gm=np.array([s['mu_b_euler'] for s in bs])[good]
    interp=PchipInterpolator(p,mu,extrapolate=False);dn=PchipInterpolator(p,n,extrapolate=False);g=interp(pg)-gm
    samples=[dict(n_bsk=float(nl),n_inner_interpolated=float(nh),P=float(pp),delta_mu=float(dm)) for nl,nh,pp,dm in zip(ng,dn(pg),pg,g)]
    results[name]={'source_path':str(path.relative_to(ROOT)),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'common_pressure_interval':[float(pg[0]),float(pg[-1])],'min_delta_mu':float(min(g)),'max_delta_mu':float(max(g)),
      'samples':samples,'status':'INTERPOLATION_DIAGNOSTIC_NOT_CROSSING_PROOF'}
    print(name,'gap',min(g),max(g),'n_inner_at_upper',dn(pg[-1]),flush=True)
(OUT/'initial_original_gap_scan.json').write_text(json.dumps({'representation':'unchanged source inner + independent outer anchored BSk22',
    'scope':'Actual committed inner samples interpolated solely to locate candidate crossings; no original-grid absence proof','models':results},indent=2)+'\n')
