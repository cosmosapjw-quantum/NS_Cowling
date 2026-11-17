"""A geometry comparison must preserve missing branches and numerical failures."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from refine_gppva_geometry import family_select,transition_evidence,FAMILIES

def row(family,energy,valid=True):
    return dict(requested_family=family,dimension=FAMILIES[family][0],seed=FAMILIES[family][1],
                n_b=.06,baryon_constraint_error=0.,charge_constraint_error=0.,
                scalar_pde_residual=0.,vector_pde_residual=0.,electron_pde_residual=0.,
                epsilon=energy,energy_per_baryon=energy,
                numerical_status='STATIONARY' if valid else 'REQUIRES_REFINEMENT',
                muon_threshold_crossed=False)

def test_missing_geometry_family_cannot_be_hidden_by_success():
    best,status,missing=family_select([row('droplet',1.),dict(requested_family='rod',execution='FAILED')],
                                      ('droplet','rod'))
    assert best is None and status=='HOLD_UNRESOLVED_GEOMETRY' and missing==['rod']

def test_recovered_family_does_not_require_every_arbitrary_start():
    trials=[row('droplet',1.),dict(requested_family='rod',execution='FAILED'),row('rod',.9)]
    best,status,missing=family_select(trials,('droplet','rod'))
    assert best is trials[-1] and status=='LOWEST_EVALUATED_STATIONARY_CANDIDATE' and missing==[]

def test_known_lower_unconverged_energy_is_retained():
    trials=[row('droplet',1.),row('rod',.9),row('droplet',.8,False)]
    best,status,_=family_select(trials,('droplet','rod'))
    assert best is None and status=='LOWER_UNCONVERGED_ATTEMPT'

def test_grid_energy_ordering_is_not_automatically_a_continuum_sign():
    # Real profile shape metadata supplies physical geometry independently of
    # requested_family; no source-table values are involved.
    def physical(g,e):
        d=3 if g=='droplet' else 2
        return dict(row(g,e),baryon_density_contrast=1.,
                    profile={'nn':[.1,.05],'np':[.02,.01]},pressure=.1,mu_B=950.)
    records=[]
    for n,coarse,fine in [(.06,-.001,.00001),(.065,.001,-.00001)]:
        for mesh,delta in [(192,coarse),(384,fine)]:
            records.append({'n_b':n,'points':mesh,'attempts':[physical('droplet',950.),physical('rod',950.+delta)]})
    result=transition_evidence([{'left':.06,'right':.065,'geometries':['droplet','rod']}],records,[192,384])[0]
    assert len(result['fine_mesh_sign_change_brackets'])==1
    assert not result['fine_mesh_sign_change_brackets'][0]['both_signs_exceed_mesh_change']

def test_collapsed_seed_does_not_create_a_nonexistent_physical_branch():
    from refine_gppva_geometry import branch_minima
    collapsed=dict(row('rod',950.),baryon_density_contrast=0.,
                   profile={'nn':[.05,.05],'np':[.002,.002]},pressure=.1,mu_B=950.)
    assert set(branch_minima([collapsed]))=={'uniform'}

def test_density_envelope_diagnostic_has_correct_pressure_units():
    from summarize_gppva_geometry import derivative_checks
    # Independently choose f(n)=m+C*n^3, so P=n^2*f'=3*C*n^4.
    # Simpson is exact for the quadratic integrand P/n^2 here.
    c=7.;rows=[{'n_b':n,'energy_per_baryon':939.+c*n**3,'pressure':3*c*n**4}
               for n in (.02,.04,.06)]
    report=derivative_checks(rows)
    assert len(report)==1
    assert report[0]['absolute_difference_MeV_per_baryon']<2e-13

def test_weak_cluster_seed_is_nonuniform_but_preserves_global_constraints():
    import numpy as np
    from recover_tf_family import weak_cluster_seed
    n=.095;yp=.02;points=16
    old={'n_b':n,'Y_p':yp,'Y_n':1-yp,'points':points,'dimension':3,'radius_fm':30.,
         'profile':{'x':((np.arange(points)+.5)/points).tolist(),
                    'nn':[n*(1-yp)]*points,'np':[n*yp]*points,'ne':[n*yp]*points}}
    for dimension in (1,2,3):
        row=weak_cluster_seed(old,dimension,.3,15.);p=row['profile'];w=np.array(p['weights'])
        nn=np.array(p['nn']);np_=np.array(p['np']);ne=np.array(p['ne'])
        assert abs(w@(nn+np_)-n)<3e-17
        assert abs(w@(np_-ne))<1e-18
        assert nn[0]+np_[0]>nn[-1]+np_[-1]
        assert min(np_)>0 and row['radius_fm']==15.
    assert old['radius_fm']==30. and len(set(old['profile']['np']))==1

def test_roundoff_size_sign_flip_is_not_resolved_by_identical_mesh_values():
    def physical(g,e):
        return dict(row(g,e),baryon_density_contrast=1.,profile={'nn':[.1,.05],'np':[.02,.01]},pressure=.1,mu_B=950.)
    records=[]
    for n,delta in [(.06,1e-13),(.065,-1e-13)]:
        for mesh in [192,384]:records.append({'n_b':n,'points':mesh,'attempts':[physical('droplet',950.),physical('rod',950.+delta)]})
    result=transition_evidence([{'left':.06,'right':.065,'geometries':['droplet','rod']}],records,[192,384])[0]
    bracket=result['fine_mesh_sign_change_brackets'][0]
    assert bracket['both_signs_exceed_mesh_change']
    assert not bracket['both_endpoints_exceed_roundoff_scale']
    assert all(r['roundoff_near_degeneracy'] for r in result['samples'])

def test_morphology_collapse_does_not_drop_family_continuation():
    from refine_gppva_geometry import choose_profile_seeds
    def seed(n,radius,uniform=False):
        return dict(row('rod',950.),n_b=n,points=16,radius_fm=radius,
                    baryon_density_contrast=0. if uniform else 1.,
                    profile={'nn':[.05,.05] if uniform else [.1,.05],
                             'np':[.002,.002] if uniform else [.02,.01]})
    first,second,continuation=seed(.05,25.),seed(.06,30.),seed(.07,35.,True)
    selected=choose_profile_seeds([first,second,continuation],'rod',[continuation])
    assert len(selected)==3 and selected[-1] is continuation
