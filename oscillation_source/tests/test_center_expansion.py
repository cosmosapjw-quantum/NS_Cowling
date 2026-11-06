"""Physical central series: selectable improvement, not a looser boundary test."""
from pathlib import Path
import pytest
from tools.precision_study import run_one, StudyError

@pytest.fixture(scope='module')
def central_runs(tmp_path_factory):
    root=tmp_path_factory.mktemp('central-series')
    common=dict(alpha=0.0,target_mass_kg=0.,pc_pa=1.5e34,rtol=1e-12,mode_rtol=1e-12,step_cap=.04)
    runs={}
    for tag,radius,order in [('reference',1.,0),('leading',100.,0),('second',100.,2)]:
        runs[tag]=run_one(root/tag,{**common,'r0_m':radius,'center_order':order},workspace='calculation')
    return runs

def test_center_series_reduces_actual_start_radius_error(central_runs):
    ref,lo,hi=(central_runs[k] for k in ['reference','leading','second'])
    assert hi['center_order']==2 and lo['center_order']==0
    for key in ['mass_kg','radius_m','f_hz','p1_hz']:
        old=abs(lo[key]-ref[key]);new=abs(hi[key]-ref[key])
        assert old>0
        assert new<old/30, (key,old,new)
    for key in ['f_residual','p1_residual']:
        assert abs(hi[key])<1e-8

def test_invalid_center_order_is_not_silently_accepted(tmp_path):
    with pytest.raises(StudyError):
        run_one(tmp_path/'bad',dict(alpha=0.,center_order=1),workspace='calculation')

def test_series_order_is_not_a_refinement_axis():
    from tools.precision_study import compare_sequence
    rows=[]
    for order in [0,0,2]:
        rows.append(dict(delta_f_hz=.1,delta_p1_hz=.2,identity_consistent=True,
            core_function_unchanged=True,reference={'center_order':order},variant={'center_order':order}))
    assert compare_sequence(rows)['status']=='controls_failed'

def test_public_model_validates_center_order():
    from tools.stellar_methods import Model
    assert Model(r0_m=100,center_order=2).parameters()['center_order']==2
    with pytest.raises(ValueError):Model(center_order=1).parameters()
