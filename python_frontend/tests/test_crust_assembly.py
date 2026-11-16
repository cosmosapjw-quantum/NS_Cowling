import numpy as np
import pytest
from dataclasses import replace
from compose_eos.eos_data import EOSData
from compose_eos.crust_assembly import rebase_header,assemble_crust,certify_interface,select_rows

def branch(n,p,mu,mn=939.,phase=1,sparse=False):
    n=np.array(n);p=np.array(p);mu=np.array(mu);z=np.zeros(len(n))
    e=n*mu-p
    d=EOSData(nb=n,mn=mn,mp=mn-1.,Q1=p/n,Q2=z,Q3=mu/mn-1,
      Q4=np.full(len(n),-20/mn),Q5=z,Q6=e/n/mn-1,Q7=e/n/mn-1,
      q_extra=z,Iphase=np.full(len(n),phase),Yn=z,Yp=z,Ye=z+.1,Ymu=z,
      Aav=z+100,Zav=z+10,Yh=z+1,composition_available=True,
      composition_file_present=True,beta_equilibrated=True)
    if sparse:
        d.Yn_present[:]=False;d.Yp_present[:]=False;d.Yh_present[:]=False
        d.Abar_present[:]=False;d.Zbar_present[:]=False
        d.C3_applicable[:]=False;d.C4_applicable[:]=False;d.row_status[:]=1
    return d

def test_header_conversion_preserves_dimensionful_values_and_masks():
    a=branch([.001,.002],[.001,.002],[940,941],sparse=True)
    b=rebase_header(a,mn=939.56542052,mp=938.27208816)
    for field in ('nb','p','eps','mub','muq','mul','f','F'):
        np.testing.assert_allclose(getattr(a,field),getattr(b,field),rtol=3e-15,atol=1e-12)
    for field in ('Ye','Yn_present','Yh_present','Iphase','row_status'):
        np.testing.assert_array_equal(getattr(a,field),getattr(b,field))
    b.Ye[0]=.2
    assert a.Ye[0]==.1

def test_assembly_retains_two_sided_density_jump_and_masks():
    a=branch([.0001,.001],[.0001,.002],[939,940],phase=2)
    b=branch([.0012,.002],[.002,.003],[940,941],mn=938.9,phase=21,sparse=True)
    out,cert=assemble_crust([a,b],labels=['analytic_lower','spatial_inner'],
                            mn=939.5654,mp=938.272)
    assert out.n==4
    np.testing.assert_array_equal(out.Iphase,[2,2,21,21])
    np.testing.assert_array_equal(out.Yn_present,[True,True,False,False])
    assert cert['interfaces'][0]['interpolation_across_density_gap']=='FORBIDDEN'
    assert cert['interfaces'][0]['left_row_zero_based']==1

@pytest.mark.parametrize('failure',['mu','pressure','euler','density','negative_pressure'])
def test_unphysical_or_unmatched_interface_refused(failure):
    a=branch([.0001,.001],[.0001,.002],[939,940])
    b=branch([.0012,.002],[.002,.003],[940,941])
    if failure=='mu':
        b.Q3[0]+=.01/b.mn;b.Q6[0]+=.01/b.mn;b.Q7[0]+=.01/b.mn
    if failure=='pressure':
        b.Q1[0]*=1.01;b.Q6[0]=b.Q7[0]=(b.mub[0]-b.Q1[0])/b.mn-1
    if failure=='euler':b.Q7[0]+=.01/b.mn
    if failure=='density':b.nb[0]=a.nb[-1]
    if failure=='negative_pressure':
        a.Q1[-1]=-1.;b.Q1[0]=-1.
    with pytest.raises(ValueError):certify_interface(a,b)

def test_slice_copies_masks_and_refuses_reordering():
    a=branch([.001,.002,.003],[.001,.002,.003],[940,941,942],sparse=True)
    b=select_rows(a,[0,2])
    np.testing.assert_array_equal(b.nb,[.001,.003])
    assert not b.Yn_present.any()
    with pytest.raises(ValueError):select_rows(a,[2,0])

def test_header_alone_does_not_erase_real_model_offset():
    a=branch([.0001,.001],[.0001,.002],[939,940],mn=938.9)
    b=branch([.0012,.002],[.002,.003],[939.4,941],mn=939.5654)
    a=rebase_header(a,mn=b.mn,mp=b.mp)
    assert b.mub[0]-a.mub[-1]==pytest.approx(-.6)
    with pytest.raises(ValueError,match='chemical potential'):certify_interface(a,b)

@pytest.mark.parametrize('mass',[0.,-1.,float('nan'),float('inf')])
def test_bad_header_mass_refused(mass):
    a=branch([.001],[.001],[940])
    with pytest.raises(ValueError):rebase_header(a,mn=mass,mp=938.)

@pytest.mark.parametrize("parameter",["pressure_rtol","mu_atol_MeV","euler_atol_MeV"])
@pytest.mark.parametrize("value",[float("nan"),float("inf"),0.,-1.])
def test_nonfinite_or_nonpositive_interface_tolerances_refused(parameter,value):
    a=branch([.001],[.002],[940])
    b=branch([.0012],[.002],[940])
    with pytest.raises(ValueError):certify_interface(a,b,**{parameter:value})
