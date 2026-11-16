"""Source-equation, fixed-pressure-minimum and discrete-composition tests."""
import math
import numpy as np
import pytest
from compose_eos import outer_crust as oc

def source_energy(ne,z,a,excess,radius):
    # Independent arbitrary-precision source energy; pressure is differentiated,
    # not copied from the analytic implementation.
    import mpmath as m
    n,z,a,excess,radius=map(m.mpf,map(str,(ne,z,a,excess,radius)))
    h,me,alpha=map(m.mpf,map(str,(oc.HBARC,oc.ME,oc.ALPHA)))
    e2=alpha*h; x=h*(3*m.pi**2*n)**(m.mpf(1)/3)/me
    atom=a*m.mpf(str(oc.U))+excess+m.mpf("1.44381e-5")*z**m.mpf("2.39")+m.mpf("1.55468e-12")*z**m.mpf("5.35")
    efree=me**4/(8*m.pi**2*h**3)*(x*(1+2*x*x)*m.sqrt(1+x*x)-m.asinh(x))
    ex=-alpha*me**4*x**4/(4*m.pi**3*h**3)*(1-m.mpf("1.5")*(m.sqrt(1+x*x)/x-m.asinh(x)/x**2)**2)
    corr=alpha**2*me**4*x**4/(12*m.pi**4*h**3)*(-m.mpf("12.51")+m.log10(x))
    r=(3*z/(4*m.pi*n))**(m.mpf(1)/3)
    hw=h*m.sqrt(4*m.pi*z*e2*n/atom)
    plasma=z*z*e2/(r*hw)
    s=1/(1+m.mpf(".01")*m.log(z)**m.mpf("1.5")+m.mpf(".097")/z**2)
    b1=1-m.mpf("1.1866")/z**m.mpf(".267")+m.mpf(".27")/z
    b2=1+m.mpf("2.25")/z**(m.mpf(1)/3)*(1+m.mpf(".684")*z**5+m.mpf(".222")*z**6)/(1+m.mpf(".222")*z**6)
    b3=m.mpf("41.5")/(1+m.log(z));b4=m.mpf(".395")*m.log(z)+m.mpf(".347")/z**m.mpf("1.5")
    f=m.mpf(54)/175*(12/m.pi)**(m.mpf(1)/3)*alpha*z**(m.mpf(2)/3)*b1*m.sqrt(1+b2/x**2)
    aa=(b3+m.mpf("17.9")*x*x)/(1+b4*x*x)
    screen=-f*(4*m.pi/3)**(m.mpf(1)/3)*e2*z**(m.mpf(2)/3)*n**(m.mpf(4)/3)*(1+aa*(m.mpf(".205")/plasma)**s)
    lattice=-m.mpf(str(oc.BCC))*z*z*e2/r*n/z*(1-5*radius**2/(9*r*r))
    zp=n/z*m.mpf("1.5")*m.mpf(str(oc.U1))*hw
    return n/z*(atom-z*me)+efree+ex+corr+screen+lattice+zp

@pytest.mark.parametrize("ne,z,a,excess,radius",[
    (5e-10,26,56,-60.607082,3.756),(4e-6,32,82,-65.415,4.12),
    (8e-5,38,128,63.45,4.72)])
def test_pressure_from_independent_source_energy(ne,z,a,excess,radius):
    import mpmath as m
    with m.workdps(60):
        n=m.mpf(str(ne))
        # Keep mpmath argument exact through the source function.
        fun=lambda q:source_energy(q,z,a,excess,radius)
        en=fun(n); pressure=n*m.diff(fun,n)-en
        actual=oc.isotope_state(ne,z,a,excess,radius)
        assert actual["epsilon"]==pytest.approx(float(en),rel=3e-13,abs=0)
        assert actual["pressure"]==pytest.approx(float(pressure),rel=2e-10,abs=0)

def test_mass_input_priority_identity():
    c=oc.masses()
    assert len(c.z)==8481
    i=int(np.where((c.z==26)&(c.a==56))[0][0])
    assert c.excess[i]==pytest.approx(-60.607082)
    assert c.source[i]=="AME2016_measured"
    i=int(np.where((c.z==29)&(c.a==79))[0][0])
    assert c.excess[i]==-42.408 and c.source[i]=="Welker2017_TableI"

def test_global_minimum_is_at_common_pressure():
    model=oc.OuterCrust()
    pressure=1e-5
    ne,ee,pp,mu,gg=model.candidates(pressure)
    assert np.max(abs(pp/pressure-1))<2e-11
    row=model.state(pressure)
    assert row["g"]==pytest.approx(np.min(gg),abs=2e-11)
    assert (row["A"],row["Z"])==(82,32)

def test_neutron_drip_and_composition():
    model=oc.OuterCrust(); row=model.neutron_drip()
    assert row["g"]==pytest.approx(oc.MN,abs=2e-9)
    assert (row["A"],row["Z"])==(128,38)
    assert 2.5e-4<row["n_b"]<2.8e-4
    with pytest.raises(ValueError,match="above neutron drip"):
        model.state(row["pressure"]*1.01)
    data=oc.eos_data([row])
    assert data.Yh[0]==1 and data.Yn[0]==0 and data.Yp[0]==0
    assert data.Ye[0]==data.Zav[0]/data.Aav[0]
    assert row["cluster_number_fraction"]==1/128

@pytest.mark.parametrize("p",[float("nan"),float("inf"),-1,0,1e-15,1.0])
def test_pressure_refusal(p):
    with pytest.raises(ValueError):oc.OuterCrust().state(p)

def test_missing_radius_is_not_silently_accepted():
    c=oc.MassCatalogue(np.array([26]),np.array([56]),np.array([-60.607082]),np.array([0.]),("test",))
    with pytest.raises(RuntimeError,match="radius"):oc.OuterCrust(c).state(1e-8)


def test_thin_intermediate_layer_is_recursively_resolved():
    c=oc.masses()
    selected=((c.z==26)&(c.a==56))|((c.z==38)&np.isin(c.a,[122,124,126,128]))
    indexes=np.flatnonzero(selected)
    sub=oc.MassCatalogue(c.z[selected],c.a[selected],c.excess[selected],
                        c.radius[selected],tuple(c.source[i] for i in indexes))
    model=oc.OuterCrust(sub)
    rows,transitions,drip=model.sequence(16)
    assert model.intermediate_layer_refinements>=1
    assert (38,124) in {(r["Z"],r["A"]) for r in rows}
    assert all(abs(t["g_difference_mev"])<2e-9 for t in transitions)
    assert all(t["after"]["n_b"]>t["before"]["n_b"] for t in transitions)
    assert all(t["after"]["pressure"]==t["before"]["pressure"] for t in transitions)
