"""Native lower-crust physics: source oracles are test-only Python."""
from pathlib import Path
import subprocess
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT/"fortran_eos/src/compose_eos_lower_crust_native.f90"

def test_native_lower_crust_implementation_exists():
    assert MODULE.is_file(), "native lower/outer-crust physics backend is missing"

@pytest.fixture(scope="module")
def executable(tmp_path_factory):
    assert MODULE.is_file(), "native lower/outer-crust physics backend is missing"
    out=tmp_path_factory.mktemp("native_lower")
    exe=out/"lower"
    subprocess.run(["gfortran","-std=f2018","-O2","-ffree-line-length-none","-fcheck=all",
        "-J"+str(out),"-I"+str(out),str(ROOT/"fortran_eos/src/compose_eos_types.f90"),
        str(MODULE),str(ROOT/"fortran_eos/tests/native_lower_crust_driver.f90"),
        "-o",str(exe)],check=True,capture_output=True,text=True)
    return exe

def run(executable,*args):
    p=subprocess.run([str(executable),*map(str,args)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+p.stderr
    return np.fromstring(p.stdout.strip().splitlines()[-1],sep=" ")

def test_native_outer_fixed_pressure_global_minimum(executable):
    from compose_eos.outer_crust import OuterCrust
    solver=OuterCrust()
    for p in (1e-10,1e-8,1e-6,1e-5,1e-4,4e-4):
        native=run(executable,"outer",p)
        row=solver.state(p)
        expected=[row[k] for k in ("n_b","epsilon","pressure","g","mu_e","Y_n","Y_p","Y_e","Y_mu","Y_h","A","Z")]
        np.testing.assert_allclose(native[:12],expected,rtol=2e-10,atol=2e-11)
        assert native[12]<2e-10

def test_native_neutron_drip(executable):
    from compose_eos.outer_crust import OuterCrust
    s=OuterCrust().neutron_drip()
    native=run(executable,"drip")
    assert tuple(native[10:12])==(128,38)
    np.testing.assert_allclose(native[:5],[s[k] for k in ("n_b","epsilon","pressure","g","mu_e")],rtol=2e-10,atol=2e-10)

@pytest.mark.parametrize("mode",["source_fit","energy_derivative"])
def test_native_bsk_analytic_and_composition(executable,mode):
    from compose_eos.bsk22_analytic import state
    for n in (2.69403e-4,.0003,.001,.002):
        s=state(n,thermodynamic_mode=mode)
        native=run(executable,"bsk",n,mode)
        expected=[s[k] for k in ("n_b","epsilon","pressure","mu_b_euler","mu_e","Y_n_free","Y_p_free","Y_e","Y_mu","Y_h","A_cluster","Z_cluster")]
        np.testing.assert_allclose(native[:12],expected,rtol=3e-13,atol=3e-12)
        assert abs(native[5]+native[6]+native[9]-1)<1e-14

def test_native_anchored_firstlaw_extension(executable):
    from compose_eos.bsk22_extended_representation import ExtendedAnchoredBSk22
    rep=ExtendedAnchoredBSk22(outer_density=.0002694320064010902,
        outer_pressure=.0004992059167689819,outer_mu_b=939.5654205200001,
        outer_epsilon=.2526497904790187)
    for n in (.0004,.001,.002,.003,.006,.01):
        s=rep.state(n);native=run(executable,"anchor",n)
        np.testing.assert_allclose(native[:12],[s[k] for k in ("n_b","epsilon","pressure","mu_b_euler","mu_e","Y_n_free","Y_p_free","Y_e","Y_mu","Y_h","A_cluster","Z_cluster")],rtol=5e-12,atol=2e-10)
        assert abs(native[14]-s["source_composition_mu_n_residual_MeV"])<2e-10
    n=.006;h=n*1e-5
    a=run(executable,"anchor",n-h);b=run(executable,"anchor",n+h);m=run(executable,"anchor",n)
    assert abs((b[1]-a[1])/(2*h)-m[3])<2e-7
    assert abs(m[14])>1e-4

@pytest.mark.parametrize("args",[("bsk",.003,"source_fit"),("anchor",.010001),("outer",0),("outer",1e-15)])
def test_native_lower_refuses_unsupported(executable,args):
    p=subprocess.run([str(executable),*map(str,args)],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode!=0

def test_native_bsk_eos_carrier(executable):
    value=run(executable,"carrier")
    assert tuple(value)==(17,1,1,1)

def test_native_outer_sequence_resolves_phase_limits(executable):
    value=run(executable,"sequence",33)
    assert value[0]>=33
    assert value[1]>=16
    assert value[2]<2e-8
    assert value[3]>0


@pytest.mark.parametrize("ne,z,a,excess,radius",[
    (5e-10,26,56,-60.607082,3.756),(4e-6,32,82,-65.415,4.12),
    (8e-5,38,128,63.45,4.72)])
def test_native_pressure_from_independent_source_energy(executable,ne,z,a,excess,radius):
    import mpmath as m
    from test_outer_crust import source_energy
    # Differentiate source energy at 60 digits; do not replay Python analytic pressure.
    with m.workdps(60):
        n=m.mpf(str(ne));fun=lambda q:source_energy(q,z,a,excess,radius)
        energy=fun(n);pressure=n*m.diff(fun,n)-energy
        native=run(executable,"isotope",ne,z,a,excess,radius)
        assert native[1]==pytest.approx(float(energy),rel=3e-13,abs=0)
        assert native[2]==pytest.approx(float(pressure),rel=2e-10,abs=0)
        assert native[3]==pytest.approx((native[1]+native[2])/native[0],rel=2e-15)
        assert native[8]==0 and native[9]==1
        assert native[7]==pytest.approx(z/a)

def test_native_source_segments_keep_drip_mismatch_visible(executable):
    value=run(executable,"source_segments",33)
    assert value[0]>=33 and value[1]==17
    assert tuple(value[2:4])==(1,0)  # overlap; independent source fits are not joined.
    assert value[4]>value[5] and value[5]==2.69403e-4
    assert value[6]<0 and abs(value[7])>1e-9 and abs(value[8])>1e-5
