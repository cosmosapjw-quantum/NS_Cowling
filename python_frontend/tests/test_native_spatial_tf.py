"""Native spatial TF verification; Python is an independent test oracle only."""
from pathlib import Path
import os, subprocess
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[2]
MODELS=("DD2","DDME2","TW","FSU2","FSU2R","FSU2H","NL3wrL55","TM1e")

def test_native_tf_sources_present():
    # Would fail if the public native solver or physical kernel is removed.
    for name in ("compose_eos_tf_models.f90","compose_eos_spatial_tf.f90"):
        assert (ROOT/"fortran_eos/src"/name).is_file(), "native spatial TF backend missing: "+name

@pytest.fixture(scope="session")
def native_tf_driver(tmp_path_factory):
    paths=[ROOT/"fortran_eos/src"/n for n in
       ("compose_eos_types.f90","compose_eos_optimize_native.f90",
        "compose_eos_tf_models.f90","compose_eos_spatial_tf.f90")]
    assert all(p.is_file() for p in paths), "native spatial TF backend has not been implemented"
    b=tmp_path_factory.mktemp("native_tf")
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    vendor=ROOT/"fortran_eos/vendor/slsqp/slsqp_optmz.f"
    subprocess.run(["gfortran","-std=legacy","-O2","-ffixed-line-length-none","-c",str(vendor),"-o",str(b/"slsqp.o")],check=True,env=env)
    exe=b/"tf_driver"
    subprocess.run(["gfortran","-std=f2018","-O2","-ffree-line-length-none","-fcheck=all",
       "-J"+str(b),"-I"+str(b),*map(str,paths),str(ROOT/"fortran_eos/tests/native_tf_driver.f90"),
       str(b/"slsqp.o"),"-llapack","-lblas","-o",str(exe)],check=True,env=env)
    return exe

def invoke(exe,args,arrays=()):
    text="\n".join(" ".join(format(float(v),".17g") for v in a) for a in arrays)+"\n"
    r=subprocess.run([str(exe),*map(str,args)],input=text,text=True,capture_output=True,
                     env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1"))
    assert r.returncode==0,r.stdout+r.stderr
    lines=r.stdout.splitlines();assert lines[0]=="ok",r.stdout+r.stderr
    vals=np.array([float(x) for x in lines[1].split()])
    profile=np.array([[float(x) for x in s.split()] for s in lines[2:]])
    return vals,profile

@pytest.mark.parametrize("dim",(1,2,3))
def test_native_electron_equilibrium(native_tf_driver,dim):
    from compose_eos.spatial_tf import RadialGrid,fermions,ELECTRON_MASS,E2
    from compose_eos.spatial_tf_equilibrium import screened_electrons
    g=RadialGrid(48,dim);p=.002*(1+.7*np.cos(np.pi*g.x));r=25.
    vals,profile=invoke(native_tf_driver,["electron",48,dim,r],[p])
    ne,v=profile.T;ref=screened_electrons(g,p,r)
    assert vals[1]<2e-10
    assert abs(g.mean(ne-p))<2e-12
    assert abs(g.mean(v))<2e-11
    assert np.max(abs(fermions(ne,ELECTRON_MASS)["ef"]-v-vals[0]))<2e-10
    assert np.max(abs(g.minus_laplacian@v/r**2-4*np.pi*E2*(p-ne)))<2e-10
    np.testing.assert_allclose(ne,ref[0],rtol=2e-9,atol=2e-13)

@pytest.mark.parametrize("model",MODELS)
def test_native_uniform_reduction(native_tf_driver,model):
    from compose_eos import gppva_rmf as rmf
    n=.03;u=rmf.beta_state(model,n,include_muons=False);p=np.full(12,n*u["Y_p"])
    vals,_=invoke(native_tf_driver,["evaluate",model,12,3,n,25.],[np.full(12,n)-p,p,p])
    assert vals[0]==pytest.approx(u["epsilon"],rel=2e-10)
    assert vals[2]==pytest.approx(u["mu_n"],rel=3e-9)
    assert abs(vals[7])<1e-10

@pytest.mark.parametrize("model",("DD2","FSU2R","NL3wrL55"))
def test_native_profile_functional_and_gradient(native_tf_driver,model):
    from compose_eos.spatial_tf import RadialGrid,initial_profile,CellFunctional
    g=RadialGrid(16,3);nn,p,ne=initial_profile(g,.03);nn+=1e-5;p+=1e-5
    ref=CellFunctional(model,g,.03).evaluate(nn,p,ne,35.)
    vals,profile=invoke(native_tf_driver,["evaluate",model,16,3,.03,35.],[nn,p,ne])
    assert vals[0]==pytest.approx(ref["epsilon"],rel=2e-10)
    np.testing.assert_allclose(profile[:,3],ref["mu_n"],rtol=3e-9,atol=2e-7)
    np.testing.assert_allclose(profile[:,4],ref["mu_p"],rtol=3e-9,atol=2e-7)
    assert vals[7]==pytest.approx(ref["d_epsilon_d_log_radius"],rel=2e-8,abs=1e-9)

@pytest.mark.parametrize("model",("DD2","FSU2R"))
def test_native_independent_energy_derivatives(native_tf_driver,model):
    from compose_eos.spatial_tf import RadialGrid,initial_profile
    g=RadialGrid(16,3);nn,p,ne=initial_profile(g,.03)
    nn=nn+2e-5;p=p+2e-5;ne=ne+2e-5
    vals,profile=invoke(native_tf_driver,["evaluate",model,16,3,.03,35.],[nn,p,ne])
    direction=np.sin(np.arange(16)+.3)*1e-3
    step=2e-5
    for species in range(3):
        a=[nn.copy(),p.copy(),ne.copy()];b=[x.copy() for x in a]
        a[species]+=step*direction;b[species]-=step*direction
        plus=invoke(native_tf_driver,["evaluate",model,16,3,.03,35.],a)[0][0]
        minus=invoke(native_tf_driver,["evaluate",model,16,3,.03,35.],b)[0][0]
        expected=g.weights@(profile[:,3+species]*direction)
        assert (plus-minus)/(2*step)==pytest.approx(expected,rel=3e-6,abs=2e-8)
    h=1e-4
    plus=invoke(native_tf_driver,["evaluate",model,16,3,.03,35*np.exp(h)],[nn,p,ne])[0][0]
    minus=invoke(native_tf_driver,["evaluate",model,16,3,.03,35*np.exp(-h)],[nn,p,ne])[0][0]
    assert (plus-minus)/(2*h)==pytest.approx(vals[7],rel=3e-6,abs=2e-8)

def test_native_preserved_n384_profile(native_tf_driver):
    import gzip,json
    path=ROOT/"data/research_finish_20261116/vector_band_prototype/FSU2R_N384_CELL.json.gz"
    state=json.loads(gzip.decompress(path.read_bytes()));p=state["profile"]
    vals,profile=invoke(native_tf_driver,["evaluate","FSU2R",384,state["dimension"],state["n_b"],state["radius_fm"]],
                        [p["nn"],p["np"],p["ne"]])
    assert vals[0]==pytest.approx(state["epsilon"],rel=2e-11,abs=2e-10)
    assert max(vals[4:6])<1e-8
    for col,key in ((3,"mu_n"),(4,"mu_p"),(5,"mu_e")):
        np.testing.assert_allclose(profile[:,col],p[key],rtol=2e-9,atol=2e-7)

@pytest.mark.parametrize("model",MODELS)
def test_native_finite_family_actual_minimization(native_tf_driver,model):
    # Actual constrained minimization, not merely evaluation of reference rows.
    vals,_=invoke(native_tf_driver,["family",model,16,3,.03,30,"droplet"])
    assert vals[15]==1
    assert vals[8]<.02 and vals[9]<2e-4
    assert max(vals[4:6])<1e-8 and vals[6]<2e-10
    assert vals[1]>0

def test_native_failed_initialization_is_retained_before_recovery(native_tf_driver):
    r=subprocess.run([str(native_tf_driver),"family","FSU2H","16","3",".075","30","droplet"],
                     capture_output=True,text=True,env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1"))
    assert r.returncode==0,r.stdout+r.stderr
    assert "EVALUATION_FAILED vector field Newton line search failed" in r.stderr
    assert "ATTEMPT 5 STATIONARY ok" in r.stderr
    # The accepted clustered solution is lower than the native uniform seed.
    records=[line.split() for line in r.stderr.splitlines() if line.startswith("ATTEMPT ")]
    good=[line for line in records if line[2] in ("STATIONARY","UNIFORM_STATIONARY")]
    assert float(good[-1][4]) < float(good[0][4])-1e-3
