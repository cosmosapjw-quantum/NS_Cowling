"""Native SPG numerical backend checks; Python is an independent test oracle."""
import json
import os
from pathlib import Path
import subprocess
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
MODELS = [f"SPG_M{i}" for i in range(1, 6)]

@pytest.fixture(scope="module")
def native_spg(tmp_path_factory):
    required = ["compose_eos_gdfm_native.f90", "compose_eos_spg_native.f90"]
    for name in required:
        assert (ROOT / "fortran_eos/src" / name).is_file(), "native SPG numerical backend is not implemented"
    build = tmp_path_factory.mktemp("native_spg")
    script = ROOT / "fortran_eos/tests/build_native_spg_test.sh"
    run = subprocess.run(["bash", str(script), str(build)], cwd=ROOT, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    return build / "native_spg_probe"

def probe(exe, *args):
    p = subprocess.run([str(exe), *map(str, args)], capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    lines = p.stdout.strip().splitlines()
    assert lines[-1] == "ok", p.stdout
    return np.fromstring(lines[0], sep=" ")

@pytest.mark.parametrize("model", MODELS)
def test_native_gdfm_preserves_author_bulk_and_rearrangement(native_spg, model):
    from compose_eos.gdfm_uniform import baryons
    for n, y in [(1e-6, 0.0), (.04, .08), (.16, .5), (.3, .12)]:
        actual = probe(native_spg, "bulk", model, n, y)
        expected = baryons(model, n, y)
        np.testing.assert_allclose(actual[:4], [expected[k] for k in
            ("epsilon_baryons", "mu_n", "mu_p", "pressure_baryons")], rtol=3e-10, atol=2e-9)
        # This also tests the zero-proton gas branch and nonlinear scalar root.

@pytest.mark.parametrize("model", MODELS)
def test_native_spg_functional_gradient_and_conserved_composition(native_spg, model):
    from compose_eos.candidate_gdfm_cld import default_surface_parameters, reduced_two_phase_cell
    s = default_surface_parameters(model)
    for nb, z in [(.001, [.6, .15, .2, .0]), (.03, [.2, .14, .06, .1]), (.05, [.12,.13,.04,.3])]:
        a = probe(native_spg, "cell", model, s.sigma0_MeV_fm2, s.b, nb, *z)
        e, g, d = reduced_two_phase_cell(model, nb, z, surface=s)
        np.testing.assert_allclose(a[:5], [e, *g], rtol=2e-8, atol=2e-7)
        assert abs(a[5] + a[6] + a[9] - 1) < 1e-13
        assert abs(a[7] - a[6] - a[9]*z[2]) < 1e-13

@pytest.mark.parametrize("model", MODELS)
def test_native_spg_reoptimizes_cells_without_python_numeric_calls(native_spg, model):
    from compose_eos.candidate_gdfm_cld import default_surface_parameters, solve_two_phase
    s = default_surface_parameters(model)
    for nb in (.001, .02, .05):
        a = probe(native_spg, "solve", model, s.sigma0_MeV_fm2, s.b, nb)
        old, _ = solve_two_phase(model, nb, surface=s)
        assert a[4] < 1e-5
        np.testing.assert_allclose(a[:3], [old["energy"], old["pressure_MeV_fm3"], old["mu_B_MeV"]], rtol=2e-6, atol=2e-7)
        assert abs(a[5]+a[6]+a[9]-1) < 1e-13

@pytest.mark.parametrize("model", MODELS)
def test_native_spg_source_crossing_and_reoptimized_first_law(native_spg, model):
    from compose_eos.candidate_gdfm_cld import default_surface_parameters
    s = default_surface_parameters(model)
    cross = probe(native_spg, "transition", model, s.sigma0_MeV_fm2, s.b)
    assert .04 < cross[0] < .12
    assert abs(cross[1]) < 3e-10
    assert cross[3] < 1e-5
    for n in (.002, .02, .05):
        h = n*2e-5
        center = probe(native_spg, "solve", model, s.sigma0_MeV_fm2, s.b, n)
        below = probe(native_spg, "solve", model, s.sigma0_MeV_fm2, s.b, n-h)
        above = probe(native_spg, "solve", model, s.sigma0_MeV_fm2, s.b, n+h)
        derivative = (above[0]-below[0])/(2*h)
        assert abs(derivative-center[2]) < 3e-5
        assert abs(n*derivative-center[0]-center[1]) < 2e-6

@pytest.mark.parametrize("model", MODELS)
def test_native_spg_complete_grid_from_drip_to_source_crossing(native_spg, model):
    from compose_eos.candidate_gdfm_cld import default_surface_parameters
    s = default_surface_parameters(model)
    p = subprocess.run([str(native_spg), "grid", model, str(s.sigma0_MeV_fm2),str(s.b),"17"],capture_output=True,text=True)
    assert p.returncode == 0, p.stdout+p.stderr
    lines = p.stdout.strip().splitlines()
    assert lines[-1] == "ok"
    rows = np.loadtxt(lines[:-1])
    assert rows.shape == (17,18)
    assert np.all(np.diff(rows[:,0]) > 0)
    assert np.all(rows[:,2] > 0)
    assert np.max(np.abs(rows[:,5]+rows[:,6]+rows[:,9]-1)) < 1e-12
    assert np.max(np.abs(rows[:,6]+rows[:,9]*rows[:,11]/rows[:,10]-rows[:,7]-rows[:,8])) < 1e-12
    assert np.max(rows[:,12]) < 1e-5

@pytest.mark.parametrize('variant',['geometric_area__nine','geometric_area__eight_measured'])
def test_source_crossing_tracks_cluster_branch_past_uniform_limit_trials(native_spg,variant):
    fit=json.loads((ROOT/'data/research_completion_20261116/theory_review/surface_fit_variants_321.json').read_text())['models']['SPG_M5']['variants'][variant]
    crossing=probe(native_spg,'transition','SPG_M5',fit['sigma0_MeV_fm2'],fit['b'])
    assert .07<crossing[0]<.08
    assert abs(crossing[1])<3e-10
    assert crossing[3]<1e-5
    # General cell selection still refuses a lower, nonstationary trial. The
    # source endpoint routine is explicitly a stationary branch comparison.
    if variant=='geometric_area__nine':
        p=subprocess.run([str(native_spg),'solve','SPG_M5',str(fit['sigma0_MeV_fm2']),str(fit['b']),'.07345312500000005'],capture_output=True,text=True)
        assert p.returncode!=0
        # GCC13 and GCC15 can stop at different rejected trial points. Both
        # must refuse the point; validate the reported reason numerically.
        lines=p.stdout.strip().splitlines()
        values=np.fromstring(lines[0],sep=' ')
        assert values.shape==(15,)
        assert np.isfinite(values).all()
        if lines[-1]=='STATIONARITY_RESIDUAL':
            assert values[4]>1e-5
        else:
            assert lines[-1]=='SEARCH_BOUND_ACTIVE'
            assert values[4]<=1e-5
            z=values[10:13]
            lo=np.array([1e-8,.07345312500000005*(1+1e-6),.0001])
            hi=np.array([1.,.25,.499999])
            assert np.any(z-lo<1e-6*(hi-lo)) or np.any(hi[1:]-z[1:]<1e-6*(hi[1:]-lo[1:]))
