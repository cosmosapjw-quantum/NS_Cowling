"""C16 native source-locked density continuation: scan, nodal classification.

TDD: all determinant samples and physical nodal counts in Fortran.
No Python physics root function during frequency scanning; separate
full-mode reconstruction required for scientific f/p1 admission.
"""
import ctypes as ct
import hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pytest
from compose_eos import stellar_cowling_native as native
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.whole_crust import WholeCrustBarotrope

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip"
SOURCE_SHA="299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04"
CASES=[
("SPG_M2","maxwell",.30385128926615257,"rapid",(5090.,5190.),5143.7197,39),
("SPG_M3","css",.192510767556007,"slow",(3690.,3800.),3736.2107,40)
]

@pytest.fixture(scope="module",params=CASES)
def scenario(request):
    model,kind,nc,conv,br,ref,knots=request.param
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
    with zipfile.ZipFile(SOURCE) as outer:
      with zipfile.ZipFile(io.BytesIO(outer.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inside:
        raw=inside.read(f"FRESH_5x2_WHOLE_STARS/{model}/{kind}/REPRESENTATION.json")
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(raw)))
    bg=solve_sharp_background_native(eos,nc,surface_h=3e-8,
      rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=2049)
    opt=CowlingSettings(surface_h=3e-8,rtol=1e-9,atol=1e-12,
                         max_step_h=.002,root_xtol_hz=1e-5)
    solver=NativeSharpCowlingDVODE(bg,opt,samples=4097,
                  coefficient_backend="native_source_cs2",
                  source_restart_relative_jump=.01)
    return model,solver,conv,br,ref,knots

def test_fortran_original_eos_scan_brackets_source_p1(scenario,monkeypatch):
    model,solver,conv,br,ref,expected_knots=scenario
    x=np.ascontiguousarray(np.linspace(br[0],br[1],9))
    def refuse(*args,**kwargs):
       raise AssertionError("Python matching_residual called by Fortran scan")
    monkeypatch.setattr(solver,"matching_residual",refuse)
    pairs,values=solver.scan_source_brackets_native(x,conversion=conv)
    assert any(a<ref<b for a,b in pairs)
    assert len(values)==len(x)
    assert all(v["status"]==0 for v in values)
    assert all(v["native_rhs_evaluations"]>0 for v in values)
    assert all(np.isfinite(v["determinant"]) for v in values)
    assert {v["source_restarts"] for v in values}=={expected_knots}


def test_fortran_original_eos_scan_rejects_out_of_range_without_overwriting(scenario):
    model,solver,conv,br,ref,expected_knots=scenario
    with pytest.raises(ValueError,match="frequency"):
        solver.scan_source_brackets_native([br[1],br[0]],conversion=conv)
    with pytest.raises(ValueError,match="frequency"):
        solver.scan_source_brackets_native([-1.,br[0]],conversion=conv)
    with pytest.raises(ValueError,match="conversion"):
        solver.scan_source_brackets_native(br,conversion="unknown")


def test_native_fortran_phase_node_counts_and_zeros():
    dll=native._native_library()
    phase=np.ascontiguousarray(np.array([0,0,0,0,1,1,1,1,1],np.int32))
    w=np.ascontiguousarray([0.,1.,-1.,2.,0.,1.,-1.,0.,2.],np.float64)
    v=np.ascontiguousarray([0.,1.,1.,-2.,0.,1.,1.,0.,-2.],np.float64)
    wcount,vcount,st=ct.c_int(-1),ct.c_int(-1),ct.c_int(-1)
    dll.cowling_native_phase_nodes(len(w),phase,w,v,1e-9,
        ct.byref(wcount),ct.byref(vcount),ct.byref(st))
    assert st.value==0
    assert (wcount.value,vcount.value)==(4,2)

def test_native_fortran_phase_nodes_fail_closed():
    dll=native._native_library()
    phase=np.ascontiguousarray([0,0,1,0],dtype=np.int32)
    x=np.ascontiguousarray([1.,-1.,1.,-1.],dtype=np.float64)
    a,b,status=ct.c_int(-1),ct.c_int(-1),ct.c_int(-1)
    dll.cowling_native_phase_nodes(4,phase,x,x,1e-9,
          ct.byref(a),ct.byref(b),ct.byref(status))
    assert status.value!=0
    assert (a.value,b.value)==(0,0)


def test_real_mode_counts_fortran_and_original_python_reference(scenario):
    model,sol,conv,br,ref,knots=scenario
    result=sol.solve_mode(br,conversion=conv,profile_points=300,
                          root_backend="fortran_source_exact")
    assert result["node_count_backend"]=="fortran_phasewise"
    assert (result["W_nodes_per_phase"],result["V_nodes_per_phase"])==(1,1)
    assert result["relative_join_mismatch"]<1e-5
    assert abs(result["frequency_hz"]-ref)<.02


def test_native_mode_branch_selector_retains_node_identity_without_ambiguity():
    dll=native._native_library()
    f=np.ascontiguousarray([1550.,1842.,3710.,5143.,6500.],dtype=float)
    nw=np.ascontiguousarray([1,0,1,1,2],dtype=np.int32)
    nv=np.ascontiguousarray([1,0,1,1,2],dtype=np.int32)
    stat=np.ascontiguousarray([0,0,0,0,0],dtype=np.int32)
    selected,status=ct.c_int(-9),ct.c_int(-9)
    dll.cowling_native_select_mode_branch(
       len(f),f,nw,nv,stat,0,1800.,.12,
       ct.byref(selected),ct.byref(status))
    assert status.value==0 and selected.value==1
    dll.cowling_native_select_mode_branch(
       len(f),f,nw,nv,stat,1,5000.,.08,
       ct.byref(selected),ct.byref(status))
    assert status.value==0 and selected.value==3


def test_native_mode_branch_selector_fails_closed_on_ambiguous_candidates():
    dll=native._native_library()
    f=np.ascontiguousarray([1800.,1820.,1850.],dtype=float)
    nw=np.ascontiguousarray([0,0,1],dtype=np.int32)
    nv=np.ascontiguousarray([0,0,1],dtype=np.int32)
    st=np.ascontiguousarray([0,0,0],dtype=np.int32)
    picked,status=ct.c_int(-1),ct.c_int(-1)
    dll.cowling_native_select_mode_branch(3,f,nw,nv,st,0,1800.,.10,
         ct.byref(picked),ct.byref(status))
    assert status.value!=0 and picked.value==-1
    st[0]=6
    dll.cowling_native_select_mode_branch(3,f,nw,nv,st,0,1800.,.10,
         ct.byref(picked),ct.byref(status))
    assert status.value==0 and picked.value==1


def test_native_fortran_phase_shape_overlap_identical_and_global_sign():
    dll=native._native_library()
    phase=np.ascontiguousarray([0,0,0,0,1,1,1,1],dtype=np.int32)
    r=np.ascontiguousarray([0.,.2,.4,.6,.6,.8,1.0,1.2],dtype=float)
    w=np.ascontiguousarray([0.,.7,1.,.4,.4,-.1,-.4,-.2],dtype=float)
    total,inner,outer,status=(ct.c_double(-1),ct.c_double(-1),
                              ct.c_double(-1),ct.c_int(-1))
    dll.cowling_native_mode_shape_overlap(
       len(r),phase,r,w,len(r),phase,r,w,
       ct.byref(total),ct.byref(inner),ct.byref(outer),ct.byref(status))
    assert status.value==0
    assert total.value==pytest.approx(1.,abs=1e-13)
    assert inner.value==pytest.approx(1.,abs=1e-13)
    assert outer.value==pytest.approx(1.,abs=1e-13)
    dll.cowling_native_mode_shape_overlap(
       len(r),phase,r,w,len(r),phase,r,np.ascontiguousarray(-w),
       ct.byref(total),ct.byref(inner),ct.byref(outer),ct.byref(status))
    assert status.value==0
    assert total.value==pytest.approx(1.,abs=1e-13)


def test_native_fortran_phase_overlap_rejects_nonmonotone_radius():
    dll=native._native_library()
    phase=np.ascontiguousarray([0,0,0,1,1,1],dtype=np.int32)
    r=np.ascontiguousarray([0.,.4,.6,.6,.7,.9],dtype=float)
    bad=np.ascontiguousarray([0.,.4,.2,.6,.7,.9],dtype=float)
    w=np.ascontiguousarray([0.,1.,1.,.5,-.2,-1.],dtype=float)
    a,b,c,status=ct.c_double(-1),ct.c_double(-1),ct.c_double(-1),ct.c_int(-1)
    dll.cowling_native_mode_shape_overlap(len(r),phase,r,w,
         len(r),phase,bad,w,ct.byref(a),ct.byref(b),ct.byref(c),
         ct.byref(status))
    assert status.value!=0
    assert a.value==b.value==c.value==0.


def test_native_fortran_original_eos_scan_independent_reference_residuals(scenario):
    model,solver,conv,br,ref,expected_knots=scenario
    points=np.ascontiguousarray([br[0],br[1]],dtype=np.float64)
    _,samples=solver.scan_source_brackets_native(points,conversion=conv)
    independent=[solver.matching_residual(float(f),conversion=conv) for f in points]
    assert all(v["status"]==0 for v in samples)
    for native_sample,oracle in zip(samples,independent):
        assert native_sample["determinant"]==pytest.approx(oracle,abs=5e-4)
        assert native_sample["source_restarts"]==expected_knots


def test_original_source_fortran_f_mode_zero_nodes_and_junction(scenario):
    model,solver,conv,p1_bracket,p1_ref,knots=scenario
    f_reference={"SPG_M2":1842.3116593617328,
                 "SPG_M3":1564.2448940763072}[model]
    mode=solver.solve_mode((f_reference-40.,f_reference+40.),
         conversion=conv,profile_points=360,root_backend="fortran_source_exact")
    assert abs(mode["frequency_hz"]-f_reference)<.02
    assert (mode["W_nodes_per_phase"],mode["V_nodes_per_phase"])==(0,0)
    assert mode["node_count_backend"]=="fortran_phasewise"
    assert mode["relative_join_mismatch"]<1e-5
    assert mode["source_callback_rhs_calls"]==0
    assert mode["root_fortran_restart_count"]==knots
