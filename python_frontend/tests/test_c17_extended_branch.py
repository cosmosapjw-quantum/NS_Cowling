"""C17: Fortran original-source Maxwell onset and overlap-aware branch checks.

A branch is permitted only when there is one matched nodal candidate,
the mode-frequency step is bounded, and *both physical phase* W-shape
diagnostics are finite and above threshold. This is a local numerical
classification, NOT a physical conserved mode-inner-product proof.
"""
from __future__ import annotations
import ctypes as ct
import numpy as np
import pytest

from compose_eos import stellar_cowling_native as native


def call_branch(freq,nw,nv,accepted,high,low, *,
                mode=1,previous=3700.,maxfraction=.12,minoverlap=.95):
    n=len(freq)
    selected=ct.c_int(-123)
    status=ct.c_int(-123)
    native._native_library().cowling_native_select_mode_branch_overlap(
       n,np.ascontiguousarray(freq,dtype=float),
       np.ascontiguousarray(nw,dtype=np.int32),
       np.ascontiguousarray(nv,dtype=np.int32),
       np.ascontiguousarray(accepted,dtype=np.int32),
       np.ascontiguousarray(high,dtype=float),
       np.ascontiguousarray(low,dtype=float),
       mode,previous,maxfraction,minoverlap,
       ct.byref(selected),ct.byref(status))
    return selected.value,status.value


def test_native_extended_branch_selects_unique_f_and_p1():
    a=call_branch([1850.,3680.,5000.],[0,1,1],[0,1,1],
                 [0,0,0],[.995,.989,.8],[.995,.999,.99],
                 mode=1,previous=3700.)
    assert a==(1,0)
    b=call_branch([1550.,3680.],[0,1],[0,1],
                 [0,0],[.985,.989],[.999,.999],mode=0,previous=1560.)
    assert b==(0,0)


def test_native_extended_branch_blocks_ambiguous_matching_modes():
    a=call_branch([3660.,3690.],[1,1],[1,1],
                  [0,0],[.99,.98],[.99,.99])
    assert a==(-1,3)


def test_native_extended_branch_fails_closed_for_small_core_overlap():
    a=call_branch([3735.],[1],[1],[0],[.60],[.999])
    assert a==(-1,5)
    b=call_branch([3735.],[1],[1],[0],[np.nan],[.999])
    assert b[0]==-1 and b[1]!=0


def test_native_extended_branch_does_not_promote_rejected_roots():
    a=call_branch([3740.],[1],[1],[33],[.99],[.99])
    assert a==(-1,4)
    b=call_branch([3740.],[1],[1],[0],[.99],[.99],
             previous=3000.)
    assert b==(-1,4)


@pytest.mark.parametrize("eps_low,eps_high,p,delta,crit,is_above",[
    (8.572084498151984,16.33826505200197,.014307398278768153,
     .9044754359602759,.5016662875921976,1),
    (76.94982713409907,89.34565739599732,.5657052249953693,
     .15991414732822803,.507297959618915,0),
])
def test_original_source_maxwell_onset_dimensionless_seidov(eps_low,eps_high,p,delta,crit,is_above):
    d,c,margin,flag,status=(ct.c_double(-1),ct.c_double(-1),
         ct.c_double(-1),ct.c_int(-1),ct.c_int(-1))
    native._native_library().cowling_native_maxwell_onset_threshold(
      eps_low,eps_high,p,ct.byref(d),ct.byref(c),
      ct.byref(margin),ct.byref(flag),ct.byref(status))
    assert status.value==0
    assert d.value==pytest.approx(delta,rel=1e-14)
    assert c.value==pytest.approx(crit,rel=1e-14)
    assert margin.value==pytest.approx(delta-crit,rel=1e-14)
    assert flag.value==is_above


def test_native_maxwell_onset_bad_eos_fails_closed():
    d,c,margin,flag,status=(ct.c_double(-1),ct.c_double(-1),
         ct.c_double(-1),ct.c_int(-1),ct.c_int(-1))
    native._native_library().cowling_native_maxwell_onset_threshold(
      8.,7.,1.,ct.byref(d),ct.byref(c),
      ct.byref(margin),ct.byref(flag),ct.byref(status))
    assert status.value!=0
    assert (d.value,c.value,margin.value)==(0.,0.,0.)
    assert flag.value==-1
