"""Regression contracts for native-density holdout, not proxy ground truth."""
import importlib
import json
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[2]
def get():
    return importlib.import_module("run_dense_tf_20261118")
def test_all_tf_cases_have_strict_interleaved_density_grids():
    mod=get()
    assert len(mod.CASES)==5
    for key,spec in mod.CASES.items():
        targets=mod.target_densities(spec)
        assert len(targets)>=8
        assert all(np.isfinite(targets))
        assert np.min(np.diff(targets))>0
        assert spec["n_low"]<=targets[0]<targets[-1]<=spec["n_high"]
        assert targets[-1] <0.05
def test_tf_h_source_provenance_and_anchor_available():
    m=get()
    for key,spec in m.CASES.items():
        states=m.existing_anchors(spec)
        assert states, f"no previously accepted N128 source for {key}"
        a=m.choose_anchor(states,float(np.mean(m.target_densities(spec))))
        assert a["accepted_stationary"]
        assert a["points"]==128
        assert a["geometry"]=="droplet"
        assert "profile" in a and len(a["profile"])==128
def test_all_uniform_models_frozen_and_grid_dense():
    from run_dense_rmf_20261118 import MODELS,grid
    assert len(MODELS)==5
    n=grid()
    assert len(n)==4097 and n[0]==0.015 and n[-1]==1.0
    assert np.all(np.diff(n)>0)
