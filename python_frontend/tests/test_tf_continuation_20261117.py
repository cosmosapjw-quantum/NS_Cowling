"""Targeted regression for the mesh continuation output identity."""
import importlib.util
from pathlib import Path
P=Path(__file__).resolve().parents[2]/'scripts/continue_tf_20261117.py'
spec=importlib.util.spec_from_file_location('continuation',P);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def test_mesh_label_changes_after_first_refinement():
    assert mod.refine_label('N032_d2_droplet_n0',64)=='N064_d2_droplet_n0'
def test_seed_interpolation_preserves_requested_means():
    import numpy as np
    old={'profile':[[.003,.001]]*16,'dimension':3,'radius':40}
    nn,pp,r=mod.initial(old,32,.002,.06988)
    w=mod.weights(32,3)
    assert abs(w@(nn+pp)-.002)<1e-17
    assert abs(w@pp/.002-.06988)<1e-14
    assert r==40
