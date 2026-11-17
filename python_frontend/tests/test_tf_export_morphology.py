"""Export labels must describe profiles, not just two endpoint values."""
from pathlib import Path
import importlib.util


def builder():
    path=Path(__file__).resolve().parents[2]/'scripts/build_gppva_tf.py'
    spec=importlib.util.spec_from_file_location('tf_morphology_builder',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_interior_dense_shell_has_no_bubble_label():
    row={'dimension':3,'baryon_density_contrast':1.,
         'profile':{'nn':[.02,.05,.06,.05,.024],
                    'np':[.01,.03,.06,.03,.011]}}
    assert builder().physical_geometry(row)=='noncanonical_spherical_profile'


def test_species_modulation_cannot_receive_uniform_export_label():
    row={'dimension':3,'baryon_density_contrast':0.,
         'profile':{'nn':[.055,.045,.035,.045,.055],
                    'np':[.005,.015,.025,.015,.005]}}
    assert builder().physical_geometry(row)=='near_uniform_profile'
