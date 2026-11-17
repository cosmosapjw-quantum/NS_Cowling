"""Counterexamples to geometry classification from endpoint orientation alone."""
import numpy as np
from compose_eos.tf_morphology import profile_morphology

def profile(protons,baryons,dimension=3):
    return {'dimension':dimension,'baryon_density_contrast':1.,
            'profile':{'np':protons,'nn':(np.array(baryons)-protons).tolist()}}

def test_strong_interior_shell_does_not_masquerade_as_bubble():
    row=profile([.01,.03,.06,.03,.011],[.03,.08,.12,.08,.035])
    assert profile_morphology(row)['classification']=='NONCANONICAL_SPHERICAL_PROFILE'
    row['dimension']=2
    assert profile_morphology(row)['classification']=='NONCANONICAL_CYLINDRICAL_PROFILE'

def test_mild_central_proton_depletion_is_not_shell_gate():
    row=profile([.045,.05,.04,.01,.001],[.13,.125,.10,.04,.02])
    assert profile_morphology(row)['classification']=='SIMPLE_RADIAL_PROFILE'

def test_baryon_proton_mismatch_stays_descriptive_unknown():
    row=profile([.01,.03,.06,.03,.011],[.13,.125,.10,.04,.02])
    assert profile_morphology(row)['classification']=='BARYON_PROTON_SHAPE_DIFFERENCE'

def test_planar_repeated_cell_retains_planar_geometry():
    row=profile([.05,.03,.001,.03,.05],[.12,.08,.03,.08,.12],1)
    assert profile_morphology(row)['classification']=='PLANAR_REPEATED_OR_TRANSLATED'

def test_constant_proton_density_has_defined_diagnostic():
    row=profile([.01]*5,[.03,.08,.12,.08,.03])
    assert profile_morphology(row)['classification']=='BARYON_PROTON_SHAPE_DIFFERENCE'


def test_constant_total_density_does_not_hide_composition_structure():
    row=profile([.005,.015,.025,.015,.005],[.06]*5)
    row['baryon_density_contrast']=0.
    assert profile_morphology(row)['classification']=='NEAR_UNIFORM_COMPOSITION_PROFILE'


def test_uniform_species_profile_is_uniform_without_dividing_by_zero():
    row=profile([0.]*5,[.06]*5)
    row['baryon_density_contrast']=0.
    assert profile_morphology(row)['classification']=='UNIFORM'
