import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from refine_dd2_geometry_quarters import strict_select


def test_exception_for_a_seed_cannot_be_hidden_by_other_successes():
    attempts=[{'initial_density':.06,'numerical_status':'STATIONARY','muon_threshold_crossed':False,'epsilon':1.,'energy_per_baryon':10.},
              {'initial_density':.065,'execution':'FAILED'},
              {'initial_density':.07,'numerical_status':'STATIONARY','muon_threshold_crossed':False,'epsilon':1.1,'energy_per_baryon':11.}]
    best,status=strict_select(attempts,[.06,.065,.07])
    assert best is None
    assert status=='HOLD_UNRESOLVED_SEED'


def test_recovered_seed_can_be_selected_and_lower_failed_energy_still_holds():
    def row(seed,energy,valid=True):
        return {'initial_density':seed,'numerical_status':'STATIONARY' if valid else 'FAILED',
                'muon_threshold_crossed':False,'epsilon':energy,'energy_per_baryon':energy}
    rows=[row(.06,10.),{'initial_density':.065,'execution':'FAILED'},row(.065,9.),row(.07,11.)]
    best,status=strict_select(rows,[.06,.065,.07])
    assert best['initial_density']==.065
    assert status=='LOWEST_EVALUATED_STATIONARY_CANDIDATE'
    rows.append(row(.06,8.,False))
    best,status=strict_select(rows,[.06,.065,.07])
    assert best is None
    assert status=='LOWER_UNCONVERGED_ATTEMPT'
