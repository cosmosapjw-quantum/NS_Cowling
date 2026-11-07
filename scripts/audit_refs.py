from pathlib import Path
ps=[
'/tmp/ns_oscillation_followup_20261106/source/Neutron_Star_Oscillations/data/normalized/compose_69_t14_si.csv',
'/tmp/ns_oscillation_followup_20261106/source/Neutron_Star_Oscillations/data/normalized/compose_69.csv',
'/tmp/ns_oscillation_followup_20261106/source/Neutron_Star_Oscillations/data/normalized/compose_134_si.csv']
Path('/tmp/ns_oscillation_followup_20261106/reference_paths.txt').write_text('\n'.join(ps)+'\n')
