import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]




def _synthetic_nb_text():
    return " 2\n 2\n 1.000000e-03\n 2.000000e-03\n 3.000000e-03\n"


def _thermo_row(inb, q1, q3, q7):
    return f"0 {inb} 1 {q1:.9e} 0.0 {q3:.9e} 0.0 0.0 0.0 {q7:.9e} 0\n"


def test_normalizer_preserves_order_duplicates_and_rest_mass_semantics(tmp_path):
    from tools import normalize_compose as nc

    nb = tmp_path / 'eos.nb'
    thermo = tmp_path / 'eos.thermo'
    nb.write_text(_synthetic_nb_text())
    thermo.write_text(
        '9.395654e+02 9.382720e+02 1\n'
        + _thermo_row(0, 1.0, 0.10, 0.02)
        + _thermo_row(1, 2.0, -0.10, 0.03)
        + _thermo_row(1, 2.1, -0.09, 0.031)
        + _thermo_row(2, 3.0, 0.20, 0.04)
    )

    result = nc.normalize_cold_table(nb, thermo)
    assert [r['i_nb'] for r in result['rows']] == [0, 1, 1, 2]
    assert result['duplicate_index_triples'] == [[0, 1, 1]]
    assert result['mu_b_nonmonotone_in_source_order'] is True
    assert result['energy_includes_rest_mass'] is True
    assert result['normalized_units'] == 'SI'
    assert result['rows'][0]['n_b_m3'] == pytest.approx(1.0e42, rel=1e-15)
    assert result['rows'][0]['pressure_pa'] == pytest.approx(1.602176634e29, rel=2e-15)
    expected_eps = 1.0e-3 * 939.5654 * 1.02 * 1.602176634e32
    assert result['rows'][0]['energy_j_m3'] == pytest.approx(expected_eps, rel=2e-15)


def test_normalizer_rejects_bad_grid_index_without_sorting(tmp_path):
    from tools import normalize_compose as nc

    nb = tmp_path / 'eos.nb'
    thermo = tmp_path / 'eos.thermo'
    nb.write_text(_synthetic_nb_text())
    thermo.write_text(
        '9.395654e+02 9.382720e+02 1\n'
        + _thermo_row(0, 1.0, 0.1, 0.02)
        + _thermo_row(9, 2.0, 0.2, 0.03)
    )
    with pytest.raises(nc.ComposeFormatError, match='density index'):
        nc.normalize_cold_table(nb, thermo)


def _sha256(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()




def test_compose_134_real_table_normalization_and_independent_samples():
    from tools import normalize_compose as nc

    raw = ROOT / 'data/raw/compose_134'
    result = nc.normalize_cold_table(raw / 'eos.nb', raw / 'eos.thermo')
    assert len(result['rows']) == 197
    assert result['duplicate_index_triples'] == []
    assert result['mu_b_nonmonotone_in_source_order'] is False
    assert result['energy_includes_rest_mass'] is True
    assert result['normalized_units'] == 'SI'

    first = result['rows'][0]
    assert first['n_b_fm3'] == pytest.approx(1.0e-7, rel=0, abs=1e-20)
    assert first['pressure_pa'] == pytest.approx(3.433863468643866e24, rel=2e-15)
    assert first['energy_j_m3'] == pytest.approx(1.4912317051883251e28, rel=2e-15)
    assert first['mu_b_mev'] == pytest.approx(930.9770468995112, rel=2e-15)

    middle = result['rows'][100]
    assert middle['n_b_fm3'] == pytest.approx(1.0e-3, rel=0, abs=1e-18)
    assert middle['pressure_pa'] == pytest.approx(1.805318211601494e29, rel=2e-15)
    assert middle['energy_j_m3'] == pytest.approx(1.5054237858103122e32, rel=2e-15)
    assert middle['mu_b_mev'] == pytest.approx(940.7833388971466, rel=2e-15)

    last = result['rows'][-1]
    assert last['n_b_fm3'] == pytest.approx(1.21, rel=0, abs=1e-14)
    assert last['pressure_pa'] == pytest.approx(1.3839708437412292e35, rel=2e-15)
    assert last['energy_j_m3'] == pytest.approx(2.574889486257928e35, rel=2e-15)
    assert last['mu_b_mev'] == pytest.approx(2042.0110390478, rel=2e-15)








def test_normalizer_cli_writes_declared_outputs(tmp_path):
    import subprocess
    import sys

    out_csv = tmp_path / 'table.csv'
    out_receipt = tmp_path / 'receipt.json'
    proc = subprocess.run(
        [
            sys.executable, str(ROOT / 'tools/normalize_compose.py'),
            '--metadata', str(ROOT / 'data/metadata/compose_134.json'),
            '--raw-dir', str(ROOT / 'data/raw/compose_134'),
            '--out-csv', str(out_csv),
            '--out-receipt', str(out_receipt),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    summary = json.loads(proc.stdout)
    assert summary['status'] == 'ok'
    assert summary['compose_id'] == 134
    assert out_csv.is_file() and out_receipt.is_file()


def test_compose_134_intake_contract_is_evidence_bound(contract):
    c = contract('compose_134_intake')
    assert c['input']['compose_id'] == 134
    assert c['raw_sha256_verified'] is True
    assert c['energy_includes_rest_mass'] is True
    assert c['normalized_units'] == 'SI'
    assert c['sample_rows_independently_recovered'] is True
    assert c['row_count'] == 197
    assert c['scientific_admission'] == 'T13_INTAKE_ONLY_NOT_T14_TABLE_EOS_ADMISSION'
    assert c['official_checksum_matches_observed'] is False
