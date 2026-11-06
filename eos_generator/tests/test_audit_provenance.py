from pathlib import Path
import importlib.util

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_eos_archive.py"
SPEC = importlib.util.spec_from_file_location("audit_eos_archive_under_test", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
SPEC.loader.exec_module(AUDIT)


def test_audit_provenance_never_serializes_machine_local_root(tmp_path):
    root = tmp_path / "compose_release" / "eos_list"
    table = root / "APR" / "APR"
    table.mkdir(parents=True)
    assert AUDIT.portable_root_label(root) == "eos_list"
    assert AUDIT.relative_source_path(table, root) == "APR/APR"
    assert str(tmp_path) not in AUDIT.portable_root_label(root)
    assert str(tmp_path) not in AUDIT.relative_source_path(table, root)


def test_audit_provenance_handles_files_outside_root(tmp_path):
    root = tmp_path / "root"
    external = tmp_path / "elsewhere" / "table"
    root.mkdir()
    assert AUDIT.relative_source_path(external, root) == "table"
