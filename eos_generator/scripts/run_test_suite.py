#!/usr/bin/env python3
"""
Zero-dependency pytest-compatible regression test runner.
Runs all test_*.py modules in tests/ using either installed pytest or a built-in
standard-library fixture shim (supporting tmp_path, monkeypatch, pytest.raises,
and pytest.approx) so tests execute in offline sandboxes without pytest installed.
"""
import contextlib
import importlib.util
import inspect
import math
import os
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "compose_eos_pkg"))


class _Approx:
    def __init__(self, expected, rel=1e-6, abs_=1e-12):
        self.expected = expected
        self.rel = rel
        self.abs = abs_

    def __eq__(self, actual):
        return bool(math.isclose(float(actual), float(self.expected), rel_tol=self.rel, abs_tol=self.abs))

    def __repr__(self):
        return f"approx({self.expected} ± rel={self.rel}, abs={self.abs})"


class _MonkeyPatch:
    def __init__(self):
        self._undo = []

    def setattr(self, target, name, value):
        old = getattr(target, name)
        setattr(target, name, value)
        self._undo.append((target, name, old))

    def undo(self):
        for target, name, old in reversed(self._undo):
            setattr(target, name, old)
        self._undo.clear()


class _PytestShim:
    @staticmethod
    @contextlib.contextmanager
    def raises(exc_type, match=None):
        try:
            yield
        except exc_type as exc:
            if match is not None and not re.search(match, str(exc)):
                raise AssertionError(f"Exception {exc!r} did not match pattern {match!r}")
            return
        raise AssertionError(f"Expected {exc_type.__name__} to be raised, but no exception was raised")

    @staticmethod
    def approx(expected, rel=1e-6, abs=1e-12):
        return _Approx(expected, rel=rel, abs_=abs)


if "pytest" not in sys.modules:
    sys.modules["pytest"] = _PytestShim()


def run_all_tests() -> int:
    t0 = time.time()
    test_dir = ROOT / "tests"
    test_files = sorted(test_dir.glob("test_*.py"))
    passed = 0
    failed = 0

    for tf in test_files:
        spec = importlib.util.spec_from_file_location(tf.stem, tf)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        for name, fn in inspect.getmembers(mod, inspect.isfunction):
            if not name.startswith("test_"):
                continue
            sig = inspect.signature(fn)
            tmp_dir = None
            mp = _MonkeyPatch()
            kwargs = {}
            try:
                if "tmp_path" in sig.parameters:
                    tmp_dir = Path(tempfile.mkdtemp(prefix=f"test_{name}_"))
                    kwargs["tmp_path"] = tmp_dir
                if "monkeypatch" in sig.parameters:
                    kwargs["monkeypatch"] = mp
                fn(**kwargs)
                print(f"  [PASS] {tf.name}::{name}")
                passed += 1
            except Exception as exc:
                print(f"  [FAIL] {tf.name}::{name} -> {type(exc).__name__}: {exc}")
                failed += 1
            finally:
                mp.undo()
                if tmp_dir is not None:
                    shutil.rmtree(tmp_dir, ignore_errors=True)

    dt = time.time() - t0
    print(f"\nTest Summary: {passed} passed, {failed} failed in {dt:.2f}s")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run_all_tests())
