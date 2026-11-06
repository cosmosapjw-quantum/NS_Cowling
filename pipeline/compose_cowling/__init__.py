"""Portable EOS-to-relativistic-Cowling pipeline."""
from .adapter import build_adapter, build_adapter_csv
__all__ = ["build_adapter", "build_adapter_csv"]
__version__ = "1.0.0"
