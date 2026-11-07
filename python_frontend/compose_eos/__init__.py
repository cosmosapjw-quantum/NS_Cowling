"""Fortran-backed Python frontend for cold CompOSE EOS handling."""
from .eos_data import EOSData
from .api import read_compose,write_compose,check_consistency,normalize_compose,find_density_gap,fill_density_gap,causal_completion,reconstruct_master,match_crust_core
from .native import NativeExecutionError,NativeFeatureUnavailable,executable
__all__=["EOSData","read_compose","write_compose","check_consistency","normalize_compose",
         "find_density_gap","fill_density_gap","causal_completion","reconstruct_master","match_crust_core",
         "NativeExecutionError","NativeFeatureUnavailable","executable"]
