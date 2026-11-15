"""Fortran-backed Python frontend for cold CompOSE EOS handling."""
from .eos_data import EOSData
from .api import read_compose,write_compose,check_consistency,normalize_compose,build_core_only,build_crust_only,find_density_gap,fill_density_gap,causal_completion,construct_maxwell_transition,reconstruct_master,match_crust_core,select_slice
from .native import NativeExecutionError,NativeFeatureUnavailable,executable
__all__=["EOSData","read_compose","write_compose","check_consistency","normalize_compose","build_core_only","build_crust_only",
         "find_density_gap","fill_density_gap","causal_completion","construct_maxwell_transition","reconstruct_master","match_crust_core","select_slice",
         "NativeExecutionError","NativeFeatureUnavailable","executable"]
