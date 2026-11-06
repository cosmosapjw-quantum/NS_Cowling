"""Explicit direct-compiler diagnostic, never called by canonical run_case.

This path distinguishes Fortran source/transport defects from missing fpm.
It validates no fpm dependency resolution or release integration.
"""
from pathlib import Path
import os
import shlex
import shutil
from .runner import (
    BuildArtifact, EnvironmentUnavailable, FLAGS, copy_sources, digest,
    execute_command, verify_spec, write_json,
)


def build_diagnostic(root: Path, work: Path, workspace: str) -> BuildArtifact:
    root=Path(root).resolve();work=Path(work).resolve()
    verify_spec(root)
    compiler=shutil.which(os.environ.get('FC','gfortran'))
    if compiler is None:
        raise EnvironmentUnavailable('gfortran diagnostic compiler is unavailable')
    work.mkdir(parents=True,exist_ok=False)
    info=copy_sources(root,work/'source',workspace)
    package=work/'source'/workspace
    build=work/'compile-work';build.mkdir()
    exe=build/'case_driver'
    support = work/'source/support/src'
    sources=[support/name for name in (
        'ns_precision.f90','ns_physical_constants.f90','ns_status.f90','ns_matching_types.f90',
        'ns_numerics_types.f90','ns_scan_types.f90','ns_background_records.f90','ns_background_access.f90','ns_surface_events.f90','ns_dopri_stages.f90','ns_safe_math.f90','ns_units.f90','ns_records.f90','ns_bootstrap_protocol.f90','ns_electron_gas.f90','ns_support_cases.f90')]
    sources += [package/'src/numerics'/name for name in ('ns_bisection.f90','ns_rk4.f90','ns_adaptive_rk.f90','ns_numerics_cases.f90')]
    sources += [package/'src/eos'/name for name in ('ns_eos_types.f90','ns_polytrope.f90','ns_two_polytrope.f90',
        'ns_direct_join.f90','ns_maxwell.f90','ns_mu_bridge.f90','ns_polytrope_cases.f90','ns_matching_cases.f90')]
    if workspace == 'solutions':
        sources += [package/'src/eos'/name for name in ('ns_nuclear_toy.f90','ns_beta_equilibrium.f90','ns_liquid_drop.f90','ns_nuclear_cases.f90')]
    else:
        sources += [package/'src/eos/ns_nuclear_cases.f90']
    # T10 event metadata depends on the public EOS side labels, so it must follow ns_eos_types in this
    # direct-compiler diagnostic ordering.  Canonical fpm builds resolve this dependency automatically.
    sources += [support/'ns_segment_events.f90']
    sources += [package/'src/structure'/name for name in ('ns_tov.f90','ns_target_mass.f90','ns_tov_cases.f90')]
    sources += [package/'src/modes'/name for name in ('ns_scan.f90','ns_shooting_toy.f90','ns_shooting_cases.f90','ns_center_state.f90','ns_surface_residual.f90','ns_cowling.f90','ns_cowling_cases.f90','ns_mode_normalization.f90','ns_mode_identity.f90','ns_cowling_rk4.f90','ns_mode_spectrum.f90','ns_mode_cases.f90','ns_interface_transfer.f90','ns_interface_cases.f90')]
    sources += [package/'src/ns_workspace.f90',package/'app/case_driver.f90']
    argv=[compiler,*shlex.split(FLAGS),'-J',str(build),'-I',str(build),
          *map(str,sources),'-o',str(exe)]
    execute_command(argv,build,work/'compile',30)
    execute_command([compiler,'--version'],build,work/'compiler-version',10)
    sha=digest(exe);receipt=work/'build.json'
    write_json(receipt,{
        'state':'built','backend':'gfortran_diagnostic','workspace':workspace,
        'canonical_fpm_verified':False,
        'source_sha256':info['sha256'],'source_files':info['files'],
        'binary':str(exe),'binary_sha256':sha,
        'command_receipt':str(work/'compile/command.json'),
        'compiler':compiler,'flags':FLAGS,
        'compiler_version':(work/'compiler-version/stdout.txt').read_text().strip(),
        'purpose':'Compiler/transport diagnostic, not fpm or physical validation',
    })
    return BuildArtifact(exe,work/'source',receipt,workspace,info['sha256'],sha,'gfortran_diagnostic')
