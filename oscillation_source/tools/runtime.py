"""T00 source-build and process transport. No stellar equations live here.

Canonical execution ALWAYS uses fpm. Compiler-control tests explicitly create a
separate gfortran_diagnostic artifact; it is never a fallback for missing fpm.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time
import tomllib
from typing import Any, Mapping, Sequence
from uuid import uuid4

SPEC_NAME = 'conventions.md'
SPEC_SHA256 = 'b025b942eff00abcb4b9487317577ba9e057c52354ad4738a51f23254e0b4617'
FLAGS = '-std=f2008 -O0 -g -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -fno-backslash'
NAME = re.compile(r'[A-Za-z][A-Za-z0-9_]*\Z')
ENV_PARAMETERS = {'probe_real':'real','probe_int':'int','probe_bool':'bool','probe_text':'text'}
CASE_SCHEMAS = {
    'environment': ENV_PARAMETERS,
    'constants': {},
    'units': {'energy_mev_fm3':'real', 'n_b_fm3':'real', 'omega_rad_s':'real',
              'q_solar':'real', 'radius_km':'real'},
    'conversion': {'operation':'text','value':'real','mass_b_kg':'real','special_input':'text'},
    'safe_norm': {'x':'real','y':'real','z':'real','special_input':'text','empty_vector':'bool'},
    'status_lifecycle': {},
}
CASE_SCHEMAS.update({'root_sqrt2': {'lo': 'real', 'hi': 'real', 'xtol': 'real', 'rtol': 'real', 'ftol': 'real', 'max_iter': 'int', 'root_at': 'real', 'scale': 'real', 'problem': 'text'}, 'root_probe': {'lo': 'real', 'hi': 'real', 'xtol': 'real', 'rtol': 'real', 'ftol': 'real', 'max_iter': 'int', 'root_at': 'real', 'scale': 'real', 'problem': 'text'}, 'ode': {'x0': 'real', 'x1': 'real', 'y_initial': 'real', 'y2_initial': 'real', 'rate': 'real', 'h_initial': 'real', 'h_min': 'real', 'h_max': 'real', 'rtol': 'real', 'atol1': 'real', 'atol2': 'real', 'event1': 'real', 'event2': 'real', 'n_steps': 'int', 'max_steps': 'int', 'record_trials': 'bool', 'use_events': 'bool', 'method': 'text', 'problem': 'text'}, 'rk4_exp': {'x0': 'real', 'x1': 'real', 'y_initial': 'real', 'y2_initial': 'real', 'rate': 'real', 'h_initial': 'real', 'h_min': 'real', 'h_max': 'real', 'rtol': 'real', 'atol1': 'real', 'atol2': 'real', 'event1': 'real', 'event2': 'real', 'n_steps': 'int', 'max_steps': 'int', 'record_trials': 'bool', 'use_events': 'bool', 'method': 'text', 'problem': 'text'}, 'adaptive_exp': {'x0': 'real', 'x1': 'real', 'y_initial': 'real', 'y2_initial': 'real', 'rate': 'real', 'h_initial': 'real', 'h_min': 'real', 'h_max': 'real', 'rtol': 'real', 'atol1': 'real', 'atol2': 'real', 'event1': 'real', 'event2': 'real', 'n_steps': 'int', 'max_steps': 'int', 'record_trials': 'bool', 'use_events': 'bool', 'method': 'text', 'problem': 'text'}, 'dopri_pair': {'x0': 'real', 'x1': 'real', 'y_initial': 'real', 'y2_initial': 'real', 'rate': 'real', 'h_initial': 'real', 'h_min': 'real', 'h_max': 'real', 'rtol': 'real', 'atol1': 'real', 'atol2': 'real', 'event1': 'real', 'event2': 'real', 'n_steps': 'int', 'max_steps': 'int', 'record_trials': 'bool', 'use_events': 'bool', 'method': 'text', 'problem': 'text'}, 'control': {'old1': 'real', 'old2': 'real', 'new1': 'real', 'new2': 'real', 'error1': 'real', 'error2': 'real', 'error_value': 'real', 'rtol': 'real', 'atol1': 'real', 'atol2': 'real', 'action': 'text', 'previous_rejected': 'bool'}})

CASE_SCHEMAS["poly_state"] = {'n_anchor_m3': 'real', 'p_anchor_pa': 'real', 'gamma': 'real', 'c_j': 'real', 'n_min_m3': 'real', 'n_max_m3': 'real', 'n_b_m3': 'real', 'value': 'real', 'k_si': 'real', 'query': 'text', 'constructor': 'text', 'action': 'text', 'special_input': 'text', 'side': 'int'}

TOV_SCHEMA = {'n_anchor_m3': 'real', 'p_anchor_pa': 'real', 'gamma': 'real', 'c_j': 'real', 'n_min_m3': 'real', 'n_max_m3': 'real', 'pc_pa': 'real', 'p_cut_pa': 'real', 'r0_m': 'real', 'r_limit_m': 'real', 'rtol': 'real', 'atol_r_m': 'real', 'atol_mass_kg': 'real', 'atol_nu': 'real', 'dq_initial': 'real', 'dq_min': 'real', 'dq_max': 'real', 'nu_start': 'real', 'radius_m': 'real', 'mass_kg': 'real', 'pressure_pa': 'real', 'energy_uniform': 'real', 'pc_lo_pa': 'real', 'pc_hi_pa': 'real', 'm_target_kg': 'real', 'mass_atol_kg': 'real', 'mass_rtol': 'real', 'root_xtol': 'real', 'stop1_pa': 'real', 'stop2_pa': 'real', 'max_steps': 'int', 'scan_n': 'int', 'root_max_iter': 'int', 'emit_profile': 'bool', 'use_stops': 'bool', 'eos_kind': 'text', 'special_input': 'text', 'surface_policy': 'text', 'branch_id': 'text'}
for _case in ("tov_baseline","tov_rhs","target_mass"):
    CASE_SCHEMAS[_case] = TOV_SCHEMA

CASE_SCHEMAS.update({'string_shot': {'length_m': 'real', 'k_per_m': 'real', 'k_min': 'real', 'k_max': 'real', 'scan_n': 'int', 'amplitude': 'real', 'method': 'text', 'n_steps': 'int', 'ode_rtol': 'real', 'atol_y_m': 'real', 'atol_slope': 'real', 'h_initial_m': 'real', 'h_max_m': 'real', 'max_steps': 'int', 'root_xtol': 'real', 'root_rtol': 'real', 'max_iter': 'int', 'zero_tol': 'real', 'residual_tol': 'real', 'emit_profile': 'bool', 'special_input': 'text'}, 'string_spectrum': {'length_m': 'real', 'k_per_m': 'real', 'k_min': 'real', 'k_max': 'real', 'scan_n': 'int', 'amplitude': 'real', 'method': 'text', 'n_steps': 'int', 'ode_rtol': 'real', 'atol_y_m': 'real', 'atol_slope': 'real', 'h_initial_m': 'real', 'h_max_m': 'real', 'max_steps': 'int', 'root_xtol': 'real', 'root_rtol': 'real', 'max_iter': 'int', 'zero_tol': 'real', 'residual_tol': 'real', 'emit_profile': 'bool', 'special_input': 'text'}, 'bracket_probe': {'scenario': 'text', 'zero_tol': 'real'}, 'scan_with_failed_trial': {'scenario': 'text', 'zero_tol': 'real'}})

CASE_SCHEMAS['string_refine'] = CASE_SCHEMAS['string_shot']

COWLING_SCHEMA = {'background_path':'text','special_input':'text','r0_m':'real','amplitude':'real',
                  'ell':'int','r_m':'real','r_probe_m':'real','omega_rad_s':'real','w_m3':'real',
                  'v_m2':'real','chi':'real','nu':'real','lambda_value':'real','dnu_dr':'real',
                  'nu_shift':'real','ode_rtol':'real','atol_w_m3':'real','atol_v_m2':'real',
                  'h_initial':'real','h_min':'real','h_max':'real','max_steps':'int','emit_profile':'bool'}
for _case in ('cowling_center','cowling_rhs_point','background_probe','surface_residual','cowling_shot'):
    CASE_SCHEMAS[_case] = COWLING_SCHEMA


# T09 matching input contract. Exact analytic fixture scenarios live in the selected Fortran workspace.
MATCHING_SCHEMA = {'action':'text','scenario':'text','pressure_pa':'real','mu_left_j':'real','mu_right_j':'real',
                   'mu_low_j':'real','mu_high_j':'real','mu_query_j':'real','root_xtol':'real','root_rtol':'real',
                   'residual_tol':'real','max_iter':'int','gamma_crust':'real','transition_n_m3':'real',
                   'crust_n_min_m3':'real','core_n_anchor_m3':'real','core_p_anchor_pa':'real',
                   'core_gamma':'real','core_c_j':'real','core_n_min_m3':'real','core_n_max_m3':'real',
                   'special_input':'text'}
for _case in ('two_poly_knot','direct_join','maxwell_match','mu_bridge'):
    CASE_SCHEMAS[_case] = MATCHING_SCHEMA


INTERFACE_SCHEMA = {'scenario':'text','reference_label':'text','r_m':'real','omega_rad_s':'real',
                    'pressure_pa':'real','eps_in':'real','eps_out':'real','nu':'real','lambda_value':'real',
                    'dnu_dr':'real','w_in_m3':'real','v_in_m2':'real','inner_side':'int','outer_side':'int',
                    'overlap':'real','candidate_interface_g':'bool','root_valid':'bool','transfer_valid':'bool'}
for _case in ('material_transfer','segment_events','segment_restart_probe','interface_identity'):
    CASE_SCHEMAS[_case] = INTERFACE_SCHEMA

# T11 nuclear-phenomenology teaching contracts.  The JSON registry below is
# checked byte-for-byte against these allowed parameter maps before execution.
ELECTRON_SCHEMA = {'n_e_m3':'real'}
NUCLEAR_SCHEMA = {'regulated':'bool','n0_fm3':'real','binding_mev':'real','k0_mev':'real',
                  's0_mev':'real','regulator_c':'real','n_b_fm3':'real','alpha':'real',
                  'scan_u_max':'real'}
CASE_SCHEMAS['electron_gas'] = ELECTRON_SCHEMA
CASE_SCHEMAS['liquid_drop_scaling'] = {}
LIQUID_DROP_SCHEMA = {'n_b_fm3': 'real', 'a': 'int', 'z': 'int', 'a_min': 'int', 'a_max': 'int', 'a_step': 'int', 'z_min': 'int', 'z_max': 'int', 'z_step': 'int', 'objective': 'text'}
for _case in ('liquid_drop_candidate','liquid_drop_scan'):
    CASE_SCHEMAS[_case] = LIQUID_DROP_SCHEMA
for _case in ('nuclear_constraints','nuclear_state','beta_equilibrium'):
    CASE_SCHEMAS[_case] = NUCLEAR_SCHEMA

TABLE_SCHEMA = {'table_path':'text','query_kind':'text','query_value':'real','side':'int'}
for _case in ('table_diagnostics','table_state'):
    CASE_SCHEMAS[_case] = TABLE_SCHEMA
CASE_SCHEMAS['real_table_admission'] = {'table_path':'text'}
TABLE_STAR_SCHEMA = {'table_path':'text','pc_pa':'real','p_cut_pa':'real','r0_m':'real','r_limit_m':'real',
    'rtol':'real','atol_r_m':'real','atol_mass_kg':'real','atol_nu':'real','dq_initial':'real','dq_min':'real',
    'dq_max':'real','nu_start':'real','max_steps':'int','emit_profile':'bool'}
CASE_SCHEMAS['real_table_star'] = TABLE_STAR_SCHEMA
CASE_SCHEMAS['real_table_star_split'] = {**TABLE_STAR_SCHEMA, 'split_pressure_pa':'real'}

STRUCTURE_DIAGNOSTICS_SCHEMA = {'background_path':'text','split_pressure_pa':'real','f_hz':'real','p1_hz':'real'}
MODE_DIFFERENCE_SCHEMA = {'f_ref_hz':'real','f_cmp_hz':'real','p1_ref_hz':'real','p1_cmp_hz':'real','resolution_hz':'real'}
CASE_SCHEMAS['structure_diagnostics'] = STRUCTURE_DIAGNOSTICS_SCHEMA
CASE_SCHEMAS['mode_difference'] = MODE_DIFFERENCE_SCHEMA
PAIRED_REFINEMENT_SCHEMA = {
    'delta_f_coarse_hz':'real','delta_f_medium_hz':'real','delta_f_fine_hz':'real',
    'delta_p1_coarse_hz':'real','delta_p1_medium_hz':'real','delta_p1_fine_hz':'real',
    'absolute_resolution_hz':'real','mass_residual_ref_kg':'real','mass_residual_cmp_kg':'real',
    'df_dmass_hz_per_kg':'real','dp1_dmass_hz_per_kg':'real','f_ref_label':'text','f_cmp_label':'text',
    'p1_ref_label':'text','p1_cmp_label':'text'}
CASE_SCHEMAS['paired_refinement'] = PAIRED_REFINEMENT_SCHEMA
CASE_SCHEMAS['controlled_table_variant'] = {'table_path':'text','alpha':'real','split_index':'int'}
CONTROLLED_TARGET_MASS_SCHEMA = {
    'table_path':'text','alpha':'real','split_index':'int','m_target_kg':'real',
    'pc_lo_pa':'real','pc_hi_pa':'real','p_cut_pa':'real','r0_m':'real','r_limit_m':'real','rtol':'real',
    'atol_r_m':'real','atol_mass_kg':'real','atol_nu':'real','dq_initial':'real','dq_min':'real','dq_max':'real',
    'nu_start':'real','max_steps':'int','emit_profile':'bool','mass_atol_kg':'real','mass_rtol':'real',
    'root_xtol':'real','root_rtol':'real','root_max_iter':'int','scan_n':'int','branch_id':'text'}
CASE_SCHEMAS['controlled_table_target_mass'] = CONTROLLED_TARGET_MASS_SCHEMA

SUPPORT_STATUSES = {'ok','unsupported_case','invalid_input','uninitialized','not_implemented',
                    'nonfinite_input','out_of_domain','range_error','no_bracket','identity_unresolved','max_iterations','stagnation','step_too_small','max_steps','nonfinite_result','matching_inconsistent'}
RESERVED = {'built_from_source','build_backend','source_sha256','binary_sha256','_evidence'}

# T08 public input contract
for _case in ['smooth_spectrum', 'mode_refine', 'mode_trial', 'mode_diagnostics', 'mode_identity_probe']:
    CASE_SCHEMAS[_case] = {'frequency_min_hz': 'real', 'frequency_max_hz': 'real', 'frequency_hz': 'real', 'ode_rtol': 'real', 'amplitude': 'real', 'nu_shift': 'real', 'root_xtol': 'real', 'root_rtol': 'real', 'residual_tol': 'real', 'node_zero': 'real', 'node_xlo': 'real', 'node_xhi': 'real', 'h_max': 'real', 'scan_n': 'int', 'cell_substeps': 'int', 'max_steps': 'int', 'max_iter': 'int', 'background_path': 'text', 'method': 'text', 'emit_profile': 'bool', 'diagnostic_action': 'text'}

class BootstrapError(RuntimeError):
    """A recorded harness failure, not a scientific zero."""
class SourceMismatch(BootstrapError): pass
class DependencyViolation(BootstrapError): pass
class EnvironmentUnavailable(BootstrapError): pass
class ProcessFailure(BootstrapError): pass
class ProcessTimeout(BootstrapError): pass
class OutputMalformed(BootstrapError): pass
class CalibrationMissing(BootstrapError): pass


from tools.atomic_io import digest, write_json


def verify_spec(root: Path) -> dict:
    path = Path(root)/'docs/specs'/SPEC_NAME
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise SourceMismatch(f'Missing authoritative specification: {path}') from exc
    sha = hashlib.sha256(raw).hexdigest()
    if sha != SPEC_SHA256:
        raise SourceMismatch(f'Wrong specification bytes: expected {SPEC_SHA256}, got {sha}')
    return {'sha256':sha, 'size_bytes':len(raw), 'path':str(path)}


def _regular_tree(directory: Path) -> None:
    if directory.is_symlink() or not directory.is_dir():
        raise DependencyViolation(f'Expected a regular source directory: {directory}')
    for item in directory.rglob('*'):
        if item.is_symlink():
            raise DependencyViolation(f'Source symlink is forbidden: {item}')
        if item.is_file() and item.suffix.lower() != '.f90':
            raise DependencyViolation(f'Undeclared source format in source package: {item}')


def validate_packages(root: Path, workspace: str) -> dict:
    root = Path(root).resolve()
    if workspace not in {'solutions','learner'}:
        raise DependencyViolation(f'Invalid workspace: {workspace!r}')
    for name in ('support',workspace):
        package=root/name; manifest=package/'fpm.toml'
        if package.is_symlink() or manifest.is_symlink():
            raise DependencyViolation('Package/manifest symlink is forbidden')
        try:
            data=tomllib.loads(manifest.read_text(encoding='utf-8'))
        except (OSError, ValueError) as exc:
            raise DependencyViolation(f'Invalid/missing fpm manifest: {manifest}') from exc
        build=data.get('build',{})
        for key in ('auto-executables','auto-tests','auto-examples'):
            if build.get(key) is not False:
                raise DependencyViolation(f'{name}: {key} must explicitly be false')
        if set(build)-{'auto-executables','auto-tests','auto-examples','module-naming'}:
            raise DependencyViolation(f'{name}: undeclared build options')
        if data.get('library') != {'source-dir':'src'}:
            raise DependencyViolation(f'{name}: library source must be local src')
        _regular_tree(package/'src')
        if name=='support':
            if data.get('dependencies',{}) or data.get('executable') or data.get('test'):
                raise DependencyViolation('Support must be a dependency-free library in this library')
        else:
            if data.get('dependencies') != {'ns-support':{'path':'../support'}}:
                raise DependencyViolation(f'{workspace}: only ../support is permitted')
            expected=[{'name':'case_driver','source-dir':'app','main':'case_driver.f90'}]
            if (package/'app/cell_solve.f90').is_file():
                expected.append({'name':'cell_solve','source-dir':'app','main':'cell_solve.f90'})
            expected_tests=[]
            if (package/'tests_fortran/cell_checks.f90').is_file():
                expected_tests=[{'name':'cell_checks','source-dir':'tests_fortran','main':'cell_checks.f90'}]
                _regular_tree(package/'tests_fortran')
            if data.get('executable') != expected or data.get('test',[]) != expected_tests or data.get('example'):
                raise DependencyViolation('Only the explicitly declared case_driver, cell_solve and cell_checks targets are permitted')
            _regular_tree(package/'app')
            if not (package/'app/case_driver.f90').is_file():
                raise DependencyViolation('Missing case driver source; binaries are not a substitute')
        # T00 has no preprocess, build scripts, platform features or transitive dependencies.
        allowed={'name','version','description','build','library','fortran','dependencies','executable','test'}
        if set(data)-allowed:
            raise DependencyViolation(f'{name}: undeclared manifest sections: {sorted(set(data)-allowed)}')
    return {'workspace':workspace,'dependencies':['support'],'source_boundary':'explicit'}


def source_fingerprint(root: Path, workspace: str) -> dict:
    validate_packages(root,workspace)
    root=Path(root).resolve(); rows=[]
    for name in ('support',workspace):
        package=root/name
        paths=[package/'fpm.toml']
        for sub in ('src','app','tests_fortran'):
            if (package/sub).exists():
                paths.extend(p for p in (package/sub).rglob('*') if p.is_file())
        for path in sorted(paths):
            raw=path.read_bytes()
            rows.append({'path':path.relative_to(root).as_posix(),'size_bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
    wire=json.dumps(rows,sort_keys=True,separators=(',',':')).encode()
    return {'algorithm':'sha256(canonical-json(file-path,size,sha256))','sha256':hashlib.sha256(wire).hexdigest(),'files':rows}


def copy_sources(root: Path, destination: Path, workspace: str) -> dict:
    info=source_fingerprint(root,workspace)
    destination.mkdir(parents=True,exist_ok=False)
    for row in info['files']:
        target=destination/row['path'];target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((Path(root)/row['path']).read_bytes())
    if source_fingerprint(destination,workspace)['sha256'] != info['sha256']:
        raise SourceMismatch('Source snapshot changed during copy')
    return info


def _unescape(text: str) -> str:
    out=[]; i=0; escape={'t':'\t','n':'\n','r':'\r','\\':'\\'}
    while i<len(text):
        if text[i]!='\\': out.append(text[i]);i+=1;continue
        i+=1
        if i==len(text) or text[i] not in escape:
            raise OutputMalformed('Invalid TSV escape')
        out.append(escape[text[i]]);i+=1
    return ''.join(out)


def parse_records(text: str) -> dict:
    lines=text.splitlines()
    if not lines or lines[0] != 'name\ttype\tvalue':
        raise OutputMalformed('Missing exact TSV header')
    result={}
    for number,line in enumerate(lines[1:],2):
        cells=line.split('\t')
        if len(cells)!=3:
            raise OutputMalformed(f'TSV row {number}: expected three fields')
        key,kind,raw=cells
        if not NAME.fullmatch(key) or key in result or key in RESERVED:
            raise OutputMalformed(f'Invalid, duplicated or reserved key: {key!r}')
        try:
            if kind=='text': value=_unescape(raw)
            elif kind=='bool':
                if raw not in {'true','false'}: raise ValueError('boolean')
                value=raw=='true'
            elif kind=='int':
                if not re.fullmatch(r'[+-]?\d+',raw):raise ValueError('integer')
                value=int(raw)
            elif kind=='real':
                value=float(raw.replace('D','e').replace('d','e'))
                if not math.isfinite(value):raise ValueError('nonfinite scalar')
            else: raise ValueError('unknown type')
        except (ValueError, OverflowError) as exc:
            raise OutputMalformed(f'Invalid {kind} value at row {number}') from exc
        result[key]=value
    if not result:
        raise OutputMalformed('Empty scalar record body')
    return result


def write_input(path: Path, token: str, parameters: Mapping[str, Any],
                schema: Mapping[str,str] | None = None) -> None:
    if not re.fullmatch(r'[0-9a-f]{32}',token): raise ValueError('Bad request token')
    if not isinstance(parameters,Mapping): raise ValueError('Parameters must be a mapping')
    if not all(isinstance(key,str) for key in parameters):
        raise ValueError('Parameter names must be strings')
    schema = ENV_PARAMETERS if schema is None else schema
    lines=['&request',f" request_token = '{token}'",'/', '&parameters']
    for key,value in sorted(parameters.items()):
        kind=schema.get(key)
        if kind is None:raise ValueError(f'Undeclared parameter: {key}')
        if kind=='real':
            if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value):
                raise ValueError('Expected finite real input')
            text=f'{float(value):.17e}'
        elif kind=='int':
            if type(value) is not int or not (-2**63<=value<2**63):raise ValueError('Expected int64 input')
            text=str(value)
        elif kind=='bool':
            if type(value) is not bool:raise ValueError('Expected boolean input')
            text='.true.' if value else '.false.'
        else:
            if not isinstance(value,str) or len(value.encode('utf-8'))>200 or any(ord(c)<32 for c in value):
                raise ValueError('Expected text <=200 UTF-8 bytes without control characters')
            if value.endswith(' '):
                raise ValueError('probe_text must not have trailing padding spaces')
            text="'"+value.replace("'","''")+"'"
        lines.append(f' {key} = {text}')
    lines.extend(['/',''])
    Path(path).write_text('\n'.join(lines),encoding='utf-8')


def load_contract(root: Path, fixture_id: str) -> dict:
    root=Path(root).resolve()
    if not isinstance(fixture_id,str) or not NAME.fullmatch(fixture_id):
        raise CalibrationMissing('Invalid fixture ID')
    path=root/'fixtures/contracts'/f'{fixture_id}.json'
    try:
        if path.is_symlink():raise ValueError('contract symlink')
        document=json.loads(path.read_text(encoding='utf-8'))
        if document.get('fixture_id')!=fixture_id or document.get('state')!='verified':raise ValueError('unverified contract')
        payload=document['payload']; rows=document['evidence']
        if not isinstance(payload,dict) or not isinstance(payload.get('input'),dict) or not isinstance(rows,list) or not rows:
            raise ValueError('Missing payload/evidence')
        for row in rows:
            rel=Path(row['path'])
            evidence=(root/rel).resolve()
            if rel.is_absolute() or '..' in rel.parts or not evidence.is_relative_to(root) or not evidence.is_file():
                raise ValueError('Evidence path outside package or missing')
            if digest(evidence)!=row['sha256']:raise ValueError('Evidence hash mismatch')
        return payload
    except (OSError,KeyError,ValueError,TypeError) as exc:
        raise CalibrationMissing(f'{fixture_id}: no validated calibration; {exc}') from exc


def execute_command(argv: Sequence[str], cwd: Path, log_dir: Path, timeout: float) -> dict:
    if isinstance(argv,(str,bytes)) or not argv or not all(isinstance(x,str) for x in argv):
        raise ValueError('Explicit argv list is required')
    if not math.isfinite(timeout) or timeout<=0:raise ValueError('Positive timeout required')
    cwd=Path(cwd).resolve(); log_dir=Path(log_dir).resolve()
    log_dir.mkdir(parents=True,exist_ok=False)
    receipt={'argv':list(argv),'cwd':str(cwd),'started_utc':datetime.now(timezone.utc).isoformat(),
             'returncode':None,'timed_out':False,'shell':False}
    start=time.monotonic(); failure=None; stdout=b''; stderr=b''
    try:
        proc=subprocess.Popen(list(argv),cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,start_new_session=True)
        try:
            stdout,stderr=proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            receipt['timed_out']=True
            os.killpg(proc.pid,signal.SIGKILL)
            stdout,stderr=proc.communicate()
            failure=ProcessTimeout(f'Process timeout; evidence: {log_dir}')
        receipt['returncode']=proc.returncode
        if failure is None and proc.returncode:
            failure=ProcessFailure(f'Process exited {proc.returncode}; evidence: {log_dir}')
    except OSError as exc:
        receipt['launch_error']=str(exc)
        failure=EnvironmentUnavailable(f'Cannot launch {argv[0]}: {exc}; evidence: {log_dir}')
    receipt['elapsed_s']=time.monotonic()-start
    (log_dir/'stdout.txt').write_bytes(stdout); (log_dir/'stderr.txt').write_bytes(stderr)
    receipt['stdout_sha256']=hashlib.sha256(stdout).hexdigest()
    receipt['stderr_sha256']=hashlib.sha256(stderr).hexdigest()
    write_json(log_dir/'command.json',receipt)
    if failure is not None:raise failure
    return receipt


@dataclass(frozen=True)
class BuildArtifact:
    executable: Path
    snapshot: Path
    receipt: Path
    workspace: str
    source_sha256: str
    binary_sha256: str
    backend: str


def verify_artifact(artifact: BuildArtifact) -> None:
    try:
        record=json.loads(artifact.receipt.read_text())
        command=json.loads(Path(record['command_receipt']).read_text())
        if record['state']!='built' or command['returncode']!=0 or command['timed_out']:
            raise ValueError('Build did not complete')
        if record['backend']!=artifact.backend or record['workspace']!=artifact.workspace:
            raise ValueError('Build workspace/backend mismatch')
        if digest(artifact.executable)!=artifact.binary_sha256 or record['binary_sha256']!=artifact.binary_sha256:
            raise ValueError('Binary changed after source build')
        if source_fingerprint(artifact.snapshot,artifact.workspace)['sha256']!=artifact.source_sha256:
            raise ValueError('Snapshot source changed after source build')
        if record['source_sha256']!=artifact.source_sha256:raise ValueError('Source receipt mismatch')
    except (OSError,KeyError,ValueError,TypeError) as exc:
        raise SourceMismatch(f'Build artifact is not verified: {exc}') from exc


class CaseRunner:
    def __init__(self, root: Path, work_root: Path, workspace: str='solutions', fpm: str='fpm'):
        self.root=Path(root).resolve();self.work_root=Path(work_root).resolve()
        self.workspace=workspace;self.fpm=fpm

    def build_fpm(self) -> BuildArtifact:
        verify_spec(self.root);validate_packages(self.root,self.workspace)
        chosen=os.environ.get('NS_FPM',self.fpm) if self.fpm=='fpm' else self.fpm
        fpm=shutil.which(chosen)
        if fpm is None and chosen=='fpm':
            local=Path.home()/'.local/bin/fpm'
            if local.is_file() and os.access(local,os.X_OK):fpm=str(local)
        if fpm is None:
            raise EnvironmentUnavailable(f'fpm not found ({chosen}); no compiler fallback is allowed')
        compiler=shutil.which(os.environ.get('FC','gfortran'))
        if compiler is None:raise EnvironmentUnavailable('Fortran compiler not found')
        build=self.work_root/('build-'+uuid4().hex);build.mkdir(parents=True)
        execute_command([fpm,'--version'],build,build/'fpm-version',10)
        version=(build/'fpm-version/stdout.txt').read_text()
        if not re.search(r'\b0\.13\.0\b',version):
            raise EnvironmentUnavailable(f'This runner requires fpm 0.13.0; actual version recorded at {build}/fpm-version')
        source=copy_sources(self.root,build/'source',self.workspace)
        package=build/'source'/self.workspace
        argv=[fpm,'build','--profile','debug','--compiler',compiler,'--flag',FLAGS,'--verbose']
        execute_command(argv,package,build/'compile',120)
        binaries=[p for p in (package/'build').rglob('case_driver') if p.is_file() and not p.is_symlink() and os.access(p,os.X_OK)]
        if len(binaries)!=1:raise SourceMismatch(f'Expected one freshly built case_driver, found {len(binaries)}')
        exe=binaries[0].resolve(); sha=digest(exe)
        receipt=build/'build.json'
        write_json(receipt,{'state':'built','backend':'fpm','workspace':self.workspace,
            'source_sha256':source['sha256'],'binary_sha256':sha,'binary':str(exe),
            'command_receipt':str(build/'compile/command.json'),'source_files':source['files'],
            'compiler':compiler,'flags':FLAGS,'fpm':fpm,'fpm_sha256':digest(Path(fpm)),
            'fpm_version':version.strip()})
        return BuildArtifact(exe,build/'source',receipt,self.workspace,source['sha256'],sha,'fpm')

    def run_case(self, case_id: str, parameters: Mapping[str,Any]) -> dict:
        return self.run_artifact(self.build_fpm(),case_id,parameters)

    def run_artifact(self, artifact: BuildArtifact, case_id: str, parameters: Mapping[str,Any], *, timeout: float=10) -> dict:
        verify_spec(self.root)
        if artifact.workspace!=self.workspace:raise SourceMismatch('Selected workspace differs from built artifact')
        verify_artifact(artifact)
        if not isinstance(case_id,str) or not NAME.fullmatch(case_id):raise ValueError('Invalid case ID')
        registry=json.loads((self.root/'manifests/cases.json').read_text())
        schema=CASE_SCHEMAS.get(case_id)
        if schema is None and parameters:
            raise ValueError('Unknown case cannot accept undeclared parameters')
        if schema is not None and registry['cases'].get(case_id,{}).get('parameters')!=schema:
            raise SourceMismatch('Case registry/schema differs from the declared contract')
        run=self.work_root/('run-'+uuid4().hex);run.mkdir(parents=True)
        token=uuid4().hex;inp=run/'input.nml'
        state={'status':'partial','case_id':case_id,'workspace':self.workspace,'build_backend':artifact.backend,
            'schema_version':2,'request_token':token,'run_id':token,
            'source_sha256':artifact.source_sha256,'spec_sha256':SPEC_SHA256,
            'binary_sha256':artifact.binary_sha256,'started_utc':datetime.now(timezone.utc).isoformat(),
            'build_receipt':str(artifact.receipt),
            'constants_source_sha256':digest(artifact.snapshot/'support/src/ns_physical_constants.f90'),
            'contract':{'units':('dimensionless_analytic_fixture' if case_id in {'root_sqrt2','root_probe','ode','rk4_exp','adaptive_exp','control','dopri_pair'} else 'SI'),'model':('analytic_numerics' if case_id in {'adaptive_exp', 'dopri_pair', 'rk4_exp', 'control', 'ode', 'root_probe', 'root_sqrt2'} else 'provided_support_only'),'eos':'not_applicable',
                        'closure':'not_applicable','interface':'not_applicable','surface':'not_applicable',
                        'interpolation':'not_applicable'}}
        if case_id == 'poly_state':
            state['contract'].update(model='cold_polytrope',eos='analytic_anchor_polytrope_v1',
                                     closure='barotropic_equilibrium_derivative',interface='smooth_single_branch')
        if case_id in {'two_poly_knot','direct_join','maxwell_match','mu_bridge'}:
            state['contract'].update(model='analytic_two_phase_polytropes',eos='analytic_anchor_polytrope_v1',
                                     closure='cold_zero_temperature_thermodynamics',
                                     interface='equilibrium_matching_only_no_mode_transfer',
                                     interpolation=('cubic_Hermite_P_of_mu' if case_id=='mu_bridge' else 'analytic_P_of_mu'))
        if case_id in {'tov_baseline','tov_rhs','target_mass'}:
            state['contract'].update(model='static_spherical_TOV',
                eos=parameters.get('eos_kind','polytrope'),
                closure=('barotropic_equilibrium_derivative' if parameters.get('eos_kind','polytrope')=='polytrope'
                         else 'formal_incompressible_benchmark_no_mode_response'),
                interface='smooth_only', interpolation='analytic_EOS_no_table',
                surface=parameters.get('surface_policy','finite_cutoff_schwarzschild_time'),
                integration_coordinate='log_pressure',
                lapse='effective_Schwarzschild_anchor_at_cutoff_not_exact_vacuum_junction')
        if case_id in {'table_diagnostics','table_state','real_table_admission','real_table_star','real_table_star_split'}:
            state['contract'].update(model='cold_tabulated_EOS',eos='normalized_table',
                closure='single_cubic_Hermite_P_of_mu_with_raw_n_slope',
                interface='native_rows_no_smoothing_no_extrapolation',
                interpolation='cubic_Hermite_P_of_mu_raw_n_v1')
            if case_id in {'real_table_star','real_table_star_split'}:
                state['contract'].update(model='static_spherical_TOV_from_cold_tabulated_EOS',
                    surface='finite_cutoff_schwarzschild_time',integration_coordinate='log_pressure',
                    lapse='effective_Schwarzschild_anchor_at_cutoff_not_exact_vacuum_junction')
        if case_id in {'structure_diagnostics','mode_difference','paired_refinement','controlled_table_variant','controlled_table_target_mass'}:
            state['contract'].update(model='postprocessed_stellar_diagnostics',
                eos='background_supplied_by_prior_verified_run',closure='read_only_no_equation_change',
                interface='no_new_interface_physics',interpolation='background_rows_trapezoid_v1')
        if case_id in {'smooth_spectrum','mode_refine','mode_trial'}:
            state['contract'].update(model='relativistic_Cowling_l2',
                eos='supplied_smooth_background',closure='barotropic_equilibrium_derivative',
                interface='smooth_only',surface='finite_cutoff_schwarzschild_time',
                interpolation='linear_q_v1',integration_coordinate='log_pressure',
                integrator=parameters.get('method','adaptive'),
                normalization='scaled_WV_L2_dx_v1',exact_vacuum_surface=False)
        write_json(run/'run.json',state)
        try:
            if case_id in {'smooth_spectrum','mode_refine','mode_trial','structure_diagnostics'}:
                bg=Path(parameters.get('background_path',''))
                if bg.is_symlink():raise SourceMismatch('Background symlink is not a regular input')
                if bg.is_file():
                    state['background_input']={'path':str(bg.resolve()),'sha256':digest(bg),
                                               'unchanged_after_run':False}
            write_input(inp,token,parameters, schema={} if schema is None else schema)
            state['input_sha256']=digest(inp)
            write_json(run/'run.json',state)
            execute_command([str(artifact.executable),case_id,str(inp)],run,run/'process',timeout)
            raw=(run/'process/stdout.txt').read_text(encoding='utf-8')
            result=parse_records(raw)
            if result.get('workspace')!=self.workspace or result.get('request_token')!=token or result.get('case_id')!=case_id:
                raise OutputMalformed('Output does not belong to this request/workspace')
            if result.get('protocol_version')!=1 or result.get('status') not in SUPPORT_STATUSES:
                raise OutputMalformed('Missing or invalid protocol/status')
            if digest(inp)!=state['input_sha256']:raise SourceMismatch('Input changed during execution')
            if 'background_input' in state:
                try:
                    same=digest(Path(state['background_input']['path']))==state['background_input']['sha256']
                except OSError:
                    same=False
                if not same:raise SourceMismatch('Background changed during execution')
                state['background_input']['unchanged_after_run']=True
            write_json(run/'result.json',result)
            state.update(status='complete' if result['status']=='ok' else 'rejected',domain_status=result['status'],
                         stdout_sha256=digest(run/'process/stdout.txt'),result_sha256=digest(run/'result.json'),
                         command_sha256=digest(run/'process/command.json'),finished_utc=datetime.now(timezone.utc).isoformat())
            write_json(run/'run.json',state)
            result.update(built_from_source=True,build_backend=artifact.backend,source_sha256=artifact.source_sha256,
                          binary_sha256=artifact.binary_sha256,
                          _evidence={'build_receipt':str(artifact.receipt),'run_receipt':str(run/'run.json'),
                                     'command_receipt':str(run/'process/command.json')})
            return result
        except (BootstrapError,ValueError,UnicodeError,OSError) as exc:
            state.update(status='partial',error_type=type(exc).__name__,detail=str(exc))
            write_json(run/'run.json',state)
            raise
