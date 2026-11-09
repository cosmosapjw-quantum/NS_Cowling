"""Python orchestration for the native Fortran TOV and Cowling executables."""
from __future__ import annotations
import csv,hashlib,json,re,shutil,subprocess,tempfile
from pathlib import Path
import numpy as np
from compose_eos import EOSData,write_compose,normalize_compose
from compose_eos.native import project_root,NativeExecutionError


def _compile_solver():
    root=project_root();source=root/'oscillation_source'
    support=sorted((source/'support/src').rglob('*.f90'))
    solutions=sorted((source/'solutions/src').rglob('*.f90'))
    files=support+solutions
    digest=hashlib.sha256()
    for p in files:
        digest.update(p.relative_to(source).as_posix().encode());digest.update(p.read_bytes())
    for name in ('case_driver.f90','cell_solve.f90'):
        digest.update((source/'solutions/app'/name).read_bytes())
    key=digest.hexdigest()[:20]
    build=Path.home()/'.cache'/'compose_eos_fortran'/'solver'/key
    case_driver=build/'case_driver';cell_solve=build/'cell_solve'
    if case_driver.is_file() and cell_solve.is_file():return case_driver,cell_solve
    fc=shutil.which('gfortran')
    if not fc:raise NativeExecutionError('gfortran is required for the TOV/Cowling backend')
    mod=build/'mod';obj=build/'obj';mod.mkdir(parents=True,exist_ok=True);obj.mkdir(exist_ok=True)
    module_owner={};uses={}
    for p in files:
        text=p.read_text(errors='replace')
        names=re.findall(r'^\s*module\s+(?!procedure\b)([a-z][a-z0-9_]*)',text,re.I|re.M)
        if names:
            for name in names:module_owner[name.lower()]=p
            uses[p]=set(x.lower() for x in re.findall(
                r'^\s*use\s*(?:,\s*(?:intrinsic|non_intrinsic)\s*)?(?:::)?\s*([a-z][a-z0-9_]*)',
                text,re.I|re.M))
    remaining=set(uses);done=set();order=[]
    while remaining:
        for p in sorted(remaining):
            deps={module_owner[x] for x in uses[p] if x in module_owner and module_owner[x]!=p}
            if deps.issubset(done):
                order.append(p);done.add(p);remaining.remove(p);break
        else:
            raise NativeExecutionError("Fortran module dependency cycle detected.")
    flags=['-std=f2018','-O2','-ffree-line-length-none','-fimplicit-none','-fcheck=all']
    objects=[]
    for p in order:
        o=obj/(p.stem+'.o')
        cmd=[fc,*flags,'-J',str(mod),'-I',str(mod),'-c',str(p),'-o',str(o)]
        r=subprocess.run(cmd,capture_output=True,text=True)
        if r.returncode:raise NativeExecutionError(f"Compile failed at {p}:\n{r.stderr}")
        objects.append(str(o))
    for app_name,target in (('case_driver.f90',case_driver),('cell_solve.f90',cell_solve)):
        appobj=obj/(target.name+'.o')
        r=subprocess.run([fc,*flags,'-J',str(mod),'-I',str(mod),'-c',
            str(source/'solutions/app'/app_name),'-o',str(appobj)],capture_output=True,text=True)
        if r.returncode:raise NativeExecutionError(f"Application compile failed {app_name}:\n{r.stderr}")
        r=subprocess.run([fc,*flags,*objects,str(appobj),'-o',str(target)],capture_output=True,text=True)
        if r.returncode:raise NativeExecutionError(f"Application link failed {app_name}:\n{r.stderr}")
    return case_driver,cell_solve



def _v5_source_paths():
    """Resolve the exact source list recorded by the v5 release manifest."""
    root=project_root()
    source=root/'oscillation_source'
    manifest=root/'oscillation_release_v1'/'workflow'/'compile_manifest.json'
    data=json.loads(manifest.read_text())
    marker='Neutron_Star_Oscillations/'
    if data.get('modules')!=35 or len(data.get('sources',[]))!=35:
        raise NativeExecutionError("Unexpected v5 compile-manifest module/source count")
    if len(set(data['sources']))!=35:
        raise NativeExecutionError("Duplicate source path in v5 compile manifest")
    paths=[]
    for item in data['sources']:
        if marker not in item:
            raise NativeExecutionError(f"Unexpected v5 source-manifest entry: {item}")
        rel=item.split(marker,1)[1]
        path=source/rel
        if not path.is_file():
            raise NativeExecutionError(f"Frozen v5 source is missing: {path}")
        paths.append(path)
    app=source/'calculation'/'app'/'cell_solve.f90'
    if not app.is_file():
        raise NativeExecutionError(f"Frozen v5 application is missing: {app}")
    return source,manifest,paths,app


def _compile_v5_cell_solve():
    """Compile the v5 calculation/methods/support source family, not solutions/."""
    source,manifest,paths,app=_v5_source_paths()
    digest=hashlib.sha256()
    digest.update(manifest.read_bytes())
    for path in [*paths,app]:
        digest.update(path.relative_to(source).as_posix().encode())
        digest.update(path.read_bytes())
    key=digest.hexdigest()[:24]
    build=Path.home()/'.cache'/'compose_eos_fortran'/'v5-release'/key
    executable=build/'cell_solve'
    if executable.is_file():
        return executable,key
    fc=shutil.which('gfortran')
    if not fc:
        raise NativeExecutionError("gfortran is required for the frozen v5 backend")
    mod=build/'mod';obj=build/'obj'
    mod.mkdir(parents=True,exist_ok=True);obj.mkdir(exist_ok=True)
    flags=['-std=f2018','-O2','-ffree-line-length-none','-fimplicit-none']
    objects=[]
    for i,path in enumerate(paths):
        target=obj/f"{i:02d}_{path.stem}.o"
        result=subprocess.run([fc,*flags,'-J',str(mod),'-I',str(mod),'-c',str(path),'-o',str(target)],
                              capture_output=True,text=True)
        if result.returncode:
            raise NativeExecutionError(f"v5 module compile failed at {path}:\n{result.stderr}")
        objects.append(str(target))
    appobj=obj/'cell_solve_app.o'
    result=subprocess.run([fc,*flags,'-J',str(mod),'-I',str(mod),'-c',str(app),'-o',str(appobj)],
                          capture_output=True,text=True)
    if result.returncode:
        raise NativeExecutionError(f"v5 cell_solve application compile failed:\n{result.stderr}")
    result=subprocess.run([fc,*objects,str(appobj),'-o',str(executable)],capture_output=True,text=True)
    if result.returncode:
        raise NativeExecutionError(f"v5 cell_solve link failed:\n{result.stderr}")
    return executable,key


def _run_v5_cell_solve(eos,pc_pa,*,task,target_mass_kg=0.0,pc_bracket_pa=None,
                       frequency_hz=2500.0,frequency_lo_hz=(2000.0,6500.0),
                       frequency_hi_hz=(3000.0,8000.0),p_cut_pa=None,
                       r0_m=10.0,r_limit_m=1e8,rtol=1e-9,mode_rtol=1e-9,
                       step_cap=0.05,mass_rtol=1e-10,residual_tol=1e-6,
                       match_index=72,center_order=2,background_output_path=None):
    executable,build_key=_compile_v5_cell_solve()
    if task not in ('star','shot','spectrum'):
        raise ValueError("v5 task must be star, shot or spectrum")
    values=[pc_pa,target_mass_kg,frequency_hz,r0_m,r_limit_m,rtol,mode_rtol,step_cap,mass_rtol,residual_tol]
    if not np.all(np.isfinite(values)) or pc_pa<=0 or r0_m<=0 or r_limit_m<=r0_m:
        raise ValueError("v5 input pressure/radius/tolerance must be finite and admissible")
    if min(rtol,mode_rtol,step_cap,mass_rtol,residual_tol)<=0:
        raise ValueError("v5 tolerances must be positive")
    if center_order not in (0,2):
        raise ValueError("center_order must be 0 or 2")
    if p_cut_pa is None:
        p_cut_pa=_surface_pressure(eos)
    if not np.isfinite(p_cut_pa) or p_cut_pa<=0 or pc_pa<=p_cut_pa:
        raise ValueError("require finite 0 < surface pressure cutoff < central pressure")
    if pc_bracket_pa is None:
        pc_lo_pa,pc_hi_pa=0.65*pc_pa,1.45*pc_pa
    else:
        if len(pc_bracket_pa)!=2:
            raise ValueError("pc_bracket_pa must have two entries")
        pc_lo_pa,pc_hi_pa=map(float,pc_bracket_pa)
        if not np.all(np.isfinite([pc_lo_pa,pc_hi_pa])) or pc_lo_pa<=0 or pc_hi_pa<=pc_lo_pa:
            raise ValueError("invalid explicit central-pressure bracket")
    lo=list(map(float,frequency_lo_hz));hi=list(map(float,frequency_hi_hz))
    if len(lo)!=2 or len(hi)!=2 or not np.all(np.isfinite(lo+hi)):
        raise ValueError("frequency brackets must be two finite values for f and p1")
    with tempfile.TemporaryDirectory(prefix='v5-cell-solve-') as td0:
        td=Path(td0);table=_normalised_file(eos,td)
        prefix=td/'v5_result';nml=td/'v5_request.nml'
        nml.write_text(
            "&precision\n"
            f" table_path='{table}'\n task='{task}'\n method='dopri'\n shooting='matched'\n"
            f" match_index={int(match_index)}\n split_index=69\n alpha=0.0\n"
            f" pc_pa={float(pc_pa):.17e}\n p_cut_pa={float(p_cut_pa):.17e}\n"
            f" r0_m={float(r0_m):.17e}\n r_limit_m={float(r_limit_m):.17e}\n"
            f" target_mass_kg={float(target_mass_kg):.17e}\n"
            f" pc_lo_pa={pc_lo_pa:.17e}\n pc_hi_pa={pc_hi_pa:.17e}\n"
            f" rtol={float(rtol):.17e}\n mode_rtol={float(mode_rtol):.17e}\n"
            f" step_cap={float(step_cap):.17e}\n mass_rtol={float(mass_rtol):.17e}\n"
            f" frequency_lo={lo[0]:.17e},{lo[1]:.17e}\n"
            f" frequency_hi={hi[0]:.17e},{hi[1]:.17e}\n"
            f" frequency_hz={float(frequency_hz):.17e}\n"
            " root_atol_hz=1.0e-8\n"
            f" residual_tol={float(residual_tol):.17e}\n nmode=2\n center_order={center_order}\n/\n"
        )
        run=subprocess.run([str(executable),str(nml),str(prefix)],capture_output=True,text=True)
        result_path=prefix.with_suffix('.json')
        if not result_path.is_file():
            raise NativeExecutionError(f"v5 cell_solve produced no JSON: {run.stdout}\n{run.stderr}")
        result=json.loads(result_path.read_text())
        if run.returncode or result.get('status')!='ok':
            raise NativeExecutionError(f"v5 cell_solve {task} failed: {result}; {run.stderr}")
        result['backend']='v5_release_calculation_methods_support'
        result['backend_build_key']=build_key
        result['center_order']=center_order
        result['r_limit_m']=float(r_limit_m)
        if background_output_path is not None:
            if task != 'star':
                raise ValueError("background export is available only for TOV task='star'")
            source_profile=Path(str(prefix)+'_background.csv')
            if not source_profile.is_file():
                raise NativeExecutionError("v5 TOV run did not emit its requested background profile")
            destination=Path(background_output_path)
            if destination.exists():
                raise FileExistsError(destination)
            destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source_profile,destination)
            with destination.open(newline='') as stream:
                reader=csv.DictReader(stream)
                expected=('cell','t','r_m','m_kg','p_pa','energy_j_m3','chi','nu','w_over_r3','v_over_r2')
                if tuple(reader.fieldnames or ()) != expected:
                    raise NativeExecutionError("v5 background profile has an unexpected schema")
                profile_rows=list(reader)
            if len(profile_rows)<2:
                raise NativeExecutionError("v5 background profile has fewer than two rows")
            profile_sha=hashlib.sha256(destination.read_bytes()).hexdigest()
            result['background_profile']={
                'path':str(destination),
                'sha256':profile_sha,
                'row_count':len(profile_rows),
                'columns':list(expected),
            }
        if task=='spectrum':
            nodes={}
            for j,mode in ((1,'f'),(2,'p1')):
                nodes[mode]={}
                for region in ('inner','outer'):
                    profile=Path(str(prefix)+f'_mode{j}_{region}.csv')
                    nodes[mode][region]={
                        'W':_node_count(profile,'w_over_r3'),
                        'V':_node_count(profile,'v_over_r2')}
            result['node_counts']=nodes
        return result


def solve_v5_tov_star(eos:EOSData | str | Path,pc_pa:float,*,p_cut_pa:float|None=None,
                      r0_m:float=10.0,r_limit_m:float=1e8,rtol:float=1e-9,
                      step_cap:float=0.05,center_order:int=2,
                      background_output_path:str | Path | None=None)->dict:
    """TOV-only run using the exact frozen v5 source family and center expansion."""
    return _run_v5_cell_solve(eos,pc_pa,task='star',p_cut_pa=p_cut_pa,
        r0_m=r0_m,r_limit_m=r_limit_m,rtol=rtol,mode_rtol=rtol,
        step_cap=step_cap,center_order=center_order,
        background_output_path=background_output_path)


def solve_v5_cowling_shot(eos:EOSData | str | Path,pc_pa:float,frequency_hz:float,*,
                          p_cut_pa:float|None=None,r0_m:float=10.0,r_limit_m:float=1e8,
                          rtol:float=1e-9,mode_rtol:float=1e-9,step_cap:float=0.05,
                          match_index:int=72,center_order:int=2)->dict:
    """Return one matched Cowling determinant value at a fixed frequency."""
    return _run_v5_cell_solve(eos,pc_pa,task='shot',frequency_hz=frequency_hz,
        p_cut_pa=p_cut_pa,r0_m=r0_m,r_limit_m=r_limit_m,rtol=rtol,
        mode_rtol=mode_rtol,step_cap=step_cap,match_index=match_index,
        center_order=center_order)


def solve_v5_cowling_spectrum(eos:EOSData | str | Path,pc_pa:float,frequency_lo_hz,frequency_hi_hz,*,
                              target_mass_kg:float=2.644606257493485e30,p_cut_pa:float|None=None,
                              pc_bracket_pa:tuple[float,float]|None=None,match_index:int=72,
                              r0_m:float=10.0,r_limit_m:float=1e8,rtol:float=1e-9,
                              mode_rtol:float=1e-9,step_cap:float=0.05,
                              mass_rtol:float=1e-10,residual_tol:float=1e-6,
                              center_order:int=2)->dict:
    """Fresh Cowling roots from the exact v5 source family and center expansion."""
    return _run_v5_cell_solve(eos,pc_pa,task='spectrum',target_mass_kg=target_mass_kg,
        pc_bracket_pa=pc_bracket_pa,frequency_lo_hz=frequency_lo_hz,
        frequency_hi_hz=frequency_hi_hz,p_cut_pa=p_cut_pa,r0_m=r0_m,
        r_limit_m=r_limit_m,rtol=rtol,mode_rtol=mode_rtol,
        step_cap=step_cap,mass_rtol=mass_rtol,residual_tol=residual_tol,
        match_index=match_index,center_order=center_order)

def _normalised_file(eos:EOSData | str | Path,td:Path)->Path:
    if isinstance(eos,(str,Path)):
        path=Path(eos).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    compose=td/'compose';compose.mkdir()
    write_compose(compose,eos)
    return normalize_compose(compose,td/'normalized_eos.csv')


def _surface_pressure(eos:EOSData | str | Path)->float:
    if isinstance(eos,(str,Path)):
        with Path(eos).open(newline='') as stream:
            row=next(csv.DictReader(stream))
        if 'pressure_pa' in row:return float(row['pressure_pa'])
        if 'pressure_mev_fm3' in row:return float(row['pressure_mev_fm3'])*1.602176634e32
        raise ValueError('normalized EOS CSV has no recognized pressure column')
    return float(np.min(eos.p))*1.602176634e32


def _parse_tsv(path_text:str)->dict:
    rows=[x.split('\t') for x in path_text.splitlines() if x.strip()]
    if not rows:return {}
    result={}
    for row in rows[1:]:
        if len(row)!=3:continue
        key,kind,value=row
        if value in ('true','false'):result[key]=value=='true'
        elif kind in ('int','real'):
            try:result[key]=int(value) if kind=='int' else float(value)
            except ValueError:result[key]=value
        else:result[key]=value
    return result


def solve_tov_star(eos:EOSData | str | Path,pc_pa:float,*,p_cut_pa:float|None=None,rtol:float=1e-8,
                   max_steps:int=30000,emit_profile:bool=False)->dict:
    case_driver,_=_compile_solver()
    with tempfile.TemporaryDirectory(prefix='compose-tov-') as td0:
        td=Path(td0);table=_normalised_file(eos,td)
        if p_cut_pa is None:p_cut_pa=_surface_pressure(eos)
        nml=td/'request.nml'
        nml.write_text(
          "&request request_token='python-frontend' /\n&parameters\n"
          f" table_path='{table}',pc_pa={pc_pa:.17e},p_cut_pa={p_cut_pa:.17e},"
          f"r0_m=1.0,r_limit_m=1.0e6,rtol={rtol:.17e},atol_r_m=1.0e-6,"
          f"atol_mass_kg=1.0e14,atol_nu=1.0e-12,dq_initial=1.0e-5,dq_min=1.0e-16,"
          f"dq_max=0.2,nu_start=0.0,max_steps={max_steps},emit_profile={'.true.' if emit_profile else '.false.'}\n/\n"
        )
        run=subprocess.run([str(case_driver),'real_table_star',str(nml)],capture_output=True,text=True)
        result=_parse_tsv(run.stdout)
        if run.returncode or result.get('status')!='ok':
            raise NativeExecutionError(f"TOV Fortran solve failed: {run.stdout}\n{run.stderr}")
        result['backend']='Fortran';return result


def solve_cowling_tov_star(eos:EOSData | str | Path,pc_pa:float,*,p_cut_pa:float|None=None,
                           r0_m:float=10.0,r_limit_m:float=1e8,rtol:float=1e-9,
                           step_cap:float=0.05)->dict:
    """TOV-only background using the exact cell_solve backend of Cowling spectra.

    This is the comparison-preserving backend for the historical spectrum cohort.
    The general real_table_star case_driver API remains available independently.
    """
    _,cell_solve=_compile_solver()
    if not np.all(np.isfinite([pc_pa,r0_m,r_limit_m,rtol,step_cap])) or pc_pa<=0:
        raise ValueError("Pc, radii, tolerances and step cap must be finite; Pc must be positive")
    if r0_m<=0 or r_limit_m<=r0_m or rtol<=0 or step_cap<=0:
        raise ValueError("require 0 < r0_m < r_limit_m and positive tolerances")
    with tempfile.TemporaryDirectory(prefix='compose-cell-tov-') as td0:
        td=Path(td0);table=_normalised_file(eos,td)
        if p_cut_pa is None:p_cut_pa=_surface_pressure(eos)
        if not np.isfinite(p_cut_pa) or p_cut_pa<=0 or pc_pa<=p_cut_pa:
            raise ValueError("require finite 0 < surface cutoff < central pressure")
        prefix=td/'background'
        nml=td/'background.nml'
        nml.write_text(
          "&precision\n"
          f" table_path='{table}'\n task='star'\n method='dopri'\n shooting='matched'\n"
          f" match_index=72\n split_index=69\n alpha=0.0\n pc_pa={pc_pa:.17e}\n"
          f" p_cut_pa={p_cut_pa:.17e}\n r0_m={r0_m:.17e}\n r_limit_m={r_limit_m:.17e}\n"
          " target_mass_kg=0.0\n pc_lo_pa=1.0e30\n pc_hi_pa=2.0e30\n"
          f" rtol={rtol:.17e}\n mode_rtol={rtol:.17e}\n step_cap={step_cap:.17e}\n"
          " mass_rtol=1.0e-10\n frequency_lo=2000.0,6500.0\n"
          " frequency_hi=3000.0,8000.0\n frequency_hz=2500.0\n"
          " root_atol_hz=1.0e-8\n residual_tol=1.0e-6\n nmode=2\n/\n"
        )
        run=subprocess.run([str(cell_solve),str(nml),str(prefix)],capture_output=True,text=True)
        result_path=prefix.with_suffix('.json')
        if not result_path.is_file():
            raise NativeExecutionError(f"cell_solve TOV-only run produced no JSON: {run.stdout}\n{run.stderr}")
        result=json.loads(result_path.read_text())
        if run.returncode or result.get('status')!='ok':
            raise NativeExecutionError(f"cell_solve TOV-only result invalid: {result}; {run.stderr}")
        result['backend']='Fortran cell_solve task=star'
        return result


def solve_tov_sequence(eos:EOSData | str | Path,central_pressures_pa,**kwargs)->list[dict]:
    # Python orchestrates the scan; each stellar integration and EOS query is Fortran.
    return [solve_tov_star(eos,float(pc),**kwargs) for pc in central_pressures_pa]


def _node_count(path:Path,column:str)->int|None:
    if not path.exists():return None
    import csv
    with path.open(newline='') as stream:rows=list(csv.DictReader(stream))
    if len(rows)<2:return 0
    values=np.asarray([float(row[column]) for row in rows])
    values=values[values!=0]
    return int(np.count_nonzero(values[1:]*values[:-1]<0))


def solve_cowling_spectrum(eos:EOSData | str | Path,pc_pa:float,frequency_lo_hz,frequency_hi_hz,*,
                           target_mass_kg:float=2.644606257493485e30,p_cut_pa:float|None=None,
                           match_index:int=72,rtol:float=1e-9,mode_rtol:float=1e-9,
                           step_cap:float=0.05,mass_rtol:float=1e-10,
                           residual_tol:float=1e-6,
                           pc_bracket_pa:tuple[float,float]|None=None,
                           r_limit_m:float=1e6)->dict:
    _,cell_solve=_compile_solver()
    with tempfile.TemporaryDirectory(prefix='compose-cowling-') as td0:
        td=Path(td0);table=_normalised_file(eos,td)
        if p_cut_pa is None:p_cut_pa=_surface_pressure(eos)
        lo=list(map(float,frequency_lo_hz));hi=list(map(float,frequency_hi_hz))
        if len(lo)!=2 or len(hi)!=2:raise ValueError("frequency brackets must be [f,p1]")
        if pc_bracket_pa is None:
            pc_lo,pc_hi=0.65*pc_pa,1.45*pc_pa
        else:
            if len(pc_bracket_pa)!=2:
                raise ValueError("pc_bracket_pa must contain exactly two pressures")
            pc_lo,pc_hi=map(float,pc_bracket_pa)
            if not np.all(np.isfinite([pc_lo,pc_hi])) or pc_lo<=0 or pc_hi<=pc_lo:
                raise ValueError("pc_bracket_pa must be finite, positive and strictly increasing")
            pc_pa=float(np.sqrt(pc_lo*pc_hi))
        nml=td/'spectrum.nml';prefix=td/'spectrum'
        nml.write_text(
          "&precision\n"
          f" table_path='{table}'\n task='spectrum'\n method='dopri'\n shooting='matched'\n"
          f" match_index={match_index}\n split_index=69\n alpha=0.0\n pc_pa={pc_pa:.17e}\n"
          f" p_cut_pa={p_cut_pa:.17e}\n r0_m=10.0\n r_limit_m={r_limit_m:.17e}\n target_mass_kg={target_mass_kg:.17e}\n"
          f" pc_lo_pa={pc_lo:.17e}\n pc_hi_pa={pc_hi:.17e}\n"
          f" rtol={rtol:.17e}\n mode_rtol={mode_rtol:.17e}\n step_cap={step_cap:.17e}\n"
          f" mass_rtol={mass_rtol:.17e}\n frequency_lo={lo[0]:.17e},{lo[1]:.17e}\n"
          f" frequency_hi={hi[0]:.17e},{hi[1]:.17e}\n root_atol_hz=1.0e-8\n"
          f" residual_tol={residual_tol:.17e}\n nmode=2\n/\n"
        )
        run=subprocess.run([str(cell_solve),str(nml),str(prefix)],capture_output=True,text=True)
        result_path=prefix.with_suffix('.json')
        if not result_path.exists():
            raise NativeExecutionError(f"Cowling Fortran solve failed: {run.stdout}\n{run.stderr}")
        result=json.loads(result_path.read_text())
        if run.returncode or result.get('status')!='ok':
            raise NativeExecutionError(f"Cowling Fortran result invalid: {result}; {run.stderr}")
        nodes={}
        for j,mode in ((1,'f'),(2,'p1')):
            nodes[mode]={}
            for region in ('inner','outer'):
                path=Path(str(prefix)+f'_mode{j}_{region}.csv')
                nodes[mode][region]={
                    'W':_node_count(path,'w_over_r3'),'V':_node_count(path,'v_over_r2')}
        result['node_counts']=nodes;result['backend']='Fortran'
        return result
