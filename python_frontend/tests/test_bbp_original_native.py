"""Independent source and thermodynamic checks for the native BBP functional.

Python only compiles, invokes, reads results, and provides verification oracles.
No production EOS construction or minimization is performed in Python.
"""
from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PARAMS = ROOT / "data/eos_original_papers_20261123/research_bbp/BBP1971_ORIGINAL_PARAMETERS.dat"

DRIVER = r"""
program driver
 use compose_eos_types,only:wp,eos_data_t
 use compose_eos_io,only:write_compose,read_compose
 use compose_eos_bbp_original
 implicit none
 type(eos_data_t)::e,reloaded
 type(bbp_original_state_t)::s,a,b
 character(1024)::cmd,arg,src,out
 character(256)::status
 real(wp)::k,x,n,z(3),w,dk,dx,wn,mun,mup,f,g(4)
 integer::m,np,phase
 logical::ok
 call get_command_argument(1,cmd)
 select case(trim(cmd))
 case('bulk')
 call get_command_argument(2,arg);read(arg,*)k
 call get_command_argument(3,arg);read(arg,*)x
 call bbp_original_bulk(k,x,w,dk,dx,wn,mun,mup,ok)
 print '(l1,6es25.16e3)',ok,w,dk,dx,wn,mun,mup
 case('point')
 call get_command_argument(2,arg);read(arg,*)n
 call get_command_argument(3,arg);read(arg,*)z(1)
 call get_command_argument(4,arg);read(arg,*)z(2)
 call get_command_argument(5,arg);read(arg,*)z(3)
 call get_command_argument(6,arg);read(arg,*)m
 call bbp_original_point(n,z,m==1,s,f,g,ok)
 print '(l1,11es25.16e3)',ok,f,g,s%a,s%u,s%k,s%x,s%kg,s%pressure
 case('solve')
 call get_command_argument(2,arg);read(arg,*)n
 call get_command_argument(3,arg);read(arg,*)m
 call get_command_argument(4,arg);read(arg,*)phase
 call bbp_original_solve(n,m==1,s,status,phase_mode=phase)
 print '(a)',trim(status)
 print '(17es25.16e3)',s%k,s%x,s%kg,s%a,s%z,s%pressure,s%mu,s%mue, &
 s%energy,s%gradient,s%true_gradient,s%hessian_min,s%uniform_delta
 case('maxwell')
 call bbp_original_maxwell(a,b,status)
 print '(a)',trim(status)
 print '(8es25.16e3)',a%nb,b%nb,a%pressure,b%pressure,a%mu,b%mu,a%energy,b%energy
 case('generate')
 call get_command_argument(2,arg);read(arg,*)m
 call get_command_argument(3,arg);read(arg,*)np
 call get_command_argument(4,src);call get_command_argument(5,out)
 select case(m)
 case(1)
 call generate_bbp_original_published_native(trim(src),np,e,trim(out),status)
 case(2)
 call generate_bbp_original_variational_native(trim(src),np,e,trim(out),status)
 case(3)
 call generate_bbp_variational_maxwell_native(trim(src),np,e,trim(out),status)
 end select
 print '(a)',trim(status)
 if(trim(status)=='ok')then
 if(any(e%q(5,:)/=0._wp))error stop 'beta-equilibrated generated lepton potential is not zero'
 call write_compose(trim(out),e,status,allow_source_physical_violations=.true.)
 if(trim(status)/='ok')error stop 'write_compose failed'
 call read_compose(trim(out),reloaded,status)
 if(trim(status)/='ok')error stop 'read_compose failed'
 if(any(reloaded%q(5,:)/=0._wp))error stop 'serialized beta-equilibrium lepton potential is not zero'
 if(any(reloaded%q(4,:)>=0._wp))error stop 'charge chemical potential sign is wrong'
 print '(i0,2es25.16e3)',e%n,maxval(abs(e%yn+e%yp+e%yh-1)), &
 maxval(abs(e%yp+e%zbar/e%abar*e%yh-e%ye-e%ymu))
 endif
 end select
end program
"""


@pytest.fixture(scope="module")
def native(tmp_path_factory):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran unavailable")
    build = tmp_path_factory.mktemp("bbp_original_native")
    driver = build / "driver.f90"
    driver.write_text(DRIVER)
    exe = build / "driver"
    subprocess.run(
        [compiler, "-O2", "-g", "-ffree-line-length-none", "-fcheck=all",
         "-ffpe-trap=invalid,zero,overflow", "-J", str(build),
         str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
         str(ROOT / "fortran_eos/src/compose_eos_io.f90"),
         str(ROOT / "fortran_eos/src/compose_eos_bbp_original.f90"),
         str(driver), "-o", str(exe)], cwd=build, check=True, capture_output=True, text=True,
    )

    def run(*args):
        done = subprocess.run([str(exe), *map(str, args)], text=True, capture_output=True, check=True)
        return done.stdout.strip().splitlines()

    return run


def bulk(native, k, x):
    tokens = native("bulk", k, x)[0].split()
    assert tokens[0] == "T"
    return list(map(float, tokens[1:]))


def point(native, n, z, published=False):
    tokens = native("point", n, *z, int(published))[0].split()
    assert tokens[0] == "T"
    return list(map(float, tokens[1:]))


def generate(native, directory, model, points=257, source=PARAMS):
    directory.mkdir()
    lines = native("generate", model, points, source, directory)
    assert lines[0] == "ok", lines
    summary = list(map(float, lines[1].split()))
    assert summary[1] < 2e-14
    assert summary[2] < 2e-14
    diagnostics = json.loads((directory / "bbp_original_diagnostics.json").read_text())
    with (directory / "bbp_original_states.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert all(None not in row for row in rows), "CSV columns must match header exactly"
    assert len(rows) == int(summary[0]) == diagnostics["rows"]
    return diagnostics, rows


def test_original_bulk_limits_symmetry_and_chemical_potentials(native):
    for k in (.15, .65, 1.01, 1.43, 1.5):
        neutron = bulk(native, k, 0)
        assert neutron[0] == pytest.approx(neutron[3], rel=2e-14, abs=2e-14)
        assert neutron[2] == pytest.approx(neutron[5] - neutron[4], rel=2e-14)
        assert bulk(native, k, .17)[0] == pytest.approx(bulk(native, k, .83)[0], abs=2e-13)
        assert bulk(native, k, .5)[2] == pytest.approx(0, abs=2e-13)
    sat = bulk(native, 1.43, .5)
    assert sat[0] == pytest.approx(-16.5, abs=1e-13)
    assert sat[1] == pytest.approx(0, abs=1e-13)
    h = 1e-4
    second = (bulk(native, 1.43+h, .5)[1]-bulk(native, 1.43-h, .5)[1])/(2*h)
    assert second == pytest.approx(143/1.43**2, rel=1e-6)


def test_primary_table7_with_printed_k_rounding(native):
    # Uploaded BBP1971, printed p259 Table7. k is printed to two decimals;
    # check the source box, without treating rounded central k as exact data.
    reference = [(.15,.35,.55),(.26,.91,1.43),(.36,1.58,2.46),(.49,2.63,4.08),
                 (.65,4.19,6.49),(.86,6.68,10.44),(1.01,8.80,13.93),
                 (1.17,11.46,18.49),(1.26,13.16,21.53),(1.29,13.77,22.64),
                 (1.40,16.19,27.13),(1.50,18.68,31.93)]
    for k, w, mu in reference:
        lo = bulk(native, k-.005, 0)
        hi = bulk(native, min(1.5,k+.005), 0)
        assert lo[0]-.005 <= w <= hi[0]+.005
        assert lo[4]-.005 <= mu <= hi[4]+.005


def test_automatic_derivatives_against_energy_differences(native):
    n, z = .012, [1.26, .245, .86]
    p = point(native,n,z)
    for i in range(3):
        h=2e-6
        zp=z.copy();zm=z.copy();zp[i]+=h;zm[i]-=h
        numerical=(point(native,n,zp)[0]-point(native,n,zm)[0])/(2*h)
        assert numerical == pytest.approx(p[1+i],rel=2e-7,abs=2e-7)
    h=n*2e-6
    derivative=(point(native,n+h,z)[0]-point(native,n-h,z)[0])/(2*h)
    assert derivative == pytest.approx(p[4],rel=2e-7)
    assert n*n*p[4] == pytest.approx(p[10],rel=2e-14)


def test_published_omission_is_only_the_documented_surface_term(native):
    n,z=.012,[1.26,.245,.86]
    exact=point(native,n,z)
    published=point(native,n,z,True)
    assert published[0] == pytest.approx(exact[0],abs=1e-14)
    assert published[1:3] == pytest.approx(exact[1:3],abs=1e-13)
    a,u,k,x,kg=exact[5:10]
    ni=k**3/(1.5*math.pi**2);ng=n*z[2]
    b=bulk(native,kg,0)
    dw_dng=b[1]*kg/(3*ng)
    correction=u*ni/n * a**(-1/3) * (21/16.5) * (1-ng/ni)**(2/3) * dw_dng
    assert exact[3]-published[3] == pytest.approx(correction*n,rel=2e-12)
    assert exact[4]-published[4] == pytest.approx(correction*z[2],rel=2e-12)


def test_published_microphysics_reproduces_primary_table_scale(native):
    # Tables1/2 p251-252: rounded input density, output, and bulk fits cannot
    # justify bitwise table agreement. This catches formula/units mistakes;
    # detailed source residuals are separately retained in the research note.
    data=[(2.79e-4,1.32,.313,127,40,5e-4),
          (5.25e-4,1.32,.307,134,41,6.42e-4),
          (8.79e-4,1.32,.302,140,42,8.73e-4),
          (1.59e-3,1.31,.295,149,44,1.47e-3),
          (3.73e-3,1.30,.280,170,48,4.04e-3),
          (8.91e-3,1.27,.256,211,54,1.38e-2),
          (.0204,1.24,.220,296,65,.0477),
          (.0475,1.20,.163,548,89,.176),
          (.0789,1.21,.121,990,120,.384)]
    for n,k,x,a,z,p in data:
        out=native("solve",n,1,0)
        assert out[0] == "ok"
        values=list(map(float,out[1].split()))
        assert values[0] == pytest.approx(k,abs=.006)
        if n < .0789:
            assert values[1] == pytest.approx(x,abs=.0015)
            assert values[3] == pytest.approx(a,abs=1.5)
            assert values[4] == pytest.approx(z,abs=.75)
            assert values[5] == pytest.approx(p,rel=.02)
        else:
            # Preserve, rather than hide, the ninth anchor's unresolved
            # mismatch. These regression envelopes describe discrepancies;
            # they are NOT scientific-admission or rounding tolerances.
            assert .044 < values[5]/p-1 < .045
            assert 37 < values[3]-a < 39
            assert -.0023 < values[1]-x < -.0021
        assert max(map(abs,values[9:12])) < 2e-9


@pytest.mark.parametrize("published",[True,False])
def test_pressure_matches_independent_equations_2_20_and_6_13(native,published):
    n=.0204
    out=native("solve",n,int(published),2)
    assert out[0] == "ok"
    k,x,kg,a,z,pressure,mu,mue,*_=map(float,out[1].split())
    ni=k**3/(1.5*math.pi**2);ng=kg**3/(1.5*math.pi**2)
    u=(n-ng)/(ni-ng);ne=u*ni*x
    wi=bulk(native,k,x)[0]
    gas=bulk(native,kg,0);wg=gas[0];mun0=gas[4]
    ws=(21/16.5)*(wg-wi)*(1-ng/ni)**(2/3)
    e2=1.43996454784255
    thickness=.74/(k**3-kg**3)**(1/3)
    wth=-(4*math.pi/9)*e2*thickness**2*k**3*x*x
    dws=-(2/3)*ws/(ni-ng)
    if not published:
        dws+=(21/16.5)*(1-ng/ni)**(2/3)*gas[1]*kg/(3*ng)
    dwth=(2/3)*wth/(ni-ng)
    mug=mun0+ni*u/(1-u)*(dwth+a**(-1/3)*dws)
    pg=ng*(mug-wg)
    pe=mue*ne/4
    wc0=(3/5)*(8/(9*math.pi))**(1/3)*e2*k
    pl=u*ni*wc0*x*x*a**(2/3)*(-.5*u**(1/3)+.5*u)
    assert pressure == pytest.approx(pg+pe+pl,rel=2e-9)
    assert mu == pytest.approx(939.56542052+mug,abs=2e-9)


@pytest.fixture(scope="module")
def generated(native,tmp_path_factory):
    work=tmp_path_factory.mktemp("bbp_original_generated")
    return {m:generate(native,work/str(m),m) for m in (1,2,3)}


def test_source_closures_retain_distinct_failures_and_phases(generated):
    pub,_=generated[1]
    var,vr=generated[2]
    maxwell,mr=generated[3]
    assert pub["max_closure_gradient"] < 2e-9
    assert pub["max_true_energy_gradient"] > .1
    assert .06 < pub["max_same_phase_firstlaw_relative"] < .07
    assert var["max_true_energy_gradient"] < 2e-9
    vp=[float(row["P_MeV_fm-3"]) for row in vr]
    assert any(b<a for a,b in zip(vp,vp[1:])), "raw nonconvex branch must not be silently repaired"
    assert maxwell["max_same_phase_firstlaw_relative"] < 2e-7
    assert maxwell["minimum_true_hessian_eigenvalue"] > 0
    assert maxwell["minimum_uniform_energy_gap_MeV"] > .1
    assert maxwell["phase_count"] == 2
    mp=[float(row["P_MeV_fm-3"]) for row in mr]
    assert all(b>=a for a,b in zip(mp,mp[1:]))
    assert sum(b==a for a,b in zip(mp,mp[1:])) == 1


def test_maxwell_common_tangent_and_one_sided_first_law(native):
    result=native("maxwell")
    assert result[0] == "ok"
    nl,nh,pl,ph,ml,mh,el,eh=map(float,result[1].split())
    assert nl<nh
    assert abs(pl-ph) < 1e-12
    assert abs(ml-mh) < 1e-9
    assert (eh-el)/(nh-nl) == pytest.approx((ml+mh)/2,abs=2e-8)
    # A common tangent alone does not establish that a coexistence segment
    # lowers the energy; compare against both locally stable source branches.
    for j in range(1,10):
        n=nl+(nh-nl)*j/10
        original=native("solve",n,0,0)
        assert original[0] == "ok"
        energy=float(original[1].split()[8])
        tangent=(ml+mh)*n/2-(pl+ph)/2
        assert energy >= tangent-2e-13
    for n,mode,p in [(nl,1,pl),(nh,2,ph)]:
        h=n*1e-4
        lo=list(map(float,native("solve",n-h,0,mode)[1].split()))
        hi=list(map(float,native("solve",n+h,0,mode)[1].split()))
        f_lo=lo[8]/(n-h);f_hi=hi[8]/(n+h)
        p_fd=n*n*(f_hi-f_lo)/(2*h)
        assert p_fd == pytest.approx(p,rel=2e-7)


@pytest.mark.parametrize("points",[3,9,33])
def test_sparse_maxwell_grid_has_three_one_sided_rows(native,tmp_path,points):
    diag,rows=generate(native,tmp_path/"out",3,points)
    assert len(rows)==points-diag["removed_grid_points"]+diag["inserted_grid_points"]
    assert sum(row["gas_free"]=="T" for row in rows)>=3
    assert sum(row["gas_free"]=="F" for row in rows)>=3


@pytest.mark.parametrize("mutation",[
    lambda t:t+"extra 0\n",
    lambda t:t.replace("w0_MeV 16.5","w0_MeV NaN"),
    lambda t:t.replace("w0_MeV 16.5","w0_MeV 1e9999"),
    lambda t:t.replace("w0_MeV 16.5","w0_MeV 1e300"),
    lambda t:t.replace("w0_MeV 16.5","w0_MeV 16.5 garbage"),
    lambda t:t.replace("include_exchange 0","include_exchange ,"),
    lambda t:t.replace("K_MeV 143","K_MeV 140"),
    lambda t:t.replace("BBP1971_ORIGINAL_PARAMETERS_V1","WRONG"),
    lambda t:"\n".join(t.splitlines()[:-1])+"\n",
])
def test_parameter_source_rejects_malformed_or_changed_records(native,tmp_path,mutation):
    source=tmp_path/"params.dat";source.write_text(mutation(PARAMS.read_text()))
    output=tmp_path/"out";output.mkdir()
    out=native("generate",1,3,source,output)
    assert out[0] != "ok"
    assert not (output/"bbp_original_diagnostics.json").exists()


@pytest.mark.parametrize("n",[0,1e-300,-.1,"NaN",.080001])
def test_native_point_rejects_out_of_domain_without_runtime_abort(native,n):
    assert native("point",n,1.3,.3,.8,0)[0].split()[0] == "F"


def test_parallel_generations_are_deterministic(native,tmp_path):
    paths=[tmp_path/"one",tmp_path/"two"]
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(lambda path:generate(native,path,3,33),paths))
    assert results[0] == results[1]
    for name in ["bbp_original_states.csv","bbp_original_diagnostics.json","bbp_original_transitions.csv"]:
        assert (paths[0]/name).read_bytes() == (paths[1]/name).read_bytes()
