"""NV author-fit native regression; Python only compiles, invokes, and reads evidence.

The independent derivative oracle is Fortran real128 finite differences of the
primary energy formula. No production EOS arithmetic runs in Python.
"""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/eos_primary_crust_generation_20261123/research_nv/nv1974_primary_source.dat"
MODULE = ROOT / "fortran_eos/src/compose_eos_nv_primary.f90"

DRIVER = r"""
program nv_test
 use compose_eos_types
 use compose_eos_io,only:write_compose,read_compose
 use compose_eos_nv_primary
 use,intrinsic::iso_fortran_env,only:real128
 implicit none
 integer,parameter::hp=real128
 type(nv_primary_fit_t)::fit
 type(nv_primary_state_t)::s
 type(eos_data_t)::eos,back
 character(1024)::mode,path,out,arg
 character(256)::status
 integer::i,n,ios
 real(wp)::nb,maxp,maxmu,maxsound,maxeuler,minsound,maxcs,maxcond
 real(hp)::nh,x,h,y0,ym2,ym1,yp1,yp2,d1,d2,pr,mur,csr
 call get_command_argument(1,mode)
 call get_command_argument(2,path)
 if(mode=='generate')then
  call get_command_argument(3,out);call get_command_argument(4,arg)
  read(arg,*,iostat=ios)n
  if(ios/=0)stop 3
  call generate_nv_primary_native(trim(path),n,eos,trim(out),status)
  call check()
  if(eos%composition_available.or.eos%composition_file_present)stop 4
  if(.not.eos%beta_equilibrated)stop 4
  call write_compose(trim(out),eos,status);call check()
  call read_compose(trim(out),back,status);call check()
  if(back%n/=n.or.maxval(abs(back%q-eos%q))>2e-14_wp)stop 5
  write(*,'(a,i0)')'GENERATED ',n
 else if(mode=='outside')then
  call load_nv_primary_source(trim(path),fit,status);call check()
  call get_command_argument(3,arg);read(arg,*)nb
  call nv_primary_state(fit,nb,s,status);call check()
 else if(mode=='audit')then
  call load_nv_primary_source(trim(path),fit,status);call check()
  maxp=0;maxmu=0;maxsound=0;maxeuler=0;minsound=1;maxcs=0;maxcond=0
  do i=1,4097
   nb=exp(log(fit%nlo)+real(i-1,wp)/4096*log(fit%nhi/fit%nlo))
   if(i==1)nb=fit%nlo
   if(i==4097)nb=fit%nhi
   call nv_primary_state(fit,nb,s,status);call check()
   minsound=min(minsound,s%sound2);maxcs=max(maxcs,s%sound2);maxcond=max(maxcond,s%cancellation)
   maxeuler=max(maxeuler,abs((s%epsilon+s%pressure)/nb-s%mu_b))
   ! A five-point, 128-bit finite-difference energy oracle, independent of
   ! native analytic derivative/Horner code. Test derivative only at interior points.
   if(mod(i,128)/=1.or.i==1.or.i==4097)cycle
   nh=real(nb,hp);x=log(nh/1e-4_hp);h=1e-5_hp
   y0=source_energy(x)
   ym2=source_energy(x-2*h);ym1=source_energy(x-h)
   yp1=source_energy(x+h);yp2=source_energy(x+2*h)
   d1=(ym2-8*ym1+8*yp1-yp2)/(12*h)
   d2=(-ym2+16*ym1-30*y0+16*yp1-yp2)/(12*h*h)
   pr=nh*d1;mur=939.56542052_hp+y0+d1;csr=(d1+d2)/mur
   maxp=max(maxp,abs(real(pr,wp)/s%pressure-1))
   maxmu=max(maxmu,abs(real(mur,wp)-s%mu_b))
   maxsound=max(maxsound,abs(real(csr,wp)/s%sound2-1))
  enddo
  if(maxp>2e-10_wp.or.maxmu>1e-9_wp.or.maxsound>2e-9_wp)stop 7
  if(minsound<=0.or.maxcs>=1.or.maxeuler>1e-10_wp)stop 8
  write(*,'(*(es25.16e3,:,","))')maxp,maxmu,maxsound,maxeuler,minsound,maxcs,maxcond
 else
  stop 9
 endif
contains
 subroutine check()
  if(status/='ok')then
   write(*,'(a)')'ERROR: '//trim(status)
   stop 2
  endif
 end subroutine
 real(hp) function source_energy(x)
  real(hp),intent(in)::x
  source_energy=-4._hp+exp(.28822899_hp+.59150523_hp*x+.090185940_hp*x**2 &
   -.11025614_hp*x**3+.029377479_hp*x**4-.0032618465_hp*x**5+.00013543555_hp*x**6)
 end function
end program
"""

@pytest.fixture(scope="module")
def executable(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran required")
    work = tmp_path_factory.mktemp("nv_primary_native")
    driver = work / "driver.f90"
    driver.write_text(DRIVER)
    binary = work / "nv_test"
    subprocess.run([
        "gfortran", "-std=f2018", "-O2", "-fcheck=all",
        "-ffpe-trap=invalid,zero,overflow", "-ffree-line-length-none",
        "-J" + str(work), "-I" + str(work),
        str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
        str(ROOT / "fortran_eos/src/compose_eos_io.f90"),
        str(MODULE), str(driver), "-o", str(binary),
    ], check=True, capture_output=True, text=True)
    return binary

def run(binary, *args):
    return subprocess.run([str(binary), *map(str, args)], cwd=ROOT,
                          capture_output=True, text=True)

def test_native_energy_derivatives_against_fortran_quad_oracle(executable):
    result = run(executable, "audit", SOURCE)
    assert result.returncode == 0, result.stdout + result.stderr
    values = [float(x) for x in result.stdout.strip().split(",")]
    assert len(values) == 7

@pytest.mark.parametrize("points", [3, 257, 513])
def test_native_generated_carrier_and_source_provenance(executable, tmp_path, points):
    before = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    out = tmp_path / "generated"
    out.mkdir()
    result = run(executable, "generate", SOURCE, out, points)
    assert result.returncode == 0, result.stdout + result.stderr
    record = json.loads((out / "NV_PRIMARY_GENERATION.json").read_text())
    assert record["rows"] == points
    assert record["model"] == "NV1974_author_fit"
    assert record["doi"] == "10.1017/S0074180900099885"
    assert record["full_hartree_fock_recomputed"] is False
    assert record["extrapolation"] is False
    assert record["electron_energy_included"] is True
    assert record["physical_matching_acceptance_asserted"] is False
    with (out / "nv_primary_states.csv").open() as stream:
        states = list(csv.DictReader(stream))
    assert len(states) == points
    assert float(states[0]["n_b_fm3"]) == 2.79e-4
    assert float(states[-1]["n_b_fm3"]) == 7.89e-2
    with (out / "nv_primary_source_comparison.csv").open() as stream:
        comparison = list(csv.DictReader(stream))
    assert len(comparison) == 11
    assert float(comparison[6]["e_table_1974_mev"]) == 1.966
    assert float(comparison[6]["e_historical_nuda1973_mev"]) == 1.996
    assert record["max_pressure_relative_euler_diagnostic"] > 0
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == before
    assert not (out / "eos.compo").exists()

@pytest.mark.parametrize("change", [
    lambda s: s.replace("2.8822899e-1", "NaN"),
    lambda s: s.replace("2.8822899e-1", "Infinity"),
    lambda s: s.replace("2.8822899e-1", "1e999"),
    lambda s: s.replace("2.8822899e-1", ".3"),
    lambda s: s.replace("2.8822899e-1", "2*.28822899"),
    lambda s: s.replace("2.8822899e-1", ".28822899,"),
    lambda s: s.replace("1.3543555e-4", "1.3543555e-4 0"),
    lambda s: s.replace("1.3543555e-4", ""),
    lambda s: s.replace("2.79e-4 7.89e-2", "2e-4 7.89e-2"),
    lambda s: s.replace("LOG_REFERENCE_FM3 1e-4", "LOG_REFERENCE_FM3 1e35"),
    lambda s: s.replace("939.56542052 938.27208816", "939 938"),
    lambda s: s.replace("1.966 2.408", "1.996 2.408"),
    lambda s: s.replace("2.79e35 140", "2.78e35 140"),
    lambda s: s.replace("2.79e35 140", "2.79e35 141"),
    lambda s: s.replace("40 0.2 -26.8", "40 0.3 -26.8"),
    lambda s: s.replace("7.89e37 950 32 15.0 -68.3 7.37e-2 0.16 8.611 8.657", ""),
    lambda s: s + "trailing garbage\n",
])
def test_native_rejects_wrong_primary_contract_before_output(executable, tmp_path, change):
    raw = tmp_path / "bad.dat"
    raw.write_text(change(SOURCE.read_text()))
    out = tmp_path / "rejected"
    out.mkdir()
    result = run(executable, "generate", raw, out, 257)
    assert result.returncode != 0
    assert "ERROR:" in result.stdout, result.stdout + result.stderr
    assert list(out.iterdir()) == []

@pytest.mark.parametrize("density", ["0", "-.1", ".000278999", ".078900001", "NaN"])
def test_native_never_extrapolates(executable, density):
    result = run(executable, "outside", SOURCE, density)
    assert result.returncode != 0
    assert "ERROR:" in result.stdout

@pytest.mark.parametrize("points", [2, 0, 100001])
def test_grid_size_guard(executable, tmp_path, points):
    out = tmp_path / "rejected"
    out.mkdir()
    result = run(executable, "generate", SOURCE, out, points)
    assert result.returncode != 0
    assert "ERROR:" in result.stdout
    assert list(out.iterdir()) == []

def test_native_refuses_evidence_overwrite(executable, tmp_path):
    out = tmp_path / "generated"
    out.mkdir()
    result = run(executable, "generate", SOURCE, out, 3)
    assert result.returncode == 0, result.stdout + result.stderr
    old = {p.name: p.read_bytes() for p in out.iterdir()}
    again = run(executable, "generate", SOURCE, out, 3)
    assert again.returncode != 0
    assert {p.name: p.read_bytes() for p in out.iterdir()} == old
