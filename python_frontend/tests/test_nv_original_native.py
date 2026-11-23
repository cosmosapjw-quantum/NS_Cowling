"""Original NV1973 source binding and native bulk-functional regression.

Python only orchestrates compilation/runs and reads evidence. All independent
value/derivative oracles and the thermodynamic verification run in Fortran.
"""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/eos_original_papers_20261123/research_nv/nv1973_original_source.dat"
OLD_SOURCE = ROOT / "data/eos_primary_crust_generation_20261123/research_nv/nv1974_primary_source.dat"
MODULE = ROOT / "fortran_eos/src/compose_eos_nv_original.f90"

DRIVER = r"""program nv_original_test
 use compose_eos_types
 use compose_eos_io,only:write_compose,read_compose
 use compose_eos_nv_primary
 use compose_eos_nv_original
 use,intrinsic::iso_fortran_env,only:real128
 implicit none
 integer,parameter::hp=real128
 type(nv_original_source_t)::src
 type(nv_primary_fit_t)::oldfit
 type(nv_primary_state_t)::s,old
 type(nv_uniform_state_t)::uni,trial
 type(eos_data_t)::eos,back
 character(1024)::mode,path,out,arg,oldpath
 character(256)::status
 integer::i,j,k,n,ios
 real(wp)::nb,nn,np,v,g(2),hh(2,2),vs,gs(2),hs(2,2)
 real(wp)::maxv,maxg,maxh,maxsym,maxbeta,maxeuler,maxfirst,mind2,maxgas
 real(hp)::vr,gr(2),hr(2,2),a,b,delta,dmn,dmp
 call get_command_argument(1,mode);call get_command_argument(2,path)
 if(mode=='generate')then
  call get_command_argument(3,out);call get_command_argument(4,arg);read(arg,*,iostat=ios)n
  if(ios/=0)stop 3
  call generate_nv_original_native(trim(path),n,eos,trim(out),status);call check()
  if(eos%composition_available.or.eos%composition_file_present.or..not.eos%beta_equilibrated)stop 4
  call write_compose(trim(out),eos,status);call check()
  call read_compose(trim(out),back,status);call check()
  if(back%n/=n.or.maxval(abs(back%q-eos%q))>2e-14_wp)stop 5
  write(*,'(a,i0)')'GENERATED ',n
 else if(mode=='identity')then
  call load_nv_original_source(trim(path),src,status);call check()
  call get_command_argument(3,oldpath)
  call load_nv_primary_source(trim(oldpath),oldfit,status);call check()
  if(src%fit%table(8,7)/=1.996_wp.or.oldfit%table(8,7)/=1.966_wp)stop 10
  if(any(src%fit%coeff/=oldfit%coeff))stop 11
  if(src%fit%nlo/=oldfit%nlo.or.src%fit%nhi/=oldfit%nhi.or.src%fit%nref/=oldfit%nref)stop 12
  do i=1,4097
   nb=exp(log(src%fit%nlo)+real(i-1,wp)/4096*log(src%fit%nhi/src%fit%nlo))
   if(i==1)nb=src%fit%nlo
   if(i==4097)nb=src%fit%nhi
   call nv_primary_state(src%fit,nb,s,status);call check()
   call nv_primary_state(oldfit,nb,old,status);call check()
   if(any([s%e_offset,s%epsilon,s%pressure,s%mu_b,s%sound2,s%gamma]/= &
          [old%e_offset,old%epsilon,old%pressure,old%mu_b,old%sound2,old%gamma]))stop 13
  enddo
  write(*,'(a)')'EXACT_FLOAT64_SEMANTIC_IDENTITY 4097'
 else if(mode=='bulk_domain')then
  call load_nv_original_source(trim(path),src,status);call check()
  call get_command_argument(3,arg);read(arg,*)nn
  call get_command_argument(4,arg);read(arg,*)np
  call nv1973_bulk_potential(src,nn,np,v,g,hh,status);call check()
 else if(mode=='uniform_domain')then
  call load_nv_original_source(trim(path),src,status);call check()
  call get_command_argument(3,arg);read(arg,*)nb
  call get_command_argument(4,arg);read(arg,*)np
  call nv1973_uniform_composition(src,nb,np,uni,status);call check()
 else if(mode=='audit')then
  call load_nv_original_source(trim(path),src,status);call check()
  maxv=0;maxg=0;maxh=0;maxsym=0;maxbeta=0;maxeuler=0;maxfirst=0;mind2=huge(1._wp);maxgas=0
  ! Independent direct-power real128 oracle, including endpoints and vacuum.
  do i=0,15
   nn=real(i,wp)*.01_wp
   do j=0,15
    np=real(j,wp)*.01_wp
    call nv1973_bulk_potential(src,nn,np,v,g,hh,status);call check()
    call direct_polynomial(real(nn,hp),real(np,hp),vr,gr,hr)
    maxv=max(maxv,abs(v-real(vr,wp))/(1+abs(real(vr,wp))))
    maxg=max(maxg,maxval(abs(g-real(gr,wp))/(1+abs(real(gr,wp)))))
    maxh=max(maxh,maxval(abs(hh-real(hr,wp))/(1+abs(real(hr,wp)))))
    call nv1973_bulk_potential(src,np,nn,vs,gs,hs,status);call check()
    maxsym=max(maxsym,abs(v-vs)/(1+abs(v)),maxval(abs(g-gs(2:1:-1))/(1+abs(g))), &
     abs(hh(1,2)-hs(1,2))/(1+abs(hh(1,2))))
   enddo
  enddo
  if(maxv>1e-10_wp.or.maxg>1e-10_wp.or.maxh>1e-10_wp.or.maxsym>2e-10_wp)stop 14
  do i=1,129
   nb=exp(log(src%fit%nlo)+real(i-1,wp)/128*log(src%fit%nhi/src%fit%nlo))
   if(i==1)nb=src%fit%nlo
   if(i==129)nb=src%fit%nhi
   call nv1973_uniform_beta_state(src,nb,uni,status);call check()
   maxbeta=max(maxbeta,abs(uni%beta_residual));mind2=min(mind2,uni%beta_curvature)
   maxeuler=max(maxeuler,abs(uni%nn*uni%mu_n+uni%np*(uni%mu_p+uni%mu_e)-uni%energy-uni%pressure))
   ! Five-point real128 derivative of separately assembled energy functional.
   a=real(uni%nn,hp);b=real(uni%np,hp);delta=a*1e-5_hp
   dmn=(quad_energy(a-2*delta,b)-8*quad_energy(a-delta,b)+8*quad_energy(a+delta,b)- &
         quad_energy(a+2*delta,b))/(12*delta)
   delta=b*1e-5_hp
   dmp=(quad_energy(a,b-2*delta)-8*quad_energy(a,b-delta)+8*quad_energy(a,b+delta)- &
         quad_energy(a,b+2*delta))/(12*delta)
   maxfirst=max(maxfirst,abs(real(dmn,wp)-uni%mu_n),abs(real(dmp,wp)-uni%mu_p-uni%mu_e))
  enddo
  if(maxbeta>2e-10_wp.or.maxeuler>2e-10_wp.or.maxfirst>2e-8_wp.or.mind2<=0)stop 15
  do i=1,11
   nb=src%fit%table(1,i)*1e-39_wp
   if(i==1)nb=src%fit%nlo
   if(i==11)nb=src%fit%nhi
   call nv1973_uniform_beta_state(src,nb,uni,status);call check()
   maxgas=max(maxgas,abs(uni%e_offset-src%fit%table(9,i)))
  enddo
  write(*,'(*(es25.16e3,:,","))')maxv,maxg,maxh,maxsym,maxbeta,maxeuler,maxfirst,mind2,maxgas
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
 subroutine direct_polynomial(x,y,value,grad,hess)
  real(hp),intent(in)::x,y
  real(hp),intent(out)::value,grad(2),hess(2,2)
  integer::ia,ib
  real(hp)::c
  value=0;grad=0;hess=0
  do ia=0,6
   do ib=0,6
    c=real(src%vnm(ia,ib),hp)
    value=value+c*x**ia*y**ib
    if(ia>0)grad(1)=grad(1)+ia*c*x**(ia-1)*y**ib
    if(ib>0)grad(2)=grad(2)+ib*c*x**ia*y**(ib-1)
    if(ia>1)hess(1,1)=hess(1,1)+ia*(ia-1)*c*x**(ia-2)*y**ib
    if(ib>1)hess(2,2)=hess(2,2)+ib*(ib-1)*c*x**ia*y**(ib-2)
    if(ia>0.and.ib>0)hess(1,2)=hess(1,2)+ia*ib*c*x**(ia-1)*y**(ib-1)
   enddo
  enddo
  hess(2,1)=hess(1,2)
 end subroutine
 real(hp) function quad_energy(x,y)
  real(hp),intent(in)::x,y
  real(hp),parameter::qmn=939.56542052_hp,qmp=938.27208816_hp,qme=.51099895_hp
  real(hp),parameter::qhc=197.3269804_hp,qalpha=1._hp/137.035999084_hp,qpi=acos(-1._hp)
  real(hp)::pot,grad(2),hess(2,2),kin,slater,xe,ee
  call direct_polynomial(x,y,pot,grad,hess)
  kin=3._hp/5*qhc*qhc/(2*qmn)*(3*qpi*qpi)**(2._hp/3)*(x**(5._hp/3)+y**(5._hp/3))
  slater=-.75_hp*(3/qpi)**(1._hp/3)*qalpha*qhc*y**(4._hp/3)
  xe=qhc*(3*qpi*qpi*y)**(1._hp/3)/qme
  ee=qme**4/(8*qpi*qpi*qhc**3)*(xe*(1+2*xe*xe)*sqrt(1+xe*xe)-asinh(xe))
  quad_energy=qmn*x+qmp*y+kin+pot+slater+ee
 end function
end program
"""

@pytest.fixture(scope="module")
def executable(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran required")
    work = tmp_path_factory.mktemp("nv_original_native")
    driver = work / "driver.f90"
    driver.write_text(DRIVER)
    binary = work / "nv_original_test"
    subprocess.run([
        "gfortran", "-std=f2018", "-O2", "-fcheck=all",
        "-ffpe-trap=invalid,zero,overflow", "-ffree-line-length-none",
        "-J" + str(work), "-I" + str(work),
        str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
        str(ROOT / "fortran_eos/src/compose_eos_io.f90"),
        str(ROOT / "fortran_eos/src/compose_eos_nv_primary.f90"),
        str(MODULE), str(driver), "-o", str(binary),
    ], check=True, capture_output=True, text=True)
    return binary

def run(binary, *args):
    return subprocess.run([str(binary), *map(str, args)], cwd=ROOT,
                          capture_output=True, text=True)

def test_primary_fit_semantic_identity_and_historical_table_distinction(executable):
    result = run(executable, "identity", SOURCE, OLD_SOURCE)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "EXACT_FLOAT64_SEMANTIC_IDENTITY 4097" in result.stdout

def test_table1_functional_native_derivatives_and_uniform_thermodynamics(executable):
    result = run(executable, "audit", SOURCE)
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(result.stdout.strip().split(",")) == 9

@pytest.mark.parametrize("points", [3, 257, 513])
def test_original_native_generated_carrier_and_provenance(executable, tmp_path, points):
    before = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    old_before = hashlib.sha256(OLD_SOURCE.read_bytes()).hexdigest()
    out = tmp_path / "generated"
    out.mkdir()
    result = run(executable, "generate", SOURCE, out, points)
    assert result.returncode == 0, result.stdout + result.stderr
    record = json.loads((out / "NV_ORIGINAL_GENERATION.json").read_text())
    assert record["rows"] == points
    assert record["model"] == "NV1973_author_fit"
    assert record["doi"] == "10.1016/0375-9474(73)90349-7"
    assert record["original_pdf_sha256"] == "5611ae66700d983e7f66881839802c2204e16f91cafbf19c670fbd178042bb5b"
    assert record["published_fit_semantically_equal_to"] == "NV1974_author_fit"
    assert record["semantic_equivalent_model"] == "NV1974_author_fit"
    assert record["distinct_new_physical_eos"] is False
    assert record["full_hartree_fock_recomputed"] is False
    assert record["extrapolation"] is False
    assert record["composition_available"] is False
    assert record["electron_energy_included"] is True
    assert record["physical_matching_acceptance_asserted"] is False
    assert record["full_hf_missing_inputs"] == ["B_IJ", "C_IJ", "D_IJ", "spin_orbit_W"]
    assert record["missing_inputs_reference_doi"] == "10.1103/PhysRevC.5.1472"
    with (out / "nv1973_states.csv").open() as stream:
        states = list(csv.DictReader(stream))
    assert len(states) == points
    assert float(states[0]["n_b_fm3"]) == 2.79e-4
    assert float(states[-1]["n_b_fm3"]) == 7.89e-2
    with (out / "nv1973_table3_comparison.csv").open() as stream:
        table = list(csv.DictReader(stream))
    assert len(table) == 11
    assert float(table[6]["e_1973_table3_mev"]) == 1.996
    assert float(table[6]["e_1974_tableII_mev"]) == 1.966
    with (out / "nv1973_uniform_reconstruction.csv").open() as stream:
        uniform = list(csv.DictReader(stream))
    assert len(uniform) == 11
    assert record["max_uniform_reconstruction_minus_table3_mev"] > 0
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == before
    assert hashlib.sha256(OLD_SOURCE.read_bytes()).hexdigest() == old_before
    assert not (out / "eos.compo").exists()

@pytest.mark.parametrize("old,new", [
    ("NV1973_ORIGINAL_SOURCE_V1", "NV1974_ORIGINAL_SOURCE_V1"),
    ("5611ae66700d983e7f66881839802c2204e16f91cafbf19c670fbd178042bb5b", "0"*64),
    ("2.8822899e-1", "NaN"),
    ("2.8822899e-1", "Infinity"),
    ("2.8822899e-1", "1e999"),
    ("2.8822899e-1", ".3"),
    ("2.8822899e-1", "2*.28822899"),
    ("2.8822899e-1", ".28822899,"),
    ("1.3543555e-4", "1.3543555e-4 0"),
    ("1.3543555e-4", ""),
    ("2.79e-4 7.89e-2", "2e-4 7.89e-2"),
    ("LOG_REFERENCE_FM3 1e-4", "LOG_REFERENCE_FM3 1e35"),
    ("939.56542052 938.27208816", "939 938"),
    ("1.996 2.408", "1.966 2.408"),
    ("2.79e35 140", "2.78e35 140"),
    ("6 5 1.29625e9", "6 5 1.29625e10"),
    ("6 5 1.29625e9", "5 6 1.29625e9"),
    ("6 5 1.29625e9", ""),
    ("6 6 -1.88806e11", "6 6 -1.88806e11 trailing"),
    ("END_NV1973_SOURCE", ""),
    ("END_NV1973_SOURCE", "END_NV1973_SOURCE\ntrailing garbage"),
])
def test_reject_malformed_or_inauthentic_source_before_output(executable, tmp_path, old, new):
    raw = tmp_path / "bad.dat"
    original = SOURCE.read_text()
    assert old in original
    raw.write_text(original.replace(old, new))
    out = tmp_path / "rejected"
    out.mkdir()
    result = run(executable, "generate", raw, out, 257)
    assert result.returncode != 0
    assert "ERROR:" in result.stdout, result.stdout + result.stderr
    assert list(out.iterdir()) == []

@pytest.mark.parametrize("nn,np", [("-1e-4", ".01"), (".15", ".150001"), ("NaN", ".01")])
def test_bulk_polynomial_diagnostic_domain(executable, nn, np):
    result = run(executable, "bulk_domain", SOURCE, nn, np)
    assert result.returncode != 0
    assert "ERROR:" in result.stdout, result.stdout + result.stderr

@pytest.mark.parametrize("nb,np", [("0", "0"), (".000278999", "0"), (".078900001", ".001"),
                                    (".01", "-.0001"), (".01", ".01"), ("NaN", ".001")])
def test_uniform_reconstruction_source_domain(executable, nb, np):
    result = run(executable, "uniform_domain", SOURCE, nb, np)
    assert result.returncode != 0
    assert "ERROR:" in result.stdout, result.stdout + result.stderr

@pytest.mark.parametrize("points", [2, 0, 100001])
def test_grid_size_guard(executable, tmp_path, points):
    out = tmp_path / "rejected"
    out.mkdir()
    result = run(executable, "generate", SOURCE, out, points)
    assert result.returncode != 0
    assert "ERROR:" in result.stdout
    assert list(out.iterdir()) == []

def test_refuse_evidence_overwrite(executable, tmp_path):
    out = tmp_path / "generated"
    out.mkdir()
    first = run(executable, "generate", SOURCE, out, 257)
    assert first.returncode == 0, first.stdout + first.stderr
    before = {path.name: path.read_bytes() for path in out.iterdir()}
    second = run(executable, "generate", SOURCE, out, 513)
    assert second.returncode != 0
    assert "already exists" in second.stdout
    assert before == {path.name: path.read_bytes() for path in out.iterdir()}

@pytest.mark.parametrize("np", ["0", "1e-300", "1e-30", "1e-12"])
def test_uniform_tiny_compositions_remain_finite_under_fpe_traps(executable, np):
    result = run(executable, "uniform_domain", SOURCE, ".01", np)
    assert result.returncode == 0, result.stdout + result.stderr
