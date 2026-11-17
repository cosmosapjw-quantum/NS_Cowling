"""Reproduce numerical evidence for the native core and lower-crust backends.

Python here is an independent *test oracle*, not a production generator.
Run: PYTHONPATH=python_frontend python3 fortran_eos/tests/record_native_core_lower.py
"""
from pathlib import Path
import datetime, hashlib, io, json, os, platform, subprocess, sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
BUILD=ROOT/"build/native_core_lower_receipt"
DATA=ROOT/"data/research_native_20261116"
MODELS=("APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4")
FLAGS=["-std=f2018","-O2","-ffree-line-length-none","-fimplicit-none","-fcheck=all"]
CORE_COLS=("n_b","epsilon","pressure","mu_n","mu_p","mu_e",
           "Y_n","Y_p","Y_e","Y_mu","chi","n_low","n_high","yp_low","yp_high","gibbs_residual","phase")
LOWER_COLS=("segment","n_b","epsilon","pressure","mu_b","mu_e","Y_n","Y_p","Y_e","Y_mu","Y_h","A","Z",
            "pressure_residual","quadrature_error","composition_mu_residual")
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def invoke(args,**kw):
    return subprocess.run(list(map(str,args)),cwd=ROOT,capture_output=True,text=True,check=True,**kw)
def hashes(paths):return {str(p.relative_to(ROOT)):sha(p) for p in paths}
def save(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")
def main():
    from compose_eos.core_regeneration import MODEL_CONTRACTS,state_solver
    for area in ("core","lower"): (DATA/area).mkdir(parents=True,exist_ok=True)
    BUILD.mkdir(parents=True,exist_ok=True)
    src=ROOT/"fortran_eos/src";tests=ROOT/"fortran_eos/tests"
    common=["gfortran",*FLAGS,"-J"+str(BUILD),"-I"+str(BUILD),src/"compose_eos_types.f90"]
    core=BUILD/"core_probe";lower=BUILD/"lower_probe";unit=BUILD/"core_unit"
    builds=[]
    for module,driver,exe in (
        ("compose_eos_core_native.f90","core_native_probe.f90",core),
        ("compose_eos_core_native.f90","test_core_native.f90",unit),
        ("compose_eos_lower_crust_native.f90","native_lower_crust_driver.f90",lower)):
        r=invoke([*common,src/module,tests/driver,"-o",exe]);builds.append(r.stdout+r.stderr)
    (DATA/"core/receipt_build.log").write_text("\n".join(builds))
    result=invoke([unit]);(DATA/"core/native_unit.log").write_text(result.stdout+result.stderr)
    metadata={"baseline_commit":invoke(["git","rev-parse","HEAD"]).stdout.strip(),
        "utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "compiler":invoke(["gfortran","--version"]).stdout.splitlines()[0],
        "flags":FLAGS,"python_for_validation_only":sys.version,"platform":platform.platform(),
        "command":"PYTHONPATH=python_frontend python3 fortran_eos/tests/record_native_core_lower.py"}
    core_report=dict(metadata,backend="fortran",status="IMPLEMENTATION_VERIFIED",
        final_CompOSE_table_as_input=False,source_fidelity_and_physical_validity_separate=True,models={})
    allrows=0
    core_inputs=[src/"compose_eos_types.f90",src/"compose_eos_core_native.f90",
                 tests/"core_native_probe.f90",tests/"test_core_native.f90",Path(__file__),
                 ROOT/"python_frontend/tests/test_core_native_backend.py",
                 ROOT/"data/research_generator_20261116/bbb2_published_tables.json"]
    oracle_files=[ROOT/"python_frontend/compose_eos"/(name+".py") for name in
        ("core_regeneration","apr_uniform","apr_gibbs","bbb_uniform","gm1_uniform","skyrme_uniform","cold_leptons")]
    core_inputs+=oracle_files
    before=hashes(core_inputs)
    for model in MODELS:
        lo,hi=MODEL_CONTRACTS[model]["default_interval_fm3"]
        p=invoke([core,model,161,lo,hi]);rows=np.loadtxt(io.StringIO(p.stdout))
        np.savetxt(DATA/f"core/{model}.csv",rows,delimiter=",",header=",".join(CORE_COLS),comments="")
        solve,_=state_solver(model)
        oracle_rows=[solve(n) for n in rows[:,0]]
        expected=np.array([[r[key] for key in CORE_COLS[:10]] for r in oracle_rows])
        diff=np.abs(rows[:,:10]-expected)
        scaled=diff/np.maximum(np.abs(expected),1e-12)
        assert len(rows)==(163 if model=="APR_APR" else 161)
        np.testing.assert_allclose(rows[:,:10],expected,rtol=5e-10,atol=2e-8)
        assert np.isfinite(rows).all() and np.all(np.diff(rows[:,0])>0)
        baryon=np.max(abs(rows[:,6]+rows[:,7]-1))
        charge=np.max(abs(rows[:,7]-rows[:,8]-rows[:,9]))
        euler=np.max(abs(rows[:,1]+rows[:,2]-rows[:,0]*rows[:,3]))
        beta=np.max(abs(rows[:,3]-rows[:,4]-rows[:,5]))
        assert baryon<1e-12 and charge<1e-10 and euler<2e-7 and beta<2e-7
        path=DATA/f"core/{model}.csv"
        core_report["models"][model]={
            "rows":len(rows),"density_interval_fm3":[float(rows[0,0]),float(rows[-1,0])],
            "max_absolute_oracle_difference":dict(zip(CORE_COLS[:10],map(float,diff.max(axis=0)))),
            "max_scaled_oracle_difference":dict(zip(CORE_COLS[:10],map(float,scaled.max(axis=0)))),
            "max_baryon_residual":float(baryon),"max_charge_residual":float(charge),
            "max_Euler_residual_MeV_fm3":float(euler),"max_beta_residual_MeV":float(beta),
            "max_Gibbs_residual":float(np.max(rows[:,15])),
            "csv_sha256":sha(path),"representation_contract":MODEL_CONTRACTS[model]}
        allrows+=len(rows)
    assert before==hashes(core_inputs),"core source changed during measurement"
    core_report.update(total_rows=allrows,source_sha256=before,
        independent_checks=["native finite-difference first law",
           "fixed-neutron/proton energy derivatives for chemical potentials",
           "APR two-charge Gibbs equations","nonrelativistic free-Fermi pressure limit"],
        scope_limitations=["BBB uses published SNM/PNM energy reconstruction, not a microscopic BHF re-solve",
           "Original rounded parameter precision retained; target CompOSE agreement is a diagnostic",
           "Generation over source-model density support does not certify causal stellar support"])
    save(DATA/"core/TEST_RECEIPT.json",core_report)
    lower_inputs=[src/"compose_eos_types.f90",src/"compose_eos_lower_crust_native.f90",
                  tests/"native_lower_crust_driver.f90",Path(__file__),
                  ROOT/"python_frontend/tests/test_native_lower_crust.py",
                  ROOT/"data/outer_crust/hfb22_ame2016_cu.csv"]
    lower_inputs+=[ROOT/"python_frontend/compose_eos"/(name+".py") for name in
        ("outer_crust","bsk22_analytic","bsk22_extended_representation","bsk22_thermodynamic_reconstruction")]
    lower_inputs.append(ROOT/"python_frontend/tests/test_outer_crust.py")
    before=hashes(lower_inputs)
    p=invoke([lower,"source_rows",65]);rows=np.loadtxt(io.StringIO(p.stdout))
    np.savetxt(DATA/"lower/source_segments.csv",rows,delimiter=",",header=",".join(LOWER_COLS),comments="")
    outer=rows[rows[:,0]==0];bsk=rows[rows[:,0]==1]
    assert len(outer)>=65 and len(bsk)==17
    assert np.isfinite(rows).all()
    assert np.all(np.diff(outer[:,1])>0) and np.all(np.diff(bsk[:,1])>0)
    baryon=np.max(abs(rows[:,6]+rows[:,7]+rows[:,10]-1))
    charge=np.max(abs(rows[:,7]+rows[:,12]/rows[:,11]*rows[:,10]-rows[:,8]-rows[:,9]))
    euler=np.max(abs(rows[:,2]+rows[:,3]-rows[:,1]*rows[:,4]))
    assert baryon<1e-12 and charge<1e-12 and euler<2e-12
    seam=bsk[0,1:5]-outer[-1,1:5]
    assert seam[0]<0 and abs(seam[2])>1e-9 and abs(seam[3])>1e-5
    samep=np.where(outer[1:,3]==outer[:-1,3])[0]
    assert len(samep)>=16
    assert np.max(abs(outer[samep+1,4]-outer[samep,4]))<2e-8
    assert before==hashes(lower_inputs),"lower source changed during measurement"
    lower_report=dict(metadata,backend="fortran",status="IMPLEMENTATION_VERIFIED",
        mass_catalogue_rows=sum(1 for _ in (ROOT/"data/outer_crust/hfb22_ame2016_cu.csv").open())-1,
        outer_rows=len(outer),lower_source_fit_rows=len(bsk),outer_seed_points=65,
        lower_interval_fm3=[float(bsk[0,1]),float(bsk[-1,1])],
        outer_interval_fm3=[float(outer[0,1]),float(outer[-1,1])],
        equal_pressure_nuclear_transition_pairs=len(samep),
        max_nuclear_transition_mu_residual_MeV=float(np.max(abs(outer[samep+1,4]-outer[samep,4]))),
        max_baryon_residual=float(baryon),max_charge_residual=float(charge),
        max_Euler_residual_MeV_fm3=float(euler),
        source_seam={"difference_convention":"lower-first minus outer-last",
            "delta_nb_fm3":float(seam[0]),"delta_epsilon_MeV_fm3":float(seam[1]),
            "delta_pressure_MeV_fm3":float(seam[2]),"delta_mu_b_MeV":float(seam[3]),
            "density_overlap":True,"thermodynamically_joined":False,
            "cross_seam_interpolation_supported":False},
        max_lower_composition_mu_residual_MeV=float(np.max(abs(bsk[:,15]))),
        source_sha256=before,output_sha256={"source_segments.csv":sha(DATA/"lower/source_segments.csv")},
        independent_checks=["60-digit differentiation of nuclear+electron+lattice energy for pressure",
            "fixed-pressure global catalogue minimization", "outer isotope transition Gibbs equality",
            "anchored-branch finite-difference first law"],
        scope_limitations=["BSk22 is the published analytic representation, not a new ETF/Strutinsky minimization",
            "source-fit pressure and energy are distinct fitted quantities; their mismatch is retained",
            "Regenerated outer and BSk analytic source drips overlap; they are exported as separate segments",
            "Outer-anchored extension through .01 is explicitly a different thermodynamic representation"])
    test_log=DATA/"lower/FINAL_COMBINED_PYTEST.log"
    if test_log.is_file():
        evidence={"path":str(test_log.relative_to(ROOT)),"sha256":sha(test_log),
                  "summary":test_log.read_text().strip().splitlines()[-1],
                  "scope":"test_core_native_backend.py and test_native_lower_crust.py"}
        core_report["pytest_evidence"]=evidence;lower_report["pytest_evidence"]=evidence
    save(DATA/"core/TEST_RECEIPT.json",core_report)
    save(DATA/"lower/TEST_RECEIPT.json",lower_report)
    print(json.dumps({"core_rows":allrows,"core_models":len(MODELS),"outer_rows":len(outer),
                      "lower_rows":len(bsk),"status":"IMPLEMENTATION_VERIFIED"},indent=2))
if __name__=="__main__":main()
