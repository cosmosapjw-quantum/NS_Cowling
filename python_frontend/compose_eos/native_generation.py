"""Input transport and provenance for the native EOS generator.

This module does not evaluate a nuclear functional, optimize a cell or choose a
phase. All physical decisions and actual table generation belong to Fortran.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path

NATIVE_FILES=("eos.nb","eos.t","eos.yq","eos.thermo","eos.compo")
NATIVE_DIAGNOSTICS=("sound_speed.csv", "attachment_diagnostics.json", "SUPPORT.json",
                    "states.csv", "first_law.csv", "attempts.csv", "profiles.csv", "selected.csv",
                    "primary_source_diagnostics.json", "bps_states.csv", "bps_transitions.csv",
                    "primary_source_values.csv", "primary_source_anchor_comparison.csv",
                    "NV_PRIMARY_GENERATION.json", "nv_primary_states.csv", "nv_primary_source_comparison.csv")

def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def _quoted(value):
    # Fortran namelist string quoting, not shell quoting.
    return "'"+str(value).replace("'","''")+"'"

def _real(value):
    value=float(value)
    if not math.isfinite(value):
        raise ValueError("nonfinite native request")
    return format(value,".17e")

def surface_input(root, model, variant):
    root=Path(root)
    from .native_model_contracts import SPG_INFERRED_VARIANT, SPG_INFERRED_SURFACE
    if variant==SPG_INFERRED_VARIANT:
        path=root/SPG_INFERRED_SURFACE["path"]
        if file_sha(path)!=SPG_INFERRED_SURFACE["sha256"]:
            raise ValueError("immutable inferred SPG surface input checksum mismatch")
        source=json.loads(path.read_text())
        fit=source["parameters"][model]
        sigma,b=float(fit["sigma0"]),float(fit["b"])
        if not math.isfinite(sigma+b) or sigma<=0 or b<=-16:
            raise ValueError("invalid immutable surface input")
        return {"sigma0":sigma,"surface_b":b,"label":variant,
                **SPG_INFERRED_SURFACE,
                "A_convention":source["A_convention"],
                "source_equivalence":"CONDITIONAL_RECONSTRUCTION_NOT_ORIGINAL_AUTHOR_CALIBRATION"}
    if variant=="AME2408_candidate":
        path=root/"data/research_generator_20261116/candidate_surface/candidate_surface_fit.json"
        fit=json.loads(path.read_text())["models"][model]["all_mass"]
    else:
        path=root/"data/research_completion_20261116/theory_review/surface_fit_variants_321.json"
        fit=json.loads(path.read_text())["models"][model]["variants"][variant]
    sigma,b=float(fit["sigma0_MeV_fm2"]),float(fit["b"])
    if not math.isfinite(sigma+b) or sigma<=0 or b<=-16:
        raise ValueError("invalid immutable surface input")
    return {"sigma0":sigma,"surface_b":b,"path":str(path.relative_to(root)),
            "sha256":file_sha(path),"label":variant,
            "source_equivalence":"CANDIDATE_NOT_ORIGINAL_AUTHOR_CALIBRATION"}

def write_request(root, output, request):
    root=Path(root);output=Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    path=output.parent/"request.nml"
    lines=["&generation", " model="+_quoted(request["model"]),
           " family="+_quoted(request["family"])]
    if request["family"]=="core":
        lines += [" n_points="+str(request["points"]),
                  " n_min="+_real(request["n_min"])," n_max="+_real(request["n_max"])]
    elif request["family"]=="primary_crust":
        source=root/request["primary_source_path"]
        if not source.is_file():
            raise ValueError("primary source input missing: "+str(source))
        lines += [" n_points="+str(request["points"]),
                  " primary_source_path="+_quoted(source.resolve())]
    elif request["family"]=="spg":
        surface=surface_input(root,request["model"],request["variant"])
        lines += [" n_points="+str(request["points"]),
                  " sigma0="+_real(surface["sigma0"])," surface_b="+_real(surface["surface_b"]),
                  " transition_policy='fixed_density_source'"]
    else:
        lines += [" n_density="+str(len(request["densities"])),
                  " densities="+",".join(_real(x) for x in request["densities"]),
                  " n_mesh="+str(len(request["resolutions"])),
                  " meshes="+",".join(str(x) for x in request["resolutions"])]
    if request.get("include_lower_crust"):
        mass=root/"data/outer_crust/hfb22_ame2016_cu.csv"
        if not mass.is_file():
            raise ValueError("nuclear mass catalogue missing")
        lines += [" include_lower_crust=.true.", " mass_catalogue="+_quoted(mass)]
    lines += ["/",""]
    path.write_text("\n".join(lines))
    return path

def source_files(root):
    root=Path(root)
    files=[root/"fortran_eos/build.sh",root/"fortran_eos/run-generator.sh"]
    files+=sorted((root/"fortran_eos/src").glob("*.f90"))
    files+=sorted((root/"fortran_eos/app").glob("*.f90"))
    files+=sorted((root/"fortran_eos/vendor").rglob("*.f"))
    return [p for p in files if p.is_file()]

def _read_native_table(directory, *, require_composition=True):
    """Use the Fortran reader again; a claimed receipt is not format evidence."""
    from .api import read_compose,check_consistency
    directory=Path(directory)
    hashes={}
    names=tuple(name for name in NATIVE_FILES if
                name!="eos.compo" or require_composition or (directory/name).exists())
    for name in names:
        path=directory/name
        if not path.is_file() or path.stat().st_size==0:
            raise ValueError("native output missing: "+str(path))
        hashes[name]=file_sha(path)
    try:
        table=read_compose(directory)
    except Exception as exc:
        raise ValueError("native output readback failed: "+str(directory)) from exc
    grid=[float(x) for x in table.nb]
    if (table.n != len(grid) or len(grid)<3 or
        any(not math.isfinite(x) or x<=0 for x in grid) or
        any(b<=a for a,b in zip(grid,grid[1:]))):
        raise ValueError("invalid native density grid")
    if table.nT>1 or table.nYq>1 or abs(table.temperature_mev)>1e-12:
        raise ValueError("native generator returned an unexpected thermodynamic slice")
    current_names=tuple(name for name in NATIVE_FILES if
                        name!="eos.compo" or require_composition or (directory/name).exists())
    if current_names!=names or hashes!={name:file_sha(directory/name) for name in names}:
        raise ValueError("native output changed during readback")
    return grid,hashes,check_consistency(table,verbose=False)

def _same_grid(actual,expected,label):
    if len(actual)!=len(expected) or any(
            not math.isclose(a,b,rel_tol=5e-13,abs_tol=0.) for a,b in zip(actual,expected)):
        raise ValueError(label+" density grid differs from the request")

def verify_native_generation(output,request):
    output=Path(output)
    receipt=json.loads((output/"RESULT.json").read_text())
    if (receipt.get("backend")!="fortran" or receipt.get("model")!=request["model"] or
        receipt.get("family")!=request["family"]):
        raise ValueError("native generator model/family/backend identity mismatch")
    if receipt.get("status")!="REQUESTED_ROWS_NATIVE_VERIFIED":
        raise ValueError("native physical generation incomplete: "+str(receipt.get("status")))
    if receipt.get("native_readback") is not True:
        raise ValueError("native readback evidence missing")
    rows=receipt.get("rows")
    if not isinstance(rows,int) or isinstance(rows,bool) or rows<3:
        raise ValueError("invalid native row count")
    needs_composition=request["family"]!="primary_crust" or request["model"]=="BPS1971_sequence"
    grid,hashes,checks=_read_native_table(output,require_composition=needs_composition)
    if rows!=len(grid):
        raise ValueError("native receipt row count differs from actual rows")
    readbacks={"generated":rows}
    consistency={"generated":checks}
    if request["family"]=="core":
        count=request["points"]
        if rows<count or rows>count+(2 if request["model"]=="APR_APR" else 0):
            raise ValueError("requested core rows missing")
        lo,hi=request["n_min"],request["n_max"]
        expected=[lo*(hi/lo)**(i/(count-1)) for i in range(count)]
        expected[0]=lo;expected[-1]=hi
        # APR inserts two in-window Gibbs endpoints. Consume each requested
        # value once without assuming where those extra rows occur.
        cursor=0
        for wanted in expected:
            while cursor<len(grid) and grid[cursor]<wanted and not math.isclose(grid[cursor],wanted,rel_tol=5e-13):
                cursor+=1
            if cursor>=len(grid) or not math.isclose(grid[cursor],wanted,rel_tol=5e-13,abs_tol=0.):
                raise ValueError("requested core density grid row missing")
            cursor+=1
        if (not math.isclose(grid[0],lo,rel_tol=5e-13,abs_tol=0.) or
            not math.isclose(grid[-1],hi,rel_tol=5e-13,abs_tol=0.)):
            raise ValueError("core output extends beyond the requested support")
    elif request["family"]=="primary_crust":
        count=request["points"]
        model=request["model"]
        if model in ("NV1974_author_fit", "BBP1971_table5") and rows!=count:
            raise ValueError("requested primary crust rows missing")
        if model=="BPS1971_sequence" and not count<=rows<=count+64:
            raise ValueError("BPS primary phase-complete row count is inconsistent")
        if model=="BBP1971_pressure_integral" and not count<=rows<=count+26:
            raise ValueError("BBP original knots or requested grid are missing")
        for name in request["contract"]["diagnostics"]:
            if not (output/name).is_file():
                raise ValueError("primary source diagnostic missing: "+name)
    elif request["family"]=="spg":
        if rows!=request["points"]:
            raise ValueError("requested SPG rows missing")
        expected=[grid[0]*(grid[-1]/grid[0])**(i/(rows-1)) for i in range(rows)]
        _same_grid(grid,expected,"SPG")
        transition=json.loads((output/"transition.json").read_text())
        if transition.get("policy")!="fixed_density_source":
            raise ValueError("unexpected SPG transition policy")
        if not math.isclose(float(transition["n_b"]),grid[-1],rel_tol=5e-9,abs_tol=0.):
            raise ValueError("SPG endpoint differs from its transition receipt")
        hashes["transition.json"]=file_sha(output/"transition.json")
    else:
        _same_grid(grid,request["densities"],"TF")
        if receipt.get("meshes")!=request["resolutions"]:
            raise ValueError("requested TF meshes not executed")
        if receipt.get("held_nodes")!=0:
            raise ValueError("unresolved TF density support")
        finest=None
        for mesh in request["resolutions"]:
            segment="mesh_"+str(mesh)
            mesh_grid,mesh_hashes,mesh_checks=_read_native_table(output/segment)
            _same_grid(mesh_grid,request["densities"],segment)
            readbacks[segment]=len(mesh_grid)
            consistency[segment]=mesh_checks
            for name,value in mesh_hashes.items():
                hashes[segment+"/"+name]=value
            finest=mesh_hashes
        if finest is None or any(hashes[name]!=finest[name] for name in NATIVE_FILES):
            raise ValueError("main TF table differs from the highest requested mesh")
    if request.get("include_lower_crust"):
        if receipt.get("lower_crust_segments") is not True:
            raise ValueError("requested lower crust not generated")
        segments=json.loads((output/"crust_segments.json").read_text())
        if segments.get("cross_segment_interpolation")!="UNAPPROVED":
            raise ValueError("missing explicit cross-segment support contract")
        hashes["crust_segments.json"]=file_sha(output/"crust_segments.json")
        for segment in ("outer","lower_source_fit"):
            segment_grid,segment_hashes,segment_checks=_read_native_table(output/segment)
            readbacks[segment]=len(segment_grid)
            consistency[segment]=segment_checks
            for name,value in segment_hashes.items():
                hashes[segment+"/"+name]=value
    # Optional native diagnostics are evidence, never a physical acceptance gate.
    # Hash them when present while keeping historical receipts readable.
    for name in NATIVE_DIAGNOSTICS:
        path=output/name
        if path.is_file():
            hashes[name]=file_sha(path)
    seals=[("RESULT.json",receipt)]
    # Public generation seals the verified inventory in the outer dispatcher
    # receipt. Use the structural relationship so archived runs remain movable.
    dispatch_path=output.parent/"DISPATCH.json"
    if output.name=="generated" and dispatch_path.is_file():
        dispatch=json.loads(dispatch_path.read_text())
        if not isinstance(dispatch,dict):
            raise ValueError("native generation dispatch receipt must be an object")
        if "evidence" in dispatch:
            seals.append(("DISPATCH.json",dispatch["evidence"]))
        # The first verification runs while DISPATCH is RUNNING, before the
        # evidence field exists. Historical unsealed receipts also stay readable.
    for seal_name,seal in seals:
        if not isinstance(seal,dict):
            raise ValueError("native generation seal must be an object: "+seal_name)
        expected_hashes=seal.get("output_sha256")
        if expected_hashes is None:
            continue
        if not isinstance(expected_hashes,dict):
            raise ValueError("native generation sha256 inventory is invalid: "+seal_name)
        for name in NATIVE_DIAGNOSTICS:
            if name in expected_hashes and name not in hashes:
                raise ValueError("sealed native diagnostic missing: "+name)
        # An older seal may predate collection of an optional diagnostic. Include
        # its present bytes in the new inventory without treating them as an old
        # verified hash. Required tables and previously sealed diagnostics must
        # match every supplied seal; neither receipt overrides the other.
        if any(expected_hashes.get(k)!=v for k,v in hashes.items()
               if k not in NATIVE_DIAGNOSTICS or k in expected_hashes):
            raise ValueError("native output checksum mismatch: "+seal_name)
    return {"rows":rows,"native_directory":str(output),"backend":"fortran",
            "native_result":receipt,"native_frontend_readback":readbacks,
            "output_sha256":hashes,"native_consistency":consistency}
