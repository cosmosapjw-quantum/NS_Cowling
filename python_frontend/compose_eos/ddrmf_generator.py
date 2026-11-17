"""Generate cold homogeneous npe-muon CompOSE tables with the native DD-RMF solver.

Densities are in fm^-3, chemical potentials in MeV, energy/pressure in MeV fm^-3.
The user supplies the sampled density domain. Physical diagnostics accompany every
table; a high-density limitation is recorded without deleting computed states.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .api import check_consistency, read_compose, write_compose
from .eos_data import EOSData
from .physical_domain import physical_domain_report
from .native import evaluate_ddrmf, executable, project_root

SUPPORTED_MODELS = ("GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate_ddrmf_core(model_id: str, densities, output: str | Path) -> dict:
    """Solve each requested uniform state and export its explicit composition.

    Cluster abundance is exactly zero in the assumed homogeneous phase. Abar=1
    and Zbar=0 supply a harmless normalization for that zero-abundance group.
    This does not locate the crust transition or establish a global ground state.
    """
    if model_id not in SUPPORTED_MODELS:
        raise ValueError(f"supported independent models are {SUPPORTED_MODELS}")
    density = np.asarray(densities, dtype=float)
    if (density.ndim != 1 or density.size < 3 or not np.isfinite(density).all()
            or np.any(density <= 0) or np.any(np.diff(density) <= 0)):
        raise ValueError("require at least three finite positive increasing densities")

    root = project_root()
    registry = root / "data/eos390_rmf_parameter_sets.v1.json"
    masses = json.loads(registry.read_text())["models"][model_id]["nucleon_masses_MeV"]
    mn, mp = float(masses["neutron"]), float(masses["proton"])
    output = Path(output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt_path = output / "generation_receipt.json"
    receipt = {"schema_version": 1, "model_id": model_id, "status": "RUNNING",
               "matter": "cold charge-neutral homogeneous npe-muon matter",
               "density_domain_fm3": [float(density[0]), float(density[-1])],
               "row_count": int(density.size)}
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    try:
        rows = []
        for n in density:
            state = evaluate_ddrmf(model_id, float(n))
            for key, tolerance in (("beta_residual_mev", 2e-8),
                                   ("charge_residual_fm3", 2e-12),
                                   ("scalar_field_residual_mev", 2e-8),
                                   ("pressure_identity_residual_mev_fm3", 2e-7)):
                value = float(state[key])
                if not np.isfinite(value) or abs(value) > tolerance:
                    raise ValueError(f"native {key} at n={n}: {value}")
            rows.append(state)

        def column(name):
            values = np.array([state[name] for state in rows], dtype=float)
            if not np.isfinite(values).all():
                raise ValueError(f"nonfinite native field {name}")
            return values

        yp, ye, ymu = (column(key) for key in ("yp", "ye", "ymu"))
        if np.any(yp < 0) or np.any(yp > 1) or np.any(ye < 0) or np.any(ymu < 0):
            raise ValueError("native composition outside nonnegative fraction support")
        if np.max(np.abs(yp - ye - ymu)) > 2e-11:
            raise ValueError("native charge fractions do not close")
        energy, pressure = column("eps_total_mev_fm3"), column("pressure_total_mev_fm3")
        mub, mue = column("mu_n_mev"), column("mu_e_mev")
        if not np.allclose(energy + pressure, density * mub, rtol=2e-10, atol=1e-10):
            raise ValueError("native total Euler identity failed")
        zero = np.zeros(density.size)
        qenergy = energy / (density * mn) - 1.0
        data = EOSData(
            nb=density, mn=mn, mp=mp, Q1=pressure / density,
            Q2=zero.copy(), Q3=mub / mn - 1.0, Q4=-mue / mn,
            Q5=zero.copy(), Q6=qenergy, Q7=qenergy.copy(), q_extra=zero.copy(),
            Iphase=np.full(density.size, 3, dtype=int),
            Yn=1.0 - yp, Yp=yp, Ye=ye, Ymu=ymu,
            Aav=np.ones(density.size), Zav=zero.copy(), Yh=zero.copy(),
            composition_file_present=True, composition_available=True,
            beta_equilibrated=True,
        )
        compose = output / "compose"
        write_compose(compose, data)
        reread = read_compose(compose)
        for name in ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7",
                     "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh"):
            if not np.allclose(getattr(data, name), getattr(reread, name),
                               rtol=2e-14, atol=2e-15):
                raise ValueError(f"native readback changed {name}")
        if not np.array_equal(data.Iphase, reread.Iphase):
            raise ValueError("native readback changed homogeneous phase labels")
        for name in ("Yn_present", "Yp_present", "Ye_present", "Ymu_present",
                     "Yh_present", "C3_applicable", "C4_applicable"):
            if not np.asarray(getattr(reread, name)).all():
                raise ValueError(f"native readback lost composition mask {name}")
        checks = check_consistency(reread, verbose=False)
        if not checks["C3"]["pass_"] or not checks["C4"]["pass_"]:
            raise ValueError("native baryon/charge closure failed after readback")
        (output / "native_states.json").write_text(json.dumps(rows, indent=2) + "\n")
        receipt.update(
            status="GENERATED_AND_NATIVE_READBACK_VERIFIED",
            reference_masses_mev={"neutron": mn, "proton": mp},
            composition={"species": ["n", "p", "e", "mu"], "phase_code": 3,
                         "phase_meaning": "homogeneous nucleonic matter",
                         "zero_cluster_normalization": {"Abar": 1, "Zbar": 0, "Yh": 0},
                         "C3": checks["C3"], "C4": checks["C4"]},
            physical_diagnostics=checks,
            sampled_domain_diagnostics=physical_domain_report(reread),
            max_abs_beta_residual_mev=max(abs(row["beta_residual_mev"]) for row in rows),
            max_abs_charge_residual_fm3=max(abs(row["charge_residual_fm3"]) for row in rows),
            muon_active_rows=sum(row["muon_branch"] == "ACTIVE" for row in rows),
            provenance={
                "native_source_sha256": _sha256(root / "fortran_eos/src/compose_eos_rmf.f90"),
                "parameter_registry_sha256": _sha256(registry),
                "native_executable_sha256": _sha256(executable()),
                "generator_sha256": _sha256(Path(__file__)),
                "physical_domain_sha256": _sha256(root / "python_frontend/compose_eos/physical_domain.py"),
                "native_cli_sha256": _sha256(root / "fortran_eos/app/compose_eos_fortran.f90"),
                "equations": "mu_n=mu_p+mu_e; n_p=n_e+n_mu; eps+P=n_B*mu_n",
                "pressure_rearrangement": "included by native solver",
                "source_comparison_status": "not_run"},
            output_sha256={p.relative_to(output).as_posix(): _sha256(p)
                           for p in sorted(output.rglob("*")) if p.is_file() and p != receipt_path},
        )
    except Exception as exc:
        receipt.update(status="FAILED", failure_type=type(exc).__name__, failure_message=str(exc))
        receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
        raise
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    return receipt
