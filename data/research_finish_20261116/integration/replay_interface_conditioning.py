from pathlib import Path
import gzip, hashlib, json

root = Path(".")
base = root / "data/research_finish_20261116/crust_assembly"
models = ["GPPVA_DD2", "GPPVA_TW", "GPPVA_FSU2R", "GPPVA_FSU2H"]
models += [f"SPG_M{i}" for i in range(1, 6)]
for model in models:
    relative = (
        f"tf_interfaces/{model}/INTERFACE.json.gz"
        if model.startswith("GPPVA")
        else f"restored_spg/{model}/RECEIPT.json"
    )
    path = base / relative
    raw = path.read_bytes()
    row = json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)
    lower = row["bsk_endpoint"]
    inner = row["tf_endpoint"] if "tf_endpoint" in row else row["spg_endpoint"]
    nl = lower["n_b"]
    ni = inner["n_b"] if "n_b" in inner else inner["n_B"]
    residual = lower["source_composition_mu_n_residual_MeV"]
    slope = 1 / ni - 1 / nl
    assert ni > nl > 0 and slope < 0
    sensitivity = 1 / abs(slope)
    diagnostic = residual / slope / lower["pressure"]
    print(model, sensitivity, residual, diagnostic, hashlib.sha256(raw).hexdigest())
