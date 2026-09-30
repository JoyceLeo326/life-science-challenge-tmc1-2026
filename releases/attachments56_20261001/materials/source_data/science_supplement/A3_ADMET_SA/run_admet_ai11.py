"""One CPU-only, eleven-molecule exploratory ADMET-AI v2 inference."""

import csv
import hashlib
import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path

import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
import torch
from admet_ai import ADMETModel
from rdkit import Chem

HERE = Path(__file__).resolve().parent
INPUT = HERE / "ADMET_SA_panel11.csv"
torch.set_num_threads(1)
torch.manual_seed(20260929)
assert not torch.cuda.is_available(), "CPU-only run requested"

with INPUT.open(encoding="utf-8-sig", newline="") as f:
    panel = list(csv.DictReader(f))
assert len(panel) == 11
for row in panel:
    mol = Chem.MolFromSmiles(row["smiles"])
    assert mol is not None and Chem.MolToInchiKey(mol) == row["inchikey"]

model = ADMETModel(include_physchem=False, drugbank_path=None, num_workers=0)
assert model.device == "cpu"
preds = model.predict(smiles=[row["smiles"] for row in panel])
assert list(preds.index) == [row["smiles"] for row in panel]
assert len(preds) == 11
assert preds.notna().all().all()

rows = []
for record, (_, predicted) in zip(panel, preds.iterrows()):
    rows.append({**{key: record[key] for key in ("id", "name", "smiles", "inchikey")},
                 **{key: float(value) for key, value in predicted.items()}})
with (HERE / "ADMET_AI_v2_11_full_predictions.csv").open("w", encoding="utf-8-sig", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

module_root = Path(__import__("admet_ai").__file__).parent
model_files = sorted((module_root / "resources" / "models").rglob("*.pt"))
manifest = {
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "evidence_class": "exploratory_pretrained_machine_learning_prediction_not_observation",
    "versions": {name: importlib.metadata.version(name) for name in
                 ("admet-ai", "chemprop", "lightning", "torch", "rdkit", "numpy", "pandas")},
    "device": model.device,
    "torch_num_threads": torch.get_num_threads(),
    "num_workers": model.num_workers,
    "include_physchem": model.include_physchem,
    "drugbank_path": None,
    "seed": 20260929,
    "num_ensembles": model.num_ensembles,
    "num_model_files": len(model_files),
    "model_files": [{"relative_path": str(p.relative_to(module_root)), "bytes": p.stat().st_size,
                     "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in model_files],
    "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
    "predictions_sha256": hashlib.sha256((HERE / "ADMET_AI_v2_11_full_predictions.csv").read_bytes()).hexdigest(),
    "model_endpoint_names": [x for xs in model.task_lists for x in xs],
    "repository": "https://github.com/swansonk14/admet_ai",
    "package": "https://pypi.org/project/admet-ai/2.0.1/",
}
(HERE / "ADMET_AI_v2_inference_manifest.json").write_text(
    json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"num_molecules": len(rows), "num_predicted_endpoints": len(preds.columns),
                  "num_model_files": len(model_files), "device": model.device,
                  "prediction_columns": list(preds.columns)}, ensure_ascii=False))


