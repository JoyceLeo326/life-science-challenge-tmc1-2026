"""Independent identity/count check and descriptive profile of a prepared snapshot."""
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import os

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.error")
ROOT = Path(os.environ.get("TMC1_LIBRARY_WORK", str(Path(__file__).resolve().parent)))
OUT = ROOT / "processed"


def read_csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def hash_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def quantiles(values):
    v = sorted(float(x) for x in values if x not in (None, ""))
    if not v:
        return {}
    def q(p):
        i = (len(v)-1)*p
        lo, hi = int(math.floor(i)), int(math.ceil(i))
        return round(v[lo] + (v[hi]-v[lo])*(i-lo), 3)
    return {"min": q(0), "p10": q(.1), "median": q(.5), "p90": q(.9), "max": q(1)}


summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
sources = read_csv(OUT / "all_source_records.csv")
parents = read_csv(OUT / "all_unique_parents.csv")
eligible = read_csv(OUT / "filtered_for_screening.csv")
errors = []
if len(sources) != summary["source_records_total"]:
    errors.append("source_count_mismatch")
if len(parents) != summary["standardized_unique_parents"]:
    errors.append("parent_count_mismatch")
if len(eligible) != summary["property_filtered_parents"]:
    errors.append("eligible_count_mismatch")
keys = [r["parent_key"] for r in parents]
if len(keys) != len(set(keys)):
    errors.append("duplicate_parent_key")
if not set(r["parent_key"] for r in eligible).issubset(set(keys)):
    errors.append("eligible_not_subset")
identity_errors = []
scaffolds = Counter()
for r in parents:
    m = Chem.MolFromSmiles(r["standardized_isomeric_smiles"])
    if m is None:
        identity_errors.append({"compound_id": r["compound_id"], "reason": "SMILES invalid"})
        continue
    full = Chem.MolToInchiKey(m)
    stripped = Chem.Mol(m)
    Chem.RemoveStereochemistry(stripped)
    parent = Chem.MolToInchiKey(stripped).split("-")[0]
    if full != r["inchi_key"] or parent != r["parent_key"]:
        identity_errors.append({"compound_id": r["compound_id"], "reason": "InChIKey mismatch",
                                "expected_full": full, "observed_full": r["inchi_key"],
                                "expected_parent": parent, "observed_parent": r["parent_key"]})
    if r["processing_status"] == "eligible":
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
        scaffolds[scaffold or "ACYCLIC"] += 1
sdf_count = 0
with (OUT / "filtered_for_screening_2d.sdf").open("rb") as stream:
    supplier = Chem.ForwardSDMolSupplier(stream, sanitize=True, removeHs=False)
    for mol in supplier:
        if mol is None:
            errors.append(f"invalid_sdf_record_at_{sdf_count+1}")
        sdf_count += 1
if sdf_count != len(eligible):
    errors.append("sdf_count_mismatch")
page_hash_mismatch = []
for page in summary["raw_pages"]:
    path = Path(page["path"])
    if not path.exists() or hash_file(path) != page["sha256"]:
        page_hash_mismatch.append(page["path"])
if page_hash_mismatch:
    errors.append("raw_page_hash_mismatch")
if identity_errors:
    errors.append("identity_errors")

profile = {"checked_at_utc": datetime.now(timezone.utc).isoformat(), "qc_pass": not errors,
    "errors": errors, "identity_error_count": len(identity_errors),
    "identity_error_examples": identity_errors[:10], "raw_page_hash_mismatch": page_hash_mismatch,
    "source_records": len(sources), "unique_parents": len(parents), "eligible_parents": len(eligible),
    "sdf_records": sdf_count, "eligible_scaffold_groups": len(scaffolds),
    "largest_scaffold_group": scaffolds.most_common(1),
    "acyclic_eligible": scaffolds["ACYCLIC"],
    "top_10_scaffolds": scaffolds.most_common(10),
    "property_quantiles_eligible": {key: quantiles(r[key] for r in eligible) for key in
                                   ("mw", "clogp", "tpsa", "hbd", "hba", "rotatable_bonds", "heavy_atoms", "rings", "qed", "formal_charge")},
    "stereo_form_count_distribution_all_parents": dict(Counter(r["stereo_form_count"] for r in parents)),
    "file_sha256": {p.name: hash_file(p) for p in OUT.glob("*.csv")}}
(OUT / "quality_profile.json").write_text(json.dumps(profile, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({k: profile[k] for k in ("qc_pass", "errors", "source_records", "unique_parents", "eligible_parents", "sdf_records", "eligible_scaffold_groups", "identity_error_count")}, ensure_ascii=False), flush=True)
if errors:
    raise SystemExit(1)
