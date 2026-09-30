"""Repair CTAB stereo where exact full InChIKey can be retained; quarantine the rest."""
import csv
import hashlib
import io
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import os
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.error")
ROOT = Path(os.environ.get("TMC1_LIBRARY_WORK", str(Path(__file__).resolve().parent)))
BASE = ROOT / "snapshots" / "v2_20000_new"
CSV = BASE / "filtered_for_screening.csv"
ORIGINAL = BASE / "filtered_for_screening_2d.sdf"
SAFE = BASE / "filtered_for_screening_exact_stereo_fullprops_2d.sdf"
EXCLUDED = BASE / "sdf_exact_stereo_exclusions.csv"
REPORT = BASE / "sdf_exact_stereo_manifest.json"
PROPS = ("compound_id", "source_id", "source_ids", "raw_smiles", "standardized_isomeric_smiles",
         "inchi_key", "parent_key", "source_url", "processing_status", "mw", "clogp", "tpsa",
         "hbd", "hba", "rotatable_bonds", "heavy_atoms", "rings", "qed", "formal_charge", "stereo_form_count")


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def blocks(stream):
    buf = []
    for line in stream:
        buf.append(line)
        if line.strip() == b"$$$$":
            yield b"".join(buf)
            buf = []
    if buf:
        raise ValueError("incomplete_sdf")


def block_key(block):
    mol = next(iter(Chem.ForwardSDMolSupplier(io.BytesIO(block), sanitize=True, removeHs=True)), None)
    return Chem.MolToInchiKey(mol) if mol is not None else "INVALID"


def with_props(molblock, row, method):
    props = b"".join((f"> <{key}>\n{row[key]}\n\n").encode("utf-8") for key in PROPS)
    props += (f"> <sdf_generation_method>\n{method}\n\n").encode("utf-8")
    return molblock.rstrip(b"\r\n") + b"\n" + props + b"$$$$\n"


excluded = []
methods = Counter()
source_count = 0
with CSV.open("r", encoding="utf-8-sig", newline="") as cf, ORIGINAL.open("rb") as sf, SAFE.open("wb") as out:
    for row, oldblock in zip(csv.DictReader(cf), blocks(sf), strict=True):
        source_count += 1
        expected = row["inchi_key"]
        original_key = block_key(oldblock)
        if original_key == expected:
            original_molblock = oldblock.split(b"M  END", 1)[0] + b"M  END\n"
            selected = with_props(original_molblock, row, "original_ctab_verified")
            method = "original_ctab_verified"
        else:
            mol = Chem.MolFromSmiles(row["standardized_isomeric_smiles"])
            if mol is None or Chem.MolToInchiKey(mol) != expected:
                raise ValueError("CSV source identity error: " + row["compound_id"])
            selected = None
            method = None
            for label, kwargs in (("no_stereo_flags_v2000", {"includeStereo": False}),
                                  ("no_stereo_flags_v3000", {"includeStereo": False, "forceV3000": True})):
                candidate = with_props(Chem.MolToMolBlock(mol, **kwargs).encode("utf-8"), row, label)
                if block_key(candidate) == expected:
                    selected = candidate
                    method = label
                    break
            if selected is None:
                excluded.append({"compound_id": row["compound_id"], "source_id": row["source_id"],
                                 "expected_inchi_key": expected, "original_sdf_inchi_key": original_key,
                                 "reason": "stereochemistry_cannot_roundtrip_as_rdkit_ctab"})
                continue
        out.write(selected)
        methods[method] += 1

with EXCLUDED.open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["compound_id", "source_id", "expected_inchi_key", "original_sdf_inchi_key", "reason"])
    w.writeheader()
    w.writerows(excluded)

expected_by_id = {}
with CSV.open("r", encoding="utf-8-sig", newline="") as f:
    for row in csv.DictReader(f):
        expected_by_id[row["compound_id"]] = row["inchi_key"]
count = 0
errors = []
seen = set()
with SAFE.open("rb") as stream:
    for mol in Chem.ForwardSDMolSupplier(stream, sanitize=True, removeHs=True):
        count += 1
        if mol is None:
            errors.append({"record": count, "reason": "invalid"})
            continue
        cid = mol.GetProp("compound_id")
        if cid in seen or Chem.MolToInchiKey(mol) != expected_by_id.get(cid) or not all(mol.HasProp(k) for k in PROPS):
            errors.append({"record": count, "compound_id": cid, "reason": "key_duplicate_or_missing_prop"})
        seen.add(cid)
if errors or count + len(excluded) != source_count:
    raise RuntimeError(f"SDF validation failed:{count}:{len(excluded)}:{errors[:3]}")
report = {"created_at_utc": datetime.now(timezone.utc).isoformat(),
          "source_csv_sha256": sha(CSV), "superseded_minimal_sdf_sha256": sha(ORIGINAL),
          "superseded_fullprops_sdf_sha256": sha(BASE / "filtered_for_screening_fullprops_2d.sdf"),
          "repaired_sdf_sha256": sha(SAFE), "excluded_csv_sha256": sha(EXCLUDED),
          "source_csv_rows": source_count, "verified_sdf_records": count,
          "excluded_unrepresentable_in_ctab": len(excluded), "methods": dict(methods),
          "qc_pass": not errors,
          "note": "CSV retains all 17645 eligible molecules. Excluded SDF records retain authoritative SMILES/InChIKey in CSV."}
REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False), flush=True)
