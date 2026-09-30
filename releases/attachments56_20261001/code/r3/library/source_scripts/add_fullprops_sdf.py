"""Write a 2D SDF with the full mandatory provenance fields from frozen CSV."""
import csv
from pathlib import Path
import os
from rdkit import Chem

BASE = Path(__file__).resolve().parent / "snapshots" / "v2_20000_new"
SOURCE = BASE / "filtered_for_screening.csv"
TARGET = BASE / "filtered_for_screening_fullprops_2d.sdf"
FIELDS = ("compound_id", "source_id", "source_ids", "raw_smiles", "standardized_isomeric_smiles",
          "inchi_key", "parent_key", "source_url", "processing_status", "mw", "clogp",
          "tpsa", "hbd", "hba", "rotatable_bonds", "heavy_atoms", "rings", "qed",
          "formal_charge", "stereo_form_count")
count = 0
with SOURCE.open("r", encoding="utf-8-sig", newline="") as inp, TARGET.open("w", encoding="utf-8") as out:
    reader = csv.DictReader(inp)
    writer = Chem.SDWriter(out)
    for row in reader:
        mol = Chem.MolFromSmiles(row["standardized_isomeric_smiles"])
        if mol is None:
            raise ValueError(f"bad_smiles:{row['compound_id']}")
        for field in FIELDS:
            mol.SetProp(field, row[field])
        writer.write(mol)
        count += 1
    writer.flush()
assert count == 17645, count
print(f"fullprops_sdf_records={count} path={TARGET}")
