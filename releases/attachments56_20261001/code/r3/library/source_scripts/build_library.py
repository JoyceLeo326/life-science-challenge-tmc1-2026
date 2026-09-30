"""Resumable ChEMBL acquisition and 2D parent library preparation.

Usage: python build_library.py fetch --pages 20
       python build_library.py prepare
Only one process; source JSON pages are never modified after download.
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski, QED
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(os.environ.get("TMC1_LIBRARY_WORK", str(Path(__file__).resolve().parent)))
RAW = Path(os.environ.get("TMC1_R2_RAW_CACHE", str(Path(__file__).resolve().parent.parent / "provenance/raw/expanded")))
OUT = ROOT / "processed"
LEGACY = Path(os.environ.get("TMC1_R1_APPROVED_CACHE", str(Path(__file__).resolve().parent.parent / "provenance/raw/legacy_approved")))
API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
LICENSE = "https://chembl.github.io/chembl-licensing/"
FIELDS = "molecule_chembl_id,molecule_structures,molecule_hierarchy,molecule_type,structure_type"
PAGE_SIZE = 1000
RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def append_jsonl(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(obj, ensure_ascii=False) + "\n")


def url_for(offset):
    return (f"{API}?molecule_type=Small+molecule&limit={PAGE_SIZE}"
            f"&offset={offset}&order_by=molecule_chembl_id&only={FIELDS}")


def fetch(pages, start):
    RAW.mkdir(parents=True, exist_ok=True)
    manifest = RAW / "manifest.jsonl"
    errors = RAW / "failures.jsonl"
    for n in range(start, start + pages):
        offset = n * PAGE_SIZE
        target = RAW / f"offset_{offset:07d}.json"
        if target.exists():
            try:
                body = target.read_bytes()
                if len(json.loads(body).get("molecules", [])) > 0:
                    print(f"cached {offset}: {len(body)} bytes", flush=True)
                    continue
            except (ValueError, OSError):
                target.rename(target.with_suffix(".corrupt"))
        url = url_for(offset)
        for attempt in range(1, 5):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "TMC1-research-library/1.0 (academic; contact via ChEMBL website)", "Accept": "application/json"})
                with urllib.request.urlopen(request, timeout=90) as response:
                    body = response.read()
                    status = response.status
                obj = json.loads(body)
                records = obj["molecules"]
                if not records:
                    raise ValueError("empty page")
                tmp = target.with_suffix(".tmp")
                tmp.write_bytes(body)
                tmp.replace(target)
                event = {"offset": offset, "url": url, "accessed_at_utc": utc(), "http_status": status,
                         "sha256": sha256(body), "bytes": len(body), "returned": len(records),
                         "total_count": obj.get("page_meta", {}).get("total_count"), "path": str(target)}
                append_jsonl(manifest, event)
                print(f"downloaded {offset}: {len(records)} records, {len(body)} bytes", flush=True)
                break
            except (urllib.error.URLError, TimeoutError, ValueError, KeyError, OSError) as exc:
                append_jsonl(errors, {"offset": offset, "url": url, "at_utc": utc(),
                                      "attempt": attempt, "error": repr(exc)})
                print(f"retry {offset} attempt {attempt}: {exc}", file=sys.stderr, flush=True)
                if attempt == 4:
                    print(f"FAILED page {offset}; future run can resume", file=sys.stderr, flush=True)
                else:
                    time.sleep(min(30, 3 * attempt))
        time.sleep(0.5)  # below ChEMBL documented 5 requests/second limit


def source_pages():
    # Existing approved pool is read-only. It is included for overlap accounting.
    for file in sorted(LEGACY.glob("molecules_offset_*.json")):
        yield "legacy_approved", file
    for file in sorted(RAW.glob("offset_*.json")):
        yield "expanded", file


def standardize(smiles, cleanup, frag, uncharger, tautomer):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("invalid_smiles")
    m = cleanup(m)
    m = frag.choose(m)
    m = uncharger.uncharge(m)
    m = tautomer.Canonicalize(m)
    Chem.SanitizeMol(m)
    if m.GetNumHeavyAtoms() == 0:
        raise ValueError("empty_parent")
    iso = Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
    stripped = Chem.Mol(m)
    Chem.RemoveStereochemistry(stripped)
    noniso = Chem.MolToSmiles(stripped, canonical=True, isomericSmiles=False)
    stereo_key = Chem.MolToInchiKey(m)
    parent_inchi_key = Chem.MolToInchiKey(stripped)
    if not stereo_key or not parent_inchi_key:
        raise ValueError("inchi_failure")
    return m, iso, noniso, stereo_key, parent_inchi_key.split("-")[0]


def properties(m):
    return {"mw": round(Descriptors.MolWt(m), 3), "clogp": round(Crippen.MolLogP(m), 3),
            "tpsa": round(rdMolDescriptors.CalcTPSA(m), 3), "hbd": Lipinski.NumHDonors(m),
            "hba": Lipinski.NumHAcceptors(m), "rotatable_bonds": Lipinski.NumRotatableBonds(m),
            "heavy_atoms": m.GetNumHeavyAtoms(), "rings": rdMolDescriptors.CalcNumRings(m),
            "qed": round(QED.qed(m), 4), "formal_charge": Chem.GetFormalCharge(m)}


def filter_reason(m, p):
    atoms = {a.GetSymbol() for a in m.GetAtoms()}
    if "C" not in atoms:
        return "no_carbon"
    if atoms - {"C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B"}:
        return "unsupported_element"
    if p["mw"] < 150 or p["mw"] > 650:
        return "mw_outside_150_650"
    if p["clogp"] < -2 or p["clogp"] > 7:
        return "clogp_outside_-2_7"
    if p["tpsa"] > 180:
        return "tpsa_over_180"
    if p["rotatable_bonds"] > 12:
        return "rotatable_bonds_over_12"
    if abs(p["formal_charge"]) > 2:
        return "absolute_charge_over_2"
    return "eligible"


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_sdf_from_csv():
    """Stream 2D SDF through a Python file handle for Unicode Windows paths."""
    source = OUT / "filtered_for_screening.csv"
    target = OUT / "filtered_for_screening_2d.sdf"
    count = 0
    with source.open("r", newline="", encoding="utf-8-sig") as input_file, target.open("w", encoding="utf-8") as output_file:
        reader = csv.DictReader(input_file)
        writer = Chem.SDWriter(output_file)
        for row in reader:
            m = Chem.MolFromSmiles(row["standardized_isomeric_smiles"])
            if m is None:
                continue
            for key in ("compound_id", "source_id", "inchi_key", "parent_key", "processing_status"):
                m.SetProp(key, str(row[key]))
            writer.write(m)
            count += 1
        writer.flush()
    return count


def summarize_csv():
    with (OUT / "all_source_records.csv").open("r", newline="", encoding="utf-8-sig") as f:
        sources = list(csv.DictReader(f))
    with (OUT / "all_unique_parents.csv").open("r", newline="", encoding="utf-8-sig") as f:
        parents = list(csv.DictReader(f))
    with (OUT / "filtered_for_screening.csv").open("r", newline="", encoding="utf-8-sig") as f:
        eligible = list(csv.DictReader(f))
    pages = []
    for group, file in source_pages():
        body = file.read_bytes()
        pages.append({"group": group, "path": str(file), "sha256": sha256(body),
                      "record_count": len(json.loads(body)["molecules"])})
    summary = {"created_at_utc": utc(), "rdkit_version": Chem.rdBase.rdkitVersion,
        "source": "ChEMBL molecule API", "license": "CC BY-SA 3.0", "license_url": LICENSE,
        "source_records_by_group": dict(Counter(r["source_group"] for r in sources)),
        "source_records_total": len(sources),
        "unique_source_ids": len({r["source_id"] for r in sources if r["source_id"]}),
        "standardized_unique_parents": len(parents),
        "property_filtered_parents": len(eligible),
        "excluded_unique_parents": len(parents)-len(eligible),
        "parent_filter_reasons": dict(Counter(r["processing_status"] for r in parents)),
        "source_processing_status_reasons": dict(Counter(r["processing_status"]+":"+r["processing_reason"] for r in sources)),
        "raw_pages": pages,
        "processing_rules": {"fragment": "RDKit LargestFragmentChooser preferOrganic",
                             "charge": "RDKit Uncharger", "tautomer": "RDKit canonical tautomer",
                             "parent_key": "first block of InChIKey after stereo removal",
                             "stereochemistry": "retained in isomeric SMILES, collapsed in parent count",
                             "3d": "none"}}
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    cleanup = rdMolStandardize.Cleanup
    frag = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    uncharger = rdMolStandardize.Uncharger()
    tautomer = rdMolStandardize.TautomerEnumerator()
    source_rows = []
    parents = {}
    seen_source = set()
    source_counts = Counter()
    statuses = Counter()
    pages = []
    for group, file in source_pages():
        body = file.read_bytes()
        obj = json.loads(body)
        pages.append({"group": group, "path": str(file), "sha256": sha256(body),
                      "record_count": len(obj["molecules"])})
        for rec in obj["molecules"]:
            source_counts[group] += 1
            sid = rec.get("molecule_chembl_id", "")
            structure = rec.get("molecule_structures") or {}
            raw = structure.get("canonical_smiles") or ""
            row = {"source_dataset": "ChEMBL", "source_group": group, "source_id": sid,
                   "source_url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{sid}/",
                   "raw_page": str(file), "raw_smiles": raw, "standardized_isomeric_smiles": "",
                   "standardized_parent_smiles": "", "inchi_key": "", "parent_key": "",
                   "compound_id": "", "processing_status": "", "processing_reason": ""}
            if not raw:
                row.update(processing_status="excluded", processing_reason="missing_smiles")
            elif sid in seen_source:
                row.update(processing_status="duplicate_source_id", processing_reason="overlap_legacy_expanded")
            else:
                seen_source.add(sid)
                try:
                    m, iso, noniso, key, parent_key = standardize(raw, cleanup, frag, uncharger, tautomer)
                    p = properties(m)
                    reason = filter_reason(m, p)
                    compound_id = "LIB_" + parent_key
                    row.update(standardized_isomeric_smiles=iso, standardized_parent_smiles=noniso,
                               inchi_key=key, parent_key=parent_key, compound_id=compound_id,
                               processing_status="standardized", processing_reason=reason)
                    if parent_key not in parents:
                        parents[parent_key] = {"compound_id": compound_id, "source_dataset": "ChEMBL",
                            "source_id": sid, "source_ids": [sid], "source_url": row["source_url"],
                            "raw_smiles": raw, "standardized_isomeric_smiles": iso,
                            "standardized_parent_smiles": noniso, "inchi_key": key,
                            "parent_key": parent_key, "processing_status": reason,
                            "stereo_forms_seen": {iso}, **p}
                    else:
                        parents[parent_key]["source_ids"].append(sid)
                        parents[parent_key]["stereo_forms_seen"].add(iso)
                except Exception as exc:
                    row.update(processing_status="excluded", processing_reason=str(exc)[:120])
            source_rows.append(row)
            statuses[row["processing_status"] + ":" + row["processing_reason"]] += 1
    parent_rows = []
    for k, row in sorted(parents.items()):
        row["source_ids"] = ";".join(sorted(set(row["source_ids"])))
        row["source_record_count"] = len(row["source_ids"].split(";"))
        row["stereo_form_count"] = len(row.pop("stereo_forms_seen"))
        parent_rows.append(row)
    fields_source = list(source_rows[0]) if source_rows else []
    fields_parent = list(parent_rows[0]) if parent_rows else []
    write_csv(OUT / "all_source_records.csv", source_rows, fields_source)
    write_csv(OUT / "all_unique_parents.csv", parent_rows, fields_parent)
    eligible = [r for r in parent_rows if r["processing_status"] == "eligible"]
    write_csv(OUT / "filtered_for_screening.csv", eligible, fields_parent)
    # 2D structures only; no conformer generation or docking here.
    write_sdf_from_csv()
    summary = summarize_csv()
    print(json.dumps({k: summary[k] for k in ("source_records_total","unique_source_ids","standardized_unique_parents","property_filtered_parents","parent_filter_reasons")}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="command", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--pages", type=int, default=20)
    f.add_argument("--start", type=int, default=0)
    sub.add_parser("prepare")
    sub.add_parser("finalize")
    a = ap.parse_args()
    if a.command == "fetch":
        fetch(a.pages, a.start)
    elif a.command == "prepare":
        prepare()
    else:
        count = write_sdf_from_csv()
        summary = summarize_csv()
        print(json.dumps({"sdf_records": count, "unique_parents": summary["standardized_unique_parents"],
                          "eligible": summary["property_filtered_parents"]}), flush=True)
