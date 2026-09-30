"""Join official exact-key active60 additions with the immutable formal278 background."""

import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
ROOT = PARENT.parents[1]
ACTIVE = ROOT / "03_算法与计算" / "ranking_active_v2_after_ai60" / "AI_active_60.csv"
CONFIRM = ROOT / "03_算法与计算" / "formal_confirmation_panel_20260929" / "confirmation_panel.csv"
PRIOR = ROOT / "01_机制与文献" / "既有药理与可购性_交付表.csv"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


plan = json.loads((HERE / "incremental_query_plan.json").read_text(encoding="utf-8"))
manifest = json.loads((HERE / "incremental_chembl_request_manifest.json").read_text(encoding="utf-8"))
assert sha(ACTIVE) == plan["frozen_csv_sha256"]
assert plan["issues"] == []
ids = set(plan["new_exact_source_ids_to_query"])
assert ids and manifest["exact_id_count"] == len(ids)

api = {}
url_by_id = {}
raw_by_id = {}
mechanisms = defaultdict(list)
for req in manifest["requests"]:
    path = HERE / req["raw_file"]
    assert sha(path) == req["raw_sha256"]
    assert set(req["request_ids"]) <= ids
    data = json.loads(path.read_text(encoding="utf-8"))
    field = "molecules" if req["resource"] == "molecule" else "mechanisms"
    records = data[field]
    assert len(records) == req["returned_record_count"]
    assert set(r["molecule_chembl_id"] for r in records) <= set(req["request_ids"])
    for sid in req["request_ids"]:
        url_by_id[req["resource"], sid] = req["url"]
        raw_by_id[req["resource"], sid] = req["raw_file"]
    for record in records:
        sid = record["molecule_chembl_id"]
        if req["resource"] == "molecule":
            assert sid not in api
            api[sid] = record
        else:
            mechanisms[sid].append(record)
assert set(api) == ids
assert not mechanisms  # no target-name request is needed for this batch

new_identity = read_csv(HERE / "incremental_new_candidate_identity.csv")
old_prefetch = read_csv(PARENT / "formal278_pharmacology_prefetch.csv")
old_ids = {r["compound_id"] for r in old_prefetch}
assert len(old_ids) == len(old_prefetch)
new_ids = {r["compound_id"] for r in new_identity}
assert not new_ids & old_ids
active = read_csv(ACTIVE)
active_by_id = {r["compound_id"]: r for r in active}
assert len(active_by_id) == 60
assert len(new_ids) == plan["new_compounds_not_in_formal278"]
assert set(active_by_id) - old_ids == new_ids

prior_matches = defaultdict(list)
for row in read_csv(PRIOR):
    prior_matches[row["chembl_id"], row["inchikey_chembl"]].append(row)

source_rows = []
new_rows = []
warnings = []
for identity in new_identity:
    cid = identity["compound_id"]
    key = identity["candidate_full_inchi_key"]
    source_ids = identity["exact_local_full_key_source_ids"].split(";")
    assert all(sid in ids for sid in source_ids)
    assert not identity["already_queried_exact_source_ids"]
    matched = []
    mismatched = []
    for sid in source_ids:
        record = api[sid]
        api_key = (record.get("molecule_structures") or {}).get("standard_inchi_key") or ""
        (matched if api_key == key else mismatched).append(sid)
        phase = record.get("max_phase")
        hierarchy = record.get("molecule_hierarchy") or {}
        source_rows.append({
            "compound_id": cid, "candidate_inchi_key": key, "chembl_id": sid,
            "source_identity_relation": "api_full_inchikey_match" if api_key == key else "source_id_parent_normalized_api_key_differs",
            "api_standard_inchi_key": api_key,
            "preferred_name_chembl": record.get("pref_name") or "",
            "max_phase_chembl_raw": "" if phase is None else str(phase),
            "first_approval_year_chembl": record.get("first_approval") or "",
            "withdrawn_flag_chembl": record.get("withdrawn_flag"),
            "black_box_warning_flag_chembl": record.get("black_box_warning"),
            "parent_chembl_id": hierarchy.get("parent_chembl_id") or "",
            "active_chembl_id": hierarchy.get("active_chembl_id") or "",
            "molecule_record_url": f"https://www.ebi.ac.uk/chembl/api/data/molecule/{sid}.json",
            "batch_response_file": raw_by_id["molecule", sid],
            "batch_request_url": url_by_id["molecule", sid],
        })
    if not matched:
        warnings.append({"compound_id": cid, "issue": "no_official_API_full_key_match",
                         "source_ids": ";".join(source_ids), "candidate_inchi_key": key})
    primary_id = active_by_id[cid]["source_id"]
    assert primary_id in source_ids
    primary = api[primary_id]
    phase_exact = [float(api[sid]["max_phase"]) for sid in matched if api[sid].get("max_phase") is not None]
    phase_all = [float(api[sid]["max_phase"]) for sid in source_ids if api[sid].get("max_phase") is not None]
    prior = [p for sid in source_ids for p in prior_matches[sid, key]]
    new_rows.append({
        "compound_id": cid, "arms": "AI_active", "computed_structure_label": cid,
        "primary_source_pref_name_chembl_not_candidate_name": primary.get("pref_name") or "",
        "primary_chembl_id": primary_id,
        "source_names_chembl": "; ".join(f"{sid}:{api[sid].get('pref_name') or '(no preferred name)'}" for sid in source_ids),
        "inchi_key": key,
        "standardized_isomeric_smiles": identity["candidate_standardized_isomeric_smiles"],
        "queried_exact_local_source_ids": ";".join(source_ids),
        "api_full_key_matching_ids": ";".join(matched),
        "api_key_different_source_form_ids": ";".join(mismatched),
        "same_parent_other_form_source_ids_not_queried": identity["other_form_source_ids_not_to_query"],
        "chembl_max_phase_for_API_key_match": max(phase_exact) if phase_exact else "",
        "chembl_max_phase_including_source_forms": max(phase_all) if phase_all else "",
        "recorded_exact_key_mechanism_count": 0,
        "recorded_verified_parent_form_mechanism_count": 0,
        "recorded_unlinked_form_mechanism_count": 0,
        "exact_key_mechanism_target_action": "",
        "verified_parent_form_mechanism_target_action": "",
        "unlinked_source_form_mechanism_not_attributed_to_candidate": "",
        "prior_18_exact_id_key_match": bool(prior),
        "prior_verified_context": "; ".join(sorted({r.get("known_pharmacology_context", "") for r in prior if r.get("known_pharmacology_context")})),
        "prior_context_source_url": ";".join(sorted({r.get("context_source_url", "") for r in prior if r.get("context_source_url")})),
        "pharmacology_coverage_status": "ChEMBL_no_mechanism_record_for_queried_source_ids" +
            (";no_API_full_key_match_to_standardized_candidate" if not matched else ""),
        "tmc1_function_status": "not_assessed_by_ChEMBL_source_record",
        "purchase_status": "not_checked",
        "clinical_indication_status": "not_assessed;max_phase_is_ChEMBL_record_not_jurisdiction_specific_label",
        "source_molecule_url": f"https://www.ebi.ac.uk/chembl/api/data/molecule/{primary_id}.json",
        "source_mechanism_query_url": url_by_id["mechanism", primary_id],
    })

columns = list(old_prefetch[0])
assert all(set(r) == set(columns) for r in new_rows)
combined = old_prefetch + new_rows
assert len(combined) == len(old_ids | new_ids) == 319
combined_by_id = {r["compound_id"]: r for r in combined}
active_rows = [combined_by_id[r["compound_id"]] for r in active]
confirmation = read_csv(CONFIRM)
assert len(confirmation) == 6 and len({r["compound_id"] for r in confirmation}) == 6
assert all(r["compound_id"] in combined_by_id for r in confirmation)
assert all(r["inchi_key"] == combined_by_id[r["compound_id"]]["inchi_key"] and
           r["standardized_isomeric_smiles"] == combined_by_id[r["compound_id"]]["standardized_isomeric_smiles"]
           for r in confirmation)
confirm_rows = []
for row in confirmation:
    cid = row["compound_id"]
    background = combined_by_id[cid]
    confirm_rows.append({
        "confirmation_selection_rank": row["confirmation_selection_rank"],
        "confirmation_tier": row["confirmation_tier"],
        "pharmacology_batch": "formal278_20260929" if cid in old_ids else "active60_20260929",
        **background,
    })

write_csv(HERE / "active60_new_chembl_source_records.csv", source_rows,
          list(read_csv(PARENT / "formal278_chembl_source_records.csv")[0]))
write_csv(HERE / "active60_new_pharmacology_prefetch.csv", new_rows, columns)
write_csv(HERE / "active60_new_identity_warnings.csv", warnings,
          ["compound_id", "issue", "source_ids", "candidate_inchi_key"])
write_csv(HERE / "formal319_pharmacology_prefetch.csv", combined, columns)
write_csv(HERE / "active60_pharmacology_prefetch.csv", active_rows, columns)
write_csv(HERE / "confirmation6_pharmacology_join.csv", confirm_rows, list(confirm_rows[0]))

summary = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "active_unique_compounds": len(active),
    "active_reused_from_formal278": len(set(active_by_id) & old_ids),
    "active_new_compounds": len(new_rows),
    "new_queried_source_ids": len(ids),
    "new_api_full_key_match_source_ids": len(source_rows) - sum(bool(r["source_identity_relation"] != "api_full_inchikey_match") for r in source_rows),
    "new_api_key_mismatch_source_ids": sum(r["source_identity_relation"] != "api_full_inchikey_match" for r in source_rows),
    "new_mechanism_records": sum(len(v) for v in mechanisms.values()),
    "new_identity_warning_compounds": len(warnings),
    "combined_unique_compounds": len(combined),
    "confirmation_compounds": len(confirmation),
    "confirmation_old_batch": sum(r["compound_id"] in old_ids for r in confirmation),
    "confirmation_new_batch": sum(r["compound_id"] in new_ids for r in confirmation),
    "sha256": {p.name: sha(p) for p in [ACTIVE, CONFIRM, HERE / "incremental_query_plan.json",
              HERE / "incremental_chembl_request_manifest.json",
              HERE / "active60_new_chembl_source_records.csv",
              HERE / "active60_new_pharmacology_prefetch.csv",
              HERE / "formal319_pharmacology_prefetch.csv",
              HERE / "active60_pharmacology_prefetch.csv",
              HERE / "confirmation6_pharmacology_join.csv"]},
}
(HERE / "active60_pharmacology_manifest.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({k: v for k, v in summary.items() if k != "sha256"}, ensure_ascii=False, indent=2))
