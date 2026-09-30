"""Hash the immutable library snapshot and its generating scripts."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import os

ROOT = Path(os.environ.get("TMC1_LIBRARY_WORK", str(Path(__file__).resolve().parent)))
SNAP = ROOT / "snapshots" / "v2_20000_new"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


files = {p.name: {"bytes": p.stat().st_size, "sha256": sha(p)} for p in SNAP.iterdir()
         if p.is_file() and p.name != "snapshot_manifest.json"}
scripts = {p.name: sha(p) for p in (ROOT / "build_library.py", ROOT / "make_overlap.py",
                                    ROOT / "verify_and_profile.py", ROOT / "add_fullprops_sdf.py")}
summary = json.loads((SNAP / "summary.json").read_text(encoding="utf-8"))
overlap = json.loads((SNAP / "historical_overlap_summary.json").read_text(encoding="utf-8"))
profile = json.loads((SNAP / "quality_profile.json").read_text(encoding="utf-8"))
manifest = {"created_at_utc": datetime.now(timezone.utc).isoformat(), "snapshot": "v2_20000_new",
            "status": "frozen_core_files_plus_append_only_fullprops_sdf", "files": files, "script_sha256": scripts,
            "counts": {"source_records": summary["source_records_total"],
                       "unique_source_ids": summary["unique_source_ids"],
                       "unique_parents": summary["standardized_unique_parents"],
                       "eligible_parents": summary["property_filtered_parents"],
                       "forward_unlabeled_pool": overlap["forward_unlabeled_pool"],
                       "sdf_records": profile["sdf_records"]},
            "qc_pass": profile["qc_pass"], "license": summary["license"],
            "license_url": summary["license_url"]}
(SNAP / "snapshot_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(manifest["counts"], ensure_ascii=False))
