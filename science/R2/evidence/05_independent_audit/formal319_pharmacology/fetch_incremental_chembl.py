"""Fetch only the frozen active60 additions from the official ChEMBL API."""

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
PLAN = json.loads((HERE / "incremental_query_plan.json").read_text(encoding="utf-8"))
IDS = PLAN["new_exact_source_ids_to_query"]
RAW = HERE / "chembl_raw"
RAW.mkdir(exist_ok=True)
BASE = "https://www.ebi.ac.uk/chembl/api/data"
LIMIT = 1000
BATCH_SIZE = 25


def fetch(url):
    req = Request(url, headers={"Accept": "application/json", "User-Agent": "TMC1-candidate-source-audit/1.0"})
    for attempt in range(4):
        try:
            with urlopen(req, timeout=60) as response:
                if response.status != 200:
                    raise RuntimeError(f"HTTP {response.status}: {url}")
                return response.read(), response.status
        except (HTTPError, URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(min(2 ** (attempt + 1), 12))


manifest_path = HERE / "incremental_chembl_request_manifest.json"
if manifest_path.exists():
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["exact_id_count"] != len(IDS):
        raise ValueError("Previously fetched manifest has a different source ID set")
else:
    manifest = {"source": BASE, "exact_id_count": len(IDS), "batch_size": BATCH_SIZE,
                "limit": LIMIT, "requests": []}

for resource, field in (("molecule", "molecules"), ("mechanism", "mechanisms")):
    for batch_idx, start in enumerate(range(0, len(IDS), BATCH_SIZE), start=1):
        ids = IDS[start:start + BATCH_SIZE]
        offset = 0
        page_idx = 1
        while True:
            url = f"{BASE}/{resource}.json?" + urlencode({
                "molecule_chembl_id__in": ",".join(ids), "limit": LIMIT, "offset": offset})
            raw_path = RAW / f"{resource}_batch_{batch_idx:02d}_page_{page_idx:02d}.json"
            existing = next((r for r in manifest["requests"] if r["resource"] == resource
                             and r["batch"] == batch_idx and r["page"] == page_idx), None)
            if existing:
                if existing["url"] != url or existing["request_ids"] != ids:
                    raise ValueError("Cached page metadata disagrees with current query")
                body = raw_path.read_bytes()
                if hashlib.sha256(body).hexdigest() != existing["raw_sha256"]:
                    raise ValueError("Cached response SHA mismatch")
            else:
                body, status = fetch(url)
                raw_path.write_bytes(body)
                time.sleep(0.5)
            data = json.loads(body)
            records = data[field]
            total = int(data["page_meta"]["total_count"])
            returned_ids = sorted({r["molecule_chembl_id"] for r in records})
            if not set(returned_ids) <= set(ids):
                raise ValueError("Official API returned an ID outside the exact request")
            if not existing:
                manifest["requests"].append({
                    "resource": resource, "batch": batch_idx, "page": page_idx,
                    "request_ids": ids, "url": url,
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                    "http_status": status, "raw_file": str(raw_path.relative_to(HERE)),
                    "raw_sha256": hashlib.sha256(body).hexdigest(),
                    "returned_record_count": len(records), "returned_ids": returned_ids,
                    "total_count_for_batch": total})
                manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
            print(f"{resource} {batch_idx} page {page_idx}: {len(records)}/{total}", flush=True)
            if offset + len(records) >= total:
                break
            if not records:
                raise ValueError("Unexpected empty API page")
            offset += len(records)
            page_idx += 1
