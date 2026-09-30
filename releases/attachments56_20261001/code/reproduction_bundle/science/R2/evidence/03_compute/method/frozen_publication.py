"""Exact publication hashes mapped to original frozen ledger identities."""
import hashlib
import json
from pathlib import Path
from portable_config import required_path

PUBLICATION_TABLE_SHA = "ea75d204141b8b6170172cb64513747729f82970f9ccdf1d2085a90c5e227226"

def frozen_matches(path, original_expected):
    path = Path(path)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual == original_expected:
        return True
    table_path = Path(__file__).with_name("publication_hashes.json")
    if hashlib.sha256(table_path.read_bytes()).hexdigest() != PUBLICATION_TABLE_SHA:
        raise ValueError("Publication hash mapping table changed")
    mapping = json.loads(table_path.read_text(encoding="utf-8"))
    relative = path.resolve().relative_to(required_path("TMC1_R2_ROOT")).as_posix()
    entry = mapping.get(relative)
    return bool(entry and entry["original_expected_sha256"] == original_expected
                and entry["publication_sha256"] == actual)
