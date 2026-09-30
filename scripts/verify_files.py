"""Verify the checked-in public file inventory without running scientific jobs."""
from pathlib import Path
import hashlib
import json
import sys

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / "PUBLIC_FILE_MANIFEST.json").read_text(encoding="utf-8"))
problems = []
for entry in manifest["files"]:
    path = (root / entry["path"]).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        problems.append({"path": entry["path"], "reason": "missing or invalid path"})
        continue
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    if path.stat().st_size != entry["bytes"] or digest.hexdigest() != entry["sha256"]:
        problems.append({"path": entry["path"], "reason": "content differs"})
print(json.dumps({"ok": not problems, "checked": len(manifest["files"]), "problems": problems}, ensure_ascii=True))
sys.exit(bool(problems))
