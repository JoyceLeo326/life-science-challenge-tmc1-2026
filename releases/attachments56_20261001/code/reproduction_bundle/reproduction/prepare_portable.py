"""Verify frozen replay inputs and prepare configured copies without calculations.

This publication entry accepts restored method copies by their own frozen hashes.
Original scientific-source hashes are provenance, not claims about edited bytes.
No model weights are loaded, and no scientific result is recomputed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


R1_INPUTS = {
    "historical_labels": ("coarse_v1/analysis/job_audit_join.csv", "65b19d87b9c0ee25e009e439217b35df1001c22c78129d3282f34b7820c36ea4"),
    "historical_library": ("coarse_v1/inputs/library/selected_candidates.csv", "484b0dddc59531abf6e0cd53414be0ee3ccf482aeb222e75797eedf2c0d9d809"),
    "parent_input": ("confirmation_v1/inputs/library/selected_candidates.csv", "8445df639801042763fe2cc6dfedc77045c59125aec50807f7fd19239c08062c"),
    "parent_scores": ("confirmation_v1/analysis/parent_setting_complete_denominators.csv", "d094ad3d1f0136f5c64bf2902c3a8a1a756b06d3894a4e819082694c7177aab8"),
    "legacy_eligible": ("screen_v2/inputs/library/eligible_pool.csv", "9add7f2a901293c7d2be2f748d17f08e35bcc1dc6f47c7573fcc5d6fbbd15916"),
    "box": ("coarse_v1/config/vina_box.txt", "3e4c7eb2c313edbe496068c9fb17790e7ca746567519836ed46131e5e7570d49"),
    "box_metadata": ("coarse_v1/config/common_box.json", "2330b07cf13cb289f845344898fcda828c2f6f6024be55af8dc75e45c3e52dca"),
}
R2_LIBRARY_INPUTS = {
    "prospective_pool": ("02_library/formal_snapshot/forward_unlabeled_pool.csv", "8fd8424797e441a3921ba3805dd8e9bf87d3421e5fe3440a32211753ec964c1e"),
    "generator_training_library": ("02_library/formal_snapshot/filtered_for_screening.csv", "439f74756bc786709f318a4379037c4df186bc4d40957bdbe245910ba48e8229"),
}
RECEPTOR_SHA = "c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1"
# Frozen identity of the completed publication table.
METHOD_TABLE_SHA = "87c3e8527aaa382b24c01d1c3944497c56ee6212a57aa7b53af79b0657679be3"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def checked(path: Path, expected: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required replay input is missing: {path}")
    actual = sha(path)
    if actual != expected:
        raise ValueError(f"SHA256 mismatch: {path}; expected {expected}, got {actual}")
    return path


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r1", type=Path, default=repo / "science/R1", help="Root containing the seven frozen historical inputs, including box metadata")
    p.add_argument("--r2", type=Path, default=repo / "science/R2/evidence", help="Root containing restored 02/03/04 evidence")
    p.add_argument("--work", required=True, type=Path, help="Separate, empty writable replay directory")
    p.add_argument("--receptor", required=True, type=Path, help="Original local M_PUB_ALL.pdbqt; never bypasses its frozen SHA")
    p.add_argument("--check-only", action="store_true", help="Read and hash configuration inputs without creating the work directory")
    a = p.parse_args()
    r1, r2, work = a.r1.resolve(), a.r2.resolve(), a.work.resolve()
    for root in (r1, r2):
        if work == root or work.is_relative_to(root) or root.is_relative_to(work):
            raise ValueError("--work must be separate from the read-only input roots")
    if a.receptor.resolve() == work or a.receptor.resolve().is_relative_to(work):
        raise ValueError("--receptor must be outside the writable work directory")
    if work.exists() and (not work.is_dir() or any(work.iterdir())):
        raise FileExistsError(f"Replay work directory must be empty: {work}")
    table = checked(Path(__file__).with_name("R2_METHOD_HASHES.json"), METHOD_TABLE_SHA)
    methods = json.loads(table.read_text(encoding="utf-8"))
    inputs = {}
    for label, (relative, expected) in R1_INPUTS.items():
        path = checked(r1 / relative, expected)
        inputs[label] = {"root": "R1", "relative_path": relative, "sha256": sha(path)}
    for label, (relative, expected) in R2_LIBRARY_INPUTS.items():
        path = checked(r2 / relative, expected)
        inputs[label] = {"root": "R2", "relative_path": relative, "sha256": sha(path)}
    checked(a.receptor, RECEPTOR_SHA)
    inputs["receptor"] = {"root": "external local input", "sha256": RECEPTOR_SHA, "redistributed": False}
    for relative, expected in methods.items():
        rel = Path(relative)
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError(f"Unsafe method table path: {relative}")
        checked(r2 / rel, expected)
    if a.check_only:
        print(json.dumps({"status": "configuration_checked", "method_and_support_files": len(methods), "inputs": inputs}, indent=2))
        return

    # Mutations begin only after every input and every publication source passes.
    configs = {
        "TMC1_R1_ROOT": str(r1), "TMC1_R2_ROOT": str(r2),
        "TMC1_WORK_ROOT": str(work), "TMC1_RULE_OUTPUT": str(work / "04_design/rule_design"),
        "TMC1_PREFLIGHT_ROOT": str(work / "ranking_replay"),
    }
    destinations = {
        "03_compute/method/": work / "portable_method/03_compute",
        "04_design/method_results/": work / "portable_method/04_design",
    }
    records = []
    for relative, expected in methods.items():
        for prefix, folder in destinations.items():
            if relative.startswith(prefix):
                target = folder / relative.removeprefix(prefix)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(r2 / relative, target)
                checked(target, expected)
                records.append({"source_relative_path": relative, "publication_source_sha256": expected,
                                "work_relative_path": target.relative_to(work).as_posix(), "work_copy_sha256": sha(target)})
                break
    for folder in destinations.values():
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "portable_paths.json").write_text(json.dumps(configs, indent=2), encoding="utf-8", newline="\n")
    (work / "04_design/rule_design").mkdir(parents=True, exist_ok=True)
    manifest = {"role": "Configured publication-method copies; no scientific computation executed",
                "source_identity": "Published source hashes verified; original scientific source hashes are documented in SOURCE_MANIFEST.json",
                "method_table_sha256": METHOD_TABLE_SHA, "inputs": inputs, "files": records, "local_configuration": configs,
                "receptor_path": str(a.receptor.resolve()), "model_weights_loaded": False}
    (work / "PORTABLE_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8", newline="\n")
    print(json.dumps({"status": "prepared", "files": len(records), "manifest": str(work / "PORTABLE_MANIFEST.json")}))


if __name__ == "__main__":
    main()
