"""Create auditable, path-only replay copies of the frozen 03/04 methods.

This tool does not run docking, training, ranking, acquisition, or evaluation.
It refuses unknown source layouts and writes original SHA256 plus unified diffs.
"""
from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import shutil
from pathlib import Path
import os


OLD_FIRST = r"<PACKAGE_ROOT>"
OLD_SCRATCH_03 = r"<PACKAGE_ROOT>\Documents\Codex\TMC1_R2_Scratch_20260929\03_compute"
OLD_SCRATCH_04 = r"<PACKAGE_ROOT>\Documents\Codex\TMC1_R2_Scratch_20260929\04_design"
LABEL_SHA = "65b19d87b9c0ee25e009e439217b35df1001c22c78129d3282f34b7820c36ea4"
OLD_LIB_SHA = "484b0dddc59531abf6e0cd53414be0ee3ccf482aeb222e75797eedf2c0d9d809"
PARENT_SHA = "8445df639801042763fe2cc6dfedc77045c59125aec50807f7fd19239c08062c"
SCORES_SHA = "d094ad3d1f0136f5c64bf2902c3a8a1a756b06d3894a4e819082694c7177aab8"
LEGACY_POOL_SHA = "9add7f2a901293c7d2be2f748d17f08e35bcc1dc6f47c7573fcc5d6fbbd15916"
BOX_SHA = "3e4c7eb2c313edbe496068c9fb17790e7ca746567519836ed46131e5e7570d49"
VINA_SHA = "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5"
RECEPTOR_SHA = "c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1"
STANDBY_SHA = "d0430853a75f7de8171dd2abdba3d0c9b4a02566efc2d048a301a2dfcbec6cf1"
INTERRUPTED_ACTIVE_SHA = "95156f2f5f100bc817df714e98dc31b7ba0a3c41e5b7c8b4741c3db83dd5aead"
CHARGE_SHA = "f1a83142364dac1a1dc2bbe58cb67c4158acf93c9ba6f9ebc2b94e6f46946f73"
PILOT_SHA = "1d19f7b911e533db8a01a78c5f44e7c669d0cff6b56719d70c6848c3e57446fc"
CONSTRAINTS_SHA = "2a46fc4c47441b6773cc95c72311d4ff1b7e6ea1399b79657d83bdf853629f9b"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def checked(path: Path, expected: str | None = None) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required replay input is missing: {path}")
    if expected and sha(path).lower() != expected.lower():
        raise ValueError(f"SHA256 mismatch for {path}; expected {expected}, got {sha(path)}")
    return path


def replace_once(source: str, old: str, new: str, label: str) -> str:
    expected = 2 if label.endswith("both timing functions") else 1
    count = source.count(old)
    if count != expected:
        raise ValueError(f"Source layout drift at {label}: expected {expected} matches, found {count}")
    return source.replace(old, new, expected)


def literal(path: Path) -> str:
    return f"Path({str(path)!r})"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r1", required=True, type=Path, help="Extracted tmc1_screening root or SHA-checked 14-dependency root")
    p.add_argument("--r2", required=True, type=Path, help="Extracted R2 package root")
    p.add_argument("--work", required=True, type=Path, help="New writable replay directory")
    a = p.parse_args()
    r1, r2, work = a.r1.resolve(), a.r2.resolve(), a.work.resolve()
    if work.exists() and any(work.iterdir()):
        raise FileExistsError(f"Replay work directory must be empty: {work}")
    if r1 == work or r2 == work or work.is_relative_to(r1) or work.is_relative_to(r2):
        raise ValueError("--work must be separate from the read-only R1 and R2 roots")
    inputs = {
        "historical_labels": (r1 / "coarse_v1/analysis/job_audit_join.csv", LABEL_SHA),
        "historical_library": (r1 / "coarse_v1/inputs/library/selected_candidates.csv", OLD_LIB_SHA),
        "parent_input": (r1 / "confirmation_v1/inputs/library/selected_candidates.csv", PARENT_SHA),
        "parent_scores": (r1 / "confirmation_v1/analysis/parent_setting_complete_denominators.csv", SCORES_SHA),
        "legacy_eligible": (r1 / "screen_v2/inputs/library/eligible_pool.csv", LEGACY_POOL_SHA),
        "box": (r1 / "coarse_v1/config/vina_box.txt", BOX_SHA),
        "standby_intervals": (r2 / "03_compute/method/standby_intervals_20260929.csv", STANDBY_SHA),
        "interrupted_active": (r2 / "03_compute/method/active_settings_interrupt_checkpoint_20260929/interrupted_inflight_attempts.csv", INTERRUPTED_ACTIVE_SHA),
        "charge_warnings": (r2 / "03_compute/method/formal_charge_warning_ids.csv", CHARGE_SHA),
        "excluded_pilot": (r2 / "03_compute/method/runtime_pilot_12.csv", PILOT_SHA),
        "receptor": (r2 / "03_compute/receptors/M_PUB_ALL.pdbqt", RECEPTOR_SHA),
        "design_constraints": (r2 / "04_design/method_results/design_constraints_all.csv", CONSTRAINTS_SHA),
    }
    for path, expected in inputs.values():
        checked(path, expected)
    # Everything above is read-only. Create the work tree only after identity checks.
    d03, d04 = work / "portable_method/03_compute", work / "portable_method/04_design"
    d03.mkdir(parents=True)
    d04.mkdir(parents=True)
    rule_output = work / "04_design/rule_design"
    rule_output.mkdir(parents=True)
    shutil.copy2(checked(r2 / "03_compute/method/chemistry_rules.py"), d03 / "chemistry_rules.py")
    for name in ("formal_charge_warning_ids.csv", "standby_intervals_20260929.csv"):
        shutil.copy2(checked(r2 / "03_compute/method" / name), d03 / name)
    interrupt_src = checked(r2 / "03_compute/method/active_settings_interrupt_checkpoint_20260929/interrupted_inflight_attempts.csv")
    interrupt_dst = d03 / "active_settings_interrupt_checkpoint_20260929" / "interrupted_inflight_attempts.csv"
    interrupt_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(interrupt_src, interrupt_dst)
    shutil.copy2(checked(r2 / "04_design/method_results/design_constraints_all.csv"), d04 / "design_constraints_all.csv")
    records = []

    def emit(group: str, name: str, patches: list[tuple[str, str, str]]) -> None:
        src = checked(r2 / group / name)
        target = (d03 if group.startswith("03") else d04) / name
        original = src.read_text(encoding="utf-8")
        edited = original
        for old, new, label in patches:
            edited = replace_once(edited, old, new, f"{name}: {label}")
        target.write_text(edited, encoding="utf-8", newline="\n")
        diff = "".join(difflib.unified_diff(original.splitlines(True), edited.splitlines(True),
                                            fromfile=f"original/{group}/{name}", tofile=f"portable/{name}"))
        diff_path = work / "diffs" / f"{name}.diff"
        diff_path.parent.mkdir(parents=True, exist_ok=True)
        diff_path.write_text(diff, encoding="utf-8")
        records.append({"source": str(src), "source_sha256": sha(src), "portable": str(target),
                        "portable_sha256": sha(target), "diff": str(diff_path)})

    first_old = f'FIRST = Path(r"{OLD_FIRST}")'
    first_new = f"FIRST = {literal(r1)}"
    emit("03_compute/method", "rank_library.py", [(first_old, first_new, "first-round root")])
    emit("03_compute/method", "active_acquire.py", [])
    emit("03_compute/method", "dock_batch.py", [
        (first_old, first_new, "first-round root"),
        (f'DEFAULT_SCRATCH = Path(r"{OLD_SCRATCH_03}")', f"DEFAULT_SCRATCH = {literal(work / '03_compute')}", "scratch"),
        ('VINA = FIRST / "tools/vina_1.2.7_win.exe"', 'VINA = None  # set from --vina-bin after argument parsing', "Vina location"),
        ('EXPORT = FIRST / ".venv/Scripts/mk_export.exe"', 'EXPORT = None  # set from --export-bin after argument parsing', "Meeko export location"),
        ('    ap.add_argument("--box", type=Path, default=BOX)',
         '    ap.add_argument("--receptor-dir", required=True, type=Path)\n'
         '    ap.add_argument("--vina-bin", required=True, type=Path)\n'
         '    ap.add_argument("--export-bin", required=True, type=Path)\n'
         '    ap.add_argument("--box", type=Path, default=BOX)', "docking path arguments"),
        ('    a = ap.parse_args()\n    root = a.scratch / a.batch',
         '    a = ap.parse_args()\n    global VINA, EXPORT\n    VINA, EXPORT = a.vina_bin, a.export_bin\n'
         f'    if sha(VINA) != "{VINA_SHA}":\n        raise ValueError("Vina binary SHA256 mismatch")\n'
         '    if not EXPORT.is_file():\n        raise FileNotFoundError(f"Meeko mk_export executable missing: {EXPORT}")\n'
         f'    if sha(a.box) != "{BOX_SHA}":\n        raise ValueError("Vina box SHA256 mismatch")\n'
         '    root = a.scratch / a.batch', "binary and box identity checks"),
        ('    receptor = FIRST / f"screen_v2/inputs/receptors/{a.receptor}.pdbqt"',
         '    receptor = a.receptor_dir / f"{a.receptor}.pdbqt"\n'
         f'    if a.receptor == "M_PUB_ALL" and sha(receptor) != "{RECEPTOR_SHA}":\n'
         '        raise ValueError("M_PUB_ALL receptor SHA256 mismatch")', "receptor path and identity"),
    ])
    emit("03_compute/method", "evaluate_forward.py", [
        ('def batch_wall_seconds(root: Path) -> float:',
         'def batch_wall_seconds(root: Path, old_override: Path | None = None, resumed_override: Path | None = None) -> float:', "initial timing function arguments"),
        ('def batch_awake_wall_seconds(root: Path) -> float:',
         'def batch_awake_wall_seconds(root: Path, old_override: Path | None = None, resumed_override: Path | None = None) -> float:', "awake timing function arguments"),
        ('old=Path(meta["source_paths"]["old"])',
         'old=old_override or Path(meta["source_paths"]["old"])', "old batch per merged root in both timing functions"),
        ('new=Path(meta["source_paths"]["resumed"])',
         'new=resumed_override or Path(meta["source_paths"]["resumed"])', "resumed batch per merged root in both timing functions"),
        ('    p.add_argument("--initial-batch",required=True,type=Path)',
         '    p.add_argument("--initial-old-batch",required=True,type=Path,help="Initial combination original 154-job physical batch")\n'
         '    p.add_argument("--initial-resumed-batch",required=True,type=Path,help="Initial combination resumed 124-job physical batch")\n'
         '    p.add_argument("--active-old-batch",required=True,type=Path,help="Active combination original 27-job physical batch")\n'
         '    p.add_argument("--active-resumed-batch",required=True,type=Path,help="Active combination resumed 14-job physical batch")\n'
         '    p.add_argument("--initial-batch",required=True,type=Path)', "batch arguments"),
        ('    a=p.parse_args()\n    meta=',
         '    a=p.parse_args()\n'
         '    for batch in (a.initial_old_batch, a.initial_resumed_batch, a.active_old_batch, a.active_resumed_batch):\n'
         '        if not (batch/"plan.json").is_file() or not (batch/"job_ledger.csv").is_file():\n'
         '            raise FileNotFoundError(f"Original physical batch plan/job ledger missing: {batch}")\n'
         '    meta=', "batch validation"),
        ('"initial_batch_wall_seconds":batch_wall_seconds(a.initial_batch),',
         '"initial_batch_wall_seconds":batch_wall_seconds(a.initial_batch, a.initial_old_batch, a.initial_resumed_batch),', "initial merged timing mapping"),
        ('"initial_batch_awake_wall_seconds":batch_awake_wall_seconds(a.initial_batch),',
         '"initial_batch_awake_wall_seconds":batch_awake_wall_seconds(a.initial_batch, a.initial_old_batch, a.initial_resumed_batch),', "initial awake timing mapping"),
        ('"active_batch_wall_seconds":batch_wall_seconds(a.active_batch) if new_active_ids else 0,',
         '"active_batch_wall_seconds":batch_wall_seconds(a.active_batch, a.active_old_batch, a.active_resumed_batch) if new_active_ids else 0,', "active merged timing mapping"),
        ('"active_batch_awake_wall_seconds":batch_awake_wall_seconds(a.active_batch) if new_active_ids else 0,',
         '"active_batch_awake_wall_seconds":batch_awake_wall_seconds(a.active_batch, a.active_old_batch, a.active_resumed_batch) if new_active_ids else 0,', "active awake timing mapping"),
    ])
    emit("04_design/method_results", "run_design.py", [
        ('HERE = Path(__file__).resolve().parent', f"HERE = {literal(rule_output)}", "isolated rule-design results"),
        (f'BASE = Path(r"{OLD_FIRST}\\confirmation_v1")', f"BASE = {literal(r1 / 'confirmation_v1')}", "parent input root"),
    ])
    emit("04_design/method_results", "train_smiles_gru.py", [
        (f'SCRATCH = Path(r"{OLD_SCRATCH_04}")', f"SCRATCH = {literal(work / '04_design')}", "model scratch"),
        (f'PARENT_INPUT = Path(r"{OLD_FIRST}\\confirmation_v1\\inputs\\library\\selected_candidates.csv")',
         f"PARENT_INPUT = {literal(inputs['parent_input'][0])}", "parent input before import"),
    ])
    emit("04_design/method_results", "train_selfies_gru.py", [
        (f'SCRATCH = Path(r"{OLD_SCRATCH_04}")', f"SCRATCH = {literal(work / '04_design')}", "model scratch"),
        (f'legacy_path = Path(r"{OLD_FIRST}\\screen_v2\\inputs\\library\\eligible_pool.csv")',
         f"legacy_path = {literal(inputs['legacy_eligible'][0])}", "legacy exclusion pool"),
        ('    enumeration_path = HERE/"design_constraints_all.csv"',
         '    enumeration_path = HERE/"design_constraints_all.csv"\n'
         f'    if hashlib.sha256(enumeration_path.read_bytes()).hexdigest() != "{CONSTRAINTS_SHA}":\n'
         '        raise ValueError("Frozen rule-enumeration SHA256 mismatch")', "frozen rule enumeration guard"),
    ])
    manifest = {"role": "portable path-only replay copies, not new scientific results",
                "r1": str(r1), "r2": str(r2), "work": str(work),
                "inputs": {k: {"path": str(path), "sha256": sha(path)} for k, (path, _) in inputs.items()},
                "source_and_diff": records}
    (work / "PORTABLE_MANIFEST.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": "prepared", "scripts": len(records), "manifest": str(work / "PORTABLE_MANIFEST.json")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
