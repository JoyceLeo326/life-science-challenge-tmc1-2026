"""Replay the release in releases/LATEST.json into temporary output paths.

Default: standard-library reconstruction of three frozen CSVs.
--with-supplement: saved-source A1-A3 reconstruction, using requirements-ci.txt.
Neither path downloads receptors or runs docking, training or molecular dynamics.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
import time

RESULTS = ("results.csv", "summary_metrics.csv", "r3_summary.csv")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(directory):
    return {p.relative_to(directory).as_posix(): digest(p)
            for p in sorted(directory.rglob("*")) if p.is_file()}


def load_latest(repository, config):
    metadata = json.loads(config.read_text(encoding="utf-8-sig"))
    relative = Path(metadata["code_path"])
    code = (repository / relative).resolve()
    if relative.is_absolute() or not code.is_relative_to(repository.resolve()):
        raise ValueError("LATEST code_path must remain within this repository")
    if not code.is_dir():
        raise ValueError("LATEST code directory is missing")
    if digest(code / "MANIFEST.json") != metadata["code_manifest_sha256"]:
        raise ValueError("LATEST code_manifest_sha256 does not match the release")
    return metadata, code


def validate(repository, config, with_supplement=False, report_dir=None):
    started = time.perf_counter()
    metadata, code = load_latest(repository, config)
    before = snapshot(code)
    expected = {name: digest(code / name) for name in RESULTS}
    jobs, logs = [], {}
    report = {
        "status": "RUNNING", "release_version": metadata["version"],
        "code_path": metadata["code_path"],
        "code_manifest_sha256": metadata["code_manifest_sha256"],
        "source_code_zip_sha256": metadata.get("code_zip_sha256"),
        "platform": platform.system(), "python": platform.python_version(),
        "scope": "frozen result readback and optional saved-source supplement rebuild",
        "new_docking_training_or_MD": False,
        "receptor_download_required": False, "jobs": jobs,
    }
    env = os.environ.copy()
    env.update(PYTHONUTF8="1", PYTHONDONTWRITEBYTECODE="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(prefix="tmc1-result-check-") as temporary:
        output = Path(temporary)

        def clean(text):
            for path, label in ((output, "<temporary-output>"),
                                (code, metadata["code_path"]),
                                (repository, "<repository>"),
                                (Path(sys.prefix), "<python-environment>")):
                text = text.replace(str(path), label).replace(path.as_posix(), label)
            return text

        def run(label, args):
            job_started = time.perf_counter()
            result = subprocess.run([sys.executable, "-I", "-B", *map(str, args)],
                                    cwd=code, env=env, text=True, encoding="utf-8",
                                    capture_output=True, timeout=300)
            logs[label] = clean(result.stdout + result.stderr)
            jobs.append({"name": label, "command": ["python", "-I", "-B", *[clean(str(x)) for x in args]],
                         "exit_code": result.returncode,
                         "seconds": round(time.perf_counter() - job_started, 3)})
            if result.returncode:
                raise RuntimeError(label + " failed with exit code " + str(result.returncode))
            return result.stdout

        try:
            run("release_manifest", ["verify_manifest.py"])
            run("formal_results", ["make_results.py", "--output", output / "results.csv"])
            comparison = [{"file": name, "expected_sha256": expected[name],
                           "actual_sha256": digest(output / name),
                           "identical": digest(output / name) == expected[name]} for name in RESULTS]
            report["result_comparisons"] = comparison
            if not all(item["identical"] for item in comparison):
                raise ValueError("Rebuilt result CSV differs from the frozen release")
            if with_supplement:
                base = Path("data/group_supplement")
                verified = json.loads(run("supplement_verify", [base / "verify_science_supplement.py"]))
                if verified["failed"] != 0 or verified["passed"] <= 0:
                    raise ValueError("Supplement verification did not pass")
                report["supplement_check_count"] = verified["passed"]
                rebuilt = output / "supplement"
                run("supplement_rebuild", [base / "rebuild_science_tables.py", "--output", rebuilt])
                receipt = json.loads((rebuilt / "OFFLINE_REBUILD_RECEIPT.json").read_text(encoding="utf-8"))
                if receipt["status"] != "PASS" or not receipt["comparisons"] or any(
                        x["differing_shared_fields"] != 0 for x in receipt["comparisons"]):
                    raise ValueError("Rebuilt supplement differs from saved source tables")
                report["supplement_rebuild"] = receipt
            after = snapshot(code)
            changed = sorted(k for k in before.keys() | after.keys() if before.get(k) != after.get(k))
            report["code_files_checked"] = len(before)
            report["source_files_changed"] = changed
            if changed:
                raise ValueError("Replay changed release source files")
            report["status"] = "PASS"
        except Exception as exc:
            report["status"] = "FAIL"
            report["error"] = clean(str(exc))
        finally:
            report["seconds"] = round(time.perf_counter() - started, 3)
            if report_dir:
                report_dir.mkdir(parents=True, exist_ok=True)
                for label, log in logs.items():
                    (report_dir / (label + ".log")).write_text(log, encoding="utf-8")
                (report_dir / "LATEST_RELEASE_CHECK.json").write_text(
                    json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    repository = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=repository / "releases/LATEST.json")
    parser.add_argument("--with-supplement", action="store_true")
    parser.add_argument("--report-dir", type=Path, help="Save JSON receipt and sanitized logs outside release/code")
    args = parser.parse_args()
    if args.report_dir:
        _, code = load_latest(repository, args.config)
        if args.report_dir.resolve().is_relative_to(code):
            parser.error("--report-dir must be outside the frozen code directory")
    report = validate(repository, args.config, args.with_supplement, args.report_dir)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
