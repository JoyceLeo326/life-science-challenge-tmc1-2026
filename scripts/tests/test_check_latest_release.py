"""Small fixtures test the checker, never stand in for scientific replay."""
from pathlib import Path
import importlib.util
import json
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("release_checker", Path(__file__).parents[1] / "check_latest_release.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class ReleaseCheckTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.code = self.root / "releases/fixture/code"
        self.code.mkdir(parents=True)
        self.config = self.root / "releases/LATEST.json"
        for name in checker.RESULTS:
            (self.code / name).write_bytes(b"fixture,value\r\nexample,1\r\n")
        (self.code / "MANIFEST.json").write_text('{"files": []}', encoding="utf-8")
        (self.code / "verify_manifest.py").write_text("print('fixture manifest pass')\n", encoding="utf-8")
        (self.code / "make_results.py").write_text(
            "from pathlib import Path\nimport argparse,shutil\n"
            "p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args()\n"
            "a.output.parent.mkdir(parents=True,exist_ok=True)\n"
            "for n in ('results.csv','summary_metrics.csv','r3_summary.csv'):shutil.copyfile(Path(__file__).parent/n,a.output.parent/n)\n",
            encoding="utf-8")
        self.metadata = {"version":"test-fixture", "code_path":"releases/fixture/code",
                         "code_manifest_sha256":checker.digest(self.code / "MANIFEST.json")}
        self.write_config()

    def write_config(self):
        self.config.write_text(json.dumps(self.metadata), encoding="utf-8")

    def test_dynamic_pointer_actual_execution_and_no_source_change(self):
        report = checker.validate(self.root, self.config)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(len(report["jobs"]), 2)
        self.assertTrue(all(c["identical"] for c in report["result_comparisons"]))
        self.assertEqual(report["source_files_changed"], [])

    def test_result_drift_fails(self):
        with (self.code / "make_results.py").open("a", encoding="utf-8") as stream:
            stream.write("a.output.write_text('changed fixture',encoding='utf-8')\n")
        report = checker.validate(self.root, self.config)
        self.assertEqual(report["status"], "FAIL")
        self.assertFalse(report["result_comparisons"][0]["identical"])

    def test_manifest_pointer_mismatch_fails(self):
        self.metadata["code_manifest_sha256"] = "0" * 64
        self.write_config()
        with self.assertRaisesRegex(ValueError, "manifest_sha256"):
            checker.validate(self.root, self.config)

    def test_pointer_cannot_leave_repository(self):
        self.metadata["code_path"] = "../../outside"
        self.write_config()
        with self.assertRaisesRegex(ValueError, "within this repository"):
            checker.validate(self.root, self.config)

    def test_replay_source_mutation_fails(self):
        with (self.code / "make_results.py").open("a", encoding="utf-8") as stream:
            stream.write("Path('unexpected-output.txt').write_text('fixture',encoding='utf-8')\n")
        report = checker.validate(self.root, self.config)
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["source_files_changed"], ["unexpected-output.txt"])


if __name__ == "__main__":
    unittest.main()
