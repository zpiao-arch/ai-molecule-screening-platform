import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "download_public_data.py"


def load_module():
    spec = importlib.util.spec_from_file_location("download_public_data", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class PublicDataDownloaderTests(unittest.TestCase):
    def test_dry_run_writes_manifest_without_network_files(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "data"
            code = module.main(["--out", str(out), "--pdb", "3TI6", "--dry-run"])
            self.assertEqual(code, 0)
            manifest = out / "manifest.json"
            self.assertTrue(manifest.exists())
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(payload["mode"], "dry_run")
            self.assertEqual(payload["pdb_ids"], ["3TI6"])
            self.assertFalse((out / "rcsb" / "3TI6.pdb").exists())

    def test_source_summary_includes_stable_public_sources(self):
        module = load_module()
        sources = {item["id"]: item for item in module.public_sources()}
        for required in ["rcsb_pdb", "pubchem", "chembl", "bindingdb", "drugcentral", "open_targets"]:
            self.assertIn(required, sources)
        self.assertEqual(sources["rcsb_pdb"]["runtime_policy"], "download_selected_cache")


if __name__ == "__main__":
    unittest.main()
