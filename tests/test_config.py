from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from hcs_flapping.config import load_config


class ConfigTests(unittest.TestCase):
    def test_missing_config_uses_empty_defaults(self) -> None:
        with TemporaryDirectory() as directory:
            config = load_config(Path(directory) / "missing.toml")
        self.assertIsNone(config.source)
        self.assertEqual(config.get("project", "random_seed", 20210117), 20210117)

    def test_toml_paths_are_loaded(self) -> None:
        with TemporaryDirectory() as directory:
            config_file = Path(directory) / "config.toml"
            config_file.write_text('[paths]\noutput = "products"\n', encoding="utf-8")
            config = load_config(config_file)
        self.assertEqual(config.path("output"), Path("products"))


if __name__ == "__main__":
    unittest.main()
