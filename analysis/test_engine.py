"""Exact fixture checks for executable modes; no timing assertions."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class EngineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.builds = []
        for features in ("", "u32", "u16", "u32,wide"):
            target = ROOT / "enum/target/fixture-tests" / (features.replace(",", "-") or "default")
            cmd = ["cargo", "build", "--quiet", "--manifest-path", str(ROOT / "enum/Cargo.toml"),
                   "--target-dir", str(target)]
            if features:
                cmd += ["--features", features]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(result.stderr)
            cls.builds.append((features, target / "debug"))
        cls.bulk = json.loads((ROOT / "analysis/m0_order15.json").read_text())
        cls.wall = json.loads((ROOT / "analysis/wall_order30.json").read_text())

    def run_engine(self, directory, binary, args, ref, profile=False):
        env = {k: v for k, v in os.environ.items() if not k.startswith("DP_")}
        if profile:
            env["DP_PROFILE"] = "1"
        exe = directory / (binary + (".exe" if os.name == "nt" else ""))
        result = subprocess.run([str(exe)] + args, env=env, capture_output=True, text=True, check=True)
        if binary == "dpenum":
            self.assertEqual("profile:" in result.stderr, profile)
            if profile:
                lines = [line for line in result.stderr.splitlines() if line.startswith("profile:")]
                # Counters reset between primes, and worker-local contributions
                # are merged. The two primes have identical state structure here.
                self.assertEqual(len(lines), 2)
                counts = [line.split("bound calls ")[1] for line in lines]
                self.assertEqual(counts[0], counts[1])
                self.assertGreater(int(counts[0].split()[0]), 0)
        records = [json.loads(line) for line in result.stdout.splitlines()]
        self.assertEqual(len(records), 2)
        self.assertEqual(len({rec["prime"] for rec in records}), 2)
        for rec in records:
            self.assertEqual(rec["n"], len(ref["S"])-1)
            for name, values in ref.items():
                if isinstance(values, list) and not (name == "mu20" and "--no-x" in args):
                    self.assertEqual(rec[name], [x % rec["prime"] for x in values], name)

    def test_bulk_coefficients_across_execution_and_diagnostic_modes(self):
        modes = ((["--threads", "1"], False),
                 (["--threads", "1", "--park", "-v"], True),
                 (["--threads", "3", "--block", "3"], False),
                 (["--threads", "3", "--block", "3", "-v"], False),
                 (["--threads", "3", "--block", "3", "-v"], True),
                 (["--threads", "2", "--no-mirror"], False),
                 (["--threads", "2", "--loose", "--no-x"], False))
        for features, directory in self.builds:
            for model in ("bond", "site"):
                for args, profile in modes:
                    with self.subTest(features=features, model=model, args=args, profile=profile):
                        self.run_engine(directory, "dpenum", ["-n", "15", "--"+model,
                                        "--primes", "2"] + args, self.bulk[model], profile)

    def test_wall_coefficients_in_memory_and_on_disk(self):
        for features, directory in self.builds:
            with tempfile.TemporaryDirectory() as td:
                for args in (["--threads", "1"], ["--threads", "3", "--block", "2", "-v"],
                             ["--threads", "3", "--block", "3", "--spill", td,
                              "--spill-above", "0", "--shards", "8"]):
                    with self.subTest(features=features, args=args):
                        self.run_engine(directory, "dpwall", ["-n", "30", "--primes", "2"]+args, self.wall)
                self.assertEqual(list(Path(td).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
