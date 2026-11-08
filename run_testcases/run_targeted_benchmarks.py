#!/usr/bin/env python3
"""Run the THREE targeted workloads using the existing validated benchmark runner.

Does NOT touch the main experiment's benchmark_runs.csv / benchmark_summary.csv.
Run from the project root:
    python3 run_testcases/run_targeted_benchmarks.py
    python3 run_testcases/run_targeted_benchmarks.py --resume
"""

import argparse
import subprocess
import sys
from pathlib import Path

ALGORITHMS = ("fmi", "ssi", "psi", "sss", "pss", "ssv", "psv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--limit", type=int, help="Run only the first N testcases for a smoke test")
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--algorithms", nargs="+", choices=ALGORITHMS, default=list(ALGORITHMS))
    args = parser.parse_args()

    project = Path(__file__).resolve().parent.parent
    generator = project / "targeted_testcases"
    if not (generator / "manifest.csv").is_file() or not list(generator.glob("*.in")):
        parser.error("Generate targeted testcases first: python3 gen_testcases/generate_targeted_benchmarks.py")
    runner = Path(__file__).resolve().with_name("run_benchmarks.py")
    if not runner.is_file():
        parser.error(f"Missing existing runner: {runner}")

    cmd = [
        sys.executable, str(runner),
        "--testcases", str(generator),
        "--exe", str(project / "build" / "main"),
        "--out", str(project / "results" / "targeted_benchmark_runs.csv"),
        "--summary", str(project / "results" / "targeted_benchmark_summary.csv"),
        "--repeats", str(args.repeats),
        "--timeout", str(args.timeout),
        "--algorithms", *args.algorithms,
    ]
    if args.limit is not None:
        cmd.extend(("--limit", str(args.limit)))
    if args.resume:
        cmd.append("--resume")

    print("Targeted benchmark results will be saved separately from the main campaign.", flush=True)
    result = subprocess.run(cmd, cwd=project, check=False)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
