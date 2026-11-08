#!/usr/bin/env python3
"""Benchmark the seven exact-DNA-search implementations.

Run from the project root after building the Release executable:
    python3 run_testcases/run_benchmarks.py
    python3 run_testcases/run_benchmarks.py --limit 2 --repeats 1   # quick check
    python3 run_testcases/run_benchmarks.py --resume                # continue interrupted run

Expected executable usage:
    ./build/main <fmi|ssi|psi|sss|pss|ssv|psv> < testcase.in
"""

import argparse
import csv
import random
import re
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ALGORITHMS = ("fmi", "ssi", "psi", "sss", "pss", "ssv", "psv")
FILE_RE = re.compile(r"(?P<genome>.+)__m(?P<m>\d+)__q(?P<q>\d+)__k(?P<k>\d+)\.in$")
PEAK_TAG = "BENCH_PEAK_RSS_KIB"
PROJECT_ROOT = Path(__file__).resolve().parent.parent

RAW_FIELDS = (
    "genome", "reference_length", "m", "q", "k", "algorithm", "trial",
    "construction_ms", "search_ms", "total_ms", "index_rss_increase_mib",
    "peak_rss_kib", "peak_rss_mib", "status", "error", "testcase",
)
SUMMARY_FIELDS = (
    "genome", "reference_length", "m", "q", "k", "algorithm", "successful_trials",
    "construction_median_ms", "search_median_ms", "total_median_ms",
    "index_rss_increase_median_mib", "peak_rss_median_mib", "status", "testcase",
)


def info_for_testcase(path):
    match = FILE_RE.fullmatch(path.name)
    if not match:
        raise ValueError(f"Unexpected testcase filename: {path.name}")
    return {"genome": match["genome"], "m": int(match["m"]),
            "q": int(match["q"]), "k": int(match["k"])}


def reference_lengths(testcase_dir):
    """Read lengths from the generator's manifest, without loading huge inputs."""
    result = {}
    manifest = testcase_dir / "manifest.csv"
    if manifest.is_file():
        with manifest.open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                result[Path(row["input_file"]).name] = row.get("reference_length", "")
    return result


def unsupported_reason(algorithm, m, k):
    if algorithm in ("psi", "pss", "psv") and (k == 0 or k > 16 or 2 * k > m):
        return "Paired seeds require 1 <= k <= 16 and 2*k <= m"
    if algorithm in ("ssi", "sss", "ssv") and (k == 0 or k > 32 or k > m):
        return "Single seeds require 1 <= k <= min(32, m)"
    if algorithm == "psv" and k > 5:
        return "Paired direct-vector index supports k <= 5"
    if algorithm == "ssv" and k > 10:
        return "Single direct-vector index supports k <= 10"
    return ""


def extract_number(output, label):
    m = re.search(r"(?im)^\s*" + re.escape(label) + r":\s*([+-]?\d+(?:\.\d+)?)\s*ms\b", output)
    return float(m.group(1)) if m else None


def extract_index_mib(output):
    m = re.search(
        r"(?im)^\s*(?:RSS increase|Est\.? index RAM|Estimated index RAM|Index RAM|Index memory)"
        r":\s*([+-]?\d+(?:\.\d+)?)\s*(KiB|KB|MiB|MB|GiB|GB)?\b", output
    )
    if not m:
        return None
    value = float(m.group(1))
    unit = (m.group(2) or "MiB").lower()
    return value / 1024 if unit in ("kib", "kb") else value * 1024 if unit in ("gib", "gb") else value


def run_one(executable, testcase, algorithm, timeout):
    command = ["/usr/bin/time", "-f", f"{PEAK_TAG}=%M", str(executable), algorithm]
    with testcase.open("rb") as input_file:
        try:
            process = subprocess.run(
                command, stdin=input_file, stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE, text=True, errors="replace", timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return None, f"Timed out after {timeout} seconds"
        except OSError as exc:
            return None, str(exc)

    output = process.stderr
    if process.returncode != 0:
        return None, f"Exit code {process.returncode}: {output[-700:].strip()}"

    build = extract_number(output, "Construction")
    search = extract_number(output, "Searching")
    total = extract_number(output, "Total")
    index_mib = extract_index_mib(output)
    peak_match = re.search(rf"(?m)^{PEAK_TAG}=(\d+)\s*$", output)

    if build is None or search is None or total is None or index_mib is None or not peak_match:
        return None, f"Could not parse benchmark metrics: {output[-700:].strip()}"

    peak_kib = int(peak_match.group(1))
    return {
        "construction_ms": build,
        "search_ms": search,
        "total_ms": total,
        "index_rss_increase_mib": index_mib,
        "peak_rss_kib": peak_kib,
        "peak_rss_mib": round(peak_kib / 1024, 5),
    }, ""


def write_summary(raw_path, summary_path):
    """Summarize successful runs; report unsupported/failing configurations."""
    groups = defaultdict(list)
    with raw_path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            # FMI ignores k, so consolidate it across representative input files.
            # This also keeps --resume correct if a --limit smoke test selected
            # a different representative from the subsequent full run.
            k = "" if row["algorithm"] == "fmi" else row["k"]
            key = (row["genome"], row["m"], row["q"], k, row["algorithm"])
            groups[key].append(row)

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for key, rows in sorted(groups.items()):
            successes = [r for r in rows if r["status"] == "ok"]
            status = "ok" if successes else "skipped" if any(r["status"] == "skipped" for r in rows) else "error"
            record = dict(zip(("genome", "m", "q", "k", "algorithm"), key))
            record["reference_length"] = next(
                (r["reference_length"] for r in rows if r["reference_length"]), ""
            )
            record["testcase"] = rows[0]["testcase"]  # Representative input filename for FMI.
            record.update(status=status, successful_trials=len(successes))
            for src, dst in (
                ("construction_ms", "construction_median_ms"),
                ("search_ms", "search_median_ms"),
                ("total_ms", "total_median_ms"),
                ("index_rss_increase_mib", "index_rss_increase_median_mib"),
                ("peak_rss_mib", "peak_rss_median_mib"),
            ):
                record[dst] = round(statistics.median(float(r[src]) for r in successes), 6) if successes else ""
            writer.writerow(record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--testcases", type=Path, default=PROJECT_ROOT / "testcases")
    parser.add_argument("--exe", type=Path, default=PROJECT_ROOT / "build" / "main")
    parser.add_argument("--out", type=Path, default=PROJECT_ROOT / "results" / "benchmark_runs.csv")
    parser.add_argument("--summary", type=Path, default=PROJECT_ROOT / "results" / "benchmark_summary.csv")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--algorithms", nargs="+", choices=ALGORITHMS, default=list(ALGORITHMS))
    parser.add_argument("--limit", type=int, help="Only the first N testcases (smoke test)")
    parser.add_argument("--timeout", type=int, default=600, help="Seconds allowed per run")
    parser.add_argument("--resume", action="store_true", help="Append only missing runs to the raw CSV")
    args = parser.parse_args()

    if args.repeats < 1 or args.timeout < 1 or (args.limit is not None and args.limit < 1):
        parser.error("--repeats, --timeout, and --limit must be positive")
    if not args.exe.is_file():
        parser.error(f"Executable not found: {args.exe}; build the Release executable first")
    if not Path("/usr/bin/time").is_file():
        parser.error("/usr/bin/time not found; install using: sudo apt install time")
    cases = sorted(args.testcases.glob("*.in"))
    if args.limit is not None:
        cases = cases[:args.limit]
    if not cases:
        parser.error(f"No .in testcases found under {args.testcases}")

    # Each (genome, m, q) uses the same reference and query list across k.
    # Pick the smallest available k as its representative FMI input file.
    case_info = {path: info_for_testcase(path) for path in cases}
    fmi_representatives = {}
    for path, data in case_info.items():
        group = (data["genome"], data["m"], data["q"])
        previous = fmi_representatives.get(group)
        if previous is None or data["k"] < case_info[previous]["k"]:
            fmi_representatives[group] = path

    executable = args.exe.resolve()
    lengths = reference_lengths(args.testcases)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    completed = set()
    completed_fmi = set()
    if args.resume and args.out.is_file():
        with args.out.open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row["status"] in ("ok", "skipped"):
                    completed.add((row["testcase"], row["algorithm"], row["trial"]))
                    if row["algorithm"] == "fmi":
                        completed_fmi.add((row["genome"], int(row["m"]),
                                           int(row["q"]), row["trial"]))

    mode = "a" if args.resume and args.out.is_file() else "w"
    rng = random.Random(42)
    n_ok = n_skip = n_error = 0
    try:
        with args.out.open(mode, newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=RAW_FIELDS)
            if mode == "w":
                writer.writeheader()
            for case_index, testcase in enumerate(cases, 1):
                data = {**case_info[testcase],
                        "reference_length": lengths.get(testcase.name, "")}
                group = (data["genome"], data["m"], data["q"])
                print(f"[{case_index}/{len(cases)}] {testcase.name}", flush=True)

                for trial in range(1, args.repeats + 1):
                    order = args.algorithms[:]
                    rng.shuffle(order)
                    for algorithm in order:
                        # FMI depends on genome, m, and q, but never on k.
                        if algorithm == "fmi" and testcase != fmi_representatives[group]:
                            continue

                        skip = unsupported_reason(algorithm, data["m"], data["k"])
                        trial_id = "0" if skip else str(trial)
                        if algorithm == "fmi":
                            if (*group, trial_id) in completed_fmi:
                                continue
                        elif (testcase.name, algorithm, trial_id) in completed:
                            continue
                        # Unsupported configurations need just one explanatory row.
                        if skip and trial != 1:
                            continue

                        row = {**data, "algorithm": algorithm, "trial": trial_id,
                               "testcase": testcase.name, "status": "", "error": ""}
                        if algorithm == "fmi":
                            row["k"] = ""  # FMI's input k is ignored; do not imply otherwise.
                        if skip:
                            row.update(status="skipped", error=skip)
                            n_skip += 1
                        else:
                            metrics, error = run_one(executable, testcase, algorithm, args.timeout)
                            if metrics is None:
                                row.update(status="error", error=error)
                                n_error += 1
                                print(f"  ERROR {algorithm} trial {trial}: {error}", file=sys.stderr, flush=True)
                            else:
                                row.update(metrics)
                                row["status"] = "ok"
                                n_ok += 1
                        writer.writerow(row)
                        file.flush()  # Preserve results if a long experiment is interrupted.
    except KeyboardInterrupt:
        print("\nInterrupted; writing summary from completed results...", file=sys.stderr)
    finally:
        if args.out.is_file():
            write_summary(args.out, args.summary)
            print(f"\nRaw trials: {args.out}")
            print(f"Median summary: {args.summary}")
            print(f"New runs: {n_ok} successful, {n_skip} skipped, {n_error} errors")


if __name__ == "__main__":
    main()
