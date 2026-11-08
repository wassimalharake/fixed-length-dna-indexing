#!/usr/bin/env python3
"""Generate reproducible exact-DNA-search benchmark inputs.

Run from anywhere (default paths resolve from the project root):
    python3 gen_testcases/gen_benchmarks.py

The existing C++ program expects:
    m q k
    reference
    query_1
    ...
"""

import argparse
import csv
import hashlib
import random
import re
from collections import Counter
from pathlib import Path

DNA = "ACGT"
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_GENOMES = [
    "ecoli_cm1",
    "scerevisiae_cm4",
    "scerevisiae_cm14",
    "scerevisiae_cm17",
]


def longest_acgt_segment(raw: str) -> tuple[str, int]:
    """Keep the longest ACGT segment without joining across ambiguous bases."""
    match = max(
        re.finditer(r"[ACGT]+", raw.upper()),
        key=lambda item: item.end() - item.start(),
        default=None,
    )
    return (match.group(), match.start()) if match else ("", 0)

def mutate(query: str, position: int, rng: random.Random) -> str:
    """Replace one DNA base with a different randomly selected base."""
    old = query[position]
    alternatives = [c for c in DNA if c != old]
    replacement = rng.choice(alternatives)
    return query[:position] + replacement + query[position + 1:]

def generate_queries(reference: str, m: int, count: int, max_k: int,
                     seed: int, common_prefix_frequency: Counter, anchor_size: int) -> list[str]:
    """Same query bank for every k, with smaller Q using prefixes of the bank."""
    rng = random.Random(seed)
    last = len(reference) - m
    if last < 0:
        return []

    # Sample regions enriched for a common 8-mer prefix. This does NOT
    # guarantee the full m-mer occurs at multiple positions.
    probes = [rng.randint(0, last) for _ in range(min(20000, last + 1))]
    common_positions = sorted(
        probes,
        key=lambda pos: common_prefix_frequency[reference[pos:pos + anchor_size]],
        reverse=True,
    )[:max(1, len(probes) // 10)]

    # 60% real exact substrings, 20% interior mutations, 10% seed mutations,
    # 10% regions with relatively common short prefixes.
    types = ["exact"] * 6 + ["middle"] * 2 + ["seed"] + ["common"]
    result = []
    while len(result) < count:
        block = types[:]
        rng.shuffle(block)
        for kind in block:
            if len(result) == count:
                break
            pos = rng.choice(common_positions) if kind == "common" else rng.randint(0, last)
            query = reference[pos:pos + m]
            if kind == "middle":
                # Mutation stays outside both seeds for EVERY tested k.
                query = mutate(query, rng.randrange(max_k, m - max_k), rng)
            elif kind == "seed":
                query = mutate(query, rng.choice([0, m - 1]), rng)
            result.append(query)
    return result


def chromosome_sort(path: Path):
    match = re.search(r"cm(\d+)$", path.stem)
    return (0 if path.stem.startswith("ecoli") else 1,
            int(match.group(1)) if match else 999, path.stem)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-dir", type=Path, default=PROJECT_ROOT / "genomes")
    parser.add_argument("--out", type=Path, default=PROJECT_ROOT / "testcases_test")
    parser.add_argument("--genomes", nargs="+", default=DEFAULT_GENOMES,
                        help="Filename stems; default: E. coli and yeast chromosomes 4, 14, 17")
    parser.add_argument("--lengths", nargs="+", type=int, default=[50, 100, 150, 250])
    parser.add_argument("--counts", nargs="+", type=int, default=[100, 1000, 10000])
    parser.add_argument("--ks", nargs="+", type=int, default=[4, 5, 6, 8, 10, 12, 16])
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if not args.counts or any(q < 1 for q in args.counts):
        parser.error("All query counts must be positive")
    if not args.lengths or any(m < 1 for m in args.lengths):
        parser.error("All read lengths must be positive")
    if not args.ks or any(k < 1 or k > 16 for k in args.ks):
        parser.error("For packed paired seeds, k must be between 1 and 16")

    paths = sorted(args.genome_dir.glob("*.txt"), key=chromosome_sort)
    if args.genomes:
        wanted = set(args.genomes)
        paths = [p for p in paths if p.stem in wanted or p.name in wanted]
        found = {p.stem for p in paths} | {p.name for p in paths}
        if wanted - found:
            parser.error(f"Unknown genomes: {', '.join(sorted(wanted - found))}")
    if not paths:
        parser.error(f"No genome files found under {args.genome_dir}")

    args.out.mkdir(parents=True, exist_ok=True)
    total = 0
    manifest_path = args.out / "manifest.csv"
    with manifest_path.open("w", newline="") as manifest_file:
        writer = csv.writer(manifest_file)
        writer.writerow(["genome", "m", "Q", "k", "reference_length", "trim_start", "input_file"])

        for path in paths:
            raw = "".join(path.read_text(encoding="ascii").split()).upper()
            if ">" in raw:
                raise ValueError(f"{path} still contains FASTA headers; use cleaned chromosome files")
            reference, start = longest_acgt_segment(raw)
            if not reference:
                print(f"SKIP {path.name}: no A/C/G/T segment")
                continue
            if len(reference) != len(raw):
                print(f"NOTE {path.name}: ambiguous bases present; using longest contiguous ACGT region "
                      f"(offset {start}, length {len(reference)} of {len(raw)}).")

            # Precompute counts once per chromosome, not for every (m, Q, k).
            anchor = min(8, min(args.lengths))
            frequency = Counter(reference[i:i + anchor]
                                for i in range(len(reference) - anchor + 1))

            for m in sorted(set(args.lengths)):
                valid_ks = [k for k in sorted(set(args.ks)) if 2 * k <= m]
                if not valid_ks or m > len(reference):
                    print(f"SKIP {path.stem}, m={m}: reference too short or no valid k")
                    continue
                max_k = max(valid_ks)
                if m <= 2 * max_k:
                    print(f"SKIP {path.stem}, m={m}: no interior base for middle-mutation workload "
                          f"with k={max_k}")
                    continue

                stable_seed = int.from_bytes(hashlib.sha256(
                    f"{args.seed}/{path.stem}/{m}".encode()).digest()[:8], "big")
                queries = generate_queries(reference, m, max(args.counts), max_k,
                                           stable_seed, frequency, anchor)
                for q in sorted(set(args.counts)):
                    for k in valid_ks:
                        filename = f"{path.stem}__m{m}__q{q}__k{k}.in"
                        outfile = args.out / filename
                        with outfile.open("w", encoding="ascii") as f:
                            f.write(f"{m} {q} {k}\n{reference}\n")
                            f.write("\n".join(queries[:q]))
                            f.write("\n")
                        writer.writerow([path.stem, m, q, k, len(reference), start,
                                         str(outfile.relative_to(PROJECT_ROOT))
                                         if outfile.is_relative_to(PROJECT_ROOT) else str(outfile)])
                        total += 1
            print(f"Finished {path.stem} ({len(reference):,} bases)")

    print(f"Generated {total} input files in {args.out}/")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
