#!/usr/bin/env python3
"""Generate three focused workloads for exact-DNA-search research.

Run from anywhere (default paths resolve from the project root):
    python3 gen_testcases/gen_targeted_benchmarks.py

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
from array import array
from pathlib import Path

DNA = "ACGT"
PROJECT_ROOT = Path(__file__).resolve().parent.parent
WORKLOADS = ("exact", "absent", "repeat_heavy")


def longest_acgt_segment(raw: str) -> tuple[str, int]:
    """Keep the longest ACGT segment without joining across ambiguous bases."""
    match = max(
        re.finditer(r"[ACGT]+", raw.upper()),
        key=lambda item: item.end() - item.start(),
        default=None,
    )
    return (match.group(), match.start()) if match else ("", 0)

def stable_rng(seed: int, genome: str, workload: str, m: int):
    digest = hashlib.sha256(f"{seed}/{genome}/{workload}/{m}".encode()).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def mutate(query: str, position: int, rng: random.Random) -> str:
    """Replace one DNA base with a different randomly selected base."""
    old = query[position]
    alternatives = [c for c in DNA if c != old]
    replacement = rng.choice(alternatives)
    return query[:position] + replacement + query[position + 1:]

def generate_exact(reference: str, m: int, count: int, rng: random.Random):
    last = len(reference) - m
    return [reference[p:p + m] for p in (rng.randrange(last + 1) for _ in range(count))]


def generate_absent(reference: str, m: int, count: int, max_k: int,
                    rng: random.Random):
    """Confirmed absent one-substitution reads, 50% middle / 50% endpoint.

    Middle mutations preserve both endpoint seeds for all tested k, whereas
    endpoint mutations change one of them. All queries are checked against
    the *entire* reference, not just the originally sampled position.
    """
    if m <= 2 * max_k:
        raise ValueError("m must be > 2 * max(k) for middle-mutation workload")
    last = len(reference) - m
    queries = []
    for j in range(count):
        kind_middle = (j % 2 == 0)
        for attempt in range(200):
            pos = rng.randrange(last + 1)
            query = reference[pos:pos + m]
            if kind_middle:
                mutation_pos = rng.randrange(max_k, m - max_k)
            else:
                mutation_pos = rng.randrange(max_k) if rng.randrange(2) == 0 \
                    else rng.randrange(m - max_k, m)
            candidate = mutate(query, mutation_pos, rng)
            if candidate not in reference:  # confirms ZERO exact matches
                queries.append(candidate)
                break
        else:
            raise RuntimeError(f"Could not find absent m={m} query after 200 attempts ")
    # Shuffle the two categories, so Q prefixes have a representative mix.
    rng.shuffle(queries)
    # Every accepted candidate was checked against the full reference above.
    return queries


def high_collision_positions(reference: str, m: int, anchor_k: int,
                             top_keys: int):
    """Find positions with frequent *paired* short-seed keys.

    Uses only a (4 * anchor_k)-bit count table, plus one byte per k-mer
    when anchor_k <= 4. This is O(n + 16**anchor_k) memory/time and avoids
    building Python dictionaries containing millions of genome substrings.
    Not all selected reads have multiple full-length exact matches: this
    tests seed-candidate pressure, which is the intended stress scenario.
    """
    if not (1 <= anchor_k <= 4):
        raise ValueError("anchor_k must be between 1 and 4 for compact byte lookup")
    if m < 2 * anchor_k or len(reference) < m:
        raise ValueError("Invalid m for paired-anchor selection")

    # Translate ASCII bases to small 2-bit numbers in one C-level pass.
    encoded = reference.encode("ascii").translate(
        bytes.maketrans(b"ACGT", bytes((0, 1, 2, 3))))
    mask = (1 << (2 * anchor_k)) - 1
    codes = bytearray(len(reference) - anchor_k + 1)
    rolling = 0
    for j, base in enumerate(encoded):
        rolling = ((rolling << 2) | base) & mask
        if j >= anchor_k - 1:
            codes[j - anchor_k + 1] = rolling

    key_count = 1 << (4 * anchor_k)
    frequencies = array("I", [0]) * key_count
    last = len(reference) - m
    second_offset = m - anchor_k
    for i in range(last + 1):
        key = (codes[i] << (2 * anchor_k)) | codes[i + second_offset]
        frequencies[key] += 1

    top = sorted(range(key_count), key=lambda x: frequencies[x], reverse=True)
    top = [key for key in top if frequencies[key] >= 2][:top_keys]
    if not top:
        raise ValueError("No repeated paired-anchor keys found")
    selected = set(top)
    positions = [i for i in range(last + 1)
                 if ((codes[i] << (2 * anchor_k)) | codes[i + second_offset]) in selected]
    if not positions:
        raise RuntimeError("No repeat-enriched candidate positions found")
    stats = (len(positions), min(int(frequencies[k]) for k in top),
             max(int(frequencies[k]) for k in top))
    return positions, stats


def generate_repeat_heavy(reference: str, m: int, count: int,
                          hotspots: list[int], rng: random.Random):
    # Prefer distinct genomic starting positions; repeat positions only when
    # there are fewer hotspots than requested queries.
    positions = rng.sample(hotspots, min(count, len(hotspots)))
    positions.extend(rng.choice(hotspots) for _ in range(count - len(positions)))
    rng.shuffle(positions)
    return [reference[p:p + m] for p in positions]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-dir", type=Path, default=PROJECT_ROOT / "genomes")
    parser.add_argument("--out", type=Path, default=PROJECT_ROOT / "targeted_testcases")
    parser.add_argument("--genomes", nargs="+", default=["ecoli_cm1", "scerevisiae_cm4"],
                        help="Input genome file stems (default: E. coli + yeast chromosome IV)")
    parser.add_argument("--m", type=int, default=150, help="One fixed read length, default 150")
    parser.add_argument("--counts", nargs="+", type=int, default=[1000, 10000])
    parser.add_argument("--ks", nargs="+", type=int, default=[4, 5, 8, 12, 16])
    parser.add_argument("--anchor-k", type=int, default=4,
                        help="Repeat-heavy workload ranks fixed-gap paired keys of this length")
    parser.add_argument("--top-keys", type=int, default=128,
                        help="Select positions from the most frequent paired anchor keys")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--force", action="store_true",
                        help="Remove previously generated .in files in the targeted directory")
    args = parser.parse_args()

    if not args.counts or any(q <= 0 for q in args.counts):
        parser.error("All query counts must be positive")
    if args.m <= 0 or not args.ks or any(k < 1 or k > 16 or 2 * k >= args.m for k in args.ks):
        parser.error("Require 1 <= k <= 16 and 2*k < m")
    if not (1 <= args.anchor_k <= 4) or args.top_keys <= 0:
        parser.error("Require 1 <= anchor-k <= 4, and top-keys > 0")
    if args.m < 2 * args.anchor_k:
        parser.error("m is too small for the repeat-heavy paired anchor")
    if len(set(args.genomes)) != len(args.genomes):
        parser.error("Duplicate genome stems")
    if args.out.resolve() == (PROJECT_ROOT / "testcases").resolve():
        parser.error("Targeted output must not overwrite the main testcases directory")

    paths = [args.genome_dir / f"{stem}.txt" for stem in args.genomes]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        parser.error(f"Genome files not found: {', '.join(missing)}")

    args.out.mkdir(parents=True, exist_ok=True)
    old_cases = list(args.out.glob("*.in"))
    if old_cases and not args.force:
        parser.error(f"{args.out} already contains {len(old_cases)} .in files; "
                     "use --force to regenerate them")
    if args.force:
        for old in old_cases:
            old.unlink()

    counts = sorted(set(args.counts))
    ks = sorted(set(args.ks))
    generated = 0
    with (args.out / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["genome", "m", "Q", "k", "reference_length", "trim_start",
                         "input_file", "source_genome", "workload", "seed",
                         "repeat_anchor_k", "repeat_top_keys"])
        for path in paths:
            raw = "".join(path.read_text(encoding="ascii").split()).upper()
            if ">" in raw:
                raise ValueError(f"{path} still contains a FASTA header")
            reference, trim_start = longest_acgt_segment(raw)
            if len(reference) < args.m:
                raise ValueError(f"{path}: no contiguous ACGT region of length m={args.m}")
            if len(reference) != len(raw):
                print(f"NOTE {path.stem}: using longest ACGT segment, "
                      f"length={len(reference):,}, offset={trim_start}", flush=True)

            hotspots, stats = high_collision_positions(
                reference, args.m, args.anchor_k, args.top_keys)
            print(f"{path.stem}: reference length={len(reference):,}; "
                  f"repeat-heavy pool={stats[0]:,} positions, "
                  f"paired-anchor key counts={stats[1]}..{stats[2]}", flush=True)

            for workload in WORKLOADS:
                rng = stable_rng(args.seed, path.stem, workload, args.m)
                if workload == "exact":
                    queries = generate_exact(reference, args.m, max(counts), rng)
                elif workload == "absent":
                    queries = generate_absent(reference, args.m, max(counts), max(ks), rng)
                else:
                    queries = generate_repeat_heavy(reference, args.m, max(counts), hotspots, rng)

                label = f"{path.stem}__{workload}"  # prevents FMI dedup across workloads
                for q in counts:
                    for k in ks:
                        filename = f"{label}__m{args.m}__q{q}__k{k}.in"
                        file_path = args.out / filename
                        with file_path.open("w", encoding="ascii") as output:
                            output.write(f"{args.m} {q} {k}\n{reference}\n")
                            output.write("\n".join(queries[:q]))
                            output.write("\n")
                        writer.writerow([label, args.m, q, k, len(reference), trim_start,
                                         (str(file_path.relative_to(PROJECT_ROOT))
                                          if file_path.is_relative_to(PROJECT_ROOT) else str(file_path)),
                                         path.stem, workload, args.seed,
                                         args.anchor_k, args.top_keys])
                        generated += 1
                print(f"  {workload}: generated {len(counts) * len(ks)} files", flush=True)

    print(f"\nGenerated {generated} targeted input files under {args.out}/", flush=True)
    print(f"Manifest: {args.out / 'manifest.csv'}", flush=True)


if __name__ == "__main__":
    main()
