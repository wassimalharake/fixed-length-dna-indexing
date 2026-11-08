#!/usr/bin/env python3
"""Split each FASTA record in genomes_fna/*.fna into genomes/<name>_cmN.txt.

Run from any directory:
    python3 genomes_fna/gen_chromosomes.py

The output contains uppercase sequence only, on one line per file.
Ambiguous bases (such as N) are preserved to avoid changing coordinates.
"""

from pathlib import Path


FNA_DIR = Path(__file__).resolve().parent
GENOMES_DIR = FNA_DIR.parent / "genomes"


def split_fna(path: Path, output_dir: Path) -> tuple[int, int]:
    chromosome = 0
    total_bases = 0
    output = None

    try:
        with path.open("r", encoding="ascii") as input_file:
            for line_number, line in enumerate(input_file, start=1):
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if output is not None:
                        output.write("\n")
                        output.close()

                    chromosome += 1
                    filename = f"{path.stem}_cm{chromosome}.txt"
                    output = (output_dir / filename).open("w", encoding="ascii")
                    continue

                if output is None:
                    raise ValueError(f"{path.name}:{line_number}: sequence before FASTA header")

                sequence = "".join(c for c in line.upper() if c in "ACGT")
                output.write(sequence)
                total_bases += len(sequence)
    finally:
        if output is not None and not output.closed:
            output.write("\n")
            output.close()

    return chromosome, total_bases


def main() -> None:
    paths = sorted(FNA_DIR.glob("*.fna"))
    if not paths:
        raise FileNotFoundError(f"No .fna files found in {FNA_DIR}")

    GENOMES_DIR.mkdir(parents=True, exist_ok=True)

    for path in paths:
        chromosome_count, base_count = split_fna(path, GENOMES_DIR)
        print(f"{path.name}: {chromosome_count} sequences, {base_count:,} bases")

    print(f"Output directory: {GENOMES_DIR}")


if __name__ == "__main__":
    main()