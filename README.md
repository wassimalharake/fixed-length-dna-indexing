# Setup and Run Benchmarks

These instructions are for **Ubuntu / WSL Ubuntu**. Run commands from the project root.

## 1. Clone and install dependencies

```bash
git clone https://github.com/wassimalharake/fixed-length-dna-indexing.git
cd fixed-length-dna-indexing

sudo apt update
sudo apt install -y build-essential cmake git python3 time

# SeqAn3 is not committed to this repository.
git clone https://github.com/seqan/seqan3.git seqan3
```

The project uses C++20 and SeqAn3. For strictly reproducible builds, check out the same SeqAn3 commit used for the original experiments, if that revision is available.

## 2. Build the executable

```bash
bash compile.sh
```

This creates `build/main`. Use the **Release** build for timing and memory measurements.

## 3. Prepare the genome references

Generate the individual chromosome files:

```bash
python3 genomes_fna/gen_chromosomes.py
```

This writes `genomes/ecoli_cm1.txt` and `genomes/scerevisiae_cm*.txt`. Sequences are converted to uppercase A/C/G/T; other bases are removed by the current preprocessing script.

## 4. Generate benchmark inputs

```bash
# Main mixed-query workloads
python3 gen_testcases/gen_benchmarks.py

# Targeted exact, absent, and repeat-heavy workloads
python3 gen_testcases/gen_targeted_benchmarks.py
```

By default, the main generator selects `ecoli_cm1`, `scerevisiae_cm4`, `scerevisiae_cm14`, and `scerevisiae_cm17` and writes inputs to `testcases/`. The targeted generator writes to `targeted_testcases/`. Each directory also gets its own `manifest.csv`.

## 5. Run the benchmarks

```bash
# Main benchmark campaign (all seven algorithms; five trials per supported configuration)
python3 run_testcases/run_benchmarks.py

# Three targeted workloads (also five trials per supported configuration)
python3 run_testcases/run_targeted_benchmarks.py
```

For a quick smoke test before a full run:

```bash
python3 run_testcases/run_benchmarks.py --limit 2 --repeats 1
python3 run_testcases/run_targeted_benchmarks.py --limit 2 --repeats 1
```

For an interrupted run, use `--resume` with the relevant command. Without `--resume`, the corresponding existing results CSVs are overwritten.

## 6. Output files

All benchmark results are written to `results/`:

```text
results/
├── benchmark_runs.csv
├── benchmark_summary.csv
├── targeted_benchmark_runs.csv
└── targeted_benchmark_summary.csv
```

The `*_runs.csv` files contain individual trials; the `*_summary.csv` files contain median measurements. The runner skips configurations unsupported by a particular algorithm (such as direct-address indexes with oversized seeds).

### Optional: Run one algorithm on one input

```bash
./build/main fmi < testcases/ecoli_cm1__m150__q1000__k4.in
./build/main ssi < testcases/ecoli_cm1__m150__q1000__k4.in
./build/main psi < testcases/ecoli_cm1__m150__q1000__k4.in
```

Available algorithm identifiers: `fmi`, `ssi`, `psi`, `sss`, `pss`, `ssv`, `psv`.