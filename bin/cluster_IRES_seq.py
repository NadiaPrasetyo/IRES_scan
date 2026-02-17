#!/usr/bin/env python3

import subprocess
import argparse
from pathlib import Path
import sys

def run(cmd):
    print(">>", " ".join(cmd))
    subprocess.run(cmd, check=True)

def collect_fastas(input_dir):
    fastas = []
    for ext in ("*.fa", "*.fasta"):
        fastas.extend(input_dir.rglob(ext))
    return fastas

def merge_fastas(fasta_files, output_fasta):
    with open(output_fasta, "w") as out:
        for fasta in fasta_files:
            with open(fasta) as f:
                out.write(f.read())
                if not f.read().endswith("\n"):
                    out.write("\n")

def main():
    parser = argparse.ArgumentParser(
        description="Cluster FASTA files recursively using MMseqs2"
    )
    parser.add_argument(
        "input_dir",
        nargs="?",
        default=".",
        help="Directory to search for FASTA files (default: current directory)",
    )
    parser.add_argument("--min-seq-id", type=float, default=0.9)
    parser.add_argument("--coverage", type=float, default=0.8)
    parser.add_argument("--cov-mode", type=int, default=1)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--outdir", default="mmseqs_cluster_out")

    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir)
    tmpdir = outdir / "tmp"

    merged_fasta = outdir / "all_sequences.fasta"
    db = outdir / "db"
    clu = outdir / "db_clu"
    rep_db = outdir / "db_rep"

    outdir.mkdir(exist_ok=True)
    tmpdir.mkdir(exist_ok=True)

    fasta_files = collect_fastas(input_dir)

    if not fasta_files:
        print("❌ No FASTA files found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(fasta_files)} FASTA files")
    merge_fastas(fasta_files, merged_fasta)

    run(["mmseqs", "createdb", str(merged_fasta), str(db)])

    run([
        "mmseqs", "cluster",
        str(db),
        str(clu),
        str(tmpdir),
        "--min-seq-id", str(args.min_seq_id),
        "-c", str(args.coverage),
        "--cov-mode", str(args.cov_mode),
        "--threads", str(args.threads),
    ])

    run([
        "mmseqs", "createtsv",
        str(db),
        str(db),
        str(clu),
        str(outdir / "clusters.tsv"),
    ])

    run([
        "mmseqs", "createsubdb",
        str(clu),
        str(db),
        str(rep_db),
    ])

    run([
        "mmseqs", "convert2fasta",
        str(rep_db),
        str(outdir / "representatives.fasta"),
    ])

    print("\n🎉 Done!")
    print(f"Clusters TSV: {outdir / 'clusters.tsv'}")
    print(f"Representatives FASTA: {outdir / 'representatives.fasta'}")

if __name__ == "__main__":
    main()