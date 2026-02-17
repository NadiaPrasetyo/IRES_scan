#!/usr/bin/env python3

import os
import subprocess
import argparse
from pathlib import Path
import sys
import logging

def setup_logging(verbose=False):
    log_file = "cluster_IRES_seq.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()] if verbose else [logging.StreamHandler()],
    )

def run(cmd):
    logging.info(">>: %s", " ".join(cmd))
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
        description="Cluster FASTA files recursively using MMseqs2: All vs All clustering of IRES sequences to find representative sequences and cluster assignments."
    )
    parser.add_argument(
        "--input_dir",
        nargs="?",
        help="Directory to search for FASTA files",
    )
    parser.add_argument("--input_fasta", help="Optional CSV fasta file with sequences to cluster (alternative to --input_dir)")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads to use for MMseqs2 (default: 16)")
    parser.add_argument("--outdir", default="data/mmseqs_cluster_out", help="Output directory for cluster results (default: data/mmseqs_cluster_out)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to console and file")
    args = parser.parse_args()

    setup_logging(args.verbose)

    input_dir = Path(args.input_dir).resolve() if args.input_dir else None
    input_fasta = Path(args.input_fasta).resolve() if args.input_fasta else None

    # check that there is at least one input:
    if not any([input_dir, input_fasta]):
        logging.error("❌ No input provided. Please specify either --input_dir or --input_fasta.")
        sys.exit(1)

    outdir = Path(args.outdir)

    db = "tmp/db"
    clu = "tmp/db_clu"
    rep_db = "tmp/db_rep"

    outdir.mkdir(exist_ok=True)
    os.makedirs("tmp", exist_ok=True)

    merged_fasta = None

    if input_dir:
        fasta_files = collect_fastas(input_dir)

        if not fasta_files:
            logging.error("❌ No FASTA files found.")
            sys.exit(1)

        logging.info(f"Found {len(fasta_files)} FASTA files")
        merged_fasta = outdir / "all_sequences.fasta"
        merge_fastas(fasta_files, merged_fasta)

    elif input_fasta:
        merged_fasta = input_fasta


    run(["mmseqs", "createdb", str(merged_fasta), str(db)])

    run([
        "mmseqs", "cluster",
        str(db),
        str(clu),
        "tmp",
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
        str(outdir / "cluster_representatives.fasta"),
    ])

    logging.info("\n🎉 Done!")
    logging.info(f"Clusters TSV: {outdir / 'clusters.tsv'}")
    logging.info(f"Representatives FASTA: {outdir / 'cluster_representatives.fasta'}")

if __name__ == "__main__":
    main()