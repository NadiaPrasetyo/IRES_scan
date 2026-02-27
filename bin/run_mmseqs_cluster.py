#!/usr/bin/env python3

import os
import shutil
import subprocess
import argparse
import sys
import logging

def setup_logging(verbose=False):
    log_file = "run_mmseqs_cluster.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()] if verbose else [logging.StreamHandler()],
    )


def run(cmd):
    logging.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit("Error while running: " + " ".join(cmd))


def main():
    parser = argparse.ArgumentParser(
        description="Simple MMseqs2 clustering pipeline"
    )

    parser.add_argument("-i", "--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--min-seq-id", default="0.9", help="Minimum sequence identity (default: 0.9)")
    parser.add_argument("-c", "--coverage", default="0.8", help="Coverage threshold (default: 0.8)")
    parser.add_argument("--cov-mode", default="0", help="Coverage mode (default: 0)")
    parser.add_argument("--threads", default="4", help="Number of threads (default: 4)")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    tmp_dir = os.path.join(args.output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    db = os.path.join(args.output_dir, "DB")
    clu_db = os.path.join(args.output_dir, "DB_clu")

    cluster_tsv = os.path.join(args.output_dir, "clusters.tsv")
    cluster_seq_db = os.path.join(args.output_dir, "DB_clu_seq")
    cluster_fasta = os.path.join(args.output_dir, "clustered_sequences.fasta")
    rep_db = os.path.join(args.output_dir, "DB_clu_rep")
    rep_fasta = os.path.join(args.output_dir, "cluster_representatives.fasta")

    # check that mmseqs is installed
    try:
        subprocess.run(["mmseqs", "--help"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        logging.error("Error: MMseqs2 is not installed or not in PATH. Please install MMseqs2 to use this script.")

    # 1️⃣ Create MMseqs database
    run([
        "mmseqs", "createdb",
        args.input_fasta,
        db
    ])

    logging.info(f"Created MMseqs database: {db}")

    # 2️⃣ Run clustering
    run([
        "mmseqs", "cluster",
        db,
        clu_db,
        tmp_dir,
        "--min-seq-id", args.min_seq_id,
        "-c", args.coverage,
        "--cov-mode", args.cov_mode,
        "--threads", args.threads
    ])

    logging.info(f"Clustering complete. Cluster database: {clu_db}")

    # 3️⃣ Export cluster TSV
    run([
        "mmseqs", "createtsv",
        db,
        db,
        clu_db,
        cluster_tsv
    ])

    # 4️⃣ Export clustered sequences FASTA
    run([
        "mmseqs", "createseqfiledb",
        db,
        clu_db,
        cluster_seq_db
    ])

    run([
        "mmseqs", "result2flat",
        db,
        db,
        cluster_seq_db,
        cluster_fasta
    ])

    # 5️⃣ Extract representative sequences
    run([
        "mmseqs", "createsubdb",
        clu_db,
        db,
        rep_db
    ])

    run([
        "mmseqs", "convert2fasta",
        rep_db,
        rep_fasta
    ])

    print("\n✅ Clustering complete!\n")
    print("Cluster TSV:", cluster_tsv)
    print("Clustered FASTA:", cluster_fasta)
    print("Representative FASTA:", rep_fasta)

    logging.info("Deleting intermediate files: Temporary directory: %s, MMseqs databases: %s, %s, %s, %s", tmp_dir, db, clu_db, cluster_seq_db, rep_db)

    # Cleanup temporary files
    shutil.rmtree(tmp_dir)  # Uncomment to remove temporary directory
    # cleanup intermediate MMseqs databases if desired
    os.remove(db)
    os.remove(clu_db)
    os.remove(cluster_seq_db)
    os.remove(rep_db)


if __name__ == "__main__":
    main()