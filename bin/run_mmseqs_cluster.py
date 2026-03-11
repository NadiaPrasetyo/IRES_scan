#!/usr/bin/env python3
'''
run_mmseqs_cluster.py

Usage:
    python run_mmseqs_cluster.py -i input.fasta -o output_dir -t threads -m mmseqs_path -l max_seq_len -s min_seq_len -d db_path -c min_cov -e min_id -p identity -r min_len -v verbose

Options:
    -i, --input-fasta FILE    Input FASTA file [required]
    -o, --output-dir DIRECTORY    Output directory [required]
    --min-seq-id FLOAT    Minimum sequence identity [default: 0.9]
    -c, --coverage FLOAT    Coverage threshold [default: 0.8]
    --cov-mode INT    Coverage mode [default: 0]
    --threads INT    Number of threads [default: 4]
    --verbose    Enable verbose logging

Output:
    A FASTA file containing all sequences from the input file/directory

The script will compile all sequences in the input file/directory into a single
FASTA file.

Author: Nadia Prasetyo
'''
import os
import shutil
import subprocess
import argparse
import sys
import logging

def setup_logging(verbose=False):
    """
    Set up logging for the script.

    If verbose is True, logs are written to both the console and the log file.
    If verbose is False, logs are only written to the log file.

    Parameters
    ----------
    verbose : bool
        If True, write logs to both the console and the log file.
        If False, write logs only to the log file.
    """
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
    """
    Simple MMseqs2 clustering pipeline

    Parameters
    ----------
    -i, --input-fasta : Input FASTA file
    -o, --output-dir : Output directory
    --min-seq-id : Minimum sequence identity (default: 0.9)
    -c, --coverage : Coverage threshold (default: 0.8)
    --cov-mode : Coverage mode (default: 0)
    --threads : Number of threads (default: 4)
    --verbose : Enable verbose logging to console and file

    Notes
    ----
    * MMseqs2 is required to be installed and in the PATH
    * The script will create a temporary directory in the output directory
    * The script will create the following files:
        * clusters.tsv: Cluster TSV file
        * clustered_sequences.fasta: Clustered sequences FASTA file
        * cluster_representatives.fasta: Representative sequences FASTA file
        * cluster_alignments.tsv: Alignment results in BLAST tab format
    * The script will remove temporary files and directories after completion
    """
    parser = argparse.ArgumentParser(
        description="Simple MMseqs2 clustering pipeline"
    )

    parser.add_argument("-i", "--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--min-seq-id", default="0.9", help="Minimum sequence identity (default: 0.9)")
    parser.add_argument("-c", "--coverage", default="0.8", help="Coverage threshold (default: 0.8)")
    parser.add_argument("--cov-mode", default="0", help="Coverage mode (default: 0)")
    parser.add_argument("--threads", default="4", help="Number of threads (default: 4)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to console and file")

    args = parser.parse_args()

    setup_logging(args.verbose)

    os.makedirs(args.output_dir, exist_ok=True)
    tmp_dir = os.path.join(args.output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    db = os.path.join(tmp_dir, "DB")
    clu_db = os.path.join(tmp_dir, "DB_clu")

    cluster_tsv = os.path.join(args.output_dir, "clusters.tsv")
    cluster_seq_db = os.path.join(tmp_dir, "DB_clu_seq")
    cluster_fasta = os.path.join(args.output_dir, "clustered_sequences.fasta")
    rep_db = os.path.join(tmp_dir, "DB_clu_rep")
    rep_fasta = os.path.join(args.output_dir, "cluster_representatives.fasta")
    align_path = os.path.join(tmp_dir, "aln")
    cluster_aln_tsv = os.path.join(args.output_dir, "cluster_alignments.tsv")

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

    logging.info(f"🧬 Computing alignment + k-mer similarity scores")
    # Compute alignments from clustering result
    run([
        "mmseqs", "align", str(db), str(db), str(clu_db), str(align_path), "-a"
    ])
    
    # Convert alignment results to BLAST tab format
    run([
    "mmseqs", "convertalis", str(db), str(db), str(align_path), f"{cluster_aln_tsv}",
          "--format-mode", "4",
          "--search-type", "3",
          "--format-output", "query,target,pident,alnlen,bits,evalue,tseq,theader"
    ])


    logging.info("Deleting intermediate files: Temporary files and directory")

    # Cleanup temporary files
    shutil.rmtree(tmp_dir)
    logging.info("Removed temporary directory")

    # remove any remaining intermediate files in the output directory
    for filename in os.listdir(args.output_dir):
        if filename.startswith("DB") or filename.startswith("DB_clu") or filename.startswith("DB_clu_seq") or filename.startswith("aln"):
            file_path = os.path.join(args.output_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
                logging.info(f"Removed intermediate file: {file_path}")
        
    logging.info("Cleanup complete")


if __name__ == "__main__":
    main()