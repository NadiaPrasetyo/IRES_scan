
# make sure mmseqs2 is installed and in your PATH
import subprocess
import logging
import os
import argparse

def setup_logging(verbose=False):
    # log into a file when needed (verbose mode)
    log_file = "find_mouse_homolog.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()] if verbose else [logging.StreamHandler()],
    )

def run(cmd):
    logging.info(">> %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description="Find mouse homologs of IRES sequences using MMseqs2")
    parser.add_argument("-i", "--input", default="data/human_IRES.fasta", help="Input FASTA file with human IRES sequences (default: data/human_IRES.*)")
    parser.add_argument("--mouse_genome", default="data/mm39.fa", help="Mouse genome file location (default: data/mm39.fa)")
    parser.add_argument("-o", "--output", default="data/mouse_IRES_homologs", help="Output file for homologous sequences (default: data/mouse_IRES_homologs.fasta)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to console")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads to use for MMseqs2 (default: 16)")
    args = parser.parse_args()
    setup_logging(verbose=args.verbose)

    # Step 1: Download mouse IRES sequences from NCBI if it does not exist anywhere in the current directory
    if not args.mouse_genome or not os.path.exists(args.mouse_genome):
        logging.info("Mouse genome not provided or does not exist. Downloading from UCSC...")
        run (["wget", "-O", "data/mm39.fa.gz", "https://hgdownload.gi.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"])
        # unzip the downloaded file
        run(["gunzip", "data/mm39.fa.gz"])

    # Step 2: Check that mmseqs2 is installed and in the PATH
    try:
        subprocess.run(["mmseqs", "version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except Exception:
        raise FileNotFoundError("MMseqs2 not found in PATH. Install via: conda install -c bioconda mmseqs2")

    # put all db and other intermediate files in a tmp directory
    os.makedirs("tmp", exist_ok=True)
    # Step 3: Create MMseqs2 databases for the input human IRES sequences and the mouse genome
    run(["mmseqs", "createdb", args.input, "tmp/human_IRES_db"])
    run(["mmseqs", "createdb", args.mouse_genome, "tmp/mouse_genome_db"])

    # Step 4: Search for homologous sequences in the mouse genome using MMseqs2
    run([
    "mmseqs", "search",
    "tmp/human_IRES_db",
    "tmp/mouse_genome_db",
    "tmp/search_results_db",
    "tmp",
    "--search-type", "3",
    "-s", "7.5",
    "--threads", str(args.threads),
    ])

    # Step 5: Convert the search results to FASTA format and m8 format
    run(["mmseqs", "convertalis", "tmp/human_IRES_db", "tmp/mouse_genome_db", "tmp/search_results_db",
          f"{args.output}.tsv",
          "--search-type", "3",
          "--format-mode", "4",
          "--format-output", "query,target,pident,alnlen,bits,qheader,evalue,qaln,taln"])

    # Step 6: Clean up temporary files
    run(["rm", "-rf", "tmp"])

if __name__ == "__main__":
    main()