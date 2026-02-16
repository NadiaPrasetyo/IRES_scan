
# make sure mmseqs2 is installed and in your PATH
import subprocess
import logging
import argparse

def setup_logging(verbose=False):
    # log into a file when needed (verbose mode)
    log_file = "find_mouse_homolog.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()] if verbose else [logging.FileHandler(log_file)],
    )

def run(cmd):
    logging.info(">>", " ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description="Find mouse homologs of IRES sequences using MMseqs2")
    parser.add_argument("-i", "--input", default="data/human_IRES.fasta", help="Input FASTA file with human IRES sequences (default: data/human_IRES.*)")
    parser.add_argument("--mouse_genome", default="data/mm39.fa", help="Mouse genome file location (default: data/mm39.fa)")
    parser.add_argument("-o", "--output", default="data/mouse_IRES_homologs", help="Output file for homologous sequences (default: data/mouse_IRES_homologs.fasta)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to console")
    args = parser.parse_args()
    setup_logging(verbose=args.verbose)

    # Step 1: Download mouse IRES sequences from NCBI if it does not exist anywhere in the current directory
    if not mouse_genome:
        logging.info("Mouse genome not provided. Downloading from UCSC...")
        run (["wget", "-O", "data/mouse_genome.fa.gz", "https://hgdownload.gi.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"])
        # unzip the downloaded file
        run(["gunzip", "data/mouse_genome.fa.gz"])

    # Step 2: Check that mmseqs2 is installed and in the PATH
    try:
        subprocess.run(["mmseqs", "version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except Exception:
        raise FileNotFoundError("MMseqs2 not found in PATH. Install via: conda install -c bioconda mmseqs2")

    # Step 3: Create MMseqs2 databases for the input human IRES sequences and the mouse genome
    run(["mmseqs", "createdb", args.input, "human_IRES_db"])
    run(["mmseqs", "createdb", "data/mouse_genome.fa", "mouse_genome_db"])

    # create index for target mouse genome database to speed up search
    mkdir("tmp", exist_ok=True)
    run(["mmseqs", "createindex", "mouse_genome_db", "tmp"])

    # Step 4: Search for homologous sequences in the mouse genome using MMseqs2
    run(["mmseqs", "search", "human_IRES_db", "mouse_genome_db", "search_results_db", "tmp"])

    # Step 5: Convert the search results to FASTA format and m8 format
    run(["mmseqs", "convertalis", "human_IRES_db", "mouse_genome_db", "search_results_db", "{args.output}.m8"])
    run(["mmseqs", "result2flat", "human_IRES_db", "mouse_genome_db", "search_results_db", "{args.output}.fasta"])


if __name__ == "__main__":
    main()