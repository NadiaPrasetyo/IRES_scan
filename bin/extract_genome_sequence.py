#!/usr/bin/env python3
'''
extract_genome_sequence.py

Usage:
    python extract_genome_sequence.py -f genome_fasta -i input_tsv -l length -o output_fasta

Options:
    -f, --genome_fasta STR: Mouse genome FASTA file
    -i, --input_tsv STR: Input TSV file
    -l, --length INT: Number of nucleotides to extract from exon start (default: 1248530, which is the max distance to the nearest downstream CCD exon found in the human IRES analysis + buffer (500))
    -o, --output_fasta STR: Output FASTA file (default: mouse_upstream_seq.fasta)

Output:
    A FASTA file containing all sequences from the input file/directory

The script will compile all sequences in the input file/directory into a single
FASTA file.

Author: Nadia Prasetyo
'''
import argparse
from Bio import SeqIO
import logging

def setup_logging(verbose=False):
    '''
    Set up logging for the script.

    If verbose is True, logs are written to both the console and the log file.
    If verbose is False, logs are only written to the log file.

    Parameters:
    verbose (bool): If True, write logs to both the console and the log file.
    '''
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler("extract_genome_sequence.log")] if verbose else [logging.StreamHandler()],
    )


def parse_arguments():
    '''
    Parse command line arguments for the script.

    Parameters:
    -f, --genome_fasta : str
        Mouse genome FASTA file
    -i, --input_tsv : str
        Input TSV file
    -l, --length : int
        Number of nucleotides to extract from exon start (default: 1248530, which is the max distance to the nearest downstream CCD exon found in the human IRES analysis + buffer (500))
    -o, --output_fasta : str
        Output FASTA file (default: mouse_upstream_seq.fasta)
    --verbose : bool
        Enable verbose logging

    Returns:
    argparse.Namespace
        Parsed command line arguments
    '''
    parser = argparse.ArgumentParser(
        description="Extract exon sequences from mouse genome using TSV coordinates."
    )
    parser.add_argument("-f", "--genome_fasta", help="Mouse genome FASTA file")
    parser.add_argument("-i", "--input_tsv", help="Input TSV file")
    parser.add_argument("-l", "--length", type=int, default=1248530, help="Number of nucleotides to extract from exon start (default: 1248530, which is the max distance to the nearest downstream CCD exon found in the human IRES analysis + buffer (500))")
    parser.add_argument("-o", "--output_fasta", default="mouse_upstream_seq.fasta", help="Output FASTA file (default: mouse_upstream_seq.fasta)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def load_genome(genome_fasta):
    '''
    Load a mouse genome FASTA file and store its sequences in a dictionary.
    
    Parameters
    ----------
    genome_fasta : str
        Path to the mouse genome FASTA file
    
    Returns
    -------
    dict
        A dictionary where the keys are the sequence IDs and the values are the sequences themselves.
    '''
    genome = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        genome[record.id] = record.seq
    return genome


def main():
    '''
    Main function of the script.

    Extracts the upstream sequences from the mouse genome for a given set of IRES coordinates in a TSV file and writes them to a new FASTA file.

    Parameters
    ----------
    --genome_fasta : str
        Mouse genome FASTA file
    --input_tsv : str
        Input TSV file
    --length : int
        Number of nucleotides to extract from exon start (default: 1248530, which is the max distance to the nearest downstream CCD exon found in the human IRES analysis + buffer (500))
    --output_fasta : str
        Output FASTA file (default: mouse_upstream_seq.fasta)
    --verbose : bool
        Enable verbose logging

    Returns
    -------
    None
    '''
    args = parse_arguments()
    setup_logging(verbose=args.verbose)
    genome = load_genome(args.genome_fasta)

    

    with open(args.input_tsv) as tsv, open(args.output_fasta, "w") as out:
        header = tsv.readline()  # skip header

        for line in tsv:
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue

            query = fields[0]
            chrom = fields[1]
            tstart = int(fields[5])
            tend = int(fields[6])

            ires_id = query.replace("|Homo", "").replace("|homo", "")

            if chrom not in genome:
                print(f"Warning: {chrom} not found in genome.")
                continue

            chrom_seq = genome[chrom]

            # Compute upstream start
            upstream_start = tstart - args.length

            # Clamp to chromosome start if needed
            if upstream_start < 1:
                upstream_start = 1

            # Convert to 0-based indexing
            start_index = upstream_start - 1
            end_index = tend  # tend is already correct for slicing

            sequence = chrom_seq[start_index:end_index]

            header_location = f"{chrom}:{upstream_start}-{tend}"

            out.write(f">{ires_id} | {header_location}\n")
            out.write(str(sequence) + "\n")


if __name__ == "__main__":
    main()