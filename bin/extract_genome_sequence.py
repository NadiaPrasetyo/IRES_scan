#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import logging

def setup_logging(verbose=False):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler("extract_genome_sequence.log")] if verbose else [logging.StreamHandler()],
    )


def parse_arguments():
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
    genome = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        genome[record.id] = record.seq
    return genome


def main():
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