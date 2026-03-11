#!/usr/bin/env python3
'''
get_mouse_homolog_seq.py

Usage:
    python get_mouse_homolog_seq.py --table table --fasta_file fasta_file --output_file output_file
    
Options:
    --table FILE Table file with IRES information (default: data/mouse_IRES/filtered_cmscan_table)
    --fasta_file FILE FASTA file with IRES sequences (default: data/mouse_IRES/mouse_IRES_sequences.fasta)
    --output_file FILE Output FASTA file for extracted IRES sequences (default: data/mouse_IRES/mouse_IRES_homologs.fasta)
    --verbose : bool    Enable verbose logging to console and file

Output:
    A FASTA file containing mouse homologs of the human IRES sequences

The script will use MMseqs2 to find mouse homologs of the human IRES sequences
and write them to a new FASTA file.

MMseqs2: https://github.com/soedinglab/MMseqs2 is required to run this script.

Author: Nadia Prasetyo
'''

from Bio import SeqIO
from Bio.Seq import Seq
import re
import argparse
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

    Returns
    -------
    None
    """
    log_file = "get_mouse_homolog_seq.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()] if verbose else [logging.StreamHandler()],
    )

def parse_fasta_headers(fasta_file):
    """
    Extract genomic coordinates from fasta headers.
    Expected format:
    >hsa_ires_00510.1 | chr18:37877999-39126684
    """
    coord_map = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        seq_id = record.id

        match = re.search(r'(chr[\w]+):(\d+)-(\d+)', header)
        if not match:
            raise ValueError(f"Could not parse genomic coordinates from header: {header}")

        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))

        coord_map[seq_id] = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "sequence": record.seq
        }

    return coord_map


def extract_ires_sequences(table_file, fasta_data):
    """
    Extracts IRES sequences from a table file and FASTA file.

    Table file expected format:
    target_name    accession    query_name    seq_from    seq_to    strand    description

    FASTA file expected format:
    >mouse_hsa_ires_00510.1 | chr18:37877999-39126684

    Returns a list of tuples, where each tuple contains:
    - A header string in the format: "mouse_<query_name> | <detected_location> | <accession> | <target_name> | <description>"
    - The extracted IRES sequence as a Seq object

    :param table_file: Table file with IRES information
    :param fasta_data: FASTA file with IRES sequences
    :return: List of tuples containing extracted IRES sequences and their corresponding headers
    """
    results = []

    with open(table_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()
            
            target_name = parts[0]
            accession = parts[1]
            query_name = parts[2]
            seq_from = int(parts[7])
            seq_to = int(parts[8])
            strand = parts[9]
            description = " ".join(parts[18:])

            if query_name not in fasta_data:
                logging.error(f"Warning: {query_name} not found in FASTA")
                continue

            fasta_entry = fasta_data[query_name]
            full_seq = fasta_entry["sequence"]
            chrom = fasta_entry["chrom"]
            region_start = fasta_entry["start"]
            region_end = fasta_entry["end"]

            # Convert to 0-based indexing for python slicing
            start_idx = min(seq_from, seq_to) - 1
            end_idx = max(seq_from, seq_to)

            extracted_seq = full_seq[start_idx:end_idx]

            # Handle strand
            if strand == "-":
                extracted_seq = extracted_seq.reverse_complement()

            # Calculate genomic location
            if strand == "+":
                genomic_start = region_start + start_idx
                genomic_end = region_start + end_idx - 1
            else:
                genomic_start = region_end - end_idx + 1
                genomic_end = region_end - start_idx

            detected_location = f"{chrom}:{genomic_start}-{genomic_end}({strand})"

            header = (
                f"mouse_{query_name} | {detected_location} | "
                f"{accession} | {target_name} | {description}"
            )

            results.append((header, extracted_seq))

    return results


def write_fasta(results, output_file):
    """
    Write a list of (header, sequence) tuples to a FASTA file.

    Parameters
    ----------
    results : list
        A list of tuples, where each tuple contains a FASTA header string and a Seq object.
    output_file : str
        The path to the output FASTA file.

    Returns
    -------
    None

    Notes
    -----
    This function writes a list of (header, sequence) tuples to a FASTA file.
    The format of the output file is a standard FASTA file, with each sequence
    preceded by a line starting with ">" and containing the header string.
    """
    with open(output_file, "w") as out:
        for header, seq in results:
            out.write(f">{header}\n")
            out.write(str(seq) + "\n")


def main():
    """
    Main function of the script.

    Extracts IRES sequences from cm_scan table and corresponding FASTA file, and writes them to a new FASTA file.

    Parameters
    ----------
    --table_file : str
        Path to the filtered cm_scan table (default: data/mouse_IRES/filtered_cmscan_table)
    --fasta_file : str
        Path to the FASTA file with original sequences (default: data/mouse_IRES/mouse_IRES_sequences.fasta)
    --output_file : str
        Path to the output FASTA file for extracted IRES sequences (default: data/mouse_IRES/mouse_IRES_homologs.fasta)
    --verbose : bool
        Enable verbose logging to console and file

    Outputs
    -------
    A FASTA file containing all the extracted IRES sequences.

    Notes
    -----
    This function will parse the FASTA file, extract the IRES sequences based on the cm_scan table,
    and write the extracted sequences to a new FASTA file.
    """
    parser = argparse.ArgumentParser(description="Extract IRES sequences from cm_scan table and corresponding FASTA file.")
    parser.add_argument("--table_file", default="data/mouse_IRES/filtered_cmscan_table", help="Path to the filtered cm_scan table (default: data/mouse_IRES/filtered_cmscan_table)")
    parser.add_argument("--fasta_file", default="data/mouse_IRES/mouse_IRES_sequences.fasta", help="Path to the FASTA file with original sequences (default: data/mouse_IRES/mouse_IRES_sequences.fasta)")
    parser.add_argument("--output_file", default="data/mouse_IRES/mouse_IRES_homologs.fasta", help="Path to the output FASTA file for extracted IRES sequences (default: data/mouse_IRES/mouse_IRES_homologs.fasta)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to console and file")
    args = parser.parse_args()
    setup_logging(verbose=args.verbose)

    fasta_data = parse_fasta_headers(args.fasta_file)
    logging.info(f"Parsed {len(fasta_data)} entries from FASTA file {args.fasta_file}")

    results = extract_ires_sequences(args.table_file, fasta_data)
    logging.info(f"Extracted {len(results)} IRES sequences from table {args.table_file}")

    write_fasta(results, args.output_file)
    logging.info(f"Done. Extracted {len(results)} IRES sequences to {args.output_file}")


if __name__ == "__main__":
    main()