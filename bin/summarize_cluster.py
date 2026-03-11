#!/usr/bin/env python3
'''
summarize_cluster.py

Usage:
    python summarize_cluster.py -i input.csv -o output.csv
    
Options:
    -i, --input_fasta FILE    Input Clustered FASTA file [required]
    -o, --output_csv FILE    Output directory [default: data/cluster_summary.csv]
    --verbose    Enable verbose logging
    
Output:
    A CSV file containing all sequences from the input file/directory
    
The script will compile all sequences in the input file/directory into a single
FASTA file.
'''
import csv
import os
import re
import sys
import time
import requests
from collections import defaultdict
import argparse
import logging

RFAM_SUBMIT_URL = "https://batch.rfam.org/submit-job"
HEADERS = {"accept": "application/json"}

rfam_cache = {}


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
    log_file = "summarize_cluster.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()] if verbose else [logging.StreamHandler()],
    )


def extract_rfam_from_header(header):
    """
    Extract an RFAM family ID from a sequence header.

    Parameters
    ----------
    header : str
        The sequence header to extract the RFAM family ID from.

    Returns
    -------
    str or None
        The extracted RFAM family ID, or None if no match is found.
    """
    match = re.search(r'(RF\d{5})', header)
    logging.info(f"Extracting RFAM from header: {header}, match: {match}")
    return match.group(1) if match else None


def parse_organism_and_location(theader):
    """
    Parse an organism and location from a sequence header.

    Parameters
    ----------
    theader : str
        The sequence header to parse the organism and location from.

    Returns
    -------
    str, str
        The parsed organism and location, respectively.
    """
    parts = theader.split("|")
    if len(parts) > 1:
        organism = parts[1].strip()
        if organism.lower() in ["homo sapiens", "human"]:
            organism = "Homo sapiens (human)"
        if organism.lower () in ["mus musculus (house mouse)", "mouse", "mus musculus"]:
            organism = "Mus musculus (mouse)"
        if organism.startswith("chr"):
            organism = "Mus musculus (mouse)"
    elif "vir" in theader:
        organism = "Virus"
    elif "homo" in theader.lower():
        organism = "Homo sapiens (human)"
    elif "mus" in theader.lower():
        organism = "Mus musculus (mouse)"
    else:
        organism = ""

    # formats it should take:
    # >CM000664.2:10445264-10445145 Homo sapiens ODC_IRES Homo sapiens chromosome 2, GRCh38 reference primary assembly.
    # >MNT_75-267|Homo sapiens
    # >hsa_ires_00076.1|Homo sapiens (human)|chr3:159,273,647-159,273,820+
    # >mouse_hsa_ires_00612.1 | chr15:78640294-78640461(-) | RF00223 | IRES_Bip | internal ribosome entry site (IRES)
    # >vir_ires_00069.1
    
    location = ""
    if "|" in theader:
        # location format: chr#:start-finish(strand)
        parts = theader.split("|")
        if len(parts) > 2:
            if (parts[2].strip().startswith("chr")):
                location = parts[2].strip()
            else:
                # If the third part doesn't start with "chr", it's not a location
                location = parts[1].strip().replace("(", "").replace(")", "")

    return organism, location


def query_rfam_from_sequence(sequence):
    """
    Query an RFAM family ID from a sequence.

    Parameters
    ----------
    sequence : str
        The sequence to query the RFAM API for.

    Returns
    -------
    str or None
        The queried RFAM family ID, or None if no match is found.
    """
    if sequence in rfam_cache:
        return rfam_cache[sequence]

    files = {"sequence_file": ("sequence.txt", sequence)}
    submit = requests.post(RFAM_SUBMIT_URL, headers=HEADERS, files=files)

    if submit.status_code != 200:
        rfam_cache[sequence] = ""
        return ""

    result_url = submit.json().get("resultURL")
    if not result_url:
        rfam_cache[sequence] = ""
        return ""

    time.sleep(5)

    while True:
        result = requests.get(result_url, headers=HEADERS)

        if result.status_code == 200:
            data = result.json()
            hits = data.get("hits", {})

            if hits:
                first_family = next(iter(hits.values()))
                accession = first_family[0].get("acc", "")
                rfam_cache[sequence] = accession
                logging.info(f"Found RFAM family {accession} for sequence")
                return accession

            rfam_cache[sequence] = ""
            logging.info(f"No RFAM family found for sequence")
            return ""

        elif result.status_code in (202, 204):
            time.sleep(3)
        else:
            rfam_cache[sequence] = ""
            logging.error(f"Error querying RFAM API: {result.status_code}")
            return ""

def load_processed_entries(output_csv):
    """
    Load already processed (Representative, Member) pairs
    and determine the last processed cluster number.

    Works with FASTA-based clustering (no dependency on input format).

    Parameters
    ----------
    output_csv : str
        The path to the output CSV file.

    Returns
    -------
    set
        A set of processed entries.
    int
        The last processed cluster number.

    Notes
    -----
    The output CSV should have the following columns:
    - Cluster Num
    - Representative
    - Member
    """

    processed = set()
    last_cluster_num = 0

    if not os.path.exists(output_csv):
        logging.info("No existing output file found. Starting fresh.")
        return processed, last_cluster_num

    with open(output_csv, newline='') as infile:
        reader = csv.DictReader(infile)

        required_fields = {"Cluster Num", "Representative", "Member"}
        if reader.fieldnames is None or not required_fields.issubset(reader.fieldnames):
            logging.warning("Output CSV missing expected columns. Resume may fail.")
            return processed, last_cluster_num

        for row in reader:
            representative = row.get("Representative", "").strip()
            member = row.get("Member", "").strip()

            if representative and member:
                processed.add((representative, member))

            try:
                cluster_num = int(row.get("Cluster Num", 0))
                last_cluster_num = max(last_cluster_num, cluster_num)
            except (ValueError, TypeError):
                continue

    logging.info(f"Loaded {len(processed)} processed entries.")
    logging.info(f"Last processed cluster number: {last_cluster_num}")

    return processed, last_cluster_num

def parse_clustered_fasta(fasta_path):
    """
    Parse a FASTA file into a list of clusters.

    A cluster is a representative sequence, followed by one or more
    member sequences. The representative sequence is the first sequence
    in the cluster, and is denoted by a line starting with ">" followed by
    the header of the sequence.

    The remaining sequences in the cluster are denoted by lines starting
    with ">" followed by the header of the sequence.

    Parameters
    ----------
    fasta_path : str
        The path to the FASTA file to parse.

    Returns
    -------
    list
        A list of clusters, where each cluster is a dictionary with the
        following keys:
        - "representative": The header of the representative sequence.
        - "members": A list of dictionaries, each with the following keys:
          - "header": The header of the sequence.
          - "sequence": The sequence data.

    Notes
    -----
    The FASTA file should have the following format:
    >cluster_header
    >representative_header
    sequence data
    >member_header
    sequence data

    Example:
    >representative_header
    ATCG
    >member_header
    GCAT
    """
    clusters = []
    current_cluster = None

    with open(fasta_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    while i < len(lines):

        line = lines[i]

        if line.startswith(">"):

            # Look ahead
            if i + 1 < len(lines) and lines[i + 1].startswith(">"):

                # This is a new cluster marker
                cluster_marker = line[1:].strip()
                representative_header = lines[i + 1][1:].strip()

                # Save previous cluster
                if current_cluster:
                    clusters.append(current_cluster)

                # Start new cluster
                current_cluster = {
                    "representative": representative_header,
                    "members": []
                }

                # Add representative as first member
                current_cluster["members"].append({
                    "header": representative_header,
                    "sequence": ""
                })

                i += 2
                continue

            else:
                # Normal member header
                if current_cluster is None:
                    i += 1
                    continue

                member_header = line[1:].strip()
                current_cluster["members"].append({
                    "header": member_header,
                    "sequence": ""
                })

                i += 1
                continue

        else:
            # Sequence line
            if current_cluster and current_cluster["members"]:
                current_cluster["members"][-1]["sequence"] += line

        i += 1

    # Append last cluster
    if current_cluster:
        clusters.append(current_cluster)

    return clusters

def summarize_clusters(input_fasta, output_csv):
    """
    Summarize clusters from a FASTA file into a CSV file.

    Parameters
    ----------
    input_fasta : str
        The path to the FASTA file to summarize.
    output_csv : str
        The path to the CSV file to write the cluster summary to.

    Notes
    -----
    The CSV file will contain the following columns:
    - "Cluster Num": The number of the cluster.
    - "Representative": The header of the representative sequence in the cluster.
    - "Member": The header of the member sequence in the cluster.
    - "RFAM family": The RFAM family of the member sequence.
    - "Organism": The organism of the member sequence.
    - "Location": The location of the member sequence.
    - "Length": The length of the member sequence.

    The CSV file will contain one row per member sequence in the FASTA file.
    """
    clusters = parse_clustered_fasta(input_fasta)

    logging.info(f"Parsed {len(clusters)} clusters from FASTA")

    processed_entries, last_cluster_num = load_processed_entries(output_csv)
    file_exists = os.path.exists(output_csv)

    cluster_number = last_cluster_num + 1

    with open(output_csv, "a", newline="") as outfile:
        writer = csv.writer(outfile)

        if not file_exists:
            writer.writerow([
                "Cluster Num",
                "Representative",
                "Member",
                "RFAM family",
                "Organism",
                "Location",
                "Length"
            ])

        for cluster in clusters:

            representative = cluster["representative"]
            members = cluster["members"]

            # Skip cluster if all members processed
            if all((representative, m["header"]) in processed_entries for m in members):
                continue

            logging.info(f"Processing cluster {cluster_number} ({representative})")

            for member in members:

                member_header = member["header"]
                sequence = member["sequence"]

                key = (representative, member_header)

                if key in processed_entries:
                    continue

                rfam = extract_rfam_from_header(member_header)

                if not rfam:
                    rfam = query_rfam_from_sequence(sequence)

                organism, location = parse_organism_and_location(member_header)

                writer.writerow([
                    cluster_number,
                    representative,
                    member_header,
                    rfam,
                    organism,
                    location,
                    len(sequence)
                ])

                outfile.flush()

            cluster_number += 1

    logging.info("Cluster summary complete.")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Summarize MMseqs2 cluster alignments with RFAM annotations (resume-safe)")
    parser.add_argument("-i", "--input_fasta", required=True)
    parser.add_argument("-o", "--output_csv", default="data/cluster_summary.csv")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)

    if not os.path.exists(args.input_fasta):
        logging.error(f"Input FASTA file '{args.input_fasta}' does not exist.")
        sys.exit(1)

    output_dir = os.path.dirname(args.output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    summarize_clusters(args.input_fasta, args.output_csv)