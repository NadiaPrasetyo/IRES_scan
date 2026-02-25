#!/usr/bin/env python3
import shutil
import sys
import re
import subprocess
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import logging
import argparse

############################################################
# CONFIG
############################################################

HG38_2BIT = "hg38.2bit"   # path to hg38 2bit genome
# BLAT_PATH = "blat"       # blat executable in PATH

def setup_logging(verbose=False):
    log_file = "find_nearest_downstream_ccd_exon.log"
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()] if verbose else [logging.StreamHandler()],
    )


############################################################
# PARSE FASTA HEADER FOR COORDINATES
############################################################

def parse_header_for_coordinates(header):
    """
    Extract chromosome, start, end, strand from header.
    Returns: (chrom, start, end, strand) or None
    """

    # Remove commas
    header_clean = header.replace(",", "")

    # Attempt to find chrX:start-end[+/-] anywhere in the string
    match = re.search(r'chr?(\w+):(\d+)-(\d+)([+-]?)', header_clean)
    if match:
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        strand = match.group(4) if match.group(4) else "+"
        return chrom, start, end, strand

    # If chr pattern not found, try generic pattern like CM000664.2:10445264-10445145
    match_generic = re.search(r'([\w\.]+):(\d+)-(\d+)([+-]?)?', header_clean)
    if match_generic:
        # Heuristic: if the first group looks like CM*, hsa*, etc., try to extract chromosome from description
        chrom_candidate = match_generic.group(1)
        start = int(match_generic.group(2))
        end = int(match_generic.group(3))
        strand = match_generic.group(4) if match_generic.group(4) else "+"

        # Check for "chromosome N" in header
        chrom_match = re.search(r'chromosome (\w+)', header_clean, re.IGNORECASE)
        if chrom_match:
            chrom = chrom_match.group(1)
        else:
            chrom = chrom_candidate  # fallback

        return chrom, start, end, strand

    return None


############################################################
# RUN BLAT IF NO COORDINATES
############################################################

def run_blat(sequence, temp_fasta="temp.fa", output_psl="temp.psl", blat_path="./blat"):
    """
    Run BLAT against hg38 to get genomic coordinates.
    Returns best hit: (chrom, start, end, strand)
    """

    with open(temp_fasta, "w") as f:
        f.write(">query\n")
        f.write(sequence + "\n")

    subprocess.run([
        blat_path,
        HG38_2BIT,
        temp_fasta,
        output_psl
    ], check=True)

    # parse PSL
    best_hit = None
    best_score = 0

    with open(output_psl) as f:
        for line in f:
            if line.startswith("match"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 21:
                continue

            matches = int(fields[0])
            strand = fields[8]
            chrom = fields[13]
            start = int(fields[15])
            end = int(fields[16])

            if matches > best_score:
                best_score = matches
                best_hit = (chrom, start, end, strand)

    return best_hit


############################################################
# LOAD CCDS GENOMIC EXON FILE
############################################################

def load_ccds_exons(ccd_file):

    df = pd.read_csv(ccd_file, sep="\t", dtype=str, low_memory=False)
    df.columns = df.columns.str.replace("#", "").str.strip()

    # Keep only Public CCDS
    df = df[
    (df["ccds_status"] == "Public") &
    (df["cds_locations"].notna()) &
    (df["cds_locations"] != "[]")
    ]

    exons_by_chr = defaultdict(list)

    for _, row in df.iterrows():

        chrom = "chr" + row["chromosome"].replace("chr", "")
        strand = row["cds_strand"]
        ccds_id = row["ccds_id"]

        # Parse exon list
        locations = str(row["cds_locations"]).strip()

        # Skip empty or malformed entries
        if not locations or locations == "[]":
            continue

        locations = locations.strip("[]")

        exon_ranges = [x.strip() for x in locations.split(",") if x.strip()]

        for exon_number, exon in enumerate(exon_ranges, start=1):

            if "-" not in exon:
                continue

            start_str, end_str = exon.split("-", 1)

            # Skip if either coordinate missing
            if not start_str or not end_str:
                continue

            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                logging.error("Bad exon entry:", exon, "in row:", row.to_dict())
                continue

            exons_by_chr[chrom].append({
                "start": start,
                "end": end,
                "strand": strand,
                "ccds_id": ccds_id,
                "exon_number": exon_number
            })

    return exons_by_chr

def load_ccds_protein_fasta(protein_fasta):

    protein_dict = {}

    for record in SeqIO.parse(protein_fasta, "fasta"):
        header = record.description

        # Example: CCDS2.2|Hs110|chr1|exon1
        parts = header.split("|")
        ccds_id = parts[0]
        exon_part = parts[-1]  # exon1

        key = f"{ccds_id}_{exon_part}"
        protein_dict[key] = str(record.seq)

    return protein_dict


############################################################
# FIND NEAREST DOWNSTREAM EXON
############################################################

def find_nearest_downstream(chrom, start, end, exons_by_chr):

    if chrom not in exons_by_chr:
        logging.info(f"No exons found for chromosome {chrom}")
        return None

    candidate = None
    min_distance = float("inf")

    for exon in exons_by_chr[chrom]:

        # Genomic downstream only (higher coordinate)
        if exon["start"] >= end:
            distance = exon["start"] - end
            if distance < min_distance:
                min_distance = distance
                candidate = exon

    return candidate


############################################################
# MAIN
############################################################

def main(fasta_file, ccd_file, output_csv, blat_path, ccds_protein_fasta=None):

    # check that blat executable exists
    if not shutil.which(blat_path):
        logging.error(f"BLAT executable not found at {blat_path}. Please check the path and try again.")
        sys.exit(1)

    exons_by_chr = load_ccds_exons(ccd_file)

    protein_dict = None
    if ccds_protein_fasta:
        protein_dict = load_ccds_protein_fasta(ccds_protein_fasta)

    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):

        ires_id = record.id
        header = record.description
        sequence = str(record.seq)

        coord = parse_header_for_coordinates(header)

        used_blat = False

        if coord:
            chrom, start, end, strand = coord
            logging.debug(f"Parsed coordinates for IRES ID={ires_id} Location= {chrom}:{start}-{end}{strand}")
        else:
            blat_hit = run_blat(sequence, blat_path=blat_path)
            if not blat_hit:
                logging.warning(f"BLAT failed for {ires_id}")
                continue
            chrom, start, end, strand = blat_hit
            logging.debug(f"BLAT coordinates for IRES ID={ires_id} Location= {chrom}:{start}-{end}{strand}")
            used_blat = True

        nearest = find_nearest_downstream(
            chrom, start, end, exons_by_chr
        )

        # 🔁 BACKUP: if nothing found, try BLAT remapping
        if nearest is None and not used_blat:
            logging.info(f"No downstream exon found using header coordinates for {ires_id}, trying BLAT remap")

            blat_hit = run_blat(sequence, blat_path=blat_path)

            if blat_hit:
                blat_chrom, blat_start, blat_end, blat_strand = blat_hit

                logging.debug(
                    f"BLAT remap for {ires_id}: "
                    f"{blat_chrom}:{blat_start}-{blat_end}{blat_strand}"
                )

                nearest = find_nearest_downstream(
                    blat_chrom,
                    blat_start,
                    blat_end,
                    exons_by_chr
                )

                # If BLAT worked, update coordinates
                if nearest:
                    chrom, start, end, strand = blat_chrom, blat_start, blat_end, blat_strand
                    logging.info(f"BLAT rescue succeeded for {ires_id}")
            else:
                logging.warning(f"BLAT rescue failed for {ires_id}")

        if nearest and protein_dict:
            logging.info(f"Found nearest exon for {ires_id}: {nearest['ccds_id']} exon {nearest['exon_number']} {nearest['strand']}")
            protein_key = f"{nearest['ccds_id']}_exon{nearest['exon_number']}"
            protein_seq = protein_dict.get(protein_key, "NA")

            results.append({
                "ires_id": ires_id,
                "ires_location": f"{chrom}:{start}-{end}{strand}",
                "nearest_exon_location":
                    f"{chrom}:{nearest['start']}-{nearest['end']}{nearest['strand']}",
                "ccds_id": nearest["ccds_id"],
                "exon_number": nearest["exon_number"],
                "protein_sequence": protein_seq
            })

        elif nearest:
            logging.info(f"Found nearest exon for {ires_id}: {nearest['ccds_id']} exon {nearest['exon_number']} without protein sequence")
            results.append({
                "ires_id": ires_id,
                "ires_location": f"{chrom}:{start}-{end}{strand}",
                "nearest_exon_location": f"{chrom}:{nearest['start']}-{nearest['end']}{nearest['strand']}",
                "ccds_id": nearest["ccds_id"],
                "exon_number": nearest["exon_number"],
                "protein_sequence": "NA"
            })
        else:
            logging.info(f"No downstream exon found for {ires_id}")
            results.append({
                "ires_id": ires_id,
                "ires_location": f"{chrom}:{start}-{end}{strand}",
                "nearest_exon_location": "NA",
                "ccds_id": "NA",
                "exon_number": "NA",
                "protein_sequence": "NA"
            })

    pd.DataFrame(results).to_csv(output_csv, index=False)


############################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find nearest downstream CCDS exon for IRES sequences.")
    parser.add_argument("--fasta", help="Input FASTA file")
    parser.add_argument("--ccd_file", help="CCDS exon file")
    parser.add_argument("--ccds_protein_fasta", help="CCDS exon protein FASTA")
    parser.add_argument("--output_csv", help="Output CSV file")
    parser.add_argument("--blat_path", default="./blat", help="Path to BLAT executable")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose)

    fasta_file = args.fasta
    ccd_file = args.ccd_file
    output_csv = args.output_csv
    blat_path = args.blat_path
    ccds_protein_fasta = args.ccds_protein_fasta

    main(fasta_file, ccd_file, output_csv, blat_path, ccds_protein_fasta)