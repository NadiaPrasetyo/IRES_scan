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
        level=logging.INFO,
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

    # Example formats:
    # chr9:87726220-87726396+
    # chr13:20,189,591-20,192,943-
    # CM000664.2:10445264-10445145

    header = header.replace(",", "")

    # chr format
    match = re.search(r'(chr[\w\.]+):(\d+)-(\d+)([+-]?)', header)
    if match:
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        strand = match.group(4) if match.group(4) else "+"
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
    """
    Expected columns:
    chromosome, nc_accession, ccds_id, cds_strand,
    exon_ordinal_position, genomic_start, genomic_end
    """

    df = pd.read_csv(ccd_file, sep="\t")

    exons_by_chr = defaultdict(list)

    for _, row in df.iterrows():
        exons_by_chr[row["chromosome"]].append({
            "start": int(row["genomic_start"]),
            "end": int(row["genomic_end"]),
            "strand": row["cds_strand"],
            "ccds_id": row["ccds_id"],
            "nc_accession": row["nc_accession"]
        })

    return exons_by_chr


############################################################
# FIND NEAREST DOWNSTREAM EXON
############################################################

def find_nearest_downstream(chrom, start, end, strand, exons_by_chr):
    """
    Strand-aware downstream logic.
    """

    if chrom not in exons_by_chr:
        return None

    candidate = None
    min_distance = float("inf")

    for exon in exons_by_chr[chrom]:

        if exon["strand"] != strand:
            continue

        if strand == "+":
            if exon["start"] >= end:
                distance = exon["start"] - end
                if distance < min_distance:
                    min_distance = distance
                    candidate = exon

        elif strand == "-":
            if exon["end"] <= start:
                distance = start - exon["end"]
                if distance < min_distance:
                    min_distance = distance
                    candidate = exon

    return candidate


############################################################
# MAIN
############################################################

def main(fasta_file, ccd_file, output_csv, blat_path):

    # check that blat executable exists
    if not shutil.which(blat_path):
        logging.error(f"BLAT executable not found at {blat_path}. Please check the path and try again.")
        sys.exit(1)

    exons_by_chr = load_ccds_exons(ccd_file)

    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):

        ires_id = record.id
        header = record.description
        sequence = str(record.seq)

        coord = parse_header_for_coordinates(header)

        if coord:
            chrom, start, end, strand = coord
        else:
            blat_hit = run_blat(sequence, blat_path=blat_path)
            if not blat_hit:
                print(f"BLAT failed for {ires_id}")
                continue
            chrom, start, end, strand = blat_hit

        nearest = find_nearest_downstream(
            chrom, start, end, strand, exons_by_chr
        )

        if nearest:
            results.append({
                "ires_id": ires_id,
                "ires_location": f"{chrom}:{start}-{end}{strand}",
                "nearest_exon_location":
                    f"{chrom}:{nearest['start']}-{nearest['end']}{strand}",
                "ccds_id": nearest["ccds_id"],
                "nc_accession": nearest["nc_accession"]
            })

    pd.DataFrame(results).to_csv(output_csv, index=False)


############################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find nearest downstream CCDS exon for IRES sequences.")
    parser.add_argument("--fasta", help="Input FASTA file")
    parser.add_argument("--ccd_file", help="CCDS exon file")
    parser.add_argument("--output_csv", help="Output CSV file")
    parser.add_argument("--blat_path", default="./blat", help="Path to BLAT executable")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose)

    fasta_file = args.fasta
    ccd_file = args.ccd_file
    output_csv = args.output_csv
    blat_path = args.blat_path

    main(fasta_file, ccd_file, output_csv, blat_path)