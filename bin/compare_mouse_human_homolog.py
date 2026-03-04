#!/usr/bin/env python3

import re
from Bio.Align import PairwiseAligner
from Bio import SeqIO
from math import exp
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def extract_accession(header):
    if "|" in header:
        header = header.split("|")[0]
    else:
        header = header.split()[0]

    header = header.replace("mouse_", "")
    return header


def calculate_identity(aln1, aln2):
    matches = sum(a == b for a, b in zip(aln1, aln2))
    alnlen = len(aln1)
    pident = (matches / alnlen) * 100 if alnlen > 0 else 0
    return pident, alnlen


def calculate_bit_score(raw_score):
    return raw_score


def calculate_evalue(bit_score, m, n):
    return m * n * exp(-bit_score)


def alignment_to_gapped_sequences(alignment, seq1, seq2):
    blocks1, blocks2 = alignment.aligned
    i = 0
    j = 0
    out1 = []
    out2 = []

    for (s1, e1), (s2, e2) in zip(blocks1, blocks2):
        if i < s1:
            out1.append(seq1[i:s1])
            out2.append("-" * (s1 - i))
            i = s1
        if j < s2:
            out1.append("-" * (s2 - j))
            out2.append(seq2[j:s2])
            j = s2

        out1.append(seq1[s1:e1])
        out2.append(seq2[s2:e2])
        i = e1
        j = e2

    if i < len(seq1):
        out1.append(seq1[i:])
        out2.append("-" * (len(seq1) - i))
    if j < len(seq2):
        out1.append("-" * (len(seq2) - j))
        out2.append(seq2[j:])

    aln1 = "".join(out1)
    aln2 = "".join(out2)

    if len(aln1) < len(aln2):
        aln1 += "-" * (len(aln2) - len(aln1))
    elif len(aln2) < len(aln1):
        aln2 += "-" * (len(aln1) - len(aln2))

    return aln1, aln2


# ---------------------------------------------------------------------
# Load sequences (UPDATED TO HANDLE DUPLICATES)
# ---------------------------------------------------------------------

def load_sequences(mouse_fasta, human_fasta):

    def load_with_duplicates(fasta_file):
        records = []
        seen_ids = defaultdict(int)

        for rec in SeqIO.parse(fasta_file, "fasta"):
            original_id = rec.id  # store original
            seen_ids[original_id] += 1

            # Rename only for output if duplicate
            if seen_ids[original_id] > 1:
                rec.id = f"{original_id}.{seen_ids[original_id]}"
            else:
                rec.id = original_id

            # Store original accession separately
            rec._accession_key = extract_accession(original_id)

            records.append(rec)

        return records

    mouse_records = load_with_duplicates(mouse_fasta)
    human_records = load_with_duplicates(human_fasta)

    mouse_dict = defaultdict(list)
    for rec in mouse_records:
        mouse_dict[rec._accession_key].append(rec)

    human_dict = defaultdict(list)
    for rec in human_records:
        human_dict[rec._accession_key].append(rec)

    return mouse_dict, human_dict


# ---------------------------------------------------------------------
# Perform comparisons (UPDATED FOR MANY-TO-MANY)
# ---------------------------------------------------------------------
def compare_homologs(mouse_dict, human_dict, output_file):

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    with open(output_file, "w") as out:
        out.write("query,target,pident,alnlen,bits,evalue,tseq,qseq,theader,qheader\n")

        for acc in mouse_dict:
            if acc not in human_dict:
                continue

            for mouse_rec in mouse_dict[acc]:
                for human_rec in human_dict[acc]:

                    mouse_seq = str(mouse_rec.seq)
                    human_seq = str(human_rec.seq)

                    alignments = aligner.align(mouse_seq, human_seq)
                    best = alignments[0]
                    score = aligner.score(mouse_seq, human_seq)

                    aln_mouse, aln_human = alignment_to_gapped_sequences(
                        best, mouse_seq, human_seq
                    )

                    pident, alnlen = calculate_identity(aln_mouse, aln_human)
                    bits = calculate_bit_score(score)
                    evalue = calculate_evalue(bits, len(mouse_rec), len(human_rec))

                    out.write(
                        f"{mouse_rec.id},{human_rec.id},"
                        f"{pident:.2f},{alnlen},{bits:.2f},{evalue:.2e},"
                        f"{aln_human},{aln_mouse},"
                        f"\"{human_rec.description}\","
                        f"\"{mouse_rec.description}\"\n"
                    )

    print(f"Alignment results written to {output_file}")

def plot_pident_vs_length(csv_file):
    df = pd.read_csv(csv_file)
    plt.figure(figsize=(8, 6))
    plt.scatter(df["alnlen"], df["pident"], alpha=0.7)
    plt.xlabel("Alignment Length")
    plt.ylabel("Percent Identity")
    plt.title("Percent Identity vs Alignment Length of Mouse-Human IRES Homologs")
    plt.grid(True)
    plt.savefig("data/pident_vs_length.png")
    print("Plot saved in data/pident_vs_length.png")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Compare mouse and human IRES homologs")
    parser.add_argument("--mouse", required=True, help="Mouse FASTA file")
    parser.add_argument("--human", required=True, help="Human FASTA file")
    parser.add_argument("--output", default="mouse_human_alignment.csv",
                        help="Output TSV file (default: mouse_human_alignment.csv)")
    parser.add_argument("--plot", action="store_true", help="Generate pident vs length plot")
    args = parser.parse_args()

    mouse_dict, human_dict = load_sequences(args.mouse, args.human)
    compare_homologs(mouse_dict, human_dict, args.output)

    if args.plot:
        plot_pident_vs_length(args.output)


if __name__ == "__main__":
    main()