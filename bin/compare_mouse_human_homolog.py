#!/usr/bin/env 
'''
compare_mouse_human_homolog.py

Usage:
    python compare_mouse_human_homolog.py --mouse mouse.fasta --human human.fasta

Options:
    --mouse FILE    Mouse FASTA file [required]
    --human FILE    Human FASTA file [required]
    --output FILE   Output TSV file [default: mouse_human_alignment.csv]
    --plot          Generate a plot of the percent identity vs alignment length

The script will compare the sequences of the mouse and human homologs,
write the alignment results to a TSV file, and optionally generate a plot of
the percent identity vs alignment length.

Author: Nadia Prasetyo
'''

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
    """
    Extract the accession number from a header string.

    Parameters
    ----------
    header : str
        The header string to extract the accession number from.

    Returns
    -------
    str
        The accession number extracted from the header string.
    """
    if "|" in header:
        header = header.split("|")[0]
    else:
        header = header.split()[0]

    header = header.replace("mouse_", "")
    return header


def calculate_identity(aln1, aln2):
    """
    Calculate the percentage identity between two aligned sequences.

    Parameters
    ----------
    aln1 : str
        The first aligned sequence.
    aln2 : str
        The second aligned sequence.

    Returns
    -------
    tuple
        A tuple containing the percentage identity and the length of the aligned sequences.
    """
    matches = sum(a == b for a, b in zip(aln1, aln2))
    alnlen = len(aln1)
    pident = (matches / alnlen) * 100 if alnlen > 0 else 0
    return pident, alnlen


def calculate_bit_score(raw_score):
    """
    Calculate the bit score from a raw score.

    Parameters
    ----------
    raw_score : int
        The raw score to calculate the bit score from.

    Returns
    -------
    int
        The bit score calculated from the raw score.
    """
    return raw_score


def calculate_evalue(bit_score, m, n):
    """
    Calculate the e-value from a bit score, m and n.

    Parameters
    ----------
    bit_score : int
        The bit score to calculate the e-value from.
    m : int
        The length of the first sequence.
    n : int
        The length of the second sequence.

    Returns
    -------
    float
        The e-value calculated from the bit score, m and n.
    """
    return m * n * exp(-bit_score)


def alignment_to_gapped_sequences(alignment, seq1, seq2):
    """
    Convert a pairwise alignment into two sequences of characters.

    Parameters
    ----------
    alignment : PairwiseAlignment
        The pairwise alignment to convert into two sequences of characters.
    seq1 : str
        The first sequence of the pairwise alignment.
    seq2 : str
        The second sequence of the pairwise alignment.

    Returns
    -------
    tuple
        A tuple containing two sequences of characters. The first sequence is the
        first sequence of the pairwise alignment, with gaps represented by
        hyphens. The second sequence is the second sequence of the pairwise
        alignment, with gaps represented by hyphens.
    """
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

    """
    Load mouse and human sequences from fasta files, handling duplicates.

    Parameters:
    mouse_fasta (str): Path to mouse fasta file.
    human_fasta (str): Path to human fasta file.

    Returns:
    tuple: A tuple of two dictionaries. Each dictionary has accessions as keys
    and lists of SeqRecord objects as values.
    """
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
# Perform comparisons 
# ---------------------------------------------------------------------
def compare_homologs(mouse_dict, human_dict, output_file):

    """
    Compare homologous sequences between mouse and human.

    Parameters
    ----------
    mouse_dict : dict
        A dictionary of mouse sequences, where the keys are accessions
        and the values are lists of SeqRecord objects.
    human_dict : dict
        A dictionary of human sequences, where the keys are accessions
        and the values are lists of SeqRecord objects.
    output_file : str
        The file to which the alignment results will be written.

    Returns
    -------
    None

    Notes
    -----
    This function will write the alignment results to the specified output file.
    The format of the output file is a CSV file with the following columns:
        query : The ID of the mouse SeqRecord.
        target : The ID of the human SeqRecord.
        pident : The percentage identity between the two sequences.
        alnlen : The length of the alignment in characters.
        bits : The bit score of the alignment.
        evalue : The e-value of the alignment.
        tseq : The aligned sequence of the target (human).
        qseq : The aligned sequence of the query (mouse).
        theader : The description of the target (human).
        qheader : The description of the query (mouse).
    """
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
    """
    Plot the percent identity vs alignment length of mouse-human IRES homologs.

    Parameters
    ----------
    csv_file : str
        The path to the CSV file containing the alignment results.

    Returns
    -------
    None

    Notes
    -----
    This function reads the CSV file, generates a scatter plot of the percent
    identity vs alignment length, and saves the plot to a PNG file.
    """
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
    """
    Compare mouse and human IRES homologs.

    Parameters
    ----------
    --mouse : str
        Path to the mouse FASTA file.
    --human : str
        Path to the human FASTA file.
    --output : str
        Path to the output TSV file (default: mouse_human_alignment.csv).
    --plot : bool
        If set, generate a plot of the percent identity vs alignment length.

    Outputs
    -------
    - TSV file containing the alignment results.
    - Optional plot of the percent identity vs alignment length.

    Notes
    -----
    This function will compare the sequences of the mouse and human IRES homologs,
    write the alignment results to a TSV file, and optionally generate a plot of
    the percent identity vs alignment length.
    """
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