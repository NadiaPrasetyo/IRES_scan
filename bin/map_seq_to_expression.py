# !/usr/bin/env python3
'''
map_seq_to_expression.py

Usage:
    python map_seq_to_expression.py -i input.tsv -o output.tsv

Options:
    --i, --input_file FILE    Input TSV file with IRES IDs [required]
    --o, --output FILE   Output TSV file with IRES sequences [default: data/filtered_IRES_with_seq.tsv]

Output:
    A TSV file containing IRES IDs and their corresponding sequences

The script will map IRES IDs to their corresponding sequences and add them to the input file.

Author: Nadia Prasetyo
'''
import os
import pandas as pd
import argparse

def main(input_file, output_file):
    """
    Maps the IRES sequences in the input file to the corresponding expression data in the literature excel files.

    Parameters:
    input_file (str): The input file containing the IRES sequences
    output_file (str): The output file with the IRES sequences mapped to the expression data

    Returns:
    None
    """
    forward_primer = "CTAGGGCGCGCCAGTCCT"
    reverse_primer = "CGACTCGGACCGATGGTGAG"

    directory = "data/literature/"

    files = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".xlsx")]

    # eGFP_expression (a.u)
    # Oligo_sequence

    data = []

    for file in files:
        df = pd.read_excel(file)
        for index, row in df.iterrows():
            oligo_sequence = row["Oligo_sequence"]
            seq_without_primers = oligo_sequence.replace(forward_primer, "").replace(reverse_primer, "")
            expression = row["eGFP_expression (a.u)"]
            data.append((seq_without_primers, expression))

    # add the expression data to the tsv file as a new field if the sequence matches
    # Representative	Member	RFAM	Organism	Location	Length	IRES_sequence	eGFP_expression (a.u)

    with open(input_file, "r") as f:
        lines = f.readlines()

    with open(output_file, "w") as f:
        f.write("Representative\tMember\tRFAM\tOrganism\tLocation\tLength\tIRES_sequence\teGFP_expression (a.u)\n")
        for line in lines[1:]:
            parts = line.strip().split("\t")
            ires_sequence = parts[6] if len(parts) > 6 else ""
            expression = ""
            for seq, expr in data:
                if seq in ires_sequence:
                    expression = expr
                    break
            f.write(line.strip() + "\t" + str(expression) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map IRES IDs to sequences and add them to the input file.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to the input TSV file containing IRES IDs.")
    parser.add_argument("-o", "--output_file", help="Path to the output TSV file with IRES sequences added.", default="data/filtered_IRES_with_seq.tsv")

    args = parser.parse_args()
    main(args.input_file, args.output_file)