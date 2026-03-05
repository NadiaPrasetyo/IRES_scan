import os
import pandas as pd


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
tsv_file = "data/filtered_IRES_with_seq.tsv"
# Representative	Member	RFAM	Organism	Location	Length	IRES_sequence	eGFP_expression (a.u)

with open(tsv_file, "r") as f:
    lines = f.readlines()

with open("data/filtered_IRES_with_expression.tsv", "w") as f:
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