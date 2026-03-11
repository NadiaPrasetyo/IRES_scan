# !/usr/bin/env python3
'''
map_IRES_ID_to_seq.py

Usage:
    python map_IRES_ID_to_seq.py --input input.tsv --output output.tsv --fasta_mode

Options:
    --input FILE    Input TSV file with IRES IDs [required]
    --output FILE   Output TSV file with IRES sequences [required]
    --fasta_mode : bool    If True, the script will extract sequences from the compiled FASTA file instead of the IRES information files [default: False]

Output:
    A TSV file containing IRES IDs and their corresponding sequences

The script will map IRES IDs to their corresponding sequences and add them to the input file.

Author: Nadia Prasetyo
'''
import csv
import argparse


human_file = "data/human_IRES/IRESbase/human_IRES_info.txt"
virus_file = "data/virus_IRES/IRESbase/viral_IRES_info.txt"
mouse_file = "data/mouse_IRES/IRESbase/mouse_IRES_info.txt"

compiled_sequences_file = "data/compiled_human_mouse_virus_IRES.fasta"

def main(input_file, output_file, fasta_mode=False):
    """
    Maps IRES IDs to their corresponding sequences and adds them to the input file.

    Parameters
    ----------
    input_file : str
        Path to the input TSV file containing IRES IDs
    output_file : str
        Path to the output TSV file with IRES sequences added
    fasta_mode : bool
        If True, the script will extract sequences from the compiled FASTA file instead of the IRES information files

    Returns
    -------
    None
    """
    
    ires_to_seq = {}

    if fasta_mode:
        with open(compiled_sequences_file, "r") as f:
            header = None
            current_seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header is not None: # Save the previous ID and sequence before starting a new one
                        ires_to_seq[header] = "".join(current_seq)
                    header = line.replace(">", "")  # Get the ID without the '>' and any description
                    current_seq = []
                else:
                    current_seq.append(line)
            if header is not None:
                ires_to_seq[header] = "".join(current_seq)
        
        with open(input_file, "r") as f, open(output_file, "w", newline="") as outfile: #add the new column to the end of the line if the IRES ID is in the ires_to_seq dictionary, otherwise add an empty string
            reader = csv.reader(f, delimiter="\t")
            writer = csv.writer(outfile, delimiter="\t")
            writer.writerow(next(reader) + ["IRES_sequence"])  # add header for the new columnd
            # check the header of the input file to determine which column contains the IRES ID and extract it accordingly
            
            for fields in reader:
                ires_id = fields[1] # Num IRES	RFAM	Organism	Location	Length
                
                for header in ires_to_seq.keys():
                    if ires_id in header:
                        ires_id = header
                        break
                    
                if ires_id in ires_to_seq:
                    fields.append(ires_to_seq[ires_id])
                else:

                    fields.append("")
                writer.writerow(fields)

    else:
        with open(human_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # skip header
            for row in reader:
                ires_id = row[0]
                sequence = row[11]
                ires_to_seq[ires_id] = sequence

        with open(virus_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # skip header
            for row in reader:
                ires_id = row[0]
                sequence = row[11]
                ires_to_seq[ires_id] = sequence

        with open(mouse_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # skip header
            for row in reader:
                ires_id = row[0]
                sequence = row[7]
                ires_to_seq[ires_id] = sequence

        with open(input_file, "r") as f, open(output_file, "w", newline="") as outfile: #add the new column to the end of the line if the IRES ID is in the ires_to_seq dictionary, otherwise add an empty string
            reader = csv.reader(f, delimiter="\t")
            writer = csv.writer(outfile, delimiter="\t")
            writer.writerow(next(reader) + ["IRES_sequence"])  # add header for the new columnd
            for fields in reader:
                ires_id = fields[1]
                if "|" in ires_id:
                    ires_id = ires_id.split("|")[0]


                if ires_id in ires_to_seq:
                    fields.append(ires_to_seq[ires_id])
                else:

                    fields.append("")
                writer.writerow(fields)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map IRES IDs to sequences and add them to the input file.")
    parser.add_argument("-i", "--input_file", help="Path to the input TSV file containing IRES IDs.")
    parser.add_argument("-o", "--output_file", help="Path to the output TSV file with IRES sequences added.", default="data/filtered_IRES_with_seq.tsv")
    parser.add_argument("--fasta_mode", action="store_true", help="If set, the input will be a FASTA file instead")

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.fasta_mode)
