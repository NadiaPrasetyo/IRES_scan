
import csv


human_file = "data/human_IRES/IRESbase/human_IRES_info.txt"
virus_file = "data/virus_IRES/IRESbase/viral_IRES_info.txt"
mouse_file = "data/mouse_IRES/IRESbase/mouse_IRES_info.txt"

input_file = "data/filtered_IRES.tsv"

ires_to_seq = {}

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

with open(input_file, "r") as f, open("data/filtered_IRES_with_seq.tsv", "w", newline="") as outfile: #add the new column to the end of the line if the IRES ID is in the ires_to_seq dictionary, otherwise add an empty string
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