import argparse
import logging
import sys
import csv
from pathlib import Path

def setup_logging(verbose=False):
    log_file = "compile_to_fasta.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()] if verbose else [logging.StreamHandler()],
    )

def collect_fastas(input_dir):
    # recursively find all FASTA files in the input directory
    fasta_files = list(input_dir.rglob("*.fa")) + list(input_dir.rglob("*.fasta"))
    return fasta_files


def main():
    parser = argparse.ArgumentParser(description="Compile all sequences in the input file/directory into a single FASTA file")
    parser.add_argument("--input_dir", nargs="?", default=".", help="Optional Directory to search for FASTA files (default: current directory)")
    parser.add_argument("--input_csv", help="Optional CSV file with sequence IDs and descriptions to include in the FASTA headers")
    parser.add_argument("--input_tsv", help="Optional TSV file with sequence IDs and descriptions to include in the FASTA headers")
    parser.add_argument("--output", default="data/all_sequences.fasta", help="Output FASTA file name (default: data/all_sequences.fasta)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to console and file")
    args = parser.parse_args()
    setup_logging(args.verbose)

    output_file = Path(args.output).resolve()

    input_dir = Path(args.input_dir).resolve() if args.input_dir else None
    input_csv = Path(args.input_csv).resolve() if args.input_csv else None
    input_tsv = Path(args.input_tsv).resolve() if args.input_tsv else None

    # check that there is at least one input:
    if not any([input_dir, input_csv, input_tsv]):
        logging.error("❌ No input provided. Please specify at least one of --input_dir, --input_csv, or --input_tsv.")
        sys.exit(1)

    if input_dir:
        fasta_files = collect_fastas(input_dir)
        if not fasta_files:
            logging.error("❌ No FASTA files found in the specified directory.")
            sys.exit(1)
        logging.info(f"Found {len(fasta_files)} FASTA files in {input_dir}")
        with open(output_file, "w") as out:
            for fasta in fasta_files:
                with open(fasta) as f:
                    out.write(f.read())
                    if not f.read().endswith("\n"):
                        out.write("\n")

    if input_csv:
        # ,iresite_id,ires_name,gene_name,virus_name,organism,ires_pos,ires_size,conclusion,boundary_determined,sequence,structure,notes
        with open(input_csv) as csvfile, open(output_file, "a") as out:
            reader = csv.DictReader(csvfile)
            for row in reader:
                header = f">{row['ires_name']}|{row['organism']}"
                sequence = row["sequence"]
                out.write(f"{header}\n{sequence}\n")

    if input_tsv:
        with open(input_tsv) as tsvfile, open(output_file, "a") as out:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            for row in reader:
                # robustly fetch fields with fallbacks and avoid KeyError
                ires_id = row.get("IRES ID") or row.get("Ires ID") or row.get("ires_id") or row.get("IRES_ID") or ""
                organism = row.get("Organism") or row.get("Virus Name") or ""
                sequence = (row.get("IRES sequence") or row.get("sequence") or "").strip()

                if not ires_id:
                    logging.warning(f"Skipping row with missing IRES ID: {row}")
                    continue
                if not sequence:
                    logging.warning(f"Skipping {ires_id} because sequence is empty.")
                    continue

                header = f">{ires_id}|{organism}" if organism else f">{ires_id}"
                out.write(f"{header}\n{sequence}\n")


    # remove empty lines at the end of the file and if a header does not have the sequence, remove it as well
    with open(output_file, "r") as f:
        lines = f.readlines()

    cleaned_lines = []
    current_sequence = []
    for line in lines:
        stripped_line = line.strip()
        if stripped_line.startswith(">"):
            # If we encounter a new header, save the previous sequence if it exists
            if current_sequence:
                cleaned_lines.append("".join(current_sequence) + "\n")
                current_sequence = []
            else:
                # Log headers with no sequences
                if cleaned_lines:
                    logging.warning(f"Header '{cleaned_lines[-1].strip()}' has no sequence and will be removed.")
                    cleaned_lines.pop()
            cleaned_lines.append(stripped_line + "\n")
        elif stripped_line:
            # Add sequence lines to the current sequence
            current_sequence.append(stripped_line)

    # Add the last sequence if it exists
    if current_sequence:
        cleaned_lines.append("".join(current_sequence) + "\n")
    else:
        if cleaned_lines and cleaned_lines[-1].startswith(">"):
            logging.warning(f"Header '{cleaned_lines[-1].strip()}' has no sequence and will be removed.")
            cleaned_lines.pop()
        elif not cleaned_lines:
            logging.warning("No headers or sequences found in the file.")
            cleaned_lines.pop()

    with open(output_file, "w") as f:
        f.writelines(cleaned_lines)

    logging.info(f"✅ Successfully compiled sequences into {output_file}, total {len(cleaned_lines) // 2} sequences.")

if __name__ == "__main__":
    main()

