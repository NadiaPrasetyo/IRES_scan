#!/usr/bin/env python3

import argparse
import csv
import logging
import os
import sys
import time
import requests
from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeoutError
import os

SEARCH_URL = "https://www.ebi.ac.uk/ebisearch/ws/rest/rfam"
FETCH_URL = "https://rfam.org/family"

FIELDS = (
    "id,description,rna_type,entry_type,rfamseq_acc,rfamseq_acc_description,scientific_name,num_seed,seq_start,seq_end,ENA,RFAM"
)

HEADERS = {
    "User-Agent": "rfam-ires-fetcher/1.0",
    "Connection": "close",
}

TIMEOUT = 30

def setup_logging(verbose=False):
    log_level =logging.INFO
    log_file = "rfam_fetch_IRES.log"
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout), logging.FileHandler(log_file) if verbose else logging.StreamHandler(sys.stdout)
        ],
    )


def safe_get(url, retries=3):
    for attempt in range(1, retries + 1):
        try:
            logging.debug(f"GET {url} (attempt {attempt})")
            r = requests.get(
                url,
                headers=HEADERS,
                timeout=(10, 30)  # (connect timeout, read timeout)
            )
            logging.debug(f"Status {r.status_code} for {url}")

            if r.status_code == 200:
                return r

        except requests.exceptions.ReadTimeout:
            logging.debug(f"READ TIMEOUT on {url}")
        except requests.exceptions.ConnectTimeout:
            logging.debug(f"CONNECT TIMEOUT on {url}")
        except requests.RequestException as e:
            logging.debug(f"REQUEST ERROR on {url}: {e}")

        time.sleep(2)

    logging.debug(f"GIVING UP on {url}")
    return None

def dedupe_entries(entries):
    seen = set()
    unique = []

    for e in entries:
        key = (
            e.get("entry_type"),
            e.get("id"),
            e.get("rfamseq_acc"),
            e.get("seq_start"),
            e.get("seq_end"),
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(e)

    return unique

def dedupe_rows(rows):
    seen = set()
    unique = []

    for r in rows:
        key = (
            r.get("entry_type"),
            r.get("rfam_id"),
            r.get("accession"),
            r.get("seq_start"),
            r.get("seq_end"),
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(r)

    return unique

def fetch_full_sequence(accession):
    """
    Fetch full FASTA sequence for an accession from ENA.
    Returns (header, sequence) or (None, None)
    """
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
    r = safe_get(url)
    if not r:
        return None, None

    lines = r.text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        return None, None

    header = lines[0]
    seq = "".join(lines[1:]).replace(" ", "").replace("\n", "")
    return header, seq

def fetch_sequence_regions(entries, organism, outdir):
    """
    Fetch subsequences ONLY for entries where:
      - entry_type == 'Sequence'
      - scientific_name matches organism (if provided)
      - Extract ONLY between seq_start and seq_end

    Saves results to: <output-dir>/<organism>_sequence_regions.fasta
    """

    fasta_out = os.path.join(outdir, f"{organism}_sequence_regions.fasta")

    count = 0

    with open(fasta_out, "w") as out_fh:
        for entry in entries:
            fields = entry.get("fields", {})

            entry_type = fields.get("entry_type", [])
            logging.info(f"Processing entry with type: {entry_type}")
            if not entry_type or entry_type[0] != "Sequence":
                continue

# "{'id': ['CM000663.2_39902015_39901795'], 'description': ['Homo sapiens IRES_L-myc'], 'rna_type': ['Cis-reg', 'IRES'], 'entry_type': ['Sequence'], 'rfamseq_acc': ['CM000663.2'], 'rfamseq_acc_description': ['Homo sapiens chromosome 1, GRCh38 reference primary assembly.'], 'scientific_name': ['Homo sapiens'], 'num_seed': [], 'seq_start': ['39902015'], 'seq_end': ['39901795'], 'ENA': ['CM000663'], 'RFAM': ['RF00261']}",CM000663.2_39902015_39901795,rfam
            
            accession_list = fields.get("rfamseq_acc", [])
            start_list = fields.get("seq_start", [])
            end_list = fields.get("seq_end", [])

            if not accession_list or not start_list or not end_list:
                continue

            accession = accession_list[0]
            start = int(start_list[0])
            end = int(end_list[0])

            header, full_seq = fetch_full_sequence(accession)
            if not full_seq:
                continue

            # Coordinates may be reversed
            if start <= end:
                subseq = full_seq[start - 1:end]
            else:
                subseq = full_seq[end - 1:start]

            fasta_header = (
                f">{accession}:{start}-{end} "
                f"{fields.get('description', [''])[0]} "
                f"{fields.get('rfamseq_acc_description', [''])[0]}"
            )

            out_fh.write(fasta_header + "\n")
            out_fh.write(subseq + "\n")

            count += 1

    logging.info(f"Saved {count} trimmed sequences to {fasta_out}")


def collect_rfam_ids(entries):
    """
    Collect RFAM family IDs from EBISearch entries.
    Works without relying on entry_type.
    """
    rfams = set()

    for entry in entries:
        fields = entry.get("fields", {})

        ids = fields.get("id", [])
        if not ids:
            continue

        entry_id = ids[0]

        # Case 1: family entry
        if entry_id.startswith("RF"):
            rfams.add(entry_id)
            continue

        # Case 2: sequence entry → extract RFAM field
        rfam_list = fields.get("RFAM", [])
        for rfam_id in rfam_list:
            if rfam_id.startswith("RF"):
                rfams.add(rfam_id)

    return sorted(rfams)

def fetch_search_results(org):
    query = f'rna_type:"IRES" AND {org}'
    params = {
        "query": query,
        "fields": FIELDS,
        "format": "json",
        "size": 1000,
    }
    r = safe_get(SEARCH_URL)
    if r is None:
        sys.exit("EBISearch unavailable")

    r = requests.get(
    SEARCH_URL,
    params=params,
    headers=HEADERS,
    timeout=(10, 30)
    )

    r.raise_for_status()
    return r.json()["entries"]

def fetch_family_regions(rfam_id, outdir):
    """
    Returns a list of dicts:
    accession, type, start, end, description, species
    """
    url = f"https://rfam.org/family/{rfam_id}/regions"
    r = safe_get(url)
    if not r:
        return []

    csv_path = os.path.join(outdir, rfam_id, f"{rfam_id}_regions.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["accession", "bits_score", "start", "end", "description", "species", "tax_id"])
        for line in r.text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 7:
                writer.writerow(parts[:7])
    logging.info(f"Saved family regions to {csv_path}")

def extract_secondary_structure(stockholm_text):
    for line in stockholm_text.splitlines():
        if line.startswith("#=GC SS_cons"):
            return line.split()[-1]
    return None


def fetch_family_data(rfam_id, outdir):
    fam_dir = os.path.join(outdir, rfam_id)
    os.makedirs(fam_dir, exist_ok=True)

    # Stockholm alignment
    r = safe_get(f"{FETCH_URL}/{rfam_id}/alignment/stockholm")
    if r:
        sto_path = os.path.join(fam_dir, f"{rfam_id}.sto")
        open(sto_path, "w").write(r.text)
        logging.info(f"Saved Stockholm alignment to {sto_path}")
        # save secondary structure in a separate file for easier access
        ss_path = os.path.join(fam_dir, f"{rfam_id}_ss.txt")
        open(ss_path, "w").write(extract_secondary_structure(r.text) or "")
        logging.info(f"Extracted secondary structure to {ss_path}")

    # Covariance model
    r = safe_get(f"{FETCH_URL}/{rfam_id}/cm")
    if r:
        cm_path = os.path.join(fam_dir, f"{rfam_id}.cm")
        open(cm_path, "w").write(r.text)
        logging.info(f"Saved covariance model to {cm_path}")

    # Structure mapping JSON
    r = safe_get(f"{FETCH_URL}/{rfam_id}/structures?content-type=application/json")
    if r:
        json_path = os.path.join(fam_dir, f"{rfam_id}_structures.json")
        open(json_path, "w").write(r.text)
        logging.info(f"Saved structure mapping JSON to {json_path}")

    url = f"{FETCH_URL}/{rfam_id}/image/cons"
    svg_path = os.path.join(fam_dir, f"{rfam_id}_structure.svg")

    r = safe_get(url)
    if r and r.headers.get("Content-Type", "").startswith("image/svg"):
        with open(svg_path, "wb") as fh:
            fh.write(r.content)
        logging.info(f"Saved structure SVG to {svg_path}")
    elif r:
        html_path = os.path.join(fam_dir, f"{rfam_id}_structure.html")
        with open(html_path, "w") as fh:
            fh.write(r.text)
        logging.info(f"Saved raw HTML to {html_path}")
    else:
        logging.info(f"No structure response for {rfam_id} at {url}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--org", help='Organism name or taxonomy code (E.g. TAXONOMY:"9606" )', required=True)
    parser.add_argument("--output-dir", help='Output directory path', required=True)
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose)

    os.makedirs(args.output_dir, exist_ok=True)

    rows = []
    entries = fetch_search_results(args.org)
    logging.info(f"Found {len(entries)} entries for {args.org}")

    entries = dedupe_entries(entries)
    logging.info(f"After deduplication: {len(entries)} entries")
    
    csv_path = os.path.join(args.output_dir, f"{args.org}_search_results.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=sorted(entries[0].keys()) if entries else [])
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry)

    
    fetch_sequence_regions(entries, args.org, args.output_dir)

    rfam_ids = collect_rfam_ids(entries)
    logging.info(f"Collected {len(rfam_ids)} unique RFAM families")


    for rfam_id in rfam_ids:
        logging.info(f"Processing family {rfam_id}")

        fetch_family_data(rfam_id, args.output_dir)
        fetch_family_regions(rfam_id, args.output_dir)

        time.sleep(0.3)


    logging.info("\nDone.")

if __name__ == "__main__":
    main()