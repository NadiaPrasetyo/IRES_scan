#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys
import time
import requests
from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeoutError
import os
from datetime import datetime
import traceback

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

DEBUG = True

def debug(msg):
    if DEBUG:
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] {msg}", flush=True)


def safe_get(url, retries=3):
    for attempt in range(1, retries + 1):
        try:
            debug(f"GET {url} (attempt {attempt})")
            r = requests.get(
                url,
                headers=HEADERS,
                timeout=(10, 30)  # (connect timeout, read timeout)
            )
            debug(f"Status {r.status_code} for {url}")

            if r.status_code == 200:
                return r

        except requests.exceptions.ReadTimeout:
            debug(f"READ TIMEOUT on {url}")
        except requests.exceptions.ConnectTimeout:
            debug(f"CONNECT TIMEOUT on {url}")
        except requests.RequestException as e:
            debug(f"REQUEST ERROR on {url}: {e}")

        time.sleep(2)

    debug(f"GIVING UP on {url}")
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

def extract_region_sequences(rfam_id, outdir):
    """
    For each region in {rfam_id}_regions.csv, fetch the full sequence,
    extract start:end, and save as FASTA.
    """
    fam_dir = os.path.join(outdir, rfam_id)
    regions_csv = os.path.join(fam_dir, f"{rfam_id}_regions.csv")

    if not os.path.exists(regions_csv):
        print(f"No regions CSV found for {rfam_id}")
        return

    with open(regions_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        for i, row in enumerate(reader, 1):
            accession = row["accession"]
            start = int(row["start"])
            end = int(row["end"])

            debug(f"[{rfam_id}] Region {i}: {accession}:{start}-{end}")

            out_fa = os.path.join(
                fam_dir, f"{accession}_{start}_{end}.fa"
            )

            # ⏭️ Skip if already processed
            if os.path.exists(out_fa):
                debug(f"[{accession}] {start}-{end} already exists, skipping")
                continue

            try:
                header, full_seq = fetch_full_sequence(accession)
                if not full_seq:
                    debug(f"[{accession}] FAILED to fetch full sequence")
                    continue

                if end > len(full_seq):
                    debug(
                        f"[{accession}] WARNING end={end} > length={len(full_seq)}"
                    )

                subseq = full_seq[start - 1:end]

                out_fa = os.path.join(
                    fam_dir, f"{accession}_{start}_{end}.fa"
                )

                with open(out_fa, "w") as out:
                    out.write(f">{accession}:{start}-{end}\n")
                    for j in range(0, len(subseq), 60):
                        out.write(subseq[j:j+60] + "\n")

                debug(f"[{accession}] saved {len(subseq)} bp")

                time.sleep(0.5) # Time out for 0.5 seconds between requests

            except Exception:
                debug(f"EXCEPTION processing {accession}")
                traceback.print_exc()

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
    print(f"Saved family regions to {csv_path}")

def extract_secondary_structure(stockholm_text):
    for line in stockholm_text.splitlines():
        if line.startswith("#=GC SS_cons"):
            return line.split()[-1]
    return None


def fetch_family_data(rfam_id, outdir):
    fam_dir = os.path.join(outdir, rfam_id)
    os.makedirs(fam_dir, exist_ok=True)

    # FASTA (ungapped)
    r = safe_get(f"{FETCH_URL}/{rfam_id}/alignment/fastau")
    if r:
        fasta_path = os.path.join(fam_dir, f"{rfam_id}.fa")
        open(fasta_path, "w").write(r.text)
        print(f"Saved ungapped FASTA to {fasta_path}")

    # Stockholm alignment
    r = safe_get(f"{FETCH_URL}/{rfam_id}/alignment/stockholm")
    if r:
        sto_path = os.path.join(fam_dir, f"{rfam_id}.sto")
        open(sto_path, "w").write(r.text)
        print(f"Saved Stockholm alignment to {sto_path}")
        # save secondary structure in a separate file for easier access
        ss_path = os.path.join(fam_dir, f"{rfam_id}_ss.txt")
        open(ss_path, "w").write(extract_secondary_structure(r.text) or "")
        print(f"Extracted secondary structure to {ss_path}")

    # Covariance model
    r = safe_get(f"{FETCH_URL}/{rfam_id}/cm")
    if r:
        cm_path = os.path.join(fam_dir, f"{rfam_id}.cm")
        open(cm_path, "w").write(r.text)
        print(f"Saved covariance model to {cm_path}")

    # Structure mapping JSON
    r = safe_get(f"{FETCH_URL}/{rfam_id}/structures?content-type=application/json")
    if r:
        json_path = os.path.join(fam_dir, f"{rfam_id}_structures.json")
        open(json_path, "w").write(r.text)
        print(f"Saved structure mapping JSON to {json_path}")

    url = f"{FETCH_URL}/{rfam_id}/image/cons"
    svg_path = os.path.join(fam_dir, f"{rfam_id}_structure.svg")

    r = safe_get(url)
    if r and r.headers.get("Content-Type", "").startswith("image/svg"):
        with open(svg_path, "wb") as fh:
            fh.write(r.content)
        print(f"Saved structure SVG to {svg_path}")
    elif r:
        html_path = os.path.join(fam_dir, f"{rfam_id}_structure.html")
        with open(html_path, "w") as fh:
            fh.write(r.text)
        print(f"Saved raw HTML to {html_path}")
    else:
        print(f"No structure response for {rfam_id} at {url}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--org", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    rows = []
    entries = fetch_search_results(args.org)
    print(f"Found {len(entries)} entries for {args.org}")

    entries = dedupe_entries(entries)
    print(f"After deduplication: {len(entries)} entries")
    
    csv_path = os.path.join(args.output_dir, f"{args.org}_search_results.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=sorted(entries[0].keys()) if entries else [])
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry)

    rfam_ids = collect_rfam_ids(entries)
    print(f"Collected {len(rfam_ids)} unique RFAM families")


    for rfam_id in rfam_ids:
        print(f"Processing family {rfam_id}")

        fetch_family_data(rfam_id, args.output_dir)
        fetch_family_regions(rfam_id, args.output_dir)
        extract_region_sequences(rfam_id, args.output_dir)

        time.sleep(0.3)


    print("\nDone.")

if __name__ == "__main__":
    main()