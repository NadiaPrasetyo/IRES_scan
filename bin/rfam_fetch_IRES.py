#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys
import time
import requests

SEARCH_URL = "https://www.ebi.ac.uk/ebisearch/ws/rest/rfam"
FETCH_URL = "https://rfam.org/family"

FIELDS = (
    "id,description,rna_type,entry_type,rfamseq_acc,rfamseq_acc_description,scientific_name,num_seed,seq_start,seq_end,ENA,RFAM"
)

HEADERS = {"User-Agent": "rfam-ires-fetcher/1.0"}

TIMEOUT = 30


def safe_get(url):
    try:
        r = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
        if r.status_code == 200:
            return r
    except requests.RequestException:
        pass
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

    r = requests.get(SEARCH_URL, params=params, headers=HEADERS)
    r.raise_for_status()
    return r.json()["entries"]

def fetch_family_regions(rfam_id):
    """
    Returns a list of dicts:
    accession, type, start, end, description, species
    """
    url = f"https://rfam.org/family/{rfam_id}/regions"
    r = safe_get(url)
    if not r:
        return []

    regions = []
    for line in r.text.splitlines():
        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < 7:
            continue

        regions.append({
            "accession": fields[0],
            "bit_score": fields[1],
            "region_type": fields[2],
            "seq_start": fields[3],
            "seq_end": fields[4],
            "description": fields[5],
            "species": fields[6],
        })

    return regions

def extract_secondary_structure(stockholm_text):
    for line in stockholm_text.splitlines():
        if line.startswith("#=GC SS_cons"):
            return line.split()[-1]
    return None


def fetch_family_data(rfam_id, outdir):
    fam_dir = os.path.join(outdir, rfam_id)
    os.makedirs(fam_dir, exist_ok=True)

    data = {
        "rfam_id": rfam_id,
        "sequence": None,
        "secondary_structure": None,
        "covariation_cm": None,
        "stockholm": None,
        "fasta": None,
        "structure_json": None,
    }

    # FASTA (ungapped)
    r = safe_get(f"{FETCH_URL}/{rfam_id}/alignment/fastau")
    if r:
        fasta_path = os.path.join(fam_dir, f"{rfam_id}.fa")
        open(fasta_path, "w").write(r.text)
        data["fasta"] = fasta_path

        # take first sequence
        seq = "".join(
            line.strip()
            for line in r.text.splitlines()
            if not line.startswith(">")
        )
        data["sequence"] = seq or None

    # Stockholm alignment
    r = safe_get(f"{FETCH_URL}/{rfam_id}/alignment/stockholm")
    if r:
        sto_path = os.path.join(fam_dir, f"{rfam_id}.sto")
        open(sto_path, "w").write(r.text)
        data["stockholm"] = sto_path
        data["secondary_structure"] = extract_secondary_structure(r.text)

    # Covariance model
    r = safe_get(f"{FETCH_URL}/{rfam_id}/cm")
    if r:
        cm_path = os.path.join(fam_dir, f"{rfam_id}.cm")
        open(cm_path, "w").write(r.text)
        data["covariation_cm"] = cm_path

    # Structure mapping JSON
    r = safe_get(f"{FETCH_URL}/{rfam_id}/structures?content-type=application/json")
    if r:
        json_path = os.path.join(fam_dir, f"{rfam_id}_structures.json")
        open(json_path, "w").write(r.text)
        data["structure_json"] = json_path

    return data


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

    rows = []

    for rfam_id in rfam_ids:
        print(f"Processing family {rfam_id}")

        fam_data = fetch_family_data(rfam_id, args.output_dir)
        regions = fetch_family_regions(rfam_id)

        if not regions:
            # still record the family even if no regions are returned
            rows.append({
                "rfam_id": rfam_id,
                **fam_data
            })
            continue

        for region in regions:
            row = {
                "rfam_id": rfam_id,
                **fam_data,
                **region,
            }
            rows.append(row)

        time.sleep(0.3)

    rows = dedupe_rows(rows)
    
    # Write CSV
    csv_path = os.path.join(args.output_dir, "rfam_ires_summary.csv")

    if not rows:
        print("\nNo valid IRES entries produced rows.")
        print("This usually means:")
        print("  - sequence entries lacked seq_start/seq_end")
        print("  - or RFAM returned metadata-only hits")
        print(f"\nEmpty CSV not written: {csv_path}")
        return

    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=sorted(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"\nDone. CSV written to {csv_path}")


if __name__ == "__main__":
    main()