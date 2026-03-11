# !/usr/bin/env python3
'''
iresite_scrape_IRES.py

Usage:
    python iresite_scrape_IRES.py --output-viral output_viral.csv --output-human output_human.csv --output-mouse output_mouse.csv

Options:
    --output-viral STR  Output CSV file for viral IRES data (default: data/iresite_viral.csv)
    --output-human STR  Output CSV file for human IRES data (default: data/iresite_human.csv)
    --output-mouse STR  Output CSV file for mouse IRES data (default: data/iresite_mouse.csv)

Output:
    A CSV file containing the extracted IRES data from IRESite.org

The script will scrape IRESite.org for viral, human, and mouse IRES data and
write it to a CSV file.

Author: Nadia Prasetyo
'''
from cProfile import label
import requests
from bs4 import BeautifulSoup
import pandas as pd
import re
import time
from urllib.parse import urljoin
import warnings
from urllib3.exceptions import InsecureRequestWarning
import argparse
import re
from html import unescape

BASE_URL = "https://iresite.org/"
VIRAL_URL = "https://iresite.org/IRESite_web.php?page=browse_viral_transcripts"
CELLULAR_URL = "https://iresite.org/IRESite_web.php?page=browse_cellular_transcripts"

HEADERS = {"User-Agent": "Mozilla/5.0"}

OUTPUT_COLUMNS = [
    "iresite_id",
    "ires_name",
    "gene_name",
    "virus_name",
    "organism",
    "ires_pos",
    "ires_size",
    "conclusion",
    "boundary_determined",
    "sequence",
    "structure",
    "notes",
]

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

warnings.simplefilter("ignore", InsecureRequestWarning)

def clean(text):
    return re.sub(r"\s+", " ", text.strip()) if text else None


def fetch(url):
    r = requests.get(url, headers=HEADERS, timeout=30, verify=False)
    r.raise_for_status()
    return r.text

# ------------------------------------------------------------
# Parse main browse tables
# ------------------------------------------------------------
def parse_browse_page(url, mode):
    """
    Parse IRESite browse pages into a list of dictionaries.

    Parameters
    ----------
        url (str): URL of the IRESite browse page
        mode (str): "viral" or "cellular" to specify the type of browse page

    Returns
    -------
        list: A list of dictionaries containing the extracted information
    """
    html = fetch(url)
    soup = BeautifulSoup(html, "lxml")

    table = soup.find("table", class_="tabular_listing")
    if table is None:
        print("ERROR: tabular_listing table not found")
        return []

    rows = table.find_all("tr")[1:]  # skip header
    records = []

    for row in rows:
        tds = row.find_all("td")
        if len(tds) != 15:
            continue

        # iresite_id
        iresite_id = tds[0].get_text(strip=True)
        if not iresite_id.isdigit():
            continue

        detail_url = f"{BASE_URL}IRESite_web.php?page=view&entry_id={iresite_id}"

        # IRES names and sizes
        ires_names = [clean(x) for x in tds[1].stripped_strings]
        ires_sizes = [clean(x) for x in tds[13].stripped_strings]

        # Virus / gene name
        src_anchor = tds[3].find("a")
        src_name = clean(src_anchor.text) if src_anchor else clean(tds[3].get_text())
        virus_name = src_name if mode == "viral" else None
        gene_name = src_name if mode == "cellular" else None

        # --------------------------------------------------
        # Organism extraction (DEBUGGABLE + SAFE)
        # --------------------------------------------------
        host_td = tds[7]
        host_span = host_td.find("span")
        organism = None

        if host_span:
            span_text = clean(host_span.get_text())
            onmouseover = host_span.get("onmouseover", "")

            print(f"\n[DEBUG] IRESite {iresite_id} ({mode})")
            print(f"  span_text   = {span_text}")
            print(f"  onmouseover = {onmouseover}")

            if onmouseover:
                unescaped = unescape(onmouseover)
                print(f"  unescaped   = {unescaped}")

                m = re.search(r"<a[^>]*>([^<]+)</a>", unescaped, flags=re.IGNORECASE)
                if m:
                    organism = clean(m.group(1))
                    print(f"  extracted   = {organism}")

                    if mode == "viral":
                        organism = re.sub(r"\s+virus$", "", organism, flags=re.IGNORECASE)
                        print(f"  stripped    = {organism}")
                else:
                    print("  regex FAILED")
            else:
                print("  no onmouseover")

            # Fallback: span text ONLY if nothing extracted
            if not organism:
                organism = span_text
                print(f"  fallback    = {organism}")

        else:
            organism = clean(host_td.get_text())
            print(f"\n[DEBUG] IRESite {iresite_id}: no span, using td text = {organism}")

        # --------------------------------------------------
        # normalize multi-IRES rows
        # --------------------------------------------------
        n = max(len(ires_names), len(ires_sizes))
        ires_names += [None] * (n - len(ires_names))
        ires_sizes += [None] * (n - len(ires_sizes))

        for name, size in zip(ires_names, ires_sizes):
            records.append({
                "iresite_id": iresite_id,
                "ires_name": name,
                "ires_size": size,
                "organism": organism,
                "virus_name": virus_name,
                "gene_name": gene_name,
                "detail_url": detail_url,
            })

    print(f"Found {len(records)} {mode} browse rows")
    return records


# ------------------------------------------------------------
# Parse individual IRES pages
# ------------------------------------------------------------
def parse_ires_page(record):
    """
    Parse an individual IRES page, given the record with the URL to fetch.

    Parameters
    ----------
    record : dict
        A dictionary containing the following keys:
        - iresite_id
        - detail_url

    Returns a list of dictionaries with the extracted information.
    """
    
    html = fetch(record["detail_url"])
    soup = BeautifulSoup(html, "lxml")

    # 1️⃣ Find the IRES section by <b>IRES:</b> marker
    ires_b_tag = soup.find("b", string="IRES:")
    if not ires_b_tag:
        print(f"No IRES section found for {record['iresite_id']}")
        return []

    # 2️⃣ The containing table with class black_border_table
    ires_table = ires_b_tag.find_next("td", class_="black_border_table")
    if not ires_table:
        print(f"No black_border_table found for {record['iresite_id']}")
        return []

    # ------------------------------
    # Extract basic fields
    # ------------------------------
    def extract_text(label):
        for td in ires_table.find_all("td"):
            b = td.find("b")
            if not b:
                continue

            label_text = b.get_text(" ", strip=True)
            if label in label_text:
                tts = td.find_all("tt")
                if tts:
                    return clean(tts[-1].get_text())  # ⬅️ key line
        return None


    ires_name = extract_text("The IRES name") or record["ires_name"]
    ires_pos = extract_text("The IRES absolute position")  # e.g., 6340-6538
    conclusion = extract_text("Conclusion")
    boundary = extract_text("How IRES boundaries were determined")

    # ------------------------------
    # Sequence and structure
    # ------------------------------
    sequence = None
    structure = None
    textarea = ires_table.find("textarea")
    if textarea:
        text = textarea.get_text()
        # Split lines and remove annotations like ^123 or spaces
        seq_lines = []
        struct_lines = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            # Sequence lines contain only ACGTU
            seq_match = re.fullmatch(r"[ACGTUacgtu\s]+", line)
            struct_match = re.fullmatch(r"[().]+", line)
            if seq_match:
                seq_lines.append(re.sub(r"[^ACGTUacgtu]", "", line).upper())
            elif struct_match:
                struct_lines.append(line)
        if seq_lines:
            sequence = "".join(seq_lines)
        if struct_lines:
            structure = "\n".join(struct_lines)

    # ------------------------------
    # Remarks / notes
    # ------------------------------
    remarks_tag = ires_table.find("b", string="Remarks:")
    notes = None
    if remarks_tag:
        pre = remarks_tag.find_next("pre")
        if pre:
            notes = clean(pre.get_text())

    return [{
        "iresite_id": record["iresite_id"],
        "ires_name": ires_name,
        "gene_name": record["gene_name"],
        "virus_name": record["virus_name"],
        "organism": record["organism"],
        "ires_pos": ires_pos,
        "ires_size": record["ires_size"],
        "conclusion": conclusion,
        "boundary_determined": boundary,
        "sequence": sequence,
        "structure": structure,
        "notes": notes,
    }]


# ------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------

def scrape_all(mode="viral"):
    """
    Scrape all IRES records from either the viral or cellular tables.

    Parameters
    ----------
    mode : str, default="viral"
        Either "viral" or "cellular", indicating which table to scrape.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the scraped IRES records, with columns
        as specified in the OUTPUT_COLUMNS constant.
    """
    if mode == "viral":
        print("Scraping viral IRES table...")
        records = parse_browse_page(VIRAL_URL, mode="virus")
    elif mode == "cellular":
        print("Scraping cellular IRES table...")
        records = parse_browse_page(CELLULAR_URL, mode="cellular")
    else:
        raise ValueError("mode must be 'viral' or 'cellular'")

    print(f"Found {len(records)} {mode} IRES rows")

    records_parsed = []

    for rec in records:
        print(f"Fetching IRES {rec['iresite_id']}")
        try:
            parsed = parse_ires_page(rec)
            if parsed:
                records_parsed.extend(parsed)
            time.sleep(0.5)  # be polite
        except Exception as e:
            print(f"Failed {rec['iresite_id']}: {e}")

    return pd.DataFrame(records_parsed, columns=OUTPUT_COLUMNS)


def main():  
    """
    Main function of the script.

    Scrapes IRESite database for viral, human, and mouse IRES data.

    Parameters
    ----------
    --output-viral : str
        Output CSV file for viral IRES data (default: data/iresite_viral.csv)
    --output-human : str
        Output CSV file for human IRES data (default: data/iresite_human.csv)
    --output-mouse : str
        Output CSV file for mouse IRES data (default: data/iresite_mouse.csv)

    Returns
    -------
    None
    """
    argparser = argparse.ArgumentParser(description="Scrape IRESite database")
    argparser.add_argument(
        "--output-viral",
        type=str,
        default="data/iresite_viral.csv",
        help="Output CSV file for viral IRES data",
    )
    argparser.add_argument(
        "--output-human",
        type=str,
        default="data/iresite_human.csv",
        help="Output CSV file for human IRES data",
    )
    argparser.add_argument(
        "--output-mouse",
        type=str,
        default="data/iresite_mouse.csv",
        help="Output CSV file for mouse IRES data",
    )
    args = argparser.parse_args()
    df_viral = scrape_all("viral")
    print(df_viral.head())

    df_cellular = scrape_all("cellular")
    print(df_cellular.head())

    # filter to only keep records with organism of interest
    organisms_of_interest = ["Mus musculus"]
    df_cellular_mouse = df_cellular[
        df_cellular["organism"].fillna("").str.contains(
            "|".join(organisms_of_interest), case=False
        )   
    ]

    df_cellular_human = df_cellular[
        df_cellular["organism"].fillna("").str.contains(
            "Homo sapiens", case=False
        )
    ]

    df_viral.to_csv(args.output_viral)
    df_cellular_mouse.to_csv(args.output_mouse)
    df_cellular_human.to_csv(args.output_human)

    print(f"Saved {len(df_viral)} records to {args.output_viral}")
    print(f"Saved {len(df_cellular_mouse)} records to {args.output_mouse}")
    print(f"Saved {len(df_cellular_human)} records to {args.output_human}")


if __name__ == "__main__":
    main()