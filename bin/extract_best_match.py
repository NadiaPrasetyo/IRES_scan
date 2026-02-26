import argparse
import logging
import os

def setup_logging(verbose=False):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler("extract_best_match.log")] if verbose else [logging.StreamHandler()],
    )

def main():
    parser = argparse.ArgumentParser(description="Extract best matching mouse homologs from MMseqs2 search results")
    parser.add_argument("-i", "--input", default="data/mouse_IRES_homologs.tsv", help="Input TSV file with MMseqs2 search results (default: data/mouse_IRES_homologs.tsv)")
    parser.add_argument("-o", "--output", default="data/best_mouse_homologs.tsv", help="Output TSV file for best matches (default: data/best_mouse_homologs.tsv)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to console")
    args = parser.parse_args()
    setup_logging(verbose=args.verbose)

    if not os.path.exists(args.input):
        logging.error("Input file does not exist: %s", args.input)
        return

    best_matches = {}
    with open(args.input, "r") as infile:
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                logging.warning("Skipping malformed line: %s", line.strip())
                continue
            query, target, pident, alnlen, bits, tstart, tend, tlen, evalue, qaln, taln = fields
            try:
                pident = float(pident)
                bits = float(bits)
            except ValueError:
                logging.warning("Skipping header or non-numeric line: %s", line.strip())
                continue
            if query not in best_matches or bits > best_matches[query][0] or (bits == best_matches[query][0] and pident > best_matches[query][1]):
                best_matches[query] = (bits, pident, evalue, target, tstart, tend, alnlen, taln)

    with open(args.output, "w") as outfile:
        outfile.write("query\ttarget\tbits\tpident\tevalue\ttstart\ttend\talnlen\ttaln\n")
        for query, (bits, pident, evalue, target, tstart, tend, alnlen, taln) in best_matches.items():
            outfile.write(f"{query}\t{target}\t{bits}\t{pident}\t{evalue}\t{tstart}\t{tend}\t{alnlen}\t{taln}\n")
    logging.info("Best matches extracted and saved to: %s", args.output)

if __name__ == "__main__":
    main()
