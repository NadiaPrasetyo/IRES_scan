def main():
    file_path = "data/mouse_IRES/mouse_IRES_cm_scan_table"
    output_path = "data/mouse_IRES/filtered_cmscan_table"

    existing_querys = set()

    with open(file_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            fields = line.strip().split()
            if len(fields) < 18:
                continue

            e_value = fields[15]
            inc = fields[16]

            try:
                e_value_float = float(e_value)
            except ValueError:
                e_value_float = float("inf")

            # if e_value_float < 0.05 or inc == "!":
            if inc == "!":
                outfile.write(line)


if __name__ == "__main__":
    main()
        