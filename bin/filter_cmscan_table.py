def main():
    """
    This function filters a CMS-Scan output table to remove any IRES without an "inc" column value of "!". It writes the filtered results to a new file.

    Files
    -----
    mouse_IRES_cm_scan_table : str
        The path to the input CMS-Scan output table file.
    filtered_cmscan_table : str
        The path to the output filtered CMS-Scan output table file.

    Outputs
    -------
    A filtered CMS-Scan output table file containing only IRES with an "inc" column value of "!".
    """
    file_path = "data/mouse_IRES/mouse_IRES_cm_scan_table"
    output_path = "data/mouse_IRES/filtered_cmscan_table"

    with open(file_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            fields = line.strip().split()
            if len(fields) < 18:
                continue

            inc = fields[16]

            # if e_value_float < 0.05 or inc == "!":
            if inc == "!":
                outfile.write(line)


if __name__ == "__main__":
    main()
        