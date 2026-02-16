Done:
- Fetch mice reference genome from UCSC: [Mouse (GRCm39/mm39)](https://hgdownload.gi.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz)
- Find homologous IRES sequences in human and mice: for all human IRES seq, try to find it in the mice genome using mmseqs2 search
- Cluster by sequence (between the different databases for human, mice, and virus IRES)
- Run covariance model in RFAM to classify and categorize the IRES sequences by structure
