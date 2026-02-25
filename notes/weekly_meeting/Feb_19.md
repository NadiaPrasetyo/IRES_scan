Done:
- Fetch mice reference genome from UCSC: [Mouse (GRCm39/mm39)](https://hgdownload.gi.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz)
- Find homologous IRES sequences in human and mice: for all human IRES seq, try to find it in the mice genome using mmseqs2 search

To Do:
- blat mmseqs for human-human matches
- human to mouse: mmseqs + covariance model
- Search mouse homologs using the downstream closest exon (do a search in a protein level)
- glocal alignment tool to find the IRES: very best match
- shorter kmer for IRES search 
- main priority: conserved between human and mouse, length 
- create a scatter plot between percent identity and length: between mice and humans
- Run covariance model in RFAM to classify and categorize the IRES sequences by structure



- interquartile range in brackets, and if it collapses enough, it can be turned into a box and whiskers
- collapse the ssl as well
- prune the essential table into one or two names
- include the interquartile ranges in the table
- instead of accessions, we report how many paralogs there are for each
- put PCA plot in supplementary