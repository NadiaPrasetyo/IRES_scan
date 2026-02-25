Done:
- interquartile range in brackets, and if it collapses enough, it can be turned into a box and whiskers
- collapse the ssl as well
- prune the essential table into one or two names
- include the interquartile ranges in the table
- instead of accessions, we report how many paralogs there are for each
- put PCA plot in supplementary

- blat mmseqs for human-human matches
- Search mouse homologs using the downstream closest exon (do a search in a protein level)
- human to mouse: mmseqs + covariance model

To Do:
- glocal alignment tool to find the IRES: very best match
- shorter kmer for IRES search
- main priority: conserved between human and mouse, length
- create a scatter plot between percent identity and length: between mice and humans
- Run covariance model in RFAM to classify and categorize the IRES sequences by structure
- Human-mouse: conserved
- for the viruses: sample for length and diversity
- run CM models on the upstream sequences (max of the human seq to exon distance and add a buffer (500 nt))

- build a maximum likelihood (proML, RAXML) tree for the SSL alignment
- put them in a tree viewer