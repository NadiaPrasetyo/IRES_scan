# IRES_scan
A pipeline to search to prioritize IRES sequences to test for circRNA application.

## Folder structure
```
project/
│
├── bin/                         # Analysis scripts
│   ├── compare_mouse_human_homolog.py
│   ├── compile_to_fasta.py
│   ├── extract_best_match.py
│   ├── extract_genome_sequence.py
│   ├── filter_cmscan_table.py
│   ├── find_mouse_homolog.py
│   ├── find_nearest_downstream_ccd_exon.py
│   ├── get_mouse_homolog_seq.py
│   ├── human_mouse_pident_length.py
│   ├── iresite_scrape_IRES.py
│   ├── map_IRES_ID_to_seq.py
│   ├── map_seq_to_expression.py
│   ├── rfam_fetch_IRES.py
│   ├── run_mmseqs_cluster.py
│   └── summarize_cluster.py
│
├── data/                        # Input data and intermediate outputs
│
│   ├── human_IRES/
│   │   ├── IRESbase/
│   │   ├── IRESite/
│   │   ├── Rfam/
│   │   ├── compiled_covariance_model.cm
│   │   ├── compiled_human_sequences.fasta
│   │   ├── nearest_exons.csv
│   │   └── nearest_exons.fasta
│
│   ├── mouse_IRES/
│   │   ├── IRESbase/
│   │   ├── IRESite/
│   │   ├── Rfam/
│   │   ├── best_mouse_homologs.tsv
│   │   ├── compiled_mouse_sequences.fasta
│   │   └── mouse_IRES_homologs.fasta

│   ├── virus_IRES/
│   │   ├── IRESbase/
│   │   ├── IRESite/
│   │   ├── Rfam/
│   │   └── compiled_virus_sequences.fasta

│   ├── literature/              # Literature-curated sequences and expression data
│
│   ├── compiled_human_mouse_virus_IRES.fasta
│
│   ├── clustering/
│   │   ├── clustered_sequences.fasta
│   │   ├── clusters.tsv
│   │   ├── cluster_alignments.tsv
│   │   ├── cluster_representatives.fasta
│   │   ├── cluster_summary.csv
│   │   └── cluster_summary_RFAM.csv
│
│   ├── filtering/
│   │   ├── filtered_IRES.tsv
│   │   ├── filtered_IRES_with_seq.tsv
│   │   └── filtered_IRES_with_expression.tsv
│
│   ├── conservation/
│   │   ├── human_mouse_conserved.tsv
│   │   ├── human_mouse_conserved_with_seq.tsv
│   │   ├── human_mouse_conserved_with_seq_exp.tsv
│   │   └── human_mouse_comparison.csv
│
│   ├── final_library/
│   │   ├── 200_IRES.tsv
│   │   ├── 200_IRES_with_seq.tsv
│   │   ├── 200_IRES_seq.fasta
│   │   └── IRES_cluster.xlsx
│
│   ├── CCDS.current.txt
│   └── CCDS_protein_exons.current.faa
│
└── notes/                       # Meeting notes and documentation
```

## Dependencies
External tools dependencies:
- [MMseqs2](github.com/soedinglab/MMseqs2/)
- [Infernal](http://eddylab.org/infernal/Userguide.pdf)
- [TCoffee](http://www.tcoffee.org/)

Python dependencies: (available in [requirements.txt](requirements.txt))
- beautifulsoup4==4.14.3
- biopython==1.86
- certifi==2026.1.4
- charset-normalizer==3.4.4
- contourpy==1.3.3
- cycler==0.12.1
- et_xmlfile==2.0.0
- fonttools==4.61.1
- greenlet==3.3.1
- idna==3.11
- kiwisolver==1.4.9
- lxml==6.0.2
- matplotlib==3.10.8
- numpy==2.4.2
- openpyxl==3.1.5
- packaging==26.0
- pandas==3.0.0
- pillow==12.1.1
- playwright==1.58.0
- pyee==13.0.0
- pyparsing==3.3.2
- python-dateutil==2.9.0.post0
- requests==2.32.5
- six==1.17.0
- soupsieve==2.8.3
- typing_extensions==4.15.0
- urllib3==2.6.3

## Aims
Bioinformatically create a library and to select which IRES elements to screen. 
Select a diverse set of IRES elements (length, conservation, structure, etc) from human, mouse, and viral sequences.

## Methods

### 1. Compile IRES sequences from public databases

Internal ribosome entry site (IRES) sequences were collected from three curated databases:

- IRESite (last update: 2019-03-18) — https://iresite.org/
- Rfam v15.1 (January 2026; 4,227 RNA families) — https://rfam.org/
- IRESbase — http://reprod.njmu.edu.cn/cgi-bin/iresbase/index.php

Sequences from these sources were downloaded, parsed, and converted into a standardized FASTA format.

---

### 2. Identify homologous IRES sequences between human and mouse

To identify conserved IRES elements between species:

- Human IRES sequences were mapped to the mouse genome using **BLAT**.
- Structural homology searches were performed using **Infernal (`cmsearch`)** with covariance models.
- Candidate matches were filtered to identify **human–mouse conserved IRES elements**.

---

### 3. Cluster IRES sequences by sequence similarity

Human, mouse, and viral IRES sequences were combined into a single dataset and clustered using **MMseqs2**.

Clustering grouped related IRES sequences based on sequence similarity, allowing redundant or highly similar sequences to be identified.

---

### 4. Classify IRES sequences by RNA structure

Structural classification was performed using covariance models from **Rfam**.

Each IRES sequence was scanned against the Rfam covariance model library to assign known RNA structural families when possible.

---

### 5. Collapse clusters into representative sequences

Cluster outputs were summarized and collapsed to generate **representative IRES sequences** for each cluster.  
This step reduced redundancy in the dataset while preserving sequence diversity.

---

### 6. Filter the IRES library

The candidate IRES library was filtered using several criteria:

- Sequence length
- Structural annotation (Rfam covariance models)
- Evolutionary conservation (human–mouse homologs)
- Expression evidence from transcriptomic datasets

This produced a refined set of **high-confidence IRES elements**.

---

### 7. Generate phylogenetic relationships

Representative IRES sequences were aligned using **T-Coffee** and used to generate a phylogenetic tree.

This analysis provides an overview of sequence relationships across human, mouse, and viral IRES elements.


## INFERNAL commands
```bash
cmcalibrate {file}.cm
cat {files}.cm >> compiled_covariance_model.cm
cmpress compiled_covariance_model.cm
cmscan compiled_covariance_model.cm mouse_upstream_seq.fasta -o data/mouse_IRES/mouse_IRES_cm_scan -tblout data/mouse_IRES/mouse_IRES_cm_scan_table
```

## References
Literature tables:
Supplementary Table 3, 4, 5, 7, 8 from Weingarten-Gabbay et al., 2016
Citation:

     Shira Weingarten-Gabbay et al.,Systematic discovery of cap-independent translation sequences in human and viral genomes.Science351,aad4939(2016).DOI:10.1126/science.aad4939

