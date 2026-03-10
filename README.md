# IRES_scan

## Dependencies
External tools dependencies:
- [MMseqs2](github.com/soedinglab/MMseqs2/)
- [Infernal](http://eddylab.org/infernal/Userguide.pdf)
- [TCoffee](http://www.tcoffee.org/)

Python dependencies:
Package            Version
------------------ -----------
beautifulsoup4     4.14.3
biopython          1.86
certifi            2026.1.4
charset-normalizer 3.4.4
contourpy          1.3.3
cycler             0.12.1
et_xmlfile         2.0.0
fonttools          4.61.1
greenlet           3.3.1
idna               3.11
kiwisolver         1.4.9
lxml               6.0.2
matplotlib         3.10.8
numpy              2.4.2
openpyxl           3.1.5
packaging          26.0
pandas             3.0.0
pillow             12.1.1
pip                26.0.1
playwright         1.58.0
pyee               13.0.0
pyparsing          3.3.2
python-dateutil    2.9.0.post0
requests           2.32.5
six                1.17.0
soupsieve          2.8.3
typing_extensions  4.15.0
urllib3            2.6.3


## Aims
Bioinformatically create a library and to select which IRES elements to screen. 
Select a diverse set of IRES elements (length, conservation, structure, etc) from human, mouse, and viral sequences.

# Methods

1. Source IRES from databases:
- [iressite (Last change to the database: 2019-03-18 09:32:49 GMT+1)](https://iresite.org/)
- [Rfam 15.1 (January 2026, 4227 families)](https://rfam.org/)
- [iresbase](http://reprod.njmu.edu.cn/cgi-bin/iresbase/index.php)

2. Find homologous IRES sequences in human and mice: for all human IRES seq, try to find it in the mice genome using blat and infernal cmsearch to find human-mouse conserved IRES

3. Run MMSeqs2 Cluster (seq based) across human, mice, and viral IRES sequences to cluster sequences based on sequence similarity
   
4. Run covariance model in RFAM to classify and categorize the IRES sequences by structure

5. Collapse and prune the essential table into representative IRES sequences

6. Filter the IRES library based on length, structure, conservation, and expression

7. Generate a phylogenetic tree for the IRES library using TCoffee






INFERNAL command 
cmcalibrate {file}.cm
cat {files}.cm >> compiled_covariance_model.cm
cmpress compiled_covariance_model.cm
cmscan compiled_covariance_model.cm mouse_upstream_seq.fasta -o data/mouse_IRES/mouse_IRES_cm_scan -tblout data/mouse_IRES/mouse_IRES_cm_scan_table



Literature tables:

     Shira Weingarten-Gabbay et al.,Systematic discovery of cap-independent translation sequences in human and viral genomes.Science351,aad4939(2016).DOI:10.1126/science.aad4939
Supplementary Table 3, 4, 5, 7, 8