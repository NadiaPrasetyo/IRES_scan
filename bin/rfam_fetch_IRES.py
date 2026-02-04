
search_url = "https://www.ebi.ac.uk/ebisearch/ws/rest/rfam"
fields = "id,gca_accession,rfamseq_acc,name,common_name,rna_type,description,num_seed,length,seq_start,seq_end,ENA"
org = "virus"  # or "human" or "mouse"
query='{search_url}?query=rna_type:"IRES" AND {org}&fields={fields}'

fetch_url = "https://rfam.org/family"
sequence_alignment_url = "{fetch_url}/{id}/alignment/stockholm"
fasta_url = "{fetch_url}/{id}/alignment/fastau" #ungapped fasta
structure_mapping_url = "{fetch_url}/{id}/structures?content-type=application/json"
sequence_region_url = "{fetch_url}/{acc}/regions" # plain text
covariation_url = "{fetch_url}/{id}/cm"
secondary_structure_url = "{fetch_url}/{id}/image/cons"
structure_mapping = "{fetch_url}/{id}/structures?content-type=application/json"