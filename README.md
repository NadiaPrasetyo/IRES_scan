# IRES_scan


External tools dependencies:
- [MMseqs2](github.com/soedinglab/MMseqs2/)
- [Infernal](http://eddylab.org/infernal/Userguide.pdf)




INFERNAL command 
cmcalibrate {file}.cm
cat {files}.cm >> compiled_covariance_model.cm
cmpress compiled_covariance_model.cm
cmscan compiled_covariance_model.cm mouse_upstream_seq.fasta -o data/mouse_IRES/mouse_IRES_cm_scan -tblout data/mouse_IRES/mouse_IRES_cm_scan_table