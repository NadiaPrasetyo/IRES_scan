Done:
- fetch from IRES databases: 
  - [iresbase](http://reprod.njmu.edu.cn/cgi-bin/iresbase/index.php): contains 554 viral IRESs, 691 human IRESs, and 83 IRESs from other eukaryotic species. 
  - [iressite (Last change to the database: 2019-03-18 09:32:49 GMT+1)](https://iresite.org/): contains 55 viral IRESs, 44 human IRESs, and 17 mouse/rat IRESs.
  - [Rfam 15.1 (January 2026, 4227 families)](https://rfam.org/) 

- Cluster databases IRES based on sequence:
  - mmseqs2 Cluster

- Have a meeting with Lydia to figure out exactly what she needs and how I can help bioinformatically for her circRNA IRES project

- swinging back and forth on what to do with the manuscript draft


To do:
MAYBE:
- maybe: FIONA: interesting observation - SSL7 definitely the most immunogenic of the three SSLs in both mice and also in human donors. Haven't looked at IgE responses in either setting if that is what this test is predicting. We use a mutated form of the protein in the vaccine - would be great to know whether these predicted allergenic peptides have been engineered out as part of this process (unlikely) and which peptides should be engineered out to reduce the likelihood of unwanted allergenicity.
  - probably a lot we can do about that, it would be useful for her

IRES:
- Try to map the human IRES to the mouse genome and pull out the homologous sequences in mice and get the downstream protein: human-mouse conservation 
- run rfam models across the IRES that I've pulled from the other databases: more sensitive
  - INFERNAL (HMMER for RNA structure): looking for conserved secondary structure and sequence - can run the Rfam covariance models on the sequence database (that you want to search)
  - make a file that contains all the cm for the IRES concatenated together and run the model across
- Cluster by sequence: see how that goes
- rank based on length and give them maybe 10 to screen 


SSL Alignment 
- deep alignment: build some trees with what I've done and highlight the ones they chose and tell them whether they have used reasonable ones or redundant
- highlight their representative strains - based on


SET up regular meetings (monthly) with Fiona and friends (SSL)
Thursday regular meetings with Paul mornings 