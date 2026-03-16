Done:
- check the Science paper (https://www-science-org.ezproxy.otago.ac.nz/doi/10.1126/science.aad4939) to see rank of expression for all the different IRESes (174nt oligonucleotides)
- trim down too small dubious sequences (length < 100nt)
- Prioritize the inclusion of mouse-human conserved IRES
- prioritize the inclusion of diverse virus sourced IRES
- Create a phylogenetic tree of the 200 IRES sequences
- Send the list of 200 IRES to Lydia to get tested



To do:
- get the conservation score (length normalized conservation score) to the final TSV


- Is there a way to evaluate the sequences they've selected? diversity of the sequences
  - aligning the sequences, identifying a non-synonymous variation (protein level)
- Machine learning paper: influenza vaccine: coverage score (how well does the sequence represent the circulating virus)
- Run immunogenictiy score or immune response score (multiple tools and compare scores)

- send email to fiona, nikki, ries to get comments and feedback on the manuscript and ask if its worth publishing. CC paul



To talk about:
- doing an honors project: can I do lab work as well as bioinformatical work AND lab work?
- is there a possibility that I can do the honors project in Auckland? is there any co-supervisor that may be suitable?: Auckland Uni Co-supervisor from PHF
- Options: Olin Silander (problem: Justin O'sulivan (boss)), Austin ganley, Ant poole (may have to pay bench fee)
- Scholarship options: Keep you as a part time (ARF) + full time honors, 
- Projects options:
  - look at ways to assess sequence randomization tool (quite aligned with Austin ganley random chromosome project): wet lab to synthesize and see activity in saccharomyces - Null models for genomes => large language models for sequence randomization (genome wide)
    => can I give a sequence and ask it to generate a good control sequence. 

  - BUSCO: assessing how complete a genome is - for a given clade, there are a whole bunch of protein coding genes that are single copies that are found only in the clade (e.g. animal). When a new animal genome comes in, you can use BUSCO to see how well it matches the clade and how much it's coverage. 
    - but there is no non-coding RNA check: and extra quality control - taking RGAM models, tRNA scans, running those on the animal, bacteria, plants, fungi, etc.

Possibly work with Malaghan: run something past Gary, William (Waikato), 