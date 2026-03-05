Done:
- Run covariance model in RFAM to classify and categorize the IRES sequences by structure
- Find Human-mouse: conserved IRES
- Run MMSeqs2 Cluster (seq based) across human, mice, and viral IRES sequences
- OKey doke. This may be the best you can do. IRES' seem to be ridiculously diverse.
- So, one measure that might be handy is the similarity between the human and mouse sequences (& length) -- not sure what to do with the viral ones, unless there are equivalent mouse viruses.
- Well some of the viruses aren't even human viruses, should we prioritise human or mouse virus IRESs?
- So our priority list at the moment:
    - length (the shorter the better)
    - human-mouse conservation
    - diversity between human, mouse, and viral IRES
    - viral IRES from virus that infect human, mice, or closely related eukarotes
    - IRES are also classified into a range of groups based on structure, so we can also select IRES that have different structural groups (type I-IV in virus and type 1 and 2 cellular IRES) - they also have different methods of recruiting and interacting with the ribosomes



To do:
- check the Science paper (https://www-science-org.ezproxy.otago.ac.nz/doi/10.1126/science.aad4939) to see rank of expression for all the different IRESes (174nt oligonucleotides)
- trim down too small dubious sequences (length < 100nt)
  