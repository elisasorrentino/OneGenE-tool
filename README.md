# OneGenE-tool
OneGenE it is part of the gene@home project, a BOINC-based project that expand Gene Regulatory Network based on public gene expression data. OneGenE pre-compute the Gene Network Expansion of every gene of the organism, using the NES2RA algorithm, producing their expansion lists. The aim of this tool is to reconstruct the expansion of a network from the expansion lists.

Given as input a Gene Regulatory Network (GRN) that we want to expand, the tool takes in input the expansion lists and aggregates them with the Borda Count aggregator. The output will be a list of genes that may interact with the GRN of input and the p-value of the interactions.

This tool is based on Vitis Vinifera expression data coming from the Vespucci compendia ([http://vespucci.colombos.fmach.it/](http://vespucci.colombos.fmach.it)).

To have more info about the gene@home project see [https://gene.disi.unitn.it/test/genehome/en/hf/index.php](https://gene.disi.unitn.it/test/genehome/en/hf/index.php)

To have more info of some previouse experiments on this organism with NES2RA algorithm see [Malacarne et al. 2018, Discovering Causal Relationships in Grapevine Expression Data to Expand Gene Networks. A Case Study: Four Networks Related to Climate Change](https://doi.org/10.3389/fpls.2018.01385)
