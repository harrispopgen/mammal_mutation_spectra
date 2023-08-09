### Phylogenetic signal of technical and biological confounders

The script `runMantelTestOnConfounders.R` was used to calculate whether numerous technical and biological variables (contig N50, scaffold N50, average sequencing coverage, Watterson's theta, age at first reproduction, reproductive lifespan) show any significant 'phylogenetic signal' (the correlation between pairwise differences in the 'trait' and the square root of phylogenetic distance between species). See SI Methods for details.

###### input files:
* `confoundersToPlotPhyloSignal.GenomeN50.etc.updatedthetavals.updatedCoverage.txt`: table of technical and biological 'trait' values for each species (provided here)
* `RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick`: RAXML tree from Upham et al. (2019) - this contains all the species from Upham et al. and gets subset to the species in our study during this script. Branch lengths represent expected substitutions (in phylogenetic_trees.tar.gz directory on Dryad)
