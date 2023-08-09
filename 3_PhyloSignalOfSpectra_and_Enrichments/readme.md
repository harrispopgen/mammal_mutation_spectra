# Phylogenetic signal and k-mer enrichments:

### Overview:
This directory of scripts was used to generate a large number of analyses for the paper, including Figures 3, 4, and 5 (and numerous supplemental figures). 

#### Step 1: Phylogenetic signal script

`step_1_GetDistances.usingProjectedSpectra.MultinomialDownsample_ALT_RescaleToHumanTargets_testedCLRIssues.NowUsingNonProjectedSpectra.USETHIS.R`: this long and heavily commented script was used to generate a large number of analyses and many intermediate plots, along with main text and SI plots related to phylogenetic signal that were presented in the manuscript (Figure 3, Figure 4, and spectrum-distance related SI figures). Analyses carried out using this script include:

* it processes mutation spectra (collapsing the 7-mer spectrum down to the 1-, 3-, and 5-mer spectra, carrying out genome target correction, multinomial downsampling, regularization and clr-transformation). 
* it outputs mutation and target counts that are ready for use in the `step_2_onSAGE_FishersExactTest.R` script
* it calculates the phylogenetic (cophenetic) distance between pairs of species based on two different phylogenetic trees (substitution-based and ultrametric), and the Aitchison distance (and cosine similarity) between each pair of species' spectra
* it uses the Mantel test to test for phylogenetic signal between the square root of phylogenetic distance and spectrum distance, including a number of variations:
	* results with and without taking the square root of phylogenetic distance
	* results based on the folded mutation spectrum
	* results based on sub-spectra (subsetting 3-, 5- and 7-mer mutation spectra by central 1-mer type or by biased gene conversion mutation categories 
	* results based on using (1-cosine similarity) as the distance metric instead of Aitchison distance
* It calculates the correlation between mutation spectrum distances and a series of technical and biological confounders, using regular and phylogenetically-aware versions of the Mantel test
* It generates 'fake' randomly shuffled 5- and 7-mer mutation spectra (randomizing counts within central 3-mer and 5-mer categories (see SI Methods for details) to determine whether the empirical 5- and 7-mer spectra contain any phylogenetic signal beyond what is contained in lower-order spectra

The script is best run in RStudio step by step to see what it's doing along the way. 

###### input files:
* `allSpectra.PopulationLevel.BasedOn5individualsPerPop.notprojected.usethis.txt`: file containing species-level 7-mer mutation spectra (on Dryad; inside the mutyperOutput_7mer.tar.gz directory. This has the 7-mer spectra counts output from mutyper, concatenated across all species and sub-populations. Includes populations and species that get excluded by this script as they are not part of main text results (including ucla_wolves and sub-populations that weren't part of main text analyses))
* genomic targets files: these contain the species-specific genomic target counts for each 7-mer mutation type (on Dryad in species-specific `mutyper_targets` directories in `mutyperOutput7mer.tar.gz`)
* `speciesCodes.txt`: file with species full Latin names for matching with phylogenetic trees and shorthand labels for plotting (provided here)
* `confoundersToPlotPhyloSignal.GenomeN50.etc.updatedthetavals.updatedCoverage`: confounders table. Quantitative variables associated with each species that could be technical or biological confounders (contig and scaffold N50, average coverage, age at first reproduction, etc.) (provided here)
* `RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick`: RAXML tree from Upham et al. (2019) - this contains all the species from Upham et al. and gets subset to the species in our study during this script. Branch lengths represent expected substitutions (in phylogenetic_trees.tar.gz directory on Dryad)
* `listOfSpecies.ForTimeTree.pluswolves.20220204.nwk`: TimeTree ultrametric phylogenetic tree from TimeTree.org with branch lengths scaled in millions of years (in phylogenetic_trees.tar.gz directory on Dryad)
* `colors_and_labels.R`: file containing plotting colors and labels (provided here)

###### Alternative versions of the script are included which have the following variations (alternateVersions/ directory):
* `step_1_GetDistances.usingProjectedSpectra.MultinomialDownsample_ALT_RescaleToHumanTargets_testedCLRIssues.NowUsingNonProjectedSpectra.ExclVaqPBGOC.R`: excludes the lowest diversity species/populations (vaquita, polar bear excluded, and low-diversity Gulf of California (GOC) fin whales swapped for higher diversity Eastern North Pacific (ENP) fin whales) to explore the impacts of data sparsity on results
* `step_1_GetDistances.usingProjectedSpectra.MultinomialDownsample_ALT_RescaleToHumanTargets_testedCLRIssues.NowUsingNonProjectedSpectra.SepCpGs1merOnly.CpTTpGOnly.R`: re-analyzes 1-mer phylogenetic signal with CpG>TpG mutations as a seventh mutation category, separated from other C>T mutations
* `step_1_GetDistances.usingProjectedSpectra.MultinomialDownsample_ALT_RescaleToHumanTargets_testedCLRIssues.NowUsingNonProjectedSpectra.SepCpGsEXCLUDE.CpTTpGOnly.excluded.R`: re-analyzes 1-mer phylogenetic signal with CpG>TpG mutations excluded
* `step_1_GetDistances.usingProjectedSpectra.MultinomialDownsample_ALT_RescaleToHumanTargets_testedCLRIssues.NowUsingNonProjectedSpectra.ILR.13merOnly.R`: reanalyzes 1-mer and 3-mer results using the isometric log ratio (ILR) transformation instead of the centered log ratio transformation 


#### Step 2: Enrichments scripts
* `step_2_onSAGE_FishersExactTest.R` and `step_2_onSAGE_FishersExactTest.Wrapper.sh`: This script took in the files output in step 1 that contain uncorrected and untransformed 3-mer, 5-mer and 7-mer spectra and target information and carried out Fisher's exact test for every k-mer of the 3-mer, 5-mer and 7-mer spectra of every species/population to determine whether the rate of that k-mer is enriched or depleted relative to the background 1-mer rate (e.g. if the rate of ATAAA>ATGAA mutations is higher or lower than the overall A>G mutation rate). Note the script is best run on a remote cluster as it is quite slow and memory intensive. The "wrapper" script was used to submit it to the cluster.
* `step_2b_RemakeEnrichmentPlots_ifneeded.R`: a script to process the results of the Fisher's exact test described above, generate summary statistics, and make the plots shown in Figure 5 and enrichment-related SI figures. The results of Fisher's exact test (`step 2a`) that are used as input to step 2b are on Dryad in the `FishersExactTestResults_Enrichments.tar.gz` directory.

#### Additional script: Kmult analysis 
* `Kmult.only.R`: a script to calculate the Kmult statistic for phylogenetic signal
