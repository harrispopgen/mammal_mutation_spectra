# Evolution of the mutation spectrum across a mammalian phylogeny

This GitHub repository contains the scripts that were used to generate the analyses and figures in Beichman, Robinson, Lin, Nigenda-Morales, Moreno Estrada, & Harris. 
"Evolution of the mutation spectrum across a mammalian phylogeny," _Molecular Biology & Evolution (in revision)_

<img src="(https://github.com/harrispopgen/mammal_mutation_spectra/assets/10618123/9d66d27c-c352-479d-86be-3ce1f8ef7bfa)" width="100" height="100">

### Overview
The main directories are as follows, and has its own detailed `readme` file:

* **`0_dataProcessingandPolarizationScripts/`**: scripts used to assign ancestral allele states for species that did not already have ancestral polarization
* **`1_GenerateMutationSpectraUsingMutyper/`**: scripts used to generate mutation spectra for every species in the dataset
* **`2_PCA_onSpectra/`**: scripts used to generate principal component analyses (PCA)
* **`3_PhyloSignalOfSpectra_and_Enrichments/`**: scripts used to carry out analysis of phylogenetic signal of the mutation spectrum, correlation of the spectra with possible confouners, and enrichments/depletions of particular k-mers
* **`4_PhyloSignalOfConfounders/`**: scripts used to measure the phylogenetic signal of technical and biological confounders
* **`5_MutationSignatureFitting_inSigfit`/**: scripts used to carry out mutational signature fitting analyses
* **`6_ComparingToMouseDNMs/`**: scripts used to carry out additional comparisons of mouse-wolf 1-mer spectra similarities using two alternative datasets


Each directory has its own detailed `readme` file containing script-specific information (input and output files, overviews of functions, meanings of parameters, etc.), and scripts are heavily commented. There are extensive details on the analyses in the **SI Methods** section of Beichman et al.

### Code used to generate main-text figures in the paper:
* **Figure 2** (Principal component analyses): https://github.com/harrispopgen/mammal_mutation_spectra/tree/main/2_PCA_onSpectra
* **Figures 3, 4 and 5** (Aitchison distances between spectra, Mantel test results, enrichments/depletions of particular motifs): https://github.com/harrispopgen/mammal_mutation_spectra/tree/main/3_PhyloSignalOfSpectra_and_Enrichments
* **Figures 6 and 7** (mutational signature fitting): https://github.com/harrispopgen/mammal_mutation_spectra/tree/main/5_MutationSignatureFitting_inSigfit

### Data files:

Data files from the paper are on Dryad (**DRYAD LINK**) with an extensive `readme` describing them. Smaller auxiliary files have been placed in each script's GitHub directory for convenience. 


### Additional notes: 
You would need to update paths to input files on your own system inside script if re-running scripts. We recommend that R scripts be run in RStudio step-by-step to see what each code chunk is doing. All code is presented as-is (exactly as it was used to generate the analyses of the paper), so there may be extraneous comments and analyses present within them.
