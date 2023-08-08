### Principal components analyses

These scripts were used to generate principal components analyses based on mutation spectra, including plots in Figure 2 and supplemental PCA figures. 

##### Overview: 
The script `subsampled.PCA.allspecies_HumanRescalingUseThis.R` reads in per-individual 7-mer spectra and genomic target information (these files are output by the mutyper pipeline and can be found on Dryad repository in the mutyperOutput7mer.tar.gz), processes per-individual spectra (collapses 7-mer spectra down to  1-mer, 3-mer, 5-mer spectra, rescales to human genome content, downsamples, regularizes and CLR-transforms) then carries out principal components analysis. 

###### input files:
* mutyper per-individual spectrum and target files summed over all genomic intervals from each species' directory in the `mutyperOutput_7mer.tar.gz` directory on Dryad (specific filenames inside script)
* `colors_and_labels.R`: contains colors and nice species labels for making manuscript-ready figures (provided here)
* `revisions.confounders.SequencingPlatform.ReadLength.20230710.PERINDIVIDUALFORPCA.txt`: metadata file (contains info to color PCA by sequencer/read length) (provided here)

