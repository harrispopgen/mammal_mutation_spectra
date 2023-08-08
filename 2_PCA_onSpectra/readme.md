### Principal components analyses

These scripts were used to generate principal components analyses based on mutation spectra. 

##### Overview: 
The script `subsampled.PCA.allspecies_HumanRescalingUseThis.R` reads in per-individual 7-mer spectra and genomic target information (these files are output by the mutyper pipeline and can be found on Dryad repository in the mutyperOutput7mer.tar.gz), processes per-individual spectra (collapses 7-mer spectra down to  1-mer, 3-mer, 5-mer spectra, rescales to human genome content, downsamples, regularizes and CLR-transforms) then carries out principal components analysis. 

##### Primary input files:
* mutyper per-individual spectrum and target files from each species' directory within the mutyperOutput_7mer directory summed over all genomic intervals (Dryad; filenames inside script)

##### Additional files read in by script:: 
* `colors_and_labels.R`: contains colors and nice species labels for making manuscript-ready figures
* `revisions.confounders.SequencingPlatform.ReadLength.20230710.PERINDIVIDUALFORPCA.txt`: metadata file (contains info to color PCA by sequencer/read length

##### Figures generated: 
* Figure 2
* Supplemental PCA figures

##### Additional notes: 
You would need to update paths to input files on your own system inside script if re-running. Recommend running in RStudio chunk by chunk to see what's happening along the way.
