### Further comparisons of mouse and wolf datasets 
This script (`LindsaySpectrumComparison.mice.PlusTargetCorrection.PlusUCLAWolves.HUMANRESCALING.NewBasedNotOnProjectionAnymore.USETHIS.R`) was used to compare the Aitchison distances between 1-mer mutation spectra of mouse and wolves using alternate datasets to ensure that mouse-wolf similarities we observed were not a data artifact.

The alternative wolf polymorphism dataset (generated at UCLA, therefore called ucla_wolves) was from:

* Mooney, Jazlyn A., et al. "Long-term small population size, deleterious variation, and altitude adaptation in the ethiopian wolf, a severely endangered canid." Molecular Biology and Evolution 40.1 (2023): msac277.

and the alternative mouse dataset made up of de novo mutations (DNMs) was from:

* Lindsay, Sarah J., et al. "Similarities and differences in patterns of germline mutation between mice and humans." Nature communications 10.1 (2019): 4053

The script is heaviliy commented and details of analysis are in **SI Methods**.



##### input files used:

* `all1merSpectra_plusTargets.notyethumanrescaled.ForLindsayAnalysis.containsUCLAwolves.UsedSubsetIndividuals.notprojection.txt`: 1-mer spectra, including ucla_wolves (provided here). Note these spectra have not yet undergone any genomic target rescaling (happens in the script). 
* `41467_2019_12023_MOESM5_ESM.txt`: Mouse DNM file from Lindsay et al. (2019) (https://www.nature.com/articles/s41467-019-12023-w) (provided here for convenience)
* `Mus_musculus.summedup.mutyper.targets.SeeLogForFilters.maskALL.7mer.txt`: mouse genomic targets (on Dryad; also provided here for convenience)
* `colors_and_labels.plusUCLAwolves.R`: a file containing plotting colors and labels 







