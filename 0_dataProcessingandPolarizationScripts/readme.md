These scripts were used to run `est-sfs` (https://sourceforge.net/projects/est-usfs/) to probabilistically assign ancestral allele states based on outgroup species sequences mapped to the same reference genome and allele frequency in the focal species. Extensive details about this procedure are in the paper's SI Methods. 

The final annotated vcfs and ancestral allele fastas output by this process are on the paper's Dryad repository inside each species' mutyper_pipeline_input_[species] directory.

Each species (Mus musculus, Mus spretus, brown bear, polar bear, gray wolf, fin whale) has its own directory containing its est-sfs config file, and species-specific scripts. Note that vaquita, humans and non-human great apes were polarized prior to this study (see SI Methods).

Code was modified from scripts provided by Jacqueline Robinson.

The general overview of the process is:
* *step 1*:  convert vcf files to tables ready for input into est-sfs. This step used each species' `step_1` wrapper scripts and `Prepare.KeightlyExtSFSInfile.py` python script for species, or `Prepare.KeightlyExtSFSInfile.ONEOUTGROUPONLY.py` for species that only had one outgroup species.

* *step 2*: run `est-sfs` and probabilistically choose the ancestral allele state for each site using a binomial draw, based on the est-sfs output probabilities for each site. Note this probabilistic assignment approach means that the assignment of a given site may not always be the most probable ancestral allele, but overall the assignments should average out to generate a more accurate mutation spectrum or site frequency spectrum. This step used each species' `step_2` wrapper scripts around `est-sfs` and the `assignAncestralAllele.Probabilistically.py` python scripts. 

* *step 3*: annotate the original vcf with the ancestral states inferred in the previous steps (`step_3` scripts).

* *step 4*: generate ancestral fasta files for input into mutyper by changing the reference allele state to the inferred ancestral state. Note that only sites in the vcf file will be changed, so these ancestral fasta files should not be used with any other datasets from the same species that may contain different SNPs. Also note as stated above that the ancestral allele of a given site may not be the most probable, but overall the assignments should yield more accurate mutation spectra because the ancestral states were assigned probabilistically. This step used each species' `step_4` wrapper scripts around the `generateAncestralFastaFromPolarizedVCF.ancAlleleVersion.py` python scripts.



