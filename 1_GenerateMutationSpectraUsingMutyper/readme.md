# Mutyper pipeline: generate mutation spectra from all species

### Overview
These scripts were used to generate 7-mer mutation spectra and genomic target counts for every species in the dataset. As part of this pipeline, 5 random individuals were subset from each species/population in the dataset, and input vcf files were masked and filtered, including masking out exonic regions +- 10kb, RepeatMasker-annotated regions, CpG Islands, and low-complexity regions. 

The 7-mer spectra were subsequently collapsed down into 1-mer, 3-mer and 5-mer spectra.

The scripts are heavily commented and full details of the pipeline are in the SI Methods of the paper. A summary of each step of the pipeline is below.

Each species has a species-specific `config` file (in the `config_files_per_species/` dir, and also on Dryad) that was used to provide paths to their input files and species-specific parameters (detailed below). 

Some species had sufficient individuals across multiple populations that we subset the species by population to generate per-population spectra. One population was subsequently used in main text analyses. 

The **input** files used to run this pipeline, including config files, vcf files, negative genome masks, and ancestral fasta files are on Dryad in each species' `mutyper_pipeline_input_[species].tar.gz` directory. 
  * Note that the mouse (Mus musuclus domesticus and Mus spretus) and fin whale Dryad directory tarballs have been split into parts due to size in order to upload to Dryad (`*tar.gz.PART.aa`, `*tar.gz.PART.ab`, `.ac`, etc). Before you can un-tarball them you must cat them back together ( for example: `cat mutyper_pipeline_input_files_fin_whale.tar.gz.PART.* > mutyper_pipeline_input_files_fin_whale.tar.gz.joined` )

  * Note that the human and ape species vcfs are not on Dryad, but links to where they can be downloaded are provided in their Dryad directories and in the SI Methods of the paper.

The **output** of the pipeline for all species is in the `mutyperOutput_7mer.tar.gz` directory on Dryad.


*! Note that if this pipeline's `step 0` (detailed below) is run, it will draw a random 5 individuals per population/species that may be different from the ones we drew at random for our paper. To replicate our results, see our Table S1 for the particular individuals we drew to replicate our results, or use the subsample_inds_list_files/ directories for each species (on Dryad in the `mutyperOutput_7mer.tar.gz` directory), or instead of re-running it, work directly with the output of this pipeline that is on Dryad (`mutyperOutput_7mer.tar.gz`).*


### Config files:
The config files for each species contain the information used to run the mutyper pipeline as a job array on an SGE cluster environment and contains paths to input files and species-specific parameters. If you were to use the config, you would need to update all input paths to be specific to your system possibly modify them to account for how job arrays are set up on your server. The files are commented inside, but in summary:
* `species`: species label
* `interval_or_chr_or_all`: specifies whether input genomes are split into intervals (for non-chromosomal assemblies, chromosomes, or have all autosomes combined)
* `prepend0`: TRUE/FALSE specifies if intervals have a pre-pended 0 (e.g. interval 01, 02, 03 etc.)
* `interval_count`: the number of intervals or chromosomes (e.g. 22 for human, 98 for fin whale, NA if all autosomes are combined, etc.). The config file will put the interval being worked on from the $SGE_TASK_ID in the job array (set up by wrapper scripts of pipeline)
* `kmersize`: the size of kmer you want to generate. We generated 7mer spectra then subset them down to 1,3, and 5-mer.

* `vcfdir`: path to input vcf files
* `vcffilename`: name of input vcf files
* `vcfNeedsToBeSubsetByChr`: TRUE/FALSE: whether or not you need to subset the vcf file by chromosome 
* `NEGATIVEMASK`: path to negative genome mask .bed file
* `maskLabel`: custom label to keep track of what mask you are using 
* `ancestralFastaDir`: path to ancestral fasta file 
* `ancestralFastafilename`: name of ancestral fasta file
* `chrom_pos`: 0-based position of chr name in fasta
* `individualsToExclude`: list of any individuals you want to exclude (due to low quality)
* `pops`: list of populations if you want to split your species into populations
* `poplistdir`: for species you want to subset by population, path to population list files. These can be found on Dryad in the population_list_files.tar.gz directory.
* `popfilesuffix`: ".txt" # suffix of population list files (e.g. .txt)

* `passOption`: whether you want to only select SNPs that have “PASS” in the info field
* `strictOption`: whether you want to treat lower-case bases in the ancestral fasta files as * missing data (only applies to human ancestral fasta files)

### Step 0: randomly subsample 5 individuals per species/population
The script `step_0_subsampleIndividualsOnlyNeedToRunOnce.sh` randomly draws 5 individuals from each species (or from each population within a species if a species is being subset by population). It will also exclude any low-quality individuals specified for removal. It gets all this information from each species' `config` file. If you don't want to rerun this, you can make sure the subsequent steps of the pipeline can fine the subsample_inds_list_files/ directory for each species inside the output directory to use those already-sampled individuals instead. 

### Step 1: mutyper variants
The script `step_1_mutyperVariants.SubsettingIndividuals.unified.sh` and its wrapper that submits the job arrays to the cluster `step_1_mutyperVariants.SubsettingIndividuals.unified.WRAPPER.sh` use the info in the `config` files and the results of `step 0` to subset the input vcf file to just the individuals chosen in `step 0`, filter and mask the vcf files depending on species-specific requirements in the `config` file, and run `mutyper variants` to determine the sequence context of every site in the filtered vcf file based on the provided ancestral fasta file. The output is a greatly reduced vcf file, with ancestral state added to the INFO field, ready for use by `mutyper spectra`.

The script also generates a masked ancestral fasta file for use in `mutyper targets`, which converts the regions in the negative mask file to Ns so that they cannot contribute to the spectrum or to the targets (in `step 2`).

Note that `mutyper variants` uses an ancestral fasta file to assign ancestral and derived allele states (so the 'REF' and 'ALT' columns of the vcf in fact correspond to 'ANCESTRAL' and 'DERIVED' mutation states). The software also reverse-complements mutations (so a G -> A is coded as C -> T). This results in the REF and ALT columns and corresponding 0s and 1s in the genotype fields being different from the input vcf at some sites.

### Step 2: mutyper targets
The purpose of mutyper targets is to count up the genomic targets in the masked ancestral reference file to determine how many of each ancestral motif exists (e.g. the number of AAAATAA 7-mers in the masked ancestral reference genome of the fin whale). This is an important step, because we will eventually need to correct each species' spectra to have the same genomic content based on these counts.

The ancestral fasta file was masked in `step 1` above, so that the regions we excluded from the spectra (exonic regions, repeat masker annotated regions, etc.) were also not included in the target counts.

The process occurs in two sub-steps:
* *Step 2a* : running `mutyper targets` using `step_2a_mutyper_targets.sh` and its job-array wrapper `step_2a_mutyper_targets.WRAPPER.sh`
* *Step 2b* : summing up the target counts across all genomic intervals to get total masked-genome counts using `step_2b_sumuptargets.R` and its wrapper `step_2b_sumuptargets.WRAPPER.sh`


### Step 3: mutyper spectra
In this step, the vcfs output by `step 1: mutyper variants` will be converted into tab-delimited mutation spectra, both at the per-individual level (one spectrum for each of the 5 sub-sampled individual per species/population), and at the species/population level (one spectrum for the whole species or population). 

A crucial step in generating the per-individual spectra is the `--randomize` parameter, which randomly assigns a shared polymorphism to a single individual that carries it, so that similarities in spectra are not driven by shared variation within the species or population. When generating population-level spectra, the equivalent is that each polymorphism is only counted toward the spectrum counts once, even if it appears in multiple individuals. The population or species-level spectra are therefore simply the sum of the individual-level spectra.

The individual-level spectra are used for PCA analyses, and the species-level spectra for pairwise distance calculations between species.

The process occurs in three sub-steps:
* *Step 3a* : running `mutyper spectra` using `step_3a_mutyper_spectrum.RandomizingAcrossPopsWithinSpecies.sh` and its job-array wrapper `step_3a_mutyper_spectrum.RandomizingAcrossPopsWithinSpecies.WRAPPER.sh`
* *Step 3b* : summing up the spectra counts across all genomic intervals to get total masked-genome spectra counts using `step_3b_sumupspectra_acrossintervals.R` and its wrapper `step_3b_sumupspectra_acrossintervals.WRAPPER.sh`
* *Step 3c* : concatenating the spectra across all species/populations to make one master file `allSpectra.PopulationLevel.BasedOn5individualsPerPop.notprojected.usethis.txt` (which is in Dryad in the mutyperOutput_7mer.tar.gz directory)





