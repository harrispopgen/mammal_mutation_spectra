# Mutyper pipeline: generate mutation spectra from all species

### Overview
These scripts were used to generate 7-mer mutation spectra and genomic target counts for every species in the dataset. As part of this pipeline, 5 random individuals are subset from each species/population in the dataset, and input vcf files are masked and filtered. 
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
The config files for each species contain the information used to run the mutyper pipeline as a job array on an SGE cluster enrivonment and contains paths to input files and species-specific parameters. If you were to use the config, you would need to update all input paths to be specific to your system possibly modify them to account for how job arrays are set up on your server. The files are commented inside, but in summary:
* `species`: species label
* `interval_or_chr_or_all`: specifies whether input genomes are split into intervals (for non-chromosomal assemblies, chromosomes, or have all autosomes combined)
* `prepend0`: TRUE/FALSE specifies if intervals have a pre-pended 0 (e.g. interval 01, 02, 03 etc)
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


