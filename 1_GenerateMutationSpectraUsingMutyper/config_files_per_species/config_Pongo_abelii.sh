#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=Pongo_abelii


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=chr  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=FALSE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=22


########## set interval value from SGE_TASK_ID (*code is the same for all species*) ########
if [[ -z ${SGE_TASK_ID}  || ${SGE_TASK_ID} = "undefined" || $interval_or_chr_or_all = "allautos" ]] # if there are no sge task id and/or all autos 
then
	interval='allautos' # no intervals 
	intervalLabel="allautos"
elif [[ -n ${SGE_TASK_ID}  && ${SGE_TASK_ID} != "undefined" &&  $interval_or_chr_or_all = "interval" ]]
then
	# prepend 01 02 if needed 
	if [ $prepend0 = "TRUE" ]
	then
		interval=`printf %02d ${SGE_TASK_ID}` # prepend 0 
	elif [ $prepend0 = "FALSE" ]
	then
		interval=${SGE_TASK_ID} # don't prepend 0 
	else
		echo "invalid prepend0 option"
		exit 1
	fi
	# set interval label: 
	intervalLabel=interval_${interval}
# if there's a task id and it's chromosomes: 
elif [[ -n ${SGE_TASK_ID}  && ${SGE_TASK_ID} != "undefined" && $interval_or_chr_or_all = "chr" ]]
then
	interval=${SGE_TASK_ID}
	intervalLabel=chr${interval} # label for output files 

fi


######### kmer size  ########

kmersize=7
label=${kmersize}_mer


########## input vcf files (should be polarized and have non -pol sites removed (except for humans/apes)) ############
vcfdir="/net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs"
vcffilename="Pongo_abelii.vcf.gz" # note this whole genome; but Michael's anc fasta are per chrom so i am subsetting this vcf into chrs when it's processed
vcfNeedsToBeSubsetByChr=TRUE # because it needs to be subset during processing at mutyper variants step (apes only)


############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/apes_mapped_to_hg18/perInterval/chr${interval}.apes_mapped_to_hg18.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed" # ONE MASK FILE PER INTERVAL
# note: apes cannot use same negative mask as humans because are mapped to diff ref genomes (humans: hg38; apes: hg18 old), so have made a different hg18 specific mask for the apes
maskLabel=maskALL


############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species)
ancestralFastaDir="/net/harris/vol1/home/beichman/apes/polarized_ref_fastas/hg18_references_with_gagp_ancestral_alleles_exclude_recurrent/"
ancestralFastafilename="Pongo_abelii_chr${interval}.fa" # made by Michael
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)

########### individuals to exclude ############
individualsToExclude=''

######## divide inds into pops : ##########
pops='' # have it blank if no pops to split into and don't need to define other terms



####### mutyper variants options ######
passOption=TRUE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales

strictOption=FALSE # options: TRUE OR FALSE. ; should only be true for humans
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans

