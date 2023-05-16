#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=vaquita


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=allautos # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=FALSE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=NA # not needed because all together in one file 


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
vcfdir="/net/harris/vol1/home/beichman/vaquita/vcfs/SNPsOnly/polarized/20210216"
vcffilename="vaquita_20_simple_PASS_autos_variants_outgroupAlleles_ancAllele.vcf.gz" # note this is autosomes only (so don't need to exclude chrX,Y etc.)
vcfNeedsToBeSubsetByChr=FALSE


############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/vaquita_mPhoSin1/vaquita_mPhoSin1.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed"
# note for vaquita it's the whole genome mask
maskLabel=maskALL

############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species)
ancestralFastaDir="/net/harris/vol1/home/beichman/vaquita/analyses/ancestralReferenceFasta/"
ancestralFastafilename="MODIFIED.ANCESTRAL.vaquitaAncestral.autosomesOnly.chr1-21.fasta" # updated on 20220507 to just have autosomes chr 1-21
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)

########### individuals to exclude ############
individualsToExclude='z0001663,z0004380,z0004393,z0004394,z0185383' # keep this empty if you don't wasnt to exclude any individuals
# for vaquita, these are relatives of inds left in the dataset

######## divide inds into pops : ##########
pops='' # have it blank if no pops to split into and don't need to define other terms



####### mutyper variants options ######
passOption=TRUE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales

strictOption=FALSE # options: TRUE OR FALSE. (doesn't matter for vaquita since is all upper case, but keep as false) 
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans
