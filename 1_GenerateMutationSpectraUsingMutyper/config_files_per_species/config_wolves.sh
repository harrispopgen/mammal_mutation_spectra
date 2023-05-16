#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=wolves


##### flags for dealing with genome splitting by chr or interval ########
# be careful that non-autosomes are not included!
interval_or_chr_or_all=chr  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=FALSE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=38 # parallelize 
########## set interval value from SGE_TASK_ID (*code is the same for all species*) ########
# if there is a SGE task id then run this code: 
# -n tests for sge task id being set
# if there is task id and it's intervals:
# if task id is set as 'undefined' or doesn't exist or all autos is specified: 
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

########## input vcf files (should be polarized and have non -pol sites removed (except for humans/apes) AND NO NON-AUTOSOMES!!!) ############
vcfdir="/net/harris/vol1/home/beichman/dogs/analyses/keightley_polarization/20220421_focalwolf_out1coyote/vcfs_with_ancestral_alleles/wolves_only"
vcffilename="AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX.excludingCoyote.wolf_coyote_notallwolves.biallelicsnpsonly.keepingmonomorphicbecausepolyindogs.vcf.gz"
vcfNeedsToBeSubsetByChr=TRUE # it's all autosomes (and no x/y/mt) and doesn't need to be subset (some species do because of mismatch between anc fasta and vcf being subset)
# **** important; on 20220531 switching to vcf that has anc alleles assigned using probabilities from extsfs, isntead of using hard cutoffs *** (is in same dir; deleted previous vcfs)
# choosing to subset to parallelize; need to subset anc fasta too
# note: this vcf is all chromosomes, but I am subsetting int he script. ancestral fasta and mask are per chromosome. 
# fixing bug on 20220429 -- previously was ref polarized.(not BUGFIX in middle of name)
############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/dog_canFam3/perInterval/${intervalLabel}.dog_canFam3.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed"
maskLabel=maskALL
# this is per-interval for canids
############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species). note that non-interval 
ancestralFastaDir="/net/harris/vol1/home/beichman/dogs/analyses/ancestralReferenceFastas/20220421_focalwolf_out1coyote_ancAllelesAssignedUsingProbs/per_chr_ancestral_fastas_autosomesonly"
ancestralFastafilename="FINAL.MODIFIED.ANCESTRAL.canFam3.0_extSFS_AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX_${intervalLabel}.fasta"
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)
# **** important; on 20220531 switching to vcf that has anc alleles assigned using probabilities from extsfs, isntead of using hard cutoffs *** (is in same dir; deleted previous vcfs)
# 20220429 fixed major bug -- this was ref polarized and not ancestrally polarized. fixed now (note "BUGFIX" in filename)
########### individuals to exclude ############
individualsToExclude='' # keep this empty if you don't wasnt to exclude any individuals
# note I already removed the coyote
######## divide inds into pops : ##########
pops=''
poplistdir=''
popfilesuffix='' # format pop.txt

####### mutyper variants options ######
passOption=TRUE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales

strictOption=FALSE # options: TRUE OR FALSE. 
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans
