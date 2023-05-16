#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=ucla_wolves


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
vcfdir="/net/harris/vol1/home/beichman/canids_fromEduardo/analyses/keightley_polarization/20220718_focalucla_wolves_out1coyote/vcfs_with_ancestral_alleles/wolves_only/" # this is now the version with anc alleles assigned using keightley est for ucla wolves specifically
vcffilename="AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX.excludingCoyote.BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing3_5_v4_mergeGaTKfiltered_varnonvar_All15AW_jointcalled_allchr.vcf.gz" # this has same name as before but has been updated to be genome wide 
vcfNeedsToBeSubsetByChr=TRUE # it's all autosomes (and no x/y/mt) and doesn't need to be subset (some species do because of mismatch between anc fasta and vcf being subset)
# these are wolves from eduardo amorim (UCLA) and kirk lohmueller; they have been filtered 
# now these have been directly polarized using cal_coy coyote from the broad dataset and keightley sfs! Use this instead!! 

############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/dog_canFam3/perInterval/${intervalLabel}.dog_canFam3.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed"
maskLabel=maskALL
# this is per-interval for canids
# same mask for broad wolves and ucla wolves
############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species). note that non-interval 
ancestralFastaDir="/net/harris/vol1/home/beichman/canids_fromEduardo/analyses/ancestralReferenceFastas/20220718_focalucla_wolves_out1coyote_ancAllelesAssignedUsingProbs/per_chr_ancestral_fastas_autosomesonly/" 
# note this should be the per_chr_ancestral_fastas_autosomesonly dir! rather than whole genome
ancestralFastafilename="FINAL.MODIFIED.ANCESTRAL.canFam3.0_extSFS_AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX_${intervalLabel}.fasta"
# name  is the same between ucla and broad wolves but enclosing dir (above) is different, be careful about that
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)
# **** important; on 20220531 switching to vcf that has anc alleles assigned using probabilities from extsfs, isntead of using hard cutoffs *** (is in same dir; deleted previous vcfs)
# 20220429 fixed major bug -- this was ref polarized and not ancestrally polarized. fixed now (note "BUGFIX" in filename)
# 20220718 now adding ucla -wolf specific ancestral fasta! used cal coy from broad as outgroup, but will keep in rarer sites that may only be in the ucla dataset as a snp but could be monomorphic in broad wolves
# makes it a more independent replicate compared to restricting to snps in the broad wolf dataset, which is great. 
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
