#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=fin_whale


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=interval  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=TRUE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=96


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
label=${kmersize}_mer # can add other notes here 

########## input vcf files (should be polarized and have non -pol sites removed (except for humans/apes)) ############
vcfdir="/net/harris/vol1/home/beichman/fin_whale/analyses/keightley_polarization/20211020_focalENP_out1MegNov_out2BalMus/vcfs_with_ancestral_alleles/perInterval/"
vcffilename="interval_${interval}.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.PASSSITESPlusWarnCpGRefSitesONLY.lowqualremoved.BUGFIX.JointCalls_f50b4_08_B_VariantFiltration_${interval}.vcf.gz"
vcfNeedsToBeSubsetByChr=FALSE
# note major bug fix on 20220429; previously was ref pol ; now using "BUGFIX" files
# note on 20220524 I am updating it to be keightley est-sfs with anc alleles assigned by probability instead of hard cutoffs *important change!*

############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/minke_whale_GCF_000493695.1_BalAcu1.0/perInterval/interval${interval}.minke_whale_GCF_000493695.1_BalAcu1.0.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed"
maskLabel=maskALL


############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species)
ancestralFastaDir="/net/harris/vol1/home/beichman/fin_whale/analyses/ancestralReferenceFastas/20211020_focalENP_out1MegNov_out2BalMus_assignedAllelesUsingProbs/"
ancestralFastafilename="FINAL.MODIFIED.ANCESTRAL.interval_${interval}_finWhale_Ancestral_extSFS_AllelesAssignedUsingProbs_NOTEthisisMinkeWhaleRefGenome.butAncCallsAreforfinWhale.includesCpGandRepeatSites.BUGFIX.fasta"
# note on 20220524 I am updating it to be keightley est-sfs with anc alleles assigned by probability instead of hard cutoffs *important change!*

chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)
# note major bug fixed on 20220429 ; previously was ref pol
########### individuals to exclude ############
individualsToExclude='BalAcu02,BalMus01,MegNov01,EubGla01,ENPOR12,ENPAK28' # keep this empty if you don't wasnt to exclude any individuals

######## divide inds into pops : ##########
pops='ENP GOC' 
poplistdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/fin_whale_spectrum/samples/
popfilesuffix=".txt" # format pop.txt


####### mutyper variants options ######
passOption=FALSE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales where I want to keep in the sites previously filtered as CPGs and Repeats (already removed all otehr bad sites)
strictOption=FALSE # options: TRUE OR FALSE. 
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans
