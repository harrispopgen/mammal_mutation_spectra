#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=Mus_spretus # changing this from mice to be mus musculus 


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=chr  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=FALSE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=19


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
vcfdir="/net/harris/vol1/home/beichman/mice/analyses/keightley_polarization/20220607_focalMs_out1Mmd_polarizing_Mus_spretus/vcfs_with_ancestral_alleles/perChromosome/"
vcffilename="chr${interval}.PolarizingSpretus.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BiallelicSNPsOnly.BUGFIX.AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz"
vcfNeedsToBeSubsetByChr=FALSE
# as of 20220429 fixed major bug that this vcf was refpolarized
# note on 20220524 I am updating it to be keightley est-sfs with anc alleles assigned by probability instead of hard cutoffs *important change!*
# note used a different config file for work on Tom's paper (in the 1mer section)
# as of 20220608 !! I am using mus spretus polarized as the focal group in extsfs intead of assuming it has same anc state as mmd 
############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/mouse_mm10/perInterval/chr${interval}.mouse_mm10.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed" 
maskLabel=maskALL
# spretus has same mask

############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species)
ancestralFastaDir="/net/harris/vol1/home/beichman/mice/analyses/ancestralReferenceFastas/20220607_focalMs_out1Mmd_polarizing_Mus_spretus_KimuraModel_AncAlleleAssignedUsingProbs_POLARIZED_MUS_SPRETUS/"
ancestralFastafilename="FINAL.MODIFIED.ANCESTRAL.chr${interval}_MusSpretusAncestral_extSFS_assignedUsingProbs.BUGFIX.fasta"
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)
# as of 20220429 this ancestral fasta has major bug fix (previously was ref-polarized)
# note on 20220524 I am updating it to be keightley est-sfs with anc alleles assigned by probability instead of hard cutoffs *important change!*
# note used a different config file for work on Tom's paper (in the 1mer section)
# as of 20220608 !! I am using mus spretus polarized as the focal group in extsfs intead of assuming it has same anc state as mmd 

########### individuals to exclude ############
individualsToExclude='' # keep this empty if you don't wasnt to exclude any individuals

######## divide inds into pops : ##########
pops='Ms' # note: now excluding all Mmx because am doing it with its own config file; this is just for Ms.
poplistdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject/samples/mouseSamples/populationFilesForMutyper
popfilesuffix=".txt" # format pop.txt


####### mutyper variants options ######
passOption=TRUE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales

strictOption=FALSE # options: TRUE OR FALSE. should only be true for humans.
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans
