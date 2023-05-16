#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=polar_bear # changing to pb


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=interval  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=FALSE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=29 # pb mapped to pb genome has only 29 intervals ; bb has 31 ; before when I did them together they were both mapped to BB genome

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

########## input vcf files (should be polarized and have non -pol sites removed (except for humans/apes)) ############
vcfdir="/net/harris/vol1/home/beichman/bears/analyses/keightley_polarization/20220607_Keightley_focalpolar_bear_out1brown_bear_out2black_bear/vcfs_with_ancestral_alleles/perInterval/"
vcffilename="interval_${interval}.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX.SNPs.filt_variants_4b.LowCovIndsREMOVED.nomissfiltapplied.PASS.WARN.Only.CustomFilt.HardFilt.AnnotVar.TrimAlt.ref_polar_bear.int_${interval}.vcf.gz"
vcfNeedsToBeSubsetByChr=FALSE
# note: as of 20220608 I am no longer using the phased+parsimony polarized files ; instead am using keightley polarized unphased files
# note: these files still contain brown bear samples (and black bear?) so that will affect missing data (both bb and pb even tho pol separately will be called at same sites, I think that's good)
# but only pulling out the PB from this file 
############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/polar_bear_GCF_000687225.1/perInterval/interval${interval}.polar_bear_GCF_000687225.1.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed"
maskLabel=maskALL
### 20220608 made negative mask for pb (previously just had one for bb because was using pb mapped to bb) and split it into intervals
############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species)
ancestralFastaDir="/net/harris/vol1/home/beichman/bears/analyses/ancestralReferenceFastas/20220607_Keightley_focalpolar_bear_out1brown_bear_out2black_bear_assignedAllelesUsingProbs_POLARIZING_polar_bear/"
ancestralFastafilename="FINAL.MODIFIED.ANCESTRAL.interval_${interval}_polar_bearAncestral_extSFS_AllelesAssignedUsingProbs.BUGFIX.fasta"
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)
# 20220608 updated to use keightley polarization 

########### individuals to exclude ############
individualsToExclude='' # keep this empty if you don't wasnt to exclude any individuals

######## divide inds into pops : ##########
pops='PB' # just pb even though files contain others
poplistdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject/samples/populationFilesForMutyper
popfilesuffix=".txt" # format pop.txt

####### mutyper variants options ######
passOption=TRUE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales

strictOption=FALSE # options: TRUE OR FALSE. 
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans
