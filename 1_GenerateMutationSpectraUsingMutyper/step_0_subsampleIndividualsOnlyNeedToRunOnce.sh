
##################### can be run in the shell; only need to run once #############
module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)
module load bedtools/2.29.2
############## this is designed so that it won't overwrite previous selections if you rerun it! this is important so doesn't change analyses
# but wil change if you do a new mutyper run with diff masks 


speciesList='humans Mus_spretus Mus_musculus fin_whale vaquita Gorilla_gorilla Pan_troglodytes Pan_paniscus Pongo_abelii Pongo_pygmaeus wolves brown_bear polar_bear ucla_wolves'


subsamplecount=5 # diploid inds

scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/multispecies_spectra/analyses/mutyper/unified_mutyper_scripts_multispecies_SubsettingIndividuals

configdir=$scriptdir/config_files_per_species


################## loop through species subsampling individuals ###############

for species in $speciesList
do
echo starting $species
configfile=$configdir/config_${species}.sh
source $configfile # load species info 



###### set up outdir/outfiles #######
wd=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/$label/$species/mutyper_results_masked_${maskLabel}
mkdir -p $wd 
mkdir -p $wd/logs

######## set up log file #########
log=$wd/logs/step_0_${species}.selectInds.log
> $log

#variantdir=$wd/mutyper_variant_files
subsampleindsdir=$wd/subsample_inds_list_files


mkdir -p $subsampleindsdir
########## check if a subsample file already exists for the species; if it does: DONT OVERWRITE IT #########

if [[ -z $pops ]]
then
	echo "no pop file, subsetting from vcf directly"
	outfile=$subsampleindsdir/${species}_AllPopsCombinedIfThereArePops.subsample.${subsamplecount}.diploids.txt

# if no pops defined then need to pull out list of individuals and subsample 5 at random using shuf
# just use first interval to get species list
# doesn't work for vaq

# 20221114: using original vcf instead of mutyper variants file because I want to subset individuals prior to running mutyper variants.
# from 
	echo "file used for sampling: $vcfdir/$vcffilename" # these are loaded from config file 

# check for individuals to exclude from the list (only need to do this for the species that don't already have pop files)
	if [ -z $individualsToExclude ]
	then
		echo "no individuals to remove" >> $log
		rm_inds_snippet=''
	else
		echo "removing the following individuals : $individualsToExclude " >> $log
		rm_inds_snippet="-s ^$individualsToExclude" # changing to double quotes (worked before though; checked logs)
	fi

	# need to exclude bad individuals prior to subsetting
	# note; want to select 5 per pop, but want to remove missing data across all pop samples so they have same target areas. so am going to write my selections out to one central file as well. 

	# exclude bad individuals if there are any using rm_inds_snippet, then get the resulting list of individuals, and write it out
	bcftools view ${rm_inds_snippet} $vcfdir/$vcffilename -Ou | bcftools query -l | shuf | head -n ${subsamplecount} > $outfile 
	### need to exclude low quality individuals from being chosen using the rm inds snippet now that I'm doing this before mutyper variants. 
	exitVal=$?
	echo "exit val: $exitVal"
	if [ ${exitVal} -ne 0 ]; then
		echo "error in mutyper variants"
	exit 1
		else
	echo "finished $species"
fi
elif [[ -n $pops ]] # if pops exists: SUBSET POP TEXT FILES: NOTE THESE popfiles don't contain any individuals that need to be removed so they don't need the rm_inds_snippet from above .
then
	echo "drawing individuals from pop files that have already had low qual individuals excluded (if any)" >> $log
	outfile_allpopstogether=$subsampleindsdir/${species}_AllPopsCombinedIfThereArePops.subsample.${subsamplecount}.diploids.txt
	> $outfile_allpopstogether # create empty file that will have all pops together 
	for pop in $pops
	do
	echo -e "starting $pop"
	popfile=$poplistdir/${pop}${popfilesuffix}  # from config file
	perpop_outfile=$subsampleindsdir/${species}_${pop}.subsample.${subsamplecount}.diploids.txt

	cat $popfile | shuf | head -n 5 > ${perpop_outfile} # make a pop-specific outfile
	# but also add it to the running file that has all the pops:
	cat ${perpop_outfile} >> ${outfile_allpopstogether}

	done

	
fi
done
