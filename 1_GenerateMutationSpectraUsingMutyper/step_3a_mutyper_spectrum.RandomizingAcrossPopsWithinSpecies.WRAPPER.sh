
######### Wrapper script to submit mutyper spectrum #############
# must be run after step_1 because you need the masked ancestral fasta


scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/multispecies_spectra/analyses/mutyper/unified_mutyper_scripts_multispecies_SubsettingIndividuals
configdir=$scriptdir/config_files_per_species



script=$scriptdir/step_3a_mutyper_spectrum.RandomizingAcrossPopsWithinSpecies.sh
speciesList='humans Mus_spretus Mus_musculus fin_whale vaquita Gorilla_gorilla Pan_troglodytes Pan_paniscus Pongo_abelii Pongo_pygmaeus wolves brown_bear polar_bear ucla_wolves'

mkdir -p /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_spectrum

#speciesList="bears vaquita"
for species in $speciesList
do

configfile=$configdir/config_${species}.sh
source $configfile # load species info 


#################  NEED TO PRE-PICK 5 inds prior to randomize/pop spectra ####################
# otherwise for humans too few variants per individual
# so going to pick 5 inds per pop for each species (step 0)



# for vaquita need to submit with no intervals 
if [ $interval_or_chr_or_all = "allautos" ]
then
	
	# if is all autosomes you don't need intervals: 
	qsub -N ${species}_spectrum -P sage $script $configfile
	

else
	qsub -N ${species}_spectrum -P sage -t 1-${interval_count} $script $configfile ### -t 1-${interval_count}
fi

done 

