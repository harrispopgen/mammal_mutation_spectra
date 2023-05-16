######### Wrapper script to submit mutyper variants #############
scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/multispecies_spectra/analyses/mutyper/unified_mutyper_scripts_multispecies_SubsettingIndividuals
configdir=$scriptdir/config_files_per_species

script=$scriptdir/step_1_mutyperVariants.SubsettingIndividuals.unified.sh

speciesList='humans Mus_spretus Mus_musculus fin_whale vaquita Gorilla_gorilla Pan_troglodytes Pan_paniscus Pongo_abelii Pongo_pygmaeus wolves brown_bear polar_bear ucla_wolves'

# make outdir:
mkdir -p /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_variants

for species in $speciesList
do

configfile=$configdir/config_${species}.sh
source $configfile # load species info 


if [ $interval_or_chr_or_all = "allautos" ]
then
	# if is all autosomes you don't need intervals: 
	qsub -N ${species}_variants $script $configfile

else
	qsub -N ${species}_variants -t 1-${interval_count} $script $configfile
fi

done 
