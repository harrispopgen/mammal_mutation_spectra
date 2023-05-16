#! /bin/bash
#$ -l h_rt=50:00:00,mfree=30G
#$ -o /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_spectrum
#$ -e /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_spectrum
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N sumup_spectra

###################### sum up spectrum ###################
module load modules modules-init modules-gs # initialize modules 
module load pcre2/10.35 hdf5/1.10.1 R/4.0.4
module load gcc/10.2.0 # necessary for reshape2


scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/multispecies_spectra/analyses/mutyper/unified_mutyper_scripts_multispecies_SubsettingIndividuals
script=$scriptdir/step_3b_sumupspectra_acrossintervals.R

configdir=$scriptdir/config_files_per_species


speciesList='humans Mus_spretus Mus_musculus fin_whale vaquita Gorilla_gorilla Pan_troglodytes Pan_paniscus Pongo_abelii Pongo_pygmaeus wolves brown_bear polar_bear ucla_wolves'

for species in $speciesList
do
echo "starting $species"

configfile=$configdir/config_${species}.sh
source $configfile # load species info l 
# these are coming from config file: $species $interval_count $prepend0; but pops are coming from the indir 

indir=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/${label}/${species}/mutyper_results_masked_${maskLabel}/mutyper_spectrum_files/
outdir=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/${label}/allspecies_summed_up_over_intervals_forTransfer/${species}/mutyper_results_masked_${maskLabel}/mutyper_spectrum_files
inputfilesuffix=".mutyper.spectrum.SeeLogForFilters.${maskLabel}.${kmersize}mer" # two files: followed by PERPOPULATION; PERINDIVIDUAL.RANDOMIZEDWITHINPOPULATION

mkdir -p $outdir

Rscript $script $indir $outdir $species $interval_count $prepend0 $inputfilesuffix


done
# a way to check that it added correctly : head -n2 *mutyper.spectrum.SeeLogForFilters.maskALL.7mer.PERPOPULATION.txt |grep -v "A" | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}' # aha this has issues because of mutyper not reporting 0s

#  (doing this across target files this should equal summed up AAAAAAA )

# 
