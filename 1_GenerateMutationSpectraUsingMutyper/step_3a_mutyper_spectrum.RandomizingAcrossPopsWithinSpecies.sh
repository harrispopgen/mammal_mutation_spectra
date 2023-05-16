#! /bin/bash
#$ -l h_rt=10:00:00,mfree=2G
#$ -o /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_spectrum
#$ -e /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_spectrum
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -P sage

############## mutyper spectrum #################
# must be run after variants because you need the masked fasta file (can be run concurrently with targets/ksfs)

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)
module load bedtools/2.29.2

######### set up script ########

set -exou pipefail

configfile=$1

todaysdate=`date +%Y%m%d` # for the log file 

subsamplecount=5 

############### NOW SOURCE CONFIG FILE ONCE INTERVALS HAVE BEEN SET ###############
source $configfile
# note that this loads all the variables you need for the species, including kmer size etc ^^^


###### set up outdir/outfiles #######
wd=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/$label/$species/mutyper_results_masked_${maskLabel}
variantdir=$wd/mutyper_variant_files
spectrumdir=$wd/mutyper_spectrum_files
ksfsdir=$wd/mutyper_ksfs_files
targetdir=$wd/mutyper_target_files
subsampleindsdir=$wd/subsample_inds_list_files


mkdir -p $wd
mkdir -p $wd/logs
mkdir -p $variantdir
mkdir -p $spectrumdir
mkdir -p $ksfsdir
mkdir -p $targetdir


log=$wd/logs/${species}.${intervalLabel}.${todaysdate}.mutyper_spectrum.log
> $log

################### get mutyper variants file #############
# note this must match step_1 output exactly! trying to think of a way to make that more efficient.
mutypervariantsoutputname=${species}.int_or_chr_${interval}.mutyper.variants.SomeRevComped.SubsetOfIndividuals.SeeLogForFilters.${maskLabel}.${kmersize}mer.vcf.gz


######### iterate through populations if pops are provided #######
# if no pops provided, all inds are part of the population 
# if there aren't any pops defined then run without subsetting inds out: 

########## AM SUBSETTING TO 5 INDS *PRIOR* TO RUNNING MUTYPER SPECTRUM ###############


if [[ -z $pops ]] # if no pops defined
then
	echo "not splitting into pops" >> $log
	# per individual (for pca):
	subsetfile=$subsampleindsdir/${species}_AllPopsCombinedIfThereArePops.subsample.${subsamplecount}.diploids.txt #for spp without pops

	echo "individuals I'm selecting: " >> $log
	cat $subsetfile  >> $log
	echo "" >> $log
	# all individuals summarized spectrum: (still calling it PERPOPULATION for consistency)
	bcftools view -S $subsetfile $variantdir/$mutypervariantsoutputname -Ou |  bcftools view -c 1:minor -Ou |  mutyper spectra --population - > $spectrumdir/${species}.int_or_chr_${interval}.mutyper.spectrum.SeeLogForFilters.${maskLabel}.${kmersize}mer.SUBSETPRIORTOMUTYPERSPECTRUM.PERPOPULATION.txt
	
	exitVal=$?
	if [ ${exitVal} -ne 0 ]; then
		echo "error in mutyper spectrum"
		exit 1
	else
		echo "done with pop level spectrum"
	fi

	
	bcftools view -S $subsetfile $variantdir/$mutypervariantsoutputname -Ou |  bcftools view -c 1:minor -Ou | mutyper spectra --randomize - > $spectrumdir/${species}.int_or_chr_${interval}.mutyper.spectrum.SeeLogForFilters.${maskLabel}.${kmersize}mer.SUBSETPRIORTOMUTYPERSPECTRUM.PERINDIVIDUAL.RANDOMIZEDWITHINPOPULATION.txt

	
	exitVal=$?
	if [ ${exitVal} -ne 0 ]; then
		echo "error in mutyper spectrum"
		exit 1
	else
		echo "finished per individual spectrum"
	fi

else
	echo "splitting into pops: $pops" >> $log
	for pop in $pops
	do
		echo -e "starting $pop"
		# pop-specific subset file
		subsetfile=$subsampleindsdir/${species}_${pop}.subsample.${subsamplecount}.diploids.txt # for spp with pops when you want to separate the pops for spectra calculation so that you can remove fixed sites. 
		# note that not randomizing across pops currently so when comparing species shared variation would be double counted for each population.
		
		echo "individuals I'm selecting: " >> $log
		cat $subsetfile  >> $log
		echo "" >> $log
		spectrumpopdir=$spectrumdir/$pop
		mkdir -p $spectrumpopdir

		######## NOTE: it is *extremely* important for mutyper that the AN and AC fields are correct. BCFTOOLS correctly updates them through filtering (I checked). ####
		# select individuals in file | exclude fixed sites within that population | mutyper variants --> 
		# NOTE this REMOVES in fixed sites per population (sites fixed across all inds have already been removed)

		# per population:
		bcftools view -S $subsetfile $variantdir/$mutypervariantsoutputname -Ou |  bcftools view -c 1:minor -Ou  | mutyper spectra --population - > $spectrumpopdir/${species}.${pop}.int_or_chr_${interval}.mutyper.spectrum.SeeLogForFilters.${maskLabel}.${kmersize}mer.SUBSETPRIORTOMUTYPERSPECTRUM.PERPOPULATION.txt
		
		exitVal=$?
		if [ ${exitVal} -ne 0 ]; then
			echo "error in mutyper spectrum"
			exit 1
		else
			echo "done with per pop spectrum"
		fi
		
		# per individual (for pca):
		bcftools view -S $subsetfile $variantdir/$mutypervariantsoutputname |  bcftools view -c 1:minor  | mutyper spectra --randomize - > $spectrumpopdir/${species}.${pop}.int_or_chr_${interval}.mutyper.spectrum.SeeLogForFilters.${maskLabel}.${kmersize}mer.SUBSETPRIORTOMUTYPERSPECTRUM.PERINDIVIDUAL.RANDOMIZEDWITHINPOPULATION.txt

		exitVal=$?
		if [ ${exitVal} -ne 0 ]; then
			echo "error in mutyper spectrum"
			exit 1
		else
			echo "finished script"
		fi
	done
			
fi

# note this is generating per-individual and per-population based on the same set of 5 individual per pop. So the first column of the per pop spectrum is the sum of the 5 rows of the per individual spectrum (I checked; it is)
