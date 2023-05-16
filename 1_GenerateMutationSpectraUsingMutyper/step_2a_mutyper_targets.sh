#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_targets
#$ -e /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_targets
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -P sage
############## mutyper targets ###########
# must be run after variants because you need the masked fasta file (so can't just copy over files from previous run)

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)
module load bedtools/2.29.2

######### set up script ########

set -exou pipefail

configfile=$1

todaysdate=`date +%Y%m%d` # for the log file 



############### NOW SOURCE CONFIG FILE ONCE INTERVALS HAVE BEEN SET ###############
source $configfile
# note that this loads all the variables you need for the species, including kmer size etc ^^^


###### set up outdir/outfiles #######
wd=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/$label/$species/mutyper_results_masked_${maskLabel}
variantdir=$wd/mutyper_variant_files
spectrumdir=$wd/mutyper_spectrum_files
ksfsdir=$wd/mutyper_ksfs_files
targetdir=$wd/mutyper_target_files
maskedfastadir=$wd/masked_ancestral_fasta_files

#mkdir -p $wd
#mkdir -p $wd/logs
#mkdir -p $variantdir
#mkdir -p $spectrumdir
#mkdir -p $ksfsdir
#mkdir -p $targetdir


log=$wd/logs/${species}.${intervalLabel}.${todaysdate}.mutyper_targets.log
> $log
echo "ancestral fasta: $ancestralFastafilename (note am using negative masked version of fasta)" >> $log
echo "" >> $log
echo "negative mask: $NEGATIVEMASK" >> $log
echo "" >> $log
echo "################## COPY OF CONTIG FILE: ##########################" >> $log
echo "" >> $log
cat $configfile >> $log
echo "#########################################################" >> $log

############## masked ancestral fasta (from step 1) ###########
maskedancestralfasta=$maskedfastadir/${ancestralFastafilename%.fasta}.HARDMASKED.Ns.${maskLabel}.fasta


#################### using --strict ######################
# need strict for humans only
if [ $strictOption = "TRUE" ]
then
	echo "using --strict for $species - are you sure you want to do this? should be humans only" >> $log
	strict_snippet='--strict'
elif [ $strictOption = "FALSE" ]
then
	echo "not using --strict" >> $log
	strict_snippet=''
else
	echo "NOT A VALID strictOption option"
	exit 1
fi

######## mutyper targets: input the masked fasta ###########
targetsoutfilename=${species}.int_or_chr_${interval}.mutyper.targets.SeeLogForFilters.${maskLabel}.${kmersize}mer.txt

mutyper targets ${strict_snippet} --chrom_pos $chrom_pos --k $kmersize $maskedancestralfasta > $targetdir/$targetsoutfilename

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper targets"
	exit 1
else
	echo "finished"
fi


######## note about non-autosome chroms/scaffs not included in genotype calling: made sue they aren't in these anc genomes : #######
# bears: ancestral fasta already in intervals so doesn't contain scaffs <1mb (done)
# fin whale: already in intervals (done)
# human ancestor: already in chromosomes (done)
# apes ancestors: already in chromosomes (done)
# vaquita : aha! this one has non autosomes. have updated anc fasta to just be autosomes 1-21. note for past target calc I was using a positve bed mask so ti was fine. but now need ot deal with it. restrict to 1-21
# mice : already in chromosomes (done)


# note one issue with targets: dont have any sort of callability filter for the all-called sites, so targets is an overestimate. but proportions should be about right.
# more of an issue for mushi when rates matter. less of an issue for spectra.



