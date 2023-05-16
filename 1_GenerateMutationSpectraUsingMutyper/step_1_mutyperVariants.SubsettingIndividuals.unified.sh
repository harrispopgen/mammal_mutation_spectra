#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_variants
#$ -e /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/mutyper_variants
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -P sage

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)
module load bedtools/2.29.2

######## generic script #########
# need either/or for a few things
# removing individuals due to quality or other issues (gorilla, chimp, vaquita, fin whale)
# using strict or not (humans only)
# whether to filter on pass or not (fin whales dont just get pass sites)

# for apes, need to go to the interval in question using -R because vcf files aren't per chromosome. 

set -exou pipefail

configfile=$1

todaysdate=`date +%Y%m%d` # for the log file 



############### NOW SOURCE CONFIG FILE ONCE INTERVALS HAVE BEEN SET ###############
source $configfile

###### set up outdir/outfiles #######
wd=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/$label/$species/mutyper_results_masked_${maskLabel}
variantdir=$wd/mutyper_variant_files
spectrumdir=$wd/mutyper_spectrum_files
#ksfsdir=$wd/mutyper_ksfs_files # no longer making the ksfs 
targetdir=$wd/mutyper_target_files
maskedfastadir=$wd/masked_ancestral_fasta_files
subsampleindsdir=$wd/subsample_inds_list_files # where your subset files are for the 5 individuals per population

# from step 0:
subsamplecount=5 
subsetfile_allpopstogether=$subsampleindsdir/${species}_AllPopsCombinedIfThereArePops.subsample.${subsamplecount}.diploids.txt
# note you want to remove missing data across all these individuals so that if you do cross pop comparisons it's all the same target area

mkdir -p $wd
mkdir -p $wd/logs
mkdir -p $variantdir
mkdir -p $spectrumdir
#mkdir -p $ksfsdir
mkdir -p $targetdir
mkdir -p $maskedfastadir

######## set up log file #########
log=$wd/logs/${species}.${intervalLabel}.${todaysdate}.mutyper_variants.log
> $log

echo "vcffile: $vcfdir/$vcffilename" >> $log
echo "" >> $log
echo "ancestral fasta: $ancestralFastaDir/$ancestralFastafilename (note that I used the masked version for mutyper variants and targets)" >> $log
echo "" >> $log
echo "negative mask: $NEGATIVEMASK" >> $log
echo "" >> $log
echo "################## COPY OF CONFIG FILE: ##########################" >> $log
echo "" >> $log
cat $configfile >> $log
echo "#########################################################" >> $log

######### set up output name ###########
mutypervariantsoutputname=${species}.int_or_chr_${interval}.mutyper.variants.SomeRevComped.SubsetOfIndividuals.SeeLogForFilters.${maskLabel}.${kmersize}mer.vcf.gz

################# restricting to only PASS sites ##################
if [ $passOption = "TRUE" ]
then
	echo "select only PASS sites" >> $log
	pass_snippet='-f PASS' 
elif [ $passOption = "FALSE" ]
then
	echo "not restricting to PASS sites" >> $log
	pass_snippet=''
else
	echo "NOT A VALID passOption option"
	exit 1
fi

#################### using --strict ######################
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

############ if you need to subset the input vcf by chromosome before processing (apes need this) ########
if [ $vcfNeedsToBeSubsetByChr = "TRUE" ]
then
	echo "need to subset the vcf by chromosome/interval first (needed for apes)" >> $log
	subset_vcf_snippet='-r $intervalLabel' # note use -r not -R because it's text not a file
elif [ $vcfNeedsToBeSubsetByChr = "FALSE" ]
then
	echo "don't need to subset vcf by chr/interval" >> $log
	subset_vcf_snippet=''
else
	echo "not a valid vcfNeedsToBeSubsetByChr option"
fi


####################### first, mask the ancestral fasta file ##############
# I'm doing this out of an abundance of caution to avoid 7mers that extend into a masked region that shouldn't get called #############
#################### mask the ancestral fasta using the negative mask file *regions you DONT want * ###########
#### note that these fastas don't contain any non-autosome chroms (made sure for each species. but you should make sure before running)
# note this is a HARD MASK (NNNNN) (soft mask wouldn't work without --strict in mutyper variants).

maskedancestralfasta=$maskedfastadir/${ancestralFastafilename%.fasta}.HARDMASKED.Ns.${maskLabel}.fasta # set name of fasta that has been masked wtih whatever negative mask you're using

# negative hard-mask your fasta file with NNNNNNs
bedtools maskfasta -fi $ancestralFastaDir/$ancestralFastafilename -bed $NEGATIVEMASK -fo $maskedancestralfasta

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in bedtools maskfasta"
	exit 1
else
	echo "done with bedtools maskfasta"
fi

############ build mutyper variants code #########

# code snippets: 
# MUST use double quotes! 
initialize_subsetifneeded_snippet="bcftools view $subset_vcf_snippet $vcfdir/$vcffilename -Ou" # initialize with this. if no subsetting needed it'll just start reading the vcf. saves issues with needing to supply vcf at diff parts of pipe
# this will select your individuals and impose filters (mask, biallelic snps only)
# this filter snippet is now SELECTING individuals instead of excluding them (20221117)
filter_snippet="bcftools view -S $subsetfile_allpopstogether -T ^$NEGATIVEMASK -m2 -M2 -v snps $pass_snippet -Ou" # if passOption=False then $pass_snippet will be ''; if no inds to remove it will be blank
# taking out $rm_inds_snippet and replacing with selectin the individuals -s $subsetfile_allpopstogether (20221114). note using capital S because it's a file rather than -s which is for a list.
no_fixed_sites_snippet="bcftools view -c 1:minor -Ou" # this will removed sites that are fixed 0/0 or 1/1 across your subset individuals across all populations (will still have fixed sites within populations)
missing_data_snippet="bcftools view -g ^miss -Ou" # removes missing data (sites where at least one individual has ./. genotype -- note that this results in all individuals across your populations within a species having the same target regions.)
# note: am using masked ancestral fasta here just in case 7mer extended sequence context overlaps with a masked region (dont want to include)
mutyper_variants_snippet="mutyper variants --k $kmersize --chrom_pos $chrom_pos $strict_snippet $maskedancestralfasta -  | bcftools convert -Oz -o $variantdir/${mutypervariantsoutputname}" # if strictoption=FALSE then it will be '' and not used

######### need to subset vcf prior to processing? ############
# requires double quotes (single don't work)
# always start with bcftools view then add in other snippets ; if you don't need to pre-subset the vcf by chr then $subset_vcf_snippet will be empty and it'll just read the vcf

lineOfCode="${initialize_subsetifneeded_snippet} | ${filter_snippet} | ${no_fixed_sites_snippet} | ${missing_data_snippet} | ${mutyper_variants_snippet}" 
echo "FINAL CODE LINE: " >> $log
echo "$lineOfCode" >> $log


# run it: (need to use eval otherwise bash splits the code strangely)
eval $lineOfCode

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished script"
fi

