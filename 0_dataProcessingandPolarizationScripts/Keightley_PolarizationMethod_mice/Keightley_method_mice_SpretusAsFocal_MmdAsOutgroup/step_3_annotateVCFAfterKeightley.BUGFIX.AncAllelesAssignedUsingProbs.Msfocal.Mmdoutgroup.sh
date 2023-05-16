#! /bin/bash
#$ -l h_rt=50:00:00,h_data=10G
#$ -o /net/harris/vol1/home/beichman/mice/reports/extsfs
#$ -e /net/harris/vol1/home/beichman/mice/reports/extsfs
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N MSannotateAfterKeightley
#$ -t 1-19

set pipefail -euo
######### this script will (skipping xy and M and UN) split vcf by chr, annotate and tabix and make ancestral fasta based on ancAllele designation 
# so far only working with outside of 90-10 bound , skipping all other sites. 

### not running as an array for now because it runs quite quickly (only should take an hour or two)
module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12  ### NOTE vERSION MUST BE >1.11 or will not work correctly! 
focal=Ms
outgroup1=Mmd # note only one outgroup!

chromosome=chr${SGE_TASK_ID}


# config-file.txt: 
flagForThisRun="20220607_focal${focal%.txt}_out1${outgroup1%.txt}_polarizing_Mus_spretus"
projectdir=/net/harris/vol1/home/beichman/mice
wd=$projectdir/analyses/keightley_polarization/$flagForThisRun
outdir=$wd/vcfs_with_ancestral_alleles/perChromosome
mkdir -p $outdir
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
scriptdir=$gitdir/data_processing/variant_polarizing/KeightleyPolarizationMethod_mice/Keightley_method_mice_SpretusAsFocal_MmdAsOutgroup

######### have two vcfs I want to annotate. Want to annotate the full dog vcf and the wolf/coy only vcf #######

concatdir=$wd/concatted_extsfs_output_withAncAllele_KimuraModel_assignedUsingProbs
ancAllelesFile=$concatdir/ThreeColumnsOnly.AncAllelesAssignedUsingProbs.noneExcluded.KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.allIntervals.USETHIS.txt.gz

#### and updating to now have alleles assigned using probabilities ! (20220531)



###### BUG FIX!!! now this has 3 columns and will behave correctly with bcftools (!!!!)

colCount=`zcat $ancAllelesFile | head -n1 |  awk '{print NF}'`
if [ ${colCount} -ne 3 ]
then
	echo "THERE are less than or more than 3 columns in your anc alleles file! this will NOT WORK with bcftools annotate and may lead to subtle errors"
	exit 1
fi


# this annotates and adds the above h1.txt header to the header and outputs a bgzipped vcf file 
#mkdir -p $outdir/perChromosome
echo '##INFO=<ID=ancAllele,Number=1,Type=String,Description="Ancestral allele inferred by est-sfs (Keightley & Jackson, 2018). Alleles assigned using probabilities rather than cutoffs.">' > $outdir/h1.txt # make a header add on file for bcftools to add

# don't want h1.txt constantly overwritten by multiple jobs -- causes errors.


###################### vcf to annotate ###########
vcfdir=$projectdir/vcfs/
vcffilename=AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz 
chrspecificoutput=$outdir/${chromosome}.PolarizingSpretus.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BiallelicSNPsOnly.BUGFIX.${vcffilename}


### bcftoools will only annotate sites which are in the -a file so others won't get an annotation; don't want to include those sites in output
bcftools view -r $chromosome $vcfdir/$vcffilename | bcftools annotate -h $outdir/h1.txt -a ${ancAllelesFile} -c CHROM,POS,ancAllele | bcftools filter -e 'INFO/ancAllele="."' -Oz -o $chrspecificoutput  # this will e

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo "error in bcftools annotate"
    exit 1
fi



tabix -p vcf $chrspecificoutput


exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo "error in tabix"
    exit 1
fi

