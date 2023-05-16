#! /bin/bash
#$ -l h_rt=05:00:00,h_data=4G
#$ -o /net/harris/vol1/home/beichman/bears/reports/extsfs
#$ -e /net/harris/vol1/home/beichman/bears/reports/extsfs
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N annotateAfterKeightley
#$ -t 1-19
# eventually 1-19


######### this script will act per chromosome 1-19 (skipping xy and M and UN) to split vcf by chr, annotate and tabix and make ancestral fasta based on ancAllele designation 
# so far only working with outside of 90-10 bound , skipping all other sites. 

### not running as an array for now because it runs quite quickly (only should take an hour or two)
module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for annotating vcf

chromosome=chr${SGE_TASK_ID}

focal=Mmd
outgroup1=Mmm
outgroup2=Ms

# config-file.txt: 
flagForThisRun="20210225_focal${focal}_out1${outgroup1}_out2${outgroup2}"
wd=/net/harris/vol1/home/beichman/mice/analyses/keightley_polarization/${flagForThisRun}
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
scriptdir=$gitdir/data_processing/variant_polarizing/KeightleyPolarizationMethod_mice
projectdir=/net/harris/vol1/home/beichman/mice

vcfdir=$projectdir/vcfs/
vcffilename=AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz 

concatdir=$wd/concatted_extsfs_output_withAncAllele_KimuraModel_assignedUsingProbs

outdir=$wd/vcfs_with_ancestral_alleles_KimuraModel_assignedUsingProbs

# WITH MORE THAN 3 COLUMSN BCFTOOLS USES REF INSTEAD OF ANC ! oh no! ancAllelesFile=$concatdir/ConfidentAncAllelesSelectedWithAwk.90-10Range.extsfs.output.KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}.pvalues.WITHPOSITIONS.usethis.txt.gz
ancAllelesFile=$concatdir/ThreeColumnsOnly.AncAllelesAssignedUsingProbs.noneExcluded.KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}.allIntervals.USETHIS.txt.gz
### critically important that you use the 3 columns file here!!! 
# going to build in a check:
# count columns: 
colCount=`zcat $ancAllelesFile | head -n1 |  awk '{print NF}'`
if [ ${colCount} -ne 3 ]
then
	echo "THERE are less than or more than 3 columns in your anc alleles file! this will NOT WORK with bcftools annotate and may lead to subtle errors"
	exit 1
fi


ancFastascript=$gitdir/analyses/mutyper/generateAncestralFastaFromPolarizedVCF.ancAlleleVersion.py
refdir=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38

# this annotates and adds the above h1.txt header to the header and outputs a bgzipped vcf file 
mkdir -p $outdir/perChromosome
echo '##INFO=<ID=ancAllele,Number=1,Type=String,Description="Ancestral allele inferred by est-sfs (Keightley & Jackson, 2018). Alleles assigned using probabilities rather than cutoffs.">' > $outdir/h1.txt # make a header add on file for bcftools to add

chrspecificoutput=$outdir/perChromosome/${chromosome}.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BiallelicSNPsOnly.BUGFIX.${vcffilename}

# note: 
#total sites: 137058981
#total sites where anc!=ref: 8700447

##################### annotate vcf with anc allele ########################
### bcftoools will only annotate sites which are in the -a file so others won't get an annotation; don't want to include those sites in output
# filtering bcftools filter -e 'INFO/ancAllele="."' which will exclude sites which are missing the ancAllele (not that it doesn't ahve to actually say ancAllele=. in the vcf it can just be  missing and site will be excluded)
echo "starting bcftools"
# note even though all snps should be annotated now, need to keep the bcftools filter step because vcf contains indels and non biallelic snps that did not get polarized and throw errors downstream, so still remove those

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
