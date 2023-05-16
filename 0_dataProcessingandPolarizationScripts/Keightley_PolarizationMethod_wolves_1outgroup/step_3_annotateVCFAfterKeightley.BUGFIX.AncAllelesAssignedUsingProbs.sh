#! /bin/bash
#$ -l h_rt=50:00:00,h_data=10G
#$ -o /net/harris/vol1/home/beichman/dogs/reports.nobackup/polarizing
#$ -e /net/harris/vol1/home/beichman/dogs/reports.nobackup/polarizing
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N annotateAfterKeightley


######### this script will (skipping xy and M and UN) split vcf by chr, annotate and tabix and make ancestral fasta based on ancAllele designation 
# so far only working with outside of 90-10 bound , skipping all other sites. 

### not running as an array for now because it runs quite quickly (only should take an hour or two)
module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12  ### NOTE vERSION MUST BE >1.11 or will not work correctly! 

focal=wolf
outgroup1=coyote

# config-file.txt: 
flagForThisRun="20220421_focal${focal}_out1${outgroup1}"
projectdir=/net/harris/vol1/home/beichman/dogs
wd=$projectdir/analyses/keightley_polarization/$flagForThisRun
outdir=$wd/extsfs_output/
mkdir -p $outdir
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DogSpectrum
scriptdir=$gitdir/keightley_polarization

######### have two vcfs I want to annotate. Want to annotate the full dog vcf and the wolf/coy only vcf #######

concatdir=$wd/concatted_extsfs_output_withAncAllele_KimuraModel_assignedUsingProbs
ancAllelesFile=$concatdir/ThreeColumnsOnly.AncAllelesAssignedUsingProbs.noneExcluded.KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.allIntervals.USETHIS.txt.gz
#### and updating to now have alleles assigned using probabilities ! (20220531)
###### BUG FIX!!! now this has 3 columns and will behave correctly with bcftools (!!!!)
# this annotates and adds the above h1.txt header to the header and outputs a bgzipped vcf file 
mkdir -p $wd/vcfs_with_ancestral_alleles

headerFile=$wd/vcfs_with_ancestral_alleles/h1.interval.txt # make a header add on file for bcftools to add
echo '##INFO=<ID=ancAllele,Number=1,Type=String,Description="Ancestral allele inferred by est-sfs (Keightley & Jackson, 2018). Alleles assigned probabilistically.">' > $headerFile # make a header add on file for bcftools to add
# don't want h1.txt constantly overwritten by multiple jobs -- causes errors.


###################### vcf to annotate 1: wolf+coyote ###########
vcfdir=$projectdir/vcfs/wolves_coyote_forpolarization
vcffilename=wolf_coyote_notallwolves.biallelicsnpsonly.keepingmonomorphicbecausepolyindogs.vcf.gz


########### annotate vcf 1 with anc allele #######
### bcftoools will only annotate sites which are in the -a file so others won't get an annotation; don't want to include those sites in output
# filtering bcftools filter -e 'INFO/ancAllele="."' which will exclude sites which are missing the ancAllele (not that it doesn't ahve to actually say ancAllele=. in the vcf it can just be  missing and site will be excluded)
output=$wd/vcfs_with_ancestral_alleles/wolves_only/AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX.excludingCoyote.${vcffilename}
mkdir -p $wd/vcfs_with_ancestral_alleles/wolves_only
echo "starting bcftools for vcf1"
# remove coyote 
# exclude sites that weren't annotated (still doing this because there may be non-snps in the file that don't get annotated)
bcftools view -s ^"cal_coy" $vcfdir/$vcffilename -Ou | bcftools annotate -h $headerFile -a ${ancAllelesFile} -c CHROM,POS,ancAllele -Ou | bcftools filter -e 'INFO/ancAllele="."' -Oz -o $output  # this will be wolf-only vcf , but should include monomorphic sites that are poly in dogs (good)

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in annotation of vcf1"
exit 1
else
echo "finished annotating vcf1"
fi

tabix -p vcf $output

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in tabix of vcf1"
exit 1
else
echo "finished tabix of vcf1"
fi
# this will be better for annotation because a smaller vcf file -- maybe use it to annotate the fasta file? 

###################### vcf to annotate 2: full dog/coyote/wolf vcf (original from broad) ###########

vcfdir=$projectdir/vcfs/
vcffilename=broad_umass_canid_variants.1.2.vcf.gz
# same anc alleles file etc. 

########### annotate vcf 1 with anc allele #######
### bcftoools will only annotate sites which are in the -a file so others won't get an annotation; don't want to include those sites in output
# filtering bcftools filter -e 'INFO/ancAllele="."' which will exclude sites which are missing the ancAllele (not that it doesn't ahve to actually say ancAllele=. in the vcf it can just be  missing and site will be excluded)
output=$wd/vcfs_with_ancestral_alleles/AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX.${vcffilename}

echo "starting bcftools for vcf2"
# remove coyote 
# exclude sites that weren't annotated 
# don't need to explicitly remove non biallelic snps, they just won't get annotated and will be removed anyway
bcftools view $vcfdir/$vcffilename -Ou | bcftools annotate -h $headerFile -a ${ancAllelesFile} -c CHROM,POS,ancAllele -Ou | bcftools filter -e 'INFO/ancAllele="."' -Oz -o $output  # this will be wolf-only vcf , but should include monomorphic sites that are poly in dogs (good)

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in annotation of vcf2"
exit 1
else
echo "finished annotating vcf2"
fi

tabix -p vcf $output

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in tabix of vcf2"
exit 1
else
echo "finished tabix of vcf2 -- finished script"
fi
