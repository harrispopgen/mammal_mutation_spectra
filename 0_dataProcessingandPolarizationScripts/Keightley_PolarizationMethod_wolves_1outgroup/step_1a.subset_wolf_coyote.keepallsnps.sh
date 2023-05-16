#! /bin/bash
#$ -l h_rt=20:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/dogs/reports.nobackup/polarizing
#$ -e /net/harris/vol1/home/beichman/dogs/reports.nobackup/polarizing
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N subset_wolf_coy

module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)
module load bedtools/2.29.2

############## polarizing wolves ###############
# treating dogs as a sub-pop of wolves . not using their frequencies or allele states.
# but am keeping in *all* SNPs some of which may be 0/0 0/0 in wolf/coyote but are segregating in dogs.

######### step 1: pull out wolf/coy from vcf; be sure to keep all sites ###########
# not pulling *all* wolves; want as close to a pop as possible so excluding red wolf, iberian, tibetan etc. trying to just do grey wolves (some may or may not be)

samplesToInclude=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DogSpectrum/keightley_polarization/greywolves.coyote.forpolarizing.notallwolves.txt

# the vcf contains *indels* and triallelic stuff. so I want to exclude those. 
vcfdir=/net/harris/vol1/home/beichman/dogs/vcfs/
vcf=broad_umass_canid_variants.1.2.vcf.gz

mkdir $vcfdir/wolves_coyote_forpolarization
# *note here* I am selecting SNPs and subsetting wolf/coy because I want to keep all snps in that are snps in dogs but may not be in wolf/coy
# I have checked (on 20220421) that doing this does NOT remove 0/0 sites in wolf/coy. they are kept. which is good because I want to annotate those! 
# so don't worry that you're losing dog snps here, you're not! you're just getting rid of indels and triallelic sites. 
# note that missing data is kept in as well
bcftools view -m2 -M2 -v snps $vcfdir/$vcf -S $samplesToInclude -Oz > $vcfdir/wolves_coyote_forpolarization/wolf_coyote_notallwolves.biallelicsnpsonly.keepingmonomorphicbecausepolyindogs.vcf.gz

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in subsetting wolf/coy"
exit 1
else
echo "finished script"
fi