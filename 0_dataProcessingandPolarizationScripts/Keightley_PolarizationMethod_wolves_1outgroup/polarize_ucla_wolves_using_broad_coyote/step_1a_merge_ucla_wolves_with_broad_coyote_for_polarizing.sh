#! /bin/bash
#$ -l h_rt=20:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/dogs/reports.nobackup/polarizing
#$ -e /net/harris/vol1/home/beichman/dogs/reports.nobackup/polarizing
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N subset_coy_mergeWithUCLAWolves

module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)
module load bedtools/2.29.2

set -euo pipefail
# okay so I want to pull out coyote and merge it with ucla wolves in some fashion
# unlike when I merged apes to do hamming distances, I don't want ot use -0 --missing-to-ref
# because in this case even though I'm still only merging snps and not all sites I don't want ot assume anything about state
# because I don't care about sites where wolves are missing but coyote has a snp.
# Oh issue -- i don't have all sites for coyote and dog files are restricted to snps only
# so could be sites where coy + all broad wolves were reference, but ucla wolves aren't. coy would be missing for those sites
# but under Keightley it's ok those sites will still get annotated based on frequency in wolves
# okay so I think this is tolerable though not ideal. 


### need to pull out coyote 
broaddir=/net/harris/vol1/home/beichman/dogs/vcfs/
mkdir -p $broaddir/coyote_by_itself
ucladir=/net/harris/vol1/home/beichman/canids_fromEduardo/vcfs/All15AW/
mkdir -p $ucladir/merged_with_broad_coyote_for_polarizing
uclawolvesvcf=$ucladir/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing3_5_v4_mergeGaTKfiltered_varnonvar_All15AW_jointcalled_allchr.vcf.gz # this is snps only 	

broadwolvesvcf=$broaddir/wolves_coyote_forpolarization/wolf_coyote_notallwolves.biallelicsnpsonly.keepingmonomorphicbecausepolyindogs.vcf.gz
# this is biallelic sites only, but also keeps in all monomorphic 0/0 sites bc those are snps in dogs 
# so this is all the sites of the broad whole dataset, just subset to a set of wolves + coyote
# and indels removed

# name coy vcf to output:
coyotevcf=$broaddir/coyote_by_itself/coyote_by_itself.biallelicsnpsonly.keepingmonomorphicbecausepolyindogs.toMergeWithUCLAWolves.vcf.gz
bcftools view -s cal_coy $broadwolvesvcf -Oz -o $coyotevcf
# note small s because it's a individual name not a file 


exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in subsetting coy"
exit 1
fi

tabix -p vcf $coyotevcf

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in tabix"
exit 1
fi

# must be bgzipped and tabixed to work with merge

# name final output file: 
wolfcoyvcfoutput=$ucladir/merged_with_broad_coyote_for_polarizing/AWplusBroadCoyoteForPolarizing_BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing3_5_v4_mergeGaTKfiltered_varnonvar_All15AW_jointcalled_allchr.vcf.gz
bcftools merge -m snps $coyotevcf $uclawolvesvcf -Oz -o $wolfcoyvcfoutput

# note this results in lots of sites were all are missing (coy was missing in original broad vcf and ucla wolves don't have a nsp there), or where coy has a call but ucla wolves don't have a call and vice versa
# could filter out but also fine to just leave in and have be part of polarization 
exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in merging"
exit 1
fi

tabix -p vcf $wolfcoyvcfoutput

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in tabix"
exit 1
else
echo "finished script"
fi