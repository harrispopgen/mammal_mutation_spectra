#! /bin/bash
#$ -l h_rt=05:00:00,h_data=4G
#$ -o /net/harris/vol1/home/beichman/fin_whale/reports.nobackup/polarizing
#$ -e /net/harris/vol1/home/beichman/fin_whale/reports.nobackup/polarizing
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N annotateAfterKeightley
#$ -t 1-96
#$ -P sage
######### this script will act per interval (skipping xy and M and UN) to split vcf by chr, annotate and tabix and make ancestral fasta based on ancAllele designation 
# so far only working with outside of 90-10 bound , skipping all other sites. 

### not running as an array for now because it runs quite quickly (only should take an hour or two)
module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12  ### NOTE vERSION MUST BE >1.11 or will not work correctly! 

# bcftools/1.9 # for annotating vcf
#### NOTE; prior to 20220401 was using bcftools 1.9 but want to use >1.11 now because of filtering rules:


# starting with bcftoosl 1.11, the FILTER column can be queried as follows:
# 
# FILTER="PASS"
# FILTER="A"          .. exact match, for example "A;B" does not pass
# FILTER!="A"         .. exact match, for example "A;B" does pass
# FILTER~"A"          .. both "A" and "A;B" pass
# FILTER!~"A"         .. neither "A" nor "A;B" pass
# this matters because after 20220401 I am keeping in FAIL_CpGREP sites (going to filter downstream)
# and if I use bcftools prior to 1.9 it would keep any sites that has FAIL_CPgRep but also fails other things (like FAIL_CPGREP;FAIL_QUAL)
# don't want to do that. s


intervalID=`printf "%02d\n" ${SGE_TASK_ID}` # prepend 0

focal=ENP
outgroup1=MegNov
outgroup2=BalMus

# config-file.txt: 
flagForThisRun="20211020_focal${focal}_out1${outgroup1}_out2${outgroup2}"
wd=/net/harris/vol1/home/beichman/fin_whale/analyses/keightley_polarization/$flagForThisRun
#outdir=$wd/extsfs_output/
#mkdir -p $outdir
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/fin_whale_spectrum
scriptdir=$gitdir/KeightleyPolarizationMethod_fin_whale
projectdir=/net/harris/vol1/home/beichman/fin_whale

vcfdir=$projectdir/vcfs.nobackup/downloaddate_20211015_whole_genomes_plusOutgroups
vcffilename=JointCalls_f50b4_08_B_VariantFiltration_${intervalID}.vcf.gz 



concatdir=$wd/concatted_extsfs_output_withAncAllele_KimuraModel_assignedUsingProbs/interval_${intervalID}
ancAllelesFile=$concatdir/ThreeColumnsOnly.AncAllelesAssignedUsingProbs.noneExcluded.KimuraModel.10MLRuns.interval_${intervalID}.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}.allIntervals.USETHIS.txt.gz
### ^^^ FIXED BIG BUG! MUST BE 3 COLUMNS ONLY!!!!!
#ancFastascript=$scriptdir/postKeightleyPolarization_MakeAncestralFastasForMutyper.finWhales.sh
#refdir=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0


# this annotates and adds the above h1.txt header to the header and outputs a bgzipped vcf file 
mkdir -p $wd/vcfs_with_ancestral_alleles/perInterval

headerFile=$wd/vcfs_with_ancestral_alleles/perInterval/h1.interval${intervalID}.txt # make a header add on file for bcftools to add; labeling with interval so they don't overwrite and interfere when run in parallel
echo '##INFO=<ID=ancAllele,Number=1,Type=String,Description="Ancestral allele inferred by est-sfs (Keightley & Jackson, 2018).">' > $headerFile # make a header add on file for bcftools to add
# don't want h1.txt constantly overwritten by multiple jobs -- causes errors.


##################### annotate vcf with anc allele ########################
### bcftoools will only annotate sites which are in the -a file so others won't get an annotation; don't want to include those sites in output
# filtering bcftools filter -e 'INFO/ancAllele="."' which will exclude sites which are missing the ancAllele (not that it doesn't ahve to actually say ancAllele=. in the vcf it can just be  missing and site will be excluded)
# note this should also annotate sites that are fixed 1/1 ?
# OLD NAME: intspecificoutput=$wd/vcfs_with_ancestral_alleles/perInterval/OnlySitesWithConfidentAncAllelesFromExtSFS.90-10Range.PASSSITESONLY.failwarnandCpGIslandSitesfilteredout.${vcffilename}
# new name indicatign that I'm keeping the masked stuff in as of 20220401
intspecificoutput=$wd/vcfs_with_ancestral_alleles/perInterval/interval_${intervalID}.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.PASSSITESPlusWarnCpGRefSitesONLY.lowqualremoved.BUGFIX.${vcffilename}

echo "starting bcftools"
# ALSO WANT TO RESTRICT TO PASS-ONLY SITES AT THIS POINT!!! -- making a change and re-doing on 20220401! 
# I want to filter repeat regions a bit less conservatively as it was previously done
# So I want to leave in sites that are cpg or repeat (will filter myself later using masks) 
# but I want to eclude things that failed due to quality etc. 
# 
# old way I did this: bcftools view $vcfdir/$vcffilename | bcftools annotate -h $headerFile -a ${ancAllelesFile} -c CHROM,POS,ancAllele | bcftools view -i "%FILTER='PASS'" | bcftools filter -e 'INFO/ancAllele="."' -Oz -o $intspecificoutput  # this will e
# no longer need the -r part because each vcf is per interval but when I was doing mice this is when I would pull out each chrososome separately and where I'd ditch the Y and x chroms
# had -h in there twice for some reason; delete

# new way: keeping in sites that are CpG/Rep but have no other fail criteria
# ONLY WORKS CORRECLTY WITH BCFTOOLS >1.11 !!!
bcftools view $vcfdir/$vcffilename | bcftools annotate -h $headerFile -a ${ancAllelesFile} -c CHROM,POS,ancAllele | bcftools view -i "%FILTER='PASS' | %FILTER='FAIL_CpGRep'" | bcftools filter -e 'INFO/ancAllele="."' -Oz -o $intspecificoutput  # this will e
### annotation file MUST be 3 columsn only !!!!! 

tabix -p vcf $intspecificoutput

# note that this vcf still contains the outgroup sequences as well.
