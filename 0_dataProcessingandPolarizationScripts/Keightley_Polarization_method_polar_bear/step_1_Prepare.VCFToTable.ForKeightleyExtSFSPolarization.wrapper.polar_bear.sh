#! /bin/bash
#$ -l h_rt=05:00:00,h_data=4G
#$ -o /net/harris/vol1/home/beichman/bears/reports.nobackup/extsfs
#$ -e /net/harris/vol1/home/beichman/bears/reports.nobackup/extsfs
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N PBprepForExtSFS
#$ -t 1-29
#$ -P sage
## polar bear has 29; bb has 31
######## wrapper that will prepare a vcf for Keightley's ext sfs method

set -euo pipefail
 


module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
interval=${SGE_TASK_ID} # don't need to prepend 0


### need focal sp list (mmd), outgroup1 (mmm), outgroup2 (ms) popfiles
#python3 script --infile [] \
#   --focalSampleList [] \
#           --outgroup1SampleList [] \
#               --outgroup2SampleList [] \
#                   --outfile []
# gives you a table with site information and gneotype counts in focal group and major allele in outgroups 
# note that outgroup1 must be the *closer phylogenetically*, then outgroup2 is further out.
focal=polar_bear # focal species
reference=polar_bear # the reference genome; want to use PB mapped to PB and BB mapped to BB
outgroup1=brown_bear
outgroup2=black_bear

genotypeDate=20200916 # date genotypes were called. this is useful to distinguish runs if you add/remove individuals and redo joint GT calling

projectdir=/net/harris/vol1/home/beichman/bears

vcfdir=$projectdir/variant_calling/mapped_to_${reference}/vcfs/vcf_${genotypeDate}_${reference}/interval_${interval}/SNPsOnly/
vcffilename=$vcfdir/SNPs.filt_variants_4b.LowCovIndsREMOVED.nomissfiltapplied.PASS.WARN.Only.CustomFilt.HardFilt.AnnotVar.TrimAlt.ref_${reference}.int_${interval}.vcf.gz
# vcf with low cov individuals removed but no missingness filter applied

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject

scriptdir=$gitdir/data_processing/variant_polarizing/Keightley_Polarization_method_${focal} ### note is species specific 
script=$scriptdir/Prepare.KeightlyExtSFSInfile.py

poplistdir=$gitdir/samples/


flagForThisRun="20220607_Keightley_focal${focal}_out1${outgroup1}_out2${outgroup2}"

outdir=$projectdir/analyses/keightley_polarization/$flagForThisRun/input_tables_separated100K/interval_${interval}
mkdir -p $outdir
outfile=$outdir/interval_${interval}.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}.PreppedForExtsfs.txt
# note using high coverage individuals 
# note bear-specific format of poplist files
python3 $script --infile $vcffilename \
--focalSampleList $poplistdir/AB.${focal}.HIGHCOVONLY.txt  \
--outgroup1SampleList $poplistdir/AB.${outgroup1}.HIGHCOVONLY.txt \
--outgroup2SampleList $poplistdir/AB.${outgroup2}.txt \
--outfile $outfile 

# must not be bgzipped for extsfs
# then want to split up into 100K lines
# splitting each interval into sub intervals:
split -l 100000 --numeric-suffixes=1 --suffix-length=5 $outfile ${outfile}_

gzip -f $outfile
     
     
