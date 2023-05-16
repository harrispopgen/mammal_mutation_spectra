#! /bin/bash
#$ -l h_rt=05:00:00,h_data=4G
#$ -o /net/harris/vol1/home/beichman/bears/reports/extsfs
#$ -e /net/harris/vol1/home/beichman/bears/reports/extsfs
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N prepForExtSFS


######## wrapper that will prepare a vcf for Keightley's ext sfs method

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip

### need focal sp list (mmd), outgroup1 (mmm), outgroup2 (ms) popfiles
#python3 script --infile [] \
#   --focalSampleList [] \
#           --outgroup1SampleList [] \
#               --outgroup2SampleList [] \
#                   --outfile []
# gives you a table with site information and gneotype counts in focal group and major allele in outgroups 
# note that outgroup1 must be the *closer phylogenetically*, then outgroup2 is further out.

projectdir=/net/harris/vol1/home/beichman/mice
vcfdir=$projectdir/vcfs/
vcffilename=$vcfdir/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz 

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
scriptdir=$gitdir/data_processing/variant_polarizing/KeightleyPolarizationMethod_mice
script=$scriptdir/Prepare.KeightlyExtSFSInfile.py

poplistdir=$gitdir/samples/mouseSamples/populationFilesForMutyper
focal=Mmd.txt
outgroup1=Mmm.txt
outgroup2=Ms.txt
flagForThisRun="20210225_focal${focal%.txt}_out1${outgroup1%.txt}_out2${outgroup2%.txt}"

outdir=/net/harris/vol1/home/beichman/mice/analyses/keightley_polarization/$flagForThisRun/input_tables_separated100K
outfile=$outdir/AllSites.focal.${focal%.txt}.out1.${outgroup1%.txt}.out2.${outgroup2%.txt}.PreppedForExtsfs.txt

python3 $script --infile $vcffilename \
--focalSampleList $poplistdir/$focal \
--outgroup1SampleList $poplistdir/$outgroup1 \
--outgroup2SampleList $poplistdir/$outgroup2 \
--outfile $outfile

# must not be bgzipped for extsfs
# then want to split up into 100K lines
split -l 100000 --numeric-suffixes=1 --suffix-length=5 $outfile ${outfile}_

gzip -f $outfile
     
     
