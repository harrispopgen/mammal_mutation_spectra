#! /bin/bash
#$ -l h_rt=05:00:00,h_data=4G
#$ -o /net/harris/vol1/home/beichman/fin_whale/reports.nobackup/polarizing
#$ -e /net/harris/vol1/home/beichman/fin_whale/reports.nobackup/polarizing
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N prepForExtSFS
#$ -t 01-96

######## wrapper that will prepare a vcf for Keightley's ext sfs method

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
intervalID=`printf "%02d\n" ${SGE_TASK_ID}` # prepend 0

### need focal sp list (mmd), outgroup1 (mmm), outgroup2 (ms) popfiles
#python3 script --infile [] \
#   --focalSampleList [] \
#           --outgroup1SampleList [] \
#               --outgroup2SampleList [] \
#                   --outfile []
# gives you a table with site information and gneotype counts in focal group and major allele in outgroups 
# note that outgroup1 must be the *closer phylogenetically*, then outgroup2 is further out.

projectdir=/net/harris/vol1/home/beichman/fin_whale/

vcfdir=$projectdir/vcfs.nobackup/downloaddate_20211015_whole_genomes_plusOutgroups
vcffilename=$vcfdir/JointCalls_f50b4_08_B_VariantFiltration_${intervalID}.vcf.gz 

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/fin_whale_spectrum

scriptdir=$gitdir/KeightleyPolarizationMethod_fin_whale
script=$scriptdir/Prepare.KeightlyExtSFSInfile.py

poplistdir=$gitdir/samples/populationFilesForKeightleyPolarization

focal=ENP.txt
outgroup1=MegNov.txt # closer outgroup -- this case is MegNov (humpback) which I didn't realize was closer than blue whale. 
outgroup2=BalMus.txt # further outgroup (blue whale)
flagForThisRun="20211020_focal${focal%.txt}_out1${outgroup1%.txt}_out2${outgroup2%.txt}"

# since genome is split into 100 intervals already I don't think I need to do the 100K thing for this. So let's get rid of that?
#outdir=/net/harris/vol1/home/beichman/fin_whale/analyses/keightley_polarization/$flagForThisRun/input_tables_separated100K
outdir=/net/harris/vol1/home/beichman/fin_whale/analyses/keightley_polarization/$flagForThisRun/input_tables_separated100K/interval_${intervalID}
mkdir -p $outdir
outfile=$outdir/interval_${intervalID}.focal.${focal%.txt}.out1.${outgroup1%.txt}.out2.${outgroup2%.txt}.PreppedForExtsfs.txt

python3 $script --infile $vcffilename \
--focalSampleList $poplistdir/$focal \
--outgroup1SampleList $poplistdir/$outgroup1 \
--outgroup2SampleList $poplistdir/$outgroup2 \
--outfile $outfile

# must not be bgzipped for extsfs
# then want to split up into 100K lines
# SKIPPING 
# splitting each interval into sub intervals:
split -l 100000 --numeric-suffixes=1 --suffix-length=5 $outfile ${outfile}_

gzip -f $outfile
     
     
