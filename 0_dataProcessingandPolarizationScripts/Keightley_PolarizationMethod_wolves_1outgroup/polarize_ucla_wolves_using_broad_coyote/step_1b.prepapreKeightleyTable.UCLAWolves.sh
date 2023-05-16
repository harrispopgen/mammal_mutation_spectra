#! /bin/bash
#$ -l h_rt=25:00:00,h_data=4G
#$ -o /net/harris/vol1/home/beichman/dogs/reports.nobackup/polarizing
#$ -e /net/harris/vol1/home/beichman/dogs/reports.nobackup/polarizing
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N prepForExtSFS_uclawolves

######## wrapper that will prepare a vcf for Keightley's ext sfs method

# nice thing I discovered that keightley does: if all are missing data except outgroup it skips (useful)
module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip

### need focal sp list (mmd), outgroup1 (mmm) [ no outgroup 2 ]
#python3 script --infile [] \
#   --focalSampleList [] \
#           --outgroup1SampleList [] \
#                   --outfile []
# gives you a table with site information and gneotype counts in focal group and major allele in outgroups 
# note that outgroup1 must be the *closer phylogenetically*; trying without outgroup2 

projectdir=/net/harris/vol1/home/beichman/canids_fromEduardo


vcfdir=$projectdir/vcfs/All15AW/merged_with_broad_coyote_for_polarizing

vcffilename=$vcfdir/AWplusBroadCoyoteForPolarizing_BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing3_5_v4_mergeGaTKfiltered_varnonvar_All15AW_jointcalled_allchr.vcf.gz

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DogSpectrum

scriptdir=$gitdir/keightley_polarization/polarize_ucla_wolves_using_broad_coyote
script=$scriptdir/Prepare.KeightlyExtSFSInfile.ONEOUTGROUPONLY.py # modified to just work with a single outgroup (hopefully will work)

poplistdir=$scriptdir # have them in keightley script dir

focal=ucla_wolves.txt # ucla wolves! 
outgroup1=coyote.txt 
flagForThisRun="20220718_focal${focal%.txt}_out1${outgroup1%.txt}"

outdir=${projectdir}/analyses/keightley_polarization/$flagForThisRun/input_tables_separated100K
mkdir -p $outdir
outfile=$outdir/focal.${focal%.txt}.out1.${outgroup1%.txt}.PreppedForExtsfs.txt

python3 $script --infile $vcffilename \
--focalSampleList $poplistdir/$focal \
--outgroup1SampleList $poplistdir/$outgroup1 \
--outfile $outfile


exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in python table script"
exit 1
else
echo "finished converting to tables"
fi


# must not be bgzipped for extsfs
# then want to split up into 100K lines
split -l 100000 --numeric-suffixes=1 --suffix-length=5 $outfile ${outfile}_


exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in splitting"
exit 1
else
echo "finished splitting"
fi

gzip -f $outfile


exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in gzip"
exit 1
else
echo "finished script"
fi
     
