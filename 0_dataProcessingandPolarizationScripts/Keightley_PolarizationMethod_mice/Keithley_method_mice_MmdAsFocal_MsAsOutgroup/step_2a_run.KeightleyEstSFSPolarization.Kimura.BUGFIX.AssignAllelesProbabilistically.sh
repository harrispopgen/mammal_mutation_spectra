#! /bin/bash
#$ -l h_rt=05:00:00,h_data=4G
#$ -o /net/harris/vol1/home/beichman/bears/reports/extsfs
#$ -e /net/harris/vol1/home/beichman/bears/reports/extsfs
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N extSFS

######### this script will run extsfs, then choose ancestral alleles (currently based on gte .9 lte .10 but that may change)
# update the output with anc allele choice

### not running as an array for now because it runs quite quickly (only should take an hour or two)
module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for annotating vcf
module load python/3.7.7
########### species you used (for labelling input/output)###### 
focal=Mmd
outgroup1=Mmm
outgroup2=Ms

# config-file.txt: 
flagForThisRun="20210225_focal${focal}_out1${outgroup1}_out2${outgroup2}"
wd=/net/harris/vol1/home/beichman/mice/analyses/keightley_polarization/$flagForThisRun
outdir=$wd/extsfs_output_KimuraModel/
mkdir -p $outdir
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
scriptdir=$gitdir/data_processing/variant_polarizing/KeightleyPolarizationMethod_mice
configfile=$scriptdir/config.2outgroups.KimuraModel.10MLRuns.txt
estsfs=/net/harris/vol1/home/beichman/bin/est-sfs-release-2.03/est-sfs
indir=/net/harris/vol1/home/beichman/mice/analyses/keightley_polarization/input_tables_separated100K
projectdir=/net/harris/vol1/home/beichman/mice

assignAncAllelesScript=$scriptdir/assignAncestralAllele.Probabilistically.py # path

echo "starting extsfs"
############## run ext -sfs #################
for i in {00001..01371}
do
#intervalID=`printf "%05d\n" ${SGE_TASK_ID}`
intervalID=$i
outprefix=KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}_${intervalID}
# this is in format:
#n_outgroup 2 # number of outgroups
#model 1 # model (jukes cantor = 0, kimura 2 param =1 , 6-way rate = 2)
#nrandom 10 # number of ML runs (should be at least 10 for 6 way)
# est-sfs configfile datafile seedfile output-sfs output-pvals
# make a random seedfile (will be overwritten after job is over)
echo $RANDOM > $wd/seedfile.txt


fullinfilewithpositions=$indir/AllSites.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}.PreppedForExtsfs.txt_${intervalID}


# okay then need just the 7,8,9 columns for ext-sfs 
# use grep to skip the header
tempoutfilename=${fullinfilewithpositions}.TEMPORARYFOREXTSFS
awk 'BEGIN{OFS="\t";} {print $7,$8,$9}' $fullinfilewithpositions > $tempoutfilename # going to delete this when done
### also need to split into 100K line files and loop over them

$estsfs $configfile $tempoutfilename $wd/seedfile.txt $outdir/extsfs.output.${outprefix}.sfs.txt $outdir/extsfs.output.${outprefix}.pvalues.txt
# care out the p values

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo "error in extsfs"
    exit 1
fi

# this works!

# then use JAR code to paste the probs from the output into the original input file so you get the chr and pos annotations
# note that lines that start with 0 in the output are info (comment) lines, whereas lines that do not start with refer to the line number of the input file
# silly formatting but it works
# this code will take the third column of the output (which has the prob of major allele being ancestral) and paste it together with original input
# note that the position of the third column being the prob of maj being ancestral doesn't change based on # of outgroups
# after that column comes other probabilities about different trees that we dont' want to include 
header="chr\tpos\tref\talt\tmajor\tminor\tfocalACGT\toutgroup1ACGT\toutgroup2ACGT\tProbMajorAlleleIsAncestral"
echo -e $header > $outdir/extsfs.output.${outprefix}.pvalues.WITHPOSITIONS.usethis.txt
cat $outdir/extsfs.output.${outprefix}.pvalues.txt | grep -v "^0" | cut -d' ' -f3 | paste $fullinfilewithpositions - >> $outdir/extsfs.output.${outprefix}.pvalues.WITHPOSITIONS.usethis.txt
# I re-checked that this works on 20220506
rm $tempoutfilename

gzip -f $outdir/extsfs.output.${outprefix}.pvalues.txt
gzip -f $outdir/extsfs.output.${outprefix}.pvalues.WITHPOSITIONS.usethis.txt

done


############## CHOOSE ANC ALLELE ##############
############ now am assigning probabilistically instead of using hard cutoffs ########
echo "starting anc allele assignment"
# count up intervals:
# for file in `ls *WITHPOSITIONS.usethis.txt.gz`
intervalCount=`ls $wd/extsfs_output_KimuraModel/*WITHPOSITIONS.usethis.txt.gz | wc -l` # count up files to go through
inprefix=KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.out2.${outgroup2} # has interval ID at the end : _${intervalID} (include in python
concatdir=$wd/concatted_extsfs_output_withAncAllele_KimuraModel_assignedUsingProbs
mkdir -p $concatdir
outprefix="AncAllelesAssignedUsingProbs.noneExcluded.${inprefix}" # two files will get output, one with all information for troubleshooting/plotting and one with just 3 columns
# that is used for bcftools annotate

# this script will read in each file in the indir (previous step's outdir) 
python $assignAncAllelesScript --indir $wd/extsfs_output_KimuraModel/  --intervalCount $intervalCount --inprefix $inprefix --outdir $concatdir --outprefix $outprefix 

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo "error in assigning anc alleles"
    exit 1
fi


outfile_full=$concatdir/${outprefix}".allIntervals.allInfoForTroubleshooting.DONOTUSEFORBCFTOOLSANNOTATE.txt" # result of python script -- all columns, DO NOT USE FOR BCFTOOLS ANNOTATE
outfile3Columns=$concatdir/ThreeColumnsOnly.${outprefix}".allIntervals.USETHIS.txt" # use this for bcftools annotate
bgzip -f ${outfile_full} # must be bgzipped not gzipped 
bgzip -f ${outfile3Columns} # must be bgzipped, not gzipped for tabix to work

#### note I did this part manually but have updated script. so if script doesn't work in future i did it right by hand and so its' just some introduced bug (haven't rerun after changes) - so stay calm. 
# IT SEEMS IT CAN ONLY HAVE 3 COLUMNS !!! ### note you may need to change the final column number if fewer outgroups
#zcat ${ancAllelesFile_INTERMEDIATE}.gz | awk 'BEGIN {OFS="\t"} {print $1,$2,$11}' | bgzip > $concatdir/ThreeColumnsOnly.ConfidentAncAllelesSelectedWithAwk.90-10Range.extsfs.output.KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}.pvalues.WITHPOSITIONS.usethis.txt.gz
# then want to try to annotate the original vcf with these 
# make a little file with the new header info
### need to do the tabix (JAR instructions/bcftools annotate insturctions)
### DO NOT INDEX FULL FILE (will stop it from being used by accident for any bcftools stuff) tabix -S1 -s1 -b2 -e2 ${outfile_full}.gz # this indexes the AA table saying skip the header (-S1), that chr name is in field 1 (-s1) and that position start/end is field 2 (-b2 -e2)

tabix -S1 -s1 -b2 -e2 ${outfile3Columns}.gz # indexes the AA table saying skip the header (-S1), that chr name is in field 1 (-s1) and that position start/end is field 2 (-b2 -e2)

 
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo "error in tabix"
    exit 1
fi