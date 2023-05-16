#! /bin/bash
#$ -l h_rt=05:00:00,h_data=30G
#$ -o /net/harris/vol1/home/beichman/bears/reports.nobackup/extsfs
#$ -e /net/harris/vol1/home/beichman/bears/reports.nobackup/extsfs
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N PBextSFS
#$ -t 1-29
#$ -P sage
## polar bear has 29; bb has 31

######### this script will run extsfs, then choose ancestral alleles (currently based on gte .9 lte .10 but that may change)
# update the output with anc allele choice
######## NEED TO MAKE THIS WORK FOR FIN WHALE WITH SGE TASK ID 
intervalID=${SGE_TASK_ID} # don't need to prepend 0

### not running as an array for now because it runs quite quickly (only should take an hour or two)
module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for annotating vcf
module load python/3.7.7

########### species you used (for labelling input/output)###### 
focal=polar_bear # focal species
reference=polar_bear # the reference genome; want to use PB mapped to PB and BB mapped to BB
outgroup1=brown_bear
outgroup2=black_bear


# config-file.txt: 
flagForThisRun="20220607_Keightley_focal${focal}_out1${outgroup1}_out2${outgroup2}"
projectdir=/net/harris/vol1/home/beichman/bears
wd=$projectdir/analyses/keightley_polarization/$flagForThisRun
outdir=$wd/extsfs_output/interval_${intervalID}
mkdir -p $outdir
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
scriptdir=$gitdir/data_processing/variant_polarizing/Keightley_Polarization_method_${focal} ### note is species specific 
configfile=$scriptdir/config.2outgroups.KimuraModel.10MLRuns.txt
estsfs=/net/harris/vol1/home/beichman/bin/est-sfs-release-2.03/est-sfs
indir=$wd/input_tables_separated100K/interval_${intervalID}

assignAncAllelesScript=$scriptdir/assignAncestralAllele.Probabilistically.py # path

# okay need to add the 0 infront of the interval. 
### want to put something in in case estsfs crashes -- this would currently skip the crashed file which isn't good.

############## run ext -sfs #################

#for i in {00001..01371}
#do
# pre-remove any old temp files

cd $indir # important 
rm $indir/*TEMPORARYFOREXTSFS # shouldn't be any of these but if there are from aborted past runs they mess stuff up so get rid of them
# this can be a problem if any TEMPORARY FOREXTSFS files are left 
for fullinfilewithpositions in `ls *PreppedForExtsfs.txt_*` # want to make sure it's the just seaprated files not the original outfile, which is why you ahve _ after name (otherwise big nonsplit file gets included )
do
#subintervalID=i
ls $indir/$fullinfilewithpositions
#intervalID=${SGE_TASK_ID}
#outprefix=KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}_${intervalID}
# this is in format:
#n_outgroup 2 # number of outgroups
#model 1 # model (jukes cantor = 0, kimura 2 param =1 , 6-way rate = 2)
#nrandom 10 # number of ML runs (should be at least 10 for 6 way)
# est-sfs configfile datafile seedfile output-sfs output-pvals
# make a random seedfile (will be overwritten after job is over)
echo $RANDOM > $indir/seedfile.txt


#fullinfilewithpositions=$indir/interval_${intervalID}.focal.${focal%.txt}.out1.${outgroup1%.txt}.out2.${outgroup2%.txt}.PreppedForExtsfs.txt_${i}

#outprefix=${fullinfilewithpositions}

# okay then need just the 7,8,9 columns for ext-sfs -- this depends on your no. of outgroups (!!) - with canids its' just 7,8
# use grep to skip the header
tempoutfilename=${fullinfilewithpositions}.TEMPORARYFOREXTSFS
awk 'BEGIN{OFS="\t";} {print $7,$8,$9}' $indir/$fullinfilewithpositions > $indir/$tempoutfilename # going to delete this when done
### also need to split into 100K line files and loop over them

$estsfs $configfile $indir/$tempoutfilename $indir/seedfile.txt $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.sfs.txt $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.txt
# care out the p values
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo "estSFS failed for $tempoutfilename"
    exit 1
fi

# this works!

# then use JAR code to paste the probs from the output into the original input file so you get the chr and pos annotations
# note that lines that start with 0 in the output are info (comment) lines, whereas lines that do not start with refer to the line number of the input file
# silly formatting but it works
# this code will take the third column of the output (which has the prob of major allele being ancestral) and paste it together with original input
header="chr\tpos\tref\talt\tmajor\tminor\tfocalACGT\toutgroup1ACGT\toutgroup2ACGT\tProbMajorAlleleIsAncestral"
echo -e $header > $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.WITHPOSITIONS.usethis.txt
cat $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.txt | grep -v "^0" | cut -d' ' -f3 | paste $fullinfilewithpositions - >> $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.WITHPOSITIONS.usethis.txt
# f3 works regardless of no of outgroups
rm $tempoutfilename

gzip -f $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.txt
gzip -f $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.WITHPOSITIONS.usethis.txt

done


############## CHOOSE ANC ALLELE PROBABILISTICALLY ##############
echo "starting anc allele assignment"
# count up intervals:
# for file in `ls *WITHPOSITIONS.usethis.txt.gz`
# per interval outdir:
intervalCount=`ls $outdir/*WITHPOSITIONS.usethis.txt.gz | wc -l` # count up files to go through
inprefix=KimuraModel.10MLRuns.interval_${intervalID}.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}.PreppedForExtsfs.txt #
# files are in fmt extsfs.output.KimuraModel.10MLRuns.interval_01.focal.ENP.out1.MegNov.out2.BalMus.PreppedForExtsfs.txt_00025.pvalues.WITHPOSITIONS.usethis.txt.gz
# has interval ID at the end : _${intervalID} (include in python
concatdir=$wd/concatted_extsfs_output_withAncAllele_KimuraModel_assignedUsingProbs/interval_${intervalID}/
mkdir -p $concatdir
# the inprefix specifies the interval
outprefix="AncAllelesAssignedUsingProbs.noneExcluded.KimuraModel.10MLRuns.interval_${intervalID}.focal.${focal}.out1.${outgroup1}.out2.${outgroup2}" #I want to cut off the end of that inprefix so am specifiying the outprefix fully here 
# two files will get output, one with all information for troubleshooting/plotting and one with just 3 columns
# that is used for bcftools annotate

# this script will read in each file in the indir (previous step's outdir) 
python $assignAncAllelesScript --indir $outdir/  --intervalCount $intervalCount --inprefix $inprefix --outdir $concatdir --outprefix $outprefix 

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
