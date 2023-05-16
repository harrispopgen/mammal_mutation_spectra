#! /bin/bash
#$ -l h_rt=05:00:00,h_data=30G
#$ -o /net/harris/vol1/home/beichman/mice/reports/extsfs
#$ -e /net/harris/vol1/home/beichman/mice/reports/extsfs
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N extSFS

######### this script will run extsfs, then choose ancestral alleles (currently based on gte .9 lte .10 but that may change)
# update the output with anc allele choice
set -euo pipefail
### not running as an array for now because it runs quite quickly (only should take an hour or two)
module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for annotating vcf
module load python/3.7.7

########### species you used (for labelling input/output)###### 
focal=Ms
outgroup1=Mmd # note only one outgroup!

# config-file.txt: 
flagForThisRun="20220607_focal${focal%.txt}_out1${outgroup1%.txt}_polarizing_Mus_spretus"
projectdir=/net/harris/vol1/home/beichman/mice
wd=$projectdir/analyses/keightley_polarization/$flagForThisRun
outdir=$wd/extsfs_output
mkdir -p $outdir
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
scriptdir=$gitdir/data_processing/variant_polarizing/KeightleyPolarizationMethod_mice/Keightley_method_mice_SpretusAsFocal_MmdAsOutgroup
configfile=$scriptdir/config.1outgroup.KimuraModel.10MLRuns.txt
estsfs=/net/harris/vol1/home/beichman/bin/est-sfs-release-2.03/est-sfs
indir=$wd/input_tables_separated100K

assignAncAllelesScript=$scriptdir/assignAncestralAllele.Probabilistically.py # path ; works with 1 or 2 outgroups because indexes things by header rather than by column index

############## run ext -sfs #################

# pre-remove any old temp files

cd $indir # important to do!
if [ -e $indir/*TEMPORARYFOREXTSFS ]
then
rm $indir/*TEMPORARYFOREXTSFS # shouldn't be any of these but if there are from aborted past runs they mess stuff up so get rid of them
fi # this can be a problem if any TEMPORARY FOREXTSFS files are left 
# this will throw a warning that it didn't ahve any to delete -- that is good.

##### get list of files to loop ovr: #######
for fullinfilewithpositions in `ls *PreppedForExtsfs.txt_*` # want to make sure it's the just seaprated files not the original outfile, which is why you ahve _ after name (otherwise big nonsplit file gets included )
do

ls $indir/$fullinfilewithpositions
# make a random seedfile (will be overwritten after job is over)
echo $RANDOM > $indir/seedfile.txt

# okay then need just the 7,8 (would need 9 if doing TWO outgroups, but just 7 and 8 if just doing 1) columns for ext-sfs 
# use grep to skip the header
tempoutfilename=${fullinfilewithpositions}.TEMPORARYFOREXTSFS
awk 'BEGIN{OFS="\t";} {print $7,$8}' $indir/$fullinfilewithpositions > $indir/$tempoutfilename # going to delete this when done

$estsfs $configfile $indir/$tempoutfilename $indir/seedfile.txt $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.sfs.txt $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.txt
# care about the p values
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo "estSFS failed for $tempoutfilename"
    exit 1
fi

# this works!
######### check this carefully ###########
# then use JAR code to paste the probs from the output into the original input file so you get the chr and pos annotations
# note that lines that start with 0 in the output are info (comment) lines, whereas lines that do not start with refer to the line number of the input file
# silly formatting but it works
# this code will take the third column of the output (which has the prob of major allele being ancestral) and paste it together with original input


header="chr\tpos\tref\talt\tmajor\tminor\tfocalACGT\toutgroup1ACGT\tProbMajorAlleleIsAncestral"
echo -e $header > $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.WITHPOSITIONS.usethis.txt
cat $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.txt | grep -v "^0" | cut -d' ' -f3 | paste $fullinfilewithpositions - >> $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.WITHPOSITIONS.usethis.txt
# cut -f3 is for the 3 column Site Code P-major-ancestral; the cut part is good no matter how many outrgroups. 
rm $tempoutfilename

gzip -f $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.txt
gzip -f $outdir/extsfs.output.KimuraModel.10MLRuns.${fullinfilewithpositions}.pvalues.WITHPOSITIONS.usethis.txt

done


############## CHOOSE ANC ALLELE : PROBABILISTICALLY ##############
echo "starting anc allele assignment"
# count up intervals:
# for file in `ls *WITHPOSITIONS.usethis.txt.gz`
# per interval outdir:
intervalCount=`ls $outdir/*WITHPOSITIONS.usethis.txt.gz | wc -l` # count up files to go through
inprefix=KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}.PreppedForExtsfs.txt # fin whale has this longer prefix 
concatdir=$wd/concatted_extsfs_output_withAncAllele_KimuraModel_assignedUsingProbs/
mkdir -p $concatdir
# the inprefix specifies the interval
outprefix="AncAllelesAssignedUsingProbs.noneExcluded.KimuraModel.10MLRuns.focal.${focal}.out1.${outgroup1}" #I want to cut off the end of that inprefix so am specifiying the outprefix fully here 
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
