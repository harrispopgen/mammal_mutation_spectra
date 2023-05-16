#! /bin/bash
#$ -l h_rt=50:00:00,h_data=20G
#$ -o /net/harris/vol1/home/beichman/dogs/reports.nobackup/make_anc_fasta
#$ -e /net/harris/vol1/home/beichman/dogs/reports.nobackup/make_anc_fasta
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N ancestralFastaMaker
#$ -P sage
######### wrapper script to generate modified fasta from a reference fasta and a vcf file with AA= #######
# note that the vcf should have excluded any sites that could not be polarized (ambiguous sites) -- but note that those sites need to not be used for anything in the spectrum since the ref genome won't be correctly polarized for them! # at some point udpate to an "X" instead of just not updating perhaps
# note: I am choosing to continue to use the wolf-only vcf file because it is smaller and easier to work with. has the same sites as the dog+wolf vcf.  
module load modules modules-init modules-gs # initialize modules 
module load bedtools/2.29.2
module load python/3.6.5
module load seqtk/1.3

refdir=/net/harris/vol1/home/beichman/reference_genomes/canFam3

fullRefFasta=$refdir/canFam3.fa

## need a virtual env
#pythonpath=/net/gs/vol3/software/modules-sw/python/3.6.5/Linux/CentOS7/x86_64/bin/python
#virtualenv  -p $pythonpath ~/install/virtualenvs/ancestralFasta # set up virtualenv with python 3.6.5
# activate it:
source ~/install/virtualenvs/ancestralFasta/bin/activate # activate tsinfer virtual environment
############## make ancestral fasta files ######################
echo "starting ancestral fasta making"

focal=ucla_wolves
outgroup1=coyote

flagForThisRun="20220718_focal${focal}_out1${outgroup1}"

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DogSpectrum
projectdir=/net/harris/vol1/home/beichman/canids_fromEduardo

wd=$projectdir/analyses/keightley_polarization/${flagForThisRun}
ancFastascript=$gitdir/keightley_polarization/polarize_ucla_wolves_using_broad_coyote/generateAncestralFastaFromPolarizedVCF.ancAlleleVersion.py

########## note here: there are TWO vcfs you could use. one has all dogs. one has wolf/coy only. but should have same ***NUMBER*** of sites because wolf-only vcf contained poly dog sites as well. 
# so just use the wolf/coy one because it's smaller (fewer inds) and faster to use 
vcf=$wd/vcfs_with_ancestral_alleles/wolves_only/AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX.excludingCoyote.BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing3_5_v4_mergeGaTKfiltered_varnonvar_All15AW_jointcalled_allchr.vcf.gz # note here i am copying in vcf name 
######### FIXING BUG with polarization on 20220429

prefix=canFam3.0_extSFS_AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX ### CHANGE THIS FOR DIFFERENT RUNS BASED ON OTHER CRITERIA
outdir=$projectdir/analyses/ancestralReferenceFastas/${flagForThisRun}_ancAllelesAssignedUsingProbs
mkdir -p $outdir
### *note* this will make a copy of the ref fasta and copy it to outdir to be mutated; will not affect original.
# (but be careful with this because pyfaidx does alter fastas in place, so if your script doesn't include a copying step you could be in trouble. my script does.)
python $ancFastascript --vcffile $vcf --originalReferenceFasta $fullRefFasta --refPrefix $prefix --outdir $outdir 

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in making anc fasta"
exit 1
else
echo "finished making anc fasta"
fi

###### want to split ancestral fasta by chromosome: (for canids; may not want to do for other species) ######## 
perintervalancrefdir=$outdir/per_chr_ancestral_fastas_autosomesonly
mkdir -p $perintervalancrefdir
# use seqtk subseq
for chr in {1..38} # dog-specific ; will only select autosomes
do
echo "separating out chr $chr"
# make a little file with the chr name for seqtk to use: 
echo chr${chr} > $perintervalancrefdir/subseq.chr${chr}.txt
intervalSpecificFasta=$perintervalancrefdir/FINAL.MODIFIED.ANCESTRAL.${prefix}_chr${chr}.fasta # name of 
seqtk subseq -l 60 $outdir/FINAL.MODIFIED.ANCESTRAL.${prefix}.fasta $perintervalancrefdir/subseq.chr${chr}.txt > $intervalSpecificFasta

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
echo "error in separating anc fasta $chr"
exit 1
fi

done



deactivate