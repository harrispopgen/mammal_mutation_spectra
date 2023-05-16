#! /bin/bash
#$ -l h_rt=50:00:00,h_data=20G
#$ -o /net/harris/vol1/home/beichman/mice/reports/ancestralFasta/
#$ -e /net/harris/vol1/home/beichman/mice/reports/ancestralFasta/
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N MSancestralFastaMaker
#$ -t 1-1
#$ -P sage
######### wrapper script to generate modified fasta from a reference fasta and a vcf file with AA= #######

module load modules modules-init modules-gs # initialize modules 
module load bedtools/2.29.2
module load python/3.6.5

## need a virtual env
#pythonpath=/net/gs/vol3/software/modules-sw/python/3.6.5/Linux/CentOS7/x86_64/bin/python
#virtualenv  -p $pythonpath ~/install/virtualenvs/ancestralFasta # set up virtualenv with python 3.6.5
# activate it:
source ~/install/virtualenvs/ancestralFasta/bin/activate # activate tsinfer virtual environment
############## make ancestral fasta files ######################
echo "starting ancestral fasta making"

chromosome=chr${SGE_TASK_ID}

focal=Ms
outgroup1=Mmd

# config-file.txt: 
flagForThisRun="20220607_focal${focal%.txt}_out1${outgroup1%.txt}_polarizing_Mus_spretus"

refdir=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
projectdir=/net/harris/vol1/home/beichman/mice
wd=/net/harris/vol1/home/beichman/mice/analyses/keightley_polarization/${flagForThisRun}
ancFastascript=$gitdir/data_processing/variant_polarizing/KeightleyPolarizationMethod_mice/Keightley_method_mice_SpretusAsFocal_MmdAsOutgroup/generateAncestralFastaFromPolarizedVCF.ancAlleleVersion.py
 
chrspecificvcf=$wd/vcfs_with_ancestral_alleles/perChromosome/${chromosome}.PolarizingSpretus.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BiallelicSNPsOnly.BUGFIX.AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz 

refFasta=$refdir/${chromosome}.fa
prefix=${chromosome}_MusSpretusAncestral_extSFS_assignedUsingProbs.BUGFIX ### CHANGE THIS FOR DIFFERENT RUNS BASED ON OTHER CRITERIA
outdir=$projectdir/analyses/ancestralReferenceFastas/${flagForThisRun}_KimuraModel_AncAlleleAssignedUsingProbs_POLARIZED_MUS_SPRETUS
mkdir -p $outdir
### *note* this will make a copy of the interval ref fasta and copy it to outdir to be mutated; will not affect original.
# (but be careful with this because pyfaidx does alter fastas in place, so if your script doesn't include a copying step you could be in trouble. my script does.)
# to run you want your interval specific vcf, an interval specific fasta and an interval specific prefix.
python $ancFastascript --vcffile $chrspecificvcf --originalReferenceFasta $refFasta --refPrefix $prefix --outdir $outdir 

deactivate
