#! /bin/bash
#$ -l h_rt=50:00:00,h_data=20G
#$ -o /net/harris/vol1/home/beichman/bears/reports.nobackup/extsfs
#$ -e /net/harris/vol1/home/beichman/bears/reports.nobackup/extsfs
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N PBancestralFastaMaker
#$ -t 1-29
#$ -P sage
## polar bear has 29; bb has 31
set pipefail -euo
######### wrapper script to generate modified fasta from a reference fasta and a vcf file with AA= #######
########## you are here -- update this script
# note that the vcf should have excluded any sites that could not be polarized (ambiguous sites) -- but note that those sites need to not be used for anything in the spectrum since the ref genome won't have a call for them
# 
module load modules modules-init modules-gs # initialize modules 
module load bedtools/2.29.2
module load python/3.6.5
module load seqtk/1.3
source ~/install/virtualenvs/ancestralFasta/bin/activate # activate tsinfer virtual environment

intervalID=${SGE_TASK_ID} # don't need to prepend 0

focal=polar_bear # focal species
reference=polar_bear # the reference genome; want to use PB mapped to PB and BB mapped to BB
outgroup1=brown_bear
outgroup2=black_bear

refdir=/net/harris/vol1/home/beichman/reference_genomes/$reference

# interval-specific fasta files
intervaldir=$refdir/genomeIntervalBedFiles_ScaffsGT1Mb/intervalFastaFiles #

intervalSpecificFasta=$intervaldir/${reference}.interval_${intervalID}.fasta  

############## make ancestral fasta files ######################
echo "starting ancestral fasta making"


# config-file.txt: 
flagForThisRun="20220607_Keightley_focal${focal}_out1${outgroup1}_out2${outgroup2}"

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
projectdir=/net/harris/vol1/home/beichman/bears

wd=$projectdir/analyses/keightley_polarization/${flagForThisRun}

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject
scriptdir=$gitdir/data_processing/variant_polarizing/Keightley_Polarization_method_${focal} ### note is species specific 

ancFastascript=$scriptdir/generateAncestralFastaFromPolarizedVCF.ancAlleleVersion.py

intervalspecificvcf=$wd/vcfs_with_ancestral_alleles/perInterval/interval_${intervalID}.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.BUGFIX.SNPs.filt_variants_4b.LowCovIndsREMOVED.nomissfiltapplied.PASS.WARN.Only.CustomFilt.HardFilt.AnnotVar.TrimAlt.ref_${reference}.int_${intervalID}.vcf.gz


prefix=interval_${intervalID}_${focal}Ancestral_extSFS_AllelesAssignedUsingProbs.BUGFIX ### CHANGE THIS FOR DIFFERENT RUNS BASED ON OTHER CRITERIA
outdir=$projectdir/analyses/ancestralReferenceFastas/${flagForThisRun}_assignedAllelesUsingProbs_POLARIZING_${focal}
mkdir -p $outdir
### *note* this will make a copy of the interval ref fasta and copy it to outdir to be mutated; will not affect original.
###### (but be careful with this because pyfaidx does alter fastas in place, so if your script doesn't include a copying step you could be in trouble. my script does.)
# to run you want your interval specific vcf, an interval specific fasta and an interval specific prefix.
python $ancFastascript --vcffile $intervalspecificvcf --originalReferenceFasta $intervalSpecificFasta --refPrefix $prefix --outdir $outdir 

deactivate
