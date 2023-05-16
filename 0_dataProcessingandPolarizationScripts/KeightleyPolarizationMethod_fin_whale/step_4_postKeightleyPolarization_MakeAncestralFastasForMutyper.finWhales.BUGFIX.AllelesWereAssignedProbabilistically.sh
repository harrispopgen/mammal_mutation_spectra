#! /bin/bash
#$ -l h_rt=50:00:00,h_data=20G
#$ -o /net/harris/vol1/home/beichman/fin_whale/reports.nobackup/ancestralFasta/
#$ -e /net/harris/vol1/home/beichman/fin_whale/reports.nobackup/ancestralFasta/
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N ancestralFastaMaker
#$ -t 1-96
#$ -P sage
######### wrapper script to generate modified fasta from a reference fasta and a vcf file with AA= #######
########## you are here -- update this script
# note that the vcf should have excluded any sites that could not be polarized (ambiguous sites) -- but note that those sites need to not be used for anything in the spectrum since the ref genome won't have a call for them
# 
module load modules modules-init modules-gs # initialize modules 
module load bedtools/2.29.2
module load python/3.6.5
module load seqtk/1.3

intervalID=`printf "%02d\n" ${SGE_TASK_ID}` # prepend 0

refdir=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0


######### need to subset genome fasta into per interval using list files #########
contigList=$refdir/contiglist/BalAcu1.0_genomic.contiglist_${intervalID}.list
fullRefFasta=$refdir/BalAcu1.0.fa
perintervalrefdir=$refdir/per_interval_fastas
mkdir -p $perintervalrefdir
intervalSpecificFasta=$perintervalrefdir/BalAcu1.0.interval_${intervalID}.fa
# use seqtk subseq
seqtk subseq -l 60 $fullRefFasta $contigList > $intervalSpecificFasta
#          -l INT   sequence line length [0] 
# need -l 60 for default 
# this prints it with no wrapping. how can I get a warp? 
## need a virtual env
#pythonpath=/net/gs/vol3/software/modules-sw/python/3.6.5/Linux/CentOS7/x86_64/bin/python
#virtualenv  -p $pythonpath ~/install/virtualenvs/ancestralFasta # set up virtualenv with python 3.6.5
# activate it:
source ~/install/virtualenvs/ancestralFasta/bin/activate # activate tsinfer virtual environment
############## make ancestral fasta files ######################
echo "starting ancestral fasta making"

focal=ENP
outgroup1=MegNov
outgroup2=BalMus

# config-file.txt: 
flagForThisRun="20211020_focal${focal}_out1${outgroup1}_out2${outgroup2}"

# need to use contig list somehow. maybe getfasta? 
gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/fin_whale_spectrum
projectdir=/net/harris/vol1/home/beichman/fin_whale

wd=$projectdir/analyses/keightley_polarization/${flagForThisRun}
ancFastascript=$gitdir/KeightleyPolarizationMethod_fin_whale/generateAncestralFastaFromPolarizedVCF.ancAlleleVersion.py

 
# OLD VCF USED prior to 20220401 intervalspecificvcf=$wd/vcfs_with_ancestral_alleles/perInterval/OnlySitesWithConfidentAncAllelesFromExtSFS.90-10Range.PASSSITESONLY.failwarnandCpGIslandSitesfilteredout.JointCalls_f50b4_08_B_VariantFiltration_${intervalID}.vcf.gz
# new vcf used after 20220401 whcih contains CpG and repeat regions
intervalspecificvcf=$wd/vcfs_with_ancestral_alleles/perInterval/interval_${intervalID}.AncAllelesFromExtSFSAssignedUsingProbs.KimuraModel.PASSSITESPlusWarnCpGRefSitesONLY.lowqualremoved.BUGFIX.JointCalls_f50b4_08_B_VariantFiltration_${intervalID}.vcf.gz


prefix=interval_${intervalID}_finWhale_Ancestral_extSFS_AllelesAssignedUsingProbs_NOTEthisisMinkeWhaleRefGenome.butAncCallsAreforfinWhale.includesCpGandRepeatSites.BUGFIX ### CHANGE THIS FOR DIFFERENT RUNS BASED ON OTHER CRITERIA
outdir=$projectdir/analyses/ancestralReferenceFastas/${flagForThisRun}_assignedAllelesUsingProbs
mkdir -p $outdir
### *note* this will make a copy of the interval ref fasta and copy it to outdir to be mutated; will not affect original.
# (but be careful with this because pyfaidx does alter fastas in place, so if your script doesn't include a copying step you could be in trouble. my script does.)
# to run you want your interval specific vcf, an interval specific fasta and an interval specific prefix.
python $ancFastascript --vcffile $intervalspecificvcf --originalReferenceFasta $intervalSpecificFasta --refPrefix $prefix --outdir $outdir 

deactivate
