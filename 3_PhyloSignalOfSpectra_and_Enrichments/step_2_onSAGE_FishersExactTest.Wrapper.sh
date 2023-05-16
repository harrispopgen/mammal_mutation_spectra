#! /bin/bash
#$ -l h_rt=150:00:00,mfree=30G
#$ -o /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/fishers_exact_test
#$ -e /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/reports.nobackup/fishers_exact_test
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N fisherstest_R

# Wrapper script to run fisher's exact test on 3,5,7mer spectra


wd=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified_SubsettingIndividuals/fishers_exact_test/20221121.plots.SubsampledIndividualsNotProjection.epsilon.1.sepCpGs.no.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode.9999999.MantelPermutations/""

## this needs to match the directory name that you got the spectra from 
# ooh maybe automate?
# mkdir of that name on sage
# use fetch to xfer the following files:
#all3merSpectra_plusTargets.ToTransferToSageForEnrichments.txt
#all5merSpectra_plusTargets.ToTransferToSageForEnrichments.txt
#all7merSpectra_plusTargets.ToTransferToSageForEnrichments.txt
module load modules modules-init modules-gs # initialize modules 
module load pcre2/10.39 R/4.1.2
# make sure you've installed the modules you'll need: ggplot2, dplyr,


# then can run the R script:
scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/multispecies_spectra/analyses/compare_spectrumDist_toPhyloHammingDist
script=$scriptdir/step_2_onSAGE_FishersExactTest.R
Rscript $script $wd # supply location as an argument