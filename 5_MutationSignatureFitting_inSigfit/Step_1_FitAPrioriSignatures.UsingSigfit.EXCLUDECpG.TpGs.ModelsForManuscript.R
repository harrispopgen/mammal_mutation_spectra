####### okay time to read these in and run sigfit ########

# want to do two things

# want to fit the a priori aging signatures to the species 1mer-CpG spectra

# fit the a prior SBS1+SBS5 to the 3mer spectrum

# choose one population per species (match what was chosen for the mantel test results)
### NOTE: you can exclude CpG>TpG mutation counts at this stage as long as they were separated in the preivoius step 0 script
# don't need a special step 0 script for that. want everything rescaled by human genome target sizes with Cpg and non cpg sites treated as diff for C.T but not for C.G and C.A
# all that is good.
# so now just exclude them from this spectrum and from aging signatures
require(dplyr)
require(tidyverse)
require(reshape2)
require(ggplot2)
require(sigfit)
require(tune)
require(RColorBrewer)
# for mantel test:
require(ape)
require(vegan)
require(compositions)
require(scales)
require(spgs) # for rev comp. don't use on 1mers
mantelTestPermutations=9999999

### setting some sigfit parameters ###
set.seed(42)


signatureColors_forManuscript=c("SBS1"=brewer.pal(name="Set1",n=8)[8],"SBS5" =brewer.pal(name="Set1",n=8)[7], "maternal_age"=brewer.pal(name="Set3",n=8)[4] ,"paternal_age"=brewer.pal(name="Set3",n=8)[5], "young_parent_13" =brewer.pal(name="Set3",n=8)[6],"novel (1-mer-CpG)"= brewer.pal(name="Set3",n=8)[3],"novel (3-mer)"=brewer.pal(name="Set3",n=8)[1])

signatureColors_forManuscript_niceNames=c("SBS1"=brewer.pal(name="Set1",n=8)[8],"SBS5" =brewer.pal(name="Set1",n=8)[7], "maternal age"=brewer.pal(name="Set3",n=8)[4] ,"paternal age"=brewer.pal(name="Set3",n=8)[5], "young parent" =brewer.pal(name="Set3",n=8)[6],"novel (1-mer-CpG)"= brewer.pal(name="Set3",n=8)[3],"novel (3-mer)"=brewer.pal(name="Set3",n=8)[1]) # this includes novel signatures

#modelColors=c("SBS1 + SBS5 (3-mer)"="#A65628","aging signatures (1-mer-CpG)"="#984EA3")



signatureColors_1mer=c("maternal_age"=brewer.pal(name="Set3",n=8)[4] ,"paternal_age"=brewer.pal(name="Set3",n=8)[5], "young_parent_13" =brewer.pal(name="Set3",n=8)[3] ) ### 20221012: want to switch these up so they better match jonsson et al: making mother colors red, father blue, young parent purple -- young parent kind of hard to see maybe make it orange?


#### SET NITER AND ADAPT DELTA ON COMMAND LINE

options(warn=1) # this has warnings print as they happen! aha! super important for differentiating between sigfit warnings ! otherwise they just print at the end and you can't tell what's what.
# this script should be run from the shell. 
# args = commandArgs(trailingOnly=TRUE)
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied", call.=FALSE)
# } 

#niter = as.numeric(args[1])


# or for testing:
niter=10000
# using default adapt de

print(paste0("niter=",niter))

sigfit_model_choice="multinomial" # (NMF) (look back on documentation and make sure we like this choice)
print(paste0("sigfit_model_choice=",sigfit_model_choice))

print("excluding CpG.TpG entirely from spectrum and aging signatures")
#### set up an output directory #######
todaysdate=format(Sys.Date(),"%Y%m%d")
flagOfAnalysis=paste0(todaysdate,".sigfitResults.",sigfit_model_choice,".niter.",niter,".notprojected.CpGTpGExcluded") # a label that will be in title of plots and outdir 

outdir=paste0("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/sigfit/SigfitRuns_notcountsplit_notprojected/",flagOfAnalysis,"/")

dir.create(outdir,showWarnings = F,recursive = T)


#### subsetting species to just the populations that were used in the mantel distance analyses ###

speciesToInclude=c("Mus_musculus_Mmd","Mus_spretus_Ms","brown_bear_ABC","polar_bear_PB","vaquita","fin_whale_GOC","humans_AFR","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","Pongo_abelii","Pongo_pygmaeus","wolves") 

#### set an order of species for plots #######
#speciesOrder = c("Mus_spretus_Ms","Mus_musculus_Mmd","Mus_musculus_Mmm","Mus_musculus_Mmc","Pongo_abelii","Pongo_pygmaeus","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","humans_AFR","humans_EUR","humans_SAS","humans_EAS","humans_AMR","wolves","ucla_wolves","brown_bear_ABC","brown_bear_EUR","polar_bear_PB","fin_whale_GOC","fin_whale_ENP","vaquita")

speciesOrder = c("Mus_spretus_Ms","Mus_musculus_Mmd","Pongo_abelii","Pongo_pygmaeus","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","humans_AFR","wolves","brown_bear_ABC","polar_bear_PB","fin_whale_GOC","vaquita")


speciesShortNames = data.frame(speciesOriginal=speciesOrder,speciesSuperShort=c("Algerian mouse","House mouse","Sum. Orang","Bor. Orang","Gorilla","Bonobo","Chimp","Human","Wolf","Brown bear","Polar bear","Fin whale","Vaquita"))

speciesOrderShortNames = c("Algerian mouse","House mouse","Sum. Orang","Bor. Orang","Gorilla","Bonobo","Chimp","Human","Wolf","Brown bear","Polar bear","Fin whale","Vaquita")

# read in phylogenetic distances #
phyloDistances = read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/compare_spectrumDist_toPhyloHammingDist/20221121.plots.SubsampledIndividualsNotProjection.epsilon.1.sepCpGs.no.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode.9999999.MantelPermutations/pairwise.phylogeneticDistance.raxmltree.txt",header=T) # these are pairwise cophenetic distances from raxml tree.


### read in a priori aging signatures ### 
agingSignatures_wCpG = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/jonsson_parental_aging/jonsson.agingsignatures.tables9.convertedtoproportionsinexcel.REVCOMP.REORDEREDTOMATCHMYDATA.YOUNGPARENTCALCAT13.notintercept.txt",header=T,sep="\t",row.names="signature")

# this removes the CpG.TpG column
agingSignatures_noCpG <- agingSignatures_wCpG %>%
  select(-CpG.TpG)

# then we need to renormalize without the CpG fraction:

agingSignatures_noCpG_renorm = agingSignatures_noCpG /rowSums(agingSignatures_noCpG)

# do some checks:
if(!dim(agingSignatures_noCpG_renorm)[2]==6){
  stop("CpGs not properly removed")
}

if(sum(rowSums(agingSignatures_noCpG_renorm)!=1)!=0){
  stop("didn't renormalize correctly")
}

agingSignatures = agingSignatures_noCpG_renorm # can now use this downstream 
####### read in cosmic signatures from sigfit #########
data(cosmic_signatures_v3.2)
cosmic_signatures_v3.2
# keeping these as human targets since I'm rescaling everything to human. but note if I return to using 'opportunities' then I need to rescale the cosmic signatures away from human targets. are fine for now though. (see github sigfit notes for details)

###### read in my spectra (not projected, not count-split, not down-sampled, YES rescaled by human target sizes) #####
# when reading in > gets converted to . so need to put back to match cosmic signatures (3mers only; 1mers don't match cosmic signatures)
# note that 1mers don't use > symbol.
# this will read in the appropriate sep by cpgs version of file
# READ IN THE FILE THAT HAS CpGs SEPARATED THEN GOING TO EXCLUDE THEM:
# DONT NEED TO MAKE A DIFFERENT INFILE 
spectra_1mer_rescaledByHumanTargets_wCpG <- read.table(paste0("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/sigfit/all_spectra_1mer_Matrix_notdownsampled.notprojected.RESCALEDBYHUMANTARGETS_READYFORSIGFIT.separateCpGs.yes.txt"),header = T,row.names = 1,sep="\t")

# EXCLUDE CPG.TPG:
spectra_1mer_rescaledByHumanTargets_noCpGTpG <- select(spectra_1mer_rescaledByHumanTargets_wCpG,-CpG.TpG)

# 3mers not affected by cpg separation choice
spectra_3mer_rescaledByHumanTargets <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/sigfit/all_spectra_3mer_Matrix_notdownsampled.notprojected.RESCALEDBYHUMANTARGETS_READYFORSIGFIT.txt",header = T,row.names = 1,sep="\t")

# when reading in > gets converted to . so need to put back to match cosmic signatures (3mers only; 1mers don't match cosmic signatures)
colnames(spectra_3mer_rescaledByHumanTargets) <- gsub("\\.",">",colnames(spectra_3mer_rescaledByHumanTargets))

head(spectra_3mer_rescaledByHumanTargets)

######### get 1mer versions of SBS5 ##############
# make1merVersionOFCosmic3merSigntureAndRevCompToMatchAgingSignatures <- function(nameofSBSSignature){
#   signature_melt <- melt(cosmic_signatures_v3.2[nameofSBSSignature,])
#   signature_melt$mutation_type <- rownames(signature_melt)
#   signature_melt$central1mer <- paste0(substr(signature_melt$mutation_type,2,2),".",substr(signature_melt$mutation_type,6,6))
#   signature_melt[signature_melt$mutation_type %in% c("TCG>TTG","ACG>ATG","GCG>GTG","CCG>CTG"),]$central1mer <- "CpG.TpG"
#   # need to specify cpgs
#   
#   signature_1MERVERSION <- signature_melt %>%
#     group_by(central1mer) %>%
#     summarise(sum_over3mers=sum(value)) %>%
#     ungroup()
#   
#   
#   # need to match names aging signatures so rev comp T.C --> A.G ; T.G --> A.C ; T.A --> A.T
#   signature_1MERVERSION$nameToMatchAging <- signature_1MERVERSION$central1mer
#   signature_1MERVERSION[signature_1MERVERSION$central1mer=="T.C",]$nameToMatchAging <- "A.G"
#   signature_1MERVERSION[signature_1MERVERSION$central1mer=="T.G",]$nameToMatchAging <- "A.C"
#   signature_1MERVERSION[signature_1MERVERSION$central1mer=="T.A",]$nameToMatchAging <- "A.T"
#   
#   signature_1MERVERSION_wide <- data.frame(pivot_wider(signature_1MERVERSION[,c("nameToMatchAging","sum_over3mers")],names_from=nameToMatchAging,values_from = sum_over3mers))
#   rownames(signature_1MERVERSION_wide) <- paste0(nameofSBSSignature," (1-mer)")
#   return(signature_1MERVERSION_wide)
# }

# 1mer SBS5 (with CpGTpG)
# SBS5_1merVersion_containsCpGTpG <- make1merVersionOFCosmic3merSigntureAndRevCompToMatchAgingSignatures("SBS5")
# 
# # exclude CpG.TpG and renormalize:
# # this removes the CpG.TpG column
# SBS5_1merVersion_noCpGTpG <- SBS5_1merVersion_containsCpGTpG %>%
#   select(-CpG.TpG)
# 
# # then we need to renormalize without the CpG fraction:
# SBS5_1merVersion_noCpGTpG_renorm = SBS5_1merVersion_noCpGTpG /rowSums(SBS5_1merVersion_noCpGTpG)

######### combine SBS5 1mer with aging signatures  ########

aging_signatures_noCpG_renorm_plusSBS5_noCpG_1mer <- bind_rows(agingSignatures_noCpG_renorm,SBS5_1merVersion_noCpGTpG_renorm)


####### subset to just chosen pop per species ####
spectra_1mer_rescaledByHumanTargets_noCpGTpG_subsetPops <- spectra_1mer_rescaledByHumanTargets_noCpGTpG[speciesToInclude,] # subset to the rows of the populations you're including (1 pop per species)

# and that 1mers are in same order as aging signatures
if(sum(colnames(spectra_1mer_rescaledByHumanTargets_noCpGTpG_subsetPops)!=colnames(agingSignatures))!=0){
  stop("test_1mer isn't in same order as aging signatures!")
}

spectra_3mer_rescaledByHumanTargets_subsetPops <- spectra_3mer_rescaledByHumanTargets[speciesToInclude,]

# check 3mers are in same order as cosmic
if(sum(colnames(spectra_3mer_rescaledByHumanTargets_subsetPops)!=colnames(cosmic_signatures_v3.2))!=0){
  stop("spectra_3mer_rescaledByHumanTargets isn't in same order as COSMIC!")
}

################## some sigfit wrapper functions #################
fitSignaturesToData <- function(signaturesToFit,dataToFitSignaturesTo,niter,sigfit_model_choice){
  # make sure are in same order
  if(sum(colnames(dataToFitSignaturesTo)!=colnames(signaturesToFit))!=0){
    stop("dataset isn't in same order as signatures!")
  }
  
  sigfit_result_fittingSignaturesToData <- fit_signatures(counts = dataToFitSignaturesTo,signatures = signaturesToFit,iter=niter,model=sigfit_model_choice)
  
  
  return(sigfit_result_fittingSignaturesToData)
}

# fit signatures then extract novel signature
fitSignaturesThenExtractNovel <- function(signaturesToFit,dataToFitSignaturesTo,niter,sigfit_model_choice,num_novel){
  # make sure are in same order
  if(sum(colnames(dataToFitSignaturesTo)!=colnames(signaturesToFit))!=0){
    stop("dataset isn't in same order as signatures!")
  }
  
  sigfit_result_fittingSignaturesToDataThenExtracting <- fit_extract_signatures(counts = dataToFitSignaturesTo,signatures = signaturesToFit,iter=niter,model=sigfit_model_choice,num_extra_sigs=num_novel,control=list(adapt_delta=0.99)) # needs adapt_delta to avoid divergent transitions
  
  return(sigfit_result_fittingSignaturesToDataThenExtracting)
  
}

# function get cosine similarities between observed and reconstructed spectra for every species (slow, based on sigfit code) for a single model
GetCosineSimilarity_allSamples <- function(sigfit_results){
  counts=sigfit_results$data$counts_int
  samples=rownames(counts)
  reconstructions=retrieve_pars(sigfit_results,"reconstructions")$mean # get mean
  NSAMP=length(samples)
  allCos=data.frame()
  for(i in seq(1,NSAMP)){
    samplename=samples[i]
    cossimilarity=cosine_sim(counts[i,],reconstructions[i,]) # get cosine similarity (this is how sigfit does it in plotting code -- not very efficient, but oh well )
    cosdf = data.frame(species=samplename,cosine_similarity=cossimilarity)
    allCos=bind_rows(allCos,cosdf)
  }
  return(allCos)
}

# calculate residuals: . also want to use proportions not counts so all species are on same scale! 
# 
# calculateSigfitResiduals_Counts <- function(sigfitResults){
#   reconstructions <- retrieve_pars(sigfitResults,par="reconstructions")$mean
#   data <- sigfitResults$data$counts_real
#   residuals <-  data - reconstructions # subtract recon from data (works bc are same dimensions)
#   residuals_melt <- melt(residuals)
#   colnames(residuals_melt) <- c("species","mutationLabel","residual_data_minus_model")
#   reconstructions_melt <- melt(reconstructions)
#   colnames(reconstructions_melt) <- c("species","mutationLabel","reconstructionValue")
#   # merge:
#   # also put data in too:
#   data_melt <- melt(data)
#   colnames(data_melt) <- c("species","mutationLabel","data_counts")
#   combo_residuals_reconstructions <- merge(residuals_melt,reconstructions_melt,by=c("species","mutationLabel"))
#   combo_residuals_reconstructions_data <- merge(combo_residuals_reconstructions,data_melt,by=c("species","mutationLabel"))
#   
#   # rescale to be  (obs - exp ) / exp aka (data - recon) / recon
#   #combo_residuals_reconstructions_data$residual_data_minus_model_RescaledByDivByModel <- (combo_residuals_reconstructions_data$residual_data_minus_model / combo_residuals_reconstructions_data$reconstructionValue)
#   
#   return(combo_residuals_reconstructions_data)
# }
# USE proportions!
calculateSigfitResiduals <- function(sigfitResults){
  reconstructions <- retrieve_pars(sigfitResults,par="reconstructions")$mean
  # get proportions:
  reconstructions_prop <- reconstructions/rowSums(reconstructions)
  
  data <- sigfitResults$data$counts_real
  data_prop <- data/rowSums(data)
  
  residuals_prop <-  data_prop - reconstructions_prop # subtract recon from data (works bc are same dimensions)
  residuals_melt <- melt(residuals_prop)
  colnames(residuals_melt) <- c("species","mutationLabel","residual_data_minus_model")
  reconstructions_melt <- melt(reconstructions_prop)
  colnames(reconstructions_melt) <- c("species","mutationLabel","reconstruction_proportions")
  # merge:
  # also put data in too:
  data_melt <- melt(data_prop)
  colnames(data_melt) <- c("species","mutationLabel","data_proportions")
  combo_residuals_reconstructions <- merge(residuals_melt,reconstructions_melt,by=c("species","mutationLabel"))
  combo_residuals_reconstructions_data <- merge(combo_residuals_reconstructions,data_melt,by=c("species","mutationLabel"))
  
  # rescale to be  (obs - exp ) / exp aka (data - recon) / recon
  #combo_residuals_reconstructions_data$residual_data_minus_model_RescaledByDivByModel <- (combo_residuals_reconstructions_data$residual_data_minus_model / combo_residuals_reconstructions_data$reconstructionValue)
  
  return(combo_residuals_reconstructions_data)
}

# plot residuals for a single model:
# changing this to now plot the unscaled residuals (both are output by prev function)
makeSigfitResidualsPlot <- function(sigfitResults,runlabel,outdir){
  
  combo_residuals_reconstructions_data = calculateSigfitResiduals(sigfitResults)
  
  dir.create(outdir,showWarnings = F,recursive = T)
  # the gsub("\\.",">") makes a . turn into a > for label
  residualsPlot <- ggplot(combo_residuals_reconstructions_data,aes(x=gsub("\\.",">",mutationLabel),y=residual_data_minus_model))+
    geom_hline(yintercept=0,linetype="dashed",color="darkgray")+
    #geom_point()+
    geom_boxplot()+
    #facet_wrap(~mutationLabel,scales="free_x",nrow=1)+
    ylab("residual\n(data proportion - reconstruction proportion)")+
    xlab("")+
    theme_bw()+
    theme(legend.title = element_blank(),text=element_text(size=12))+
    ggtitle("Residuals comparison")+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("residuals (data proportion - reconstruction proportion)\n",runlabel))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8))
  residualsPlot
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".pdf"),residualsPlot,width=12,height=5)
  return(residualsPlot)
}

# plot residuals for a single model faceted by species 
makeSigfitResidualsPlot_perSpp <- function(sigfitResults,runlabel,outdir){
  combo_residuals_reconstructions_data = calculateSigfitResiduals(sigfitResults)
  
  dir.create(outdir,showWarnings = F,recursive = T)
  # the gsub("\\.",">") makes a . turn into a > for label
  residualsPlot <- ggplot(combo_residuals_reconstructions_data,aes(x=gsub("\\.",">",mutationLabel),y=residual_data_minus_model))+
    #geom_hline(yintercept=0,linetype="dashed",color="darkgray")+
    geom_point()+
    #geom_boxplot()+
    facet_wrap(~species,scales="free_x",nrow=1)+
    ylab("residual\n(data proportion - reconstruction proportion)")+
    xlab("")+
    theme_bw()+
    theme(legend.title = element_blank(),text=element_text(size=12))+
    ggtitle("Residuals comparison")+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("residuals (data proportion - reconstruction proportion)/(reconstruction)\n",runlabel))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8))
  residualsPlot
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".PERSPP.pdf"),residualsPlot,width=12,height=5)
  return(residualsPlot)
}

makeSigfitResidualsPlot_perSpp_BoxPlotOverTypes <- function(sigfitResults,runlabel,outdir,speciesOrder){
  combo_residuals_reconstructions_data = calculateSigfitResiduals(sigfitResults)
  
  dir.create(outdir,showWarnings = F,recursive = T)
  # the gsub("\\.",">") makes a . turn into a > for label
  
  combo_residuals_reconstructions_data$species <- factor(combo_residuals_reconstructions_data$species, levels=speciesOrder)
  
  residualsPlot <- ggplot(combo_residuals_reconstructions_data,aes(x=species,y=residual_data_minus_model))+
    geom_boxplot()+
    ylab("residual\n(data proportion - reconstruction proportion)")+
    xlab("")+
    theme_bw()+
    theme(legend.title = element_blank(),text=element_text(size=12))+
    ggtitle("Residuals comparison")+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("residuals (data proportion - reconstruction proportion)\n",runlabel))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8))
  residualsPlot
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".PERSPP.BOXPLOTOVERMUTATIONTYPES.pdf"),residualsPlot,width=12,height=5)
  return(residualsPlot)
}

makeSigfitResidualsPlot_perSpp_PointsPerType <- function(sigfitResults,runlabel,outdir,speciesOrder,sizeOfText){
  combo_residuals_reconstructions_data = calculateSigfitResiduals(sigfitResults)
  
  dir.create(outdir,showWarnings = F,recursive = T)
  # the gsub("\\.",">") makes a . turn into a > for label
  
  combo_residuals_reconstructions_data$species <- factor(combo_residuals_reconstructions_data$species, levels=speciesOrder)
  
  residualsPlot <- ggplot(combo_residuals_reconstructions_data,aes(x=species,y=residual_data_minus_model))+
    #geom_point()+
    geom_text(aes(label=mutationLabel),size=sizeOfText)+
    ylab("residual\n(data proportion - reconstruction proportion)")+
    xlab("")+
    theme_bw()+
    theme(legend.title = element_blank(),text=element_text(size=12))+
    ggtitle("Residuals comparison")+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("residuals (data proportion - reconstruction proportion)\n",runlabel))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8))
  residualsPlot
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".PERSPP.POINTSFORMUTATIONTYPES.pdf"),residualsPlot,width=12,height=5)
  return(residualsPlot)
}

## try new residual boxplot

ResidualsBoxPlot_WithOutlierLabels <- function(sigfitResultsFull,runlabel,outdir,speciesOrder,textSize){
  sigfitResultsFull_RESIDUALS = calculateSigfitResiduals(sigfitResultsFull)
  
  sigfitResultsFull_RESIDUALS$species <- factor(sigfitResultsFull_RESIDUALS$species, levels=speciesOrder)
  
  # identify outliers:
  # https://www.r-bloggers.com/2022/08/how-to-label-outliers-in-boxplots-in-ggplot2/
  findoutlier <- function(x) {
    return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x)) # going to try making this even mroe extreme to only label furthest outliers ? 
  }
  
  sigfitResultsFull_RESIDUALS <- sigfitResultsFull_RESIDUALS %>%
    group_by(species) %>%
    mutate(outlier = ifelse(findoutlier(residual_data_minus_model), "outlier", "")) %>%
    ungroup()
  
  
  residualsPlot <- ggplot(sigfitResultsFull_RESIDUALS,aes(x=species,y=residual_data_minus_model))+
    #geom_hline(yintercept = 0,color="darkgray")+
    geom_boxplot(color="darkgray")+
    geom_text(data=sigfitResultsFull_RESIDUALS[sigfitResultsFull_RESIDUALS$outlier=="outlier",],aes(label=gsub("\\.",">",mutationLabel)),size=textSize)+
    ylab("residual\n(data proportion - reconstruction proportion)")+
    xlab("")+
    theme_bw()+
    theme(legend.title = element_blank(),text=element_text(size=12))+
    ggtitle("Residuals comparison")+
    ggtitle(paste0("residuals (data proportion - reconstruction proportion)\n",runlabel))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8))
  
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".PERSPP.BOXPLOT.LabelOutliers.pdf"),residualsPlot,width=12,height=5)

    return(residualsPlot)
}

# side by side models:
ResidualsBoxPlot_WithOutlierLabels_MultipleModels <- function(ListOfsigfitResultsFull,runlabel,outdir,speciesOrderShortNames,textSize,modelColors,speciesShortNames){
  sigfitResultsFull_RESIDUALS_List = lapply(ListOfsigfitResultsFull,calculateSigfitResiduals)
  
  # collapse
  sigfitResultsFull_RESIDUALS_df <- bind_rows(sigfitResultsFull_RESIDUALS_List,.id="id")
  
  sigfitResultsFull_RESIDUALS_df$species <- factor(sigfitResultsFull_RESIDUALS_df$species, levels=speciesOrder)
  
  # identify outliers:
  # https://www.r-bloggers.com/2022/08/how-to-label-outliers-in-boxplots-in-ggplot2/
  findoutlier <- function(x) {
    return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
  }
  
  sigfitResultsFull_RESIDUALS_df <- sigfitResultsFull_RESIDUALS_df %>%
    group_by(species,id) %>%
    mutate(outlier = ifelse(findoutlier(residual_data_minus_model), "outlier", "")) %>%
    ungroup()
  
  sigfitResultsFull_RESIDUALS_df <- merge(sigfitResultsFull_RESIDUALS_df,speciesShortNames, by.x="species",by.y = "speciesOriginal",all = T)
  
  sigfitResultsFull_RESIDUALS_df$speciesSuperShort <- factor(sigfitResultsFull_RESIDUALS_df$speciesSuperShort,levels=speciesOrderShortNames)
  
  residualsPlot <- ggplot(sigfitResultsFull_RESIDUALS_df,aes(x=id,y=residual_data_minus_model,color=id))+
    geom_boxplot(position=position_dodge(width=1))+
    geom_text(data=sigfitResultsFull_RESIDUALS_df[sigfitResultsFull_RESIDUALS_df$outlier=="outlier",],aes(label=gsub("\\.",">",mutationLabel),group=id),size=textSize,position=position_dodge(width=1),color="black",vjust=-0.5)+
    ylab("residual\n(data proportion - reconstruction proportion)")+
    facet_wrap(~speciesSuperShort,nrow=1)+
    xlab("")+
    theme_bw()+
    theme(legend.title = element_blank(),text=element_text(size=12))+
    ggtitle("Residuals comparison")+
    #ggtitle(paste0("residuals (data proportion - reconstruction proportion)\n",runlabel))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8))+
    geom_hline(yintercept = 0,color="darkgray",linetype="dashed")+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_rect(fill="lightblue"),legend.position="bottom",legend.title = element_blank())+
    scale_color_manual(values=modelColors)
  
    
  residualsPlot
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".PERSPP.BOXPLOT.LabelOutliers.SideBySide.pdf"),residualsPlot,width=14,height=6)
  
  return(residualsPlot)
}


ResidualsBoxPlotPERMUTATIONTYPE_MultipleModels <- function(ListOfsigfitResultsFull,runlabel,outdir,speciesOrderShortNames,textSize,modelColors,speciesShortNames){
  sigfitResultsFull_RESIDUALS_List = lapply(ListOfsigfitResultsFull,calculateSigfitResiduals)
  
  # collapse
  sigfitResultsFull_RESIDUALS_df <- bind_rows(sigfitResultsFull_RESIDUALS_List,.id="id")
  
  sigfitResultsFull_RESIDUALS_df$species <- factor(sigfitResultsFull_RESIDUALS_df$species, levels=speciesOrder)
  
  sigfitResultsFull_RESIDUALS_df <- merge(sigfitResultsFull_RESIDUALS_df,speciesShortNames, by.x="species",by.y = "speciesOriginal",all = T)
  
  sigfitResultsFull_RESIDUALS_df$speciesSuperShort <- factor(sigfitResultsFull_RESIDUALS_df$speciesSuperShort,levels=speciesOrderShortNames)
  
  residualsPlot <- ggplot(sigfitResultsFull_RESIDUALS_df,aes(x=gsub("\\.",">",mutationLabel),y=residual_data_minus_model))+
    geom_boxplot(aes(fill=id),alpha=0.8,outlier.size = 0) + # setting outlier size= 0 so that don't get plotted twice 
    geom_point(position=position_jitterdodge(jitter.width = .1),aes(group=id,fill=id),size=.5)+
    ylab("residual\n(data proportion - reconstruction proportion)")+
    xlab("")+
    theme_bw()+
    theme(legend.title = element_blank(),text=element_text(size=12))+
    ggtitle("Residuals comparison")+
    #ggtitle(paste0("residuals (data proportion - reconstruction proportion)\n",runlabel))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=16))+
    geom_hline(yintercept = 0,color="darkgray",linetype="dashed")+
    theme(legend.position="bottom",legend.title = element_blank())+
    scale_color_manual(values=modelColors)+
    scale_fill_manual(values=modelColors)
  
  
  
  residualsPlot
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".PERMUTATIONTYPE.BOXPLOT.SideBySide.pdf"),residualsPlot,width=5,height=5)
  
  return(residualsPlot)
}


# for 3mer residuals, want to facet by 1mer spectrum 
ResidualsBoxPlotPERMUTATIONTYPE_WithOutlierLabels_MultipleModels_SEP3merBy1Mer <- function(ListOfsigfitResultsFull,runlabel,outdir,speciesOrderShortNames,textSize,modelColors,speciesShortNames){
  sigfitResultsFull_RESIDUALS_List = lapply(ListOfsigfitResultsFull,calculateSigfitResiduals)
  
  # collapse
  sigfitResultsFull_RESIDUALS_df <- bind_rows(sigfitResultsFull_RESIDUALS_List,.id="id")
  
  sigfitResultsFull_RESIDUALS_df$species <- factor(sigfitResultsFull_RESIDUALS_df$species, levels=speciesOrder)
  
  sigfitResultsFull_RESIDUALS_df <- merge(sigfitResultsFull_RESIDUALS_df,speciesShortNames, by.x="species",by.y = "speciesOriginal",all = T)
  
  sigfitResultsFull_RESIDUALS_df$speciesSuperShort <- factor(sigfitResultsFull_RESIDUALS_df$speciesSuperShort,levels=speciesOrderShortNames)
  
  ### get 1mer information: 
  sigfitResultsFull_RESIDUALS_df$mutationLabel_1mer <- paste0(substr(sigfitResultsFull_RESIDUALS_df$mutationLabel,2,2),">",substr(sigfitResultsFull_RESIDUALS_df$mutationLabel,6,6))
  
  sigfitResultsFull_RESIDUALS_df$mutationLabelRevComp <- as.character(sigfitResultsFull_RESIDUALS_df$mutationLabel)
  
  sigfitResultsFull_RESIDUALS_df$ancestral1mer <- substr(sigfitResultsFull_RESIDUALS_df$mutationLabel,2,2)
  
  sigfitResultsFull_RESIDUALS_df[sigfitResultsFull_RESIDUALS_df$ancestral1mer=="T",]$mutationLabelRevComp <- toupper(paste0(reverseComplement(substr(sigfitResultsFull_RESIDUALS_df[sigfitResultsFull_RESIDUALS_df$ancestral1mer=="T",]$mutationLabel,1,3)),">",reverseComplement(substr(sigfitResultsFull_RESIDUALS_df[sigfitResultsFull_RESIDUALS_df$ancestral1mer=="T",]$mutationLabel,5,7))))
  
  # flip rev comp to match other plot:
  sigfitResultsFull_RESIDUALS_df[sigfitResultsFull_RESIDUALS_df$mutationLabel_1mer=="T>A",]$mutationLabel_1mer <- "A>T"
  sigfitResultsFull_RESIDUALS_df[sigfitResultsFull_RESIDUALS_df$mutationLabel_1mer=="T>C",]$mutationLabel_1mer <- "A>G"
  sigfitResultsFull_RESIDUALS_df[sigfitResultsFull_RESIDUALS_df$mutationLabel_1mer=="T>G",]$mutationLabel_1mer <- "A>C"
  

  ###### find extreme outliers: 
  # https://www.r-bloggers.com/2022/08/how-to-label-outliers-in-boxplots-in-ggplot2/
  findoutlier <- function(x) {
    return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) +1.5*IQR(x)) # going to try making this even mroe extreme to only label furthest outliers ? original was 1.5*iqr
  }

  sigfitResultsFull_RESIDUALS_df <- sigfitResultsFull_RESIDUALS_df %>%
    group_by(species,mutationLabel_1mer,id) %>%
    mutate(outlier = ifelse(findoutlier(residual_data_minus_model), "outlier", "")) %>%
    ungroup()
 
  
  # for the 1mers that are T>A T>C and T>G , maybe rev comp so they match? maybe just do that in AI? but want to get in same order as other plots. 
  
  
  residualsPlot <- ggplot(sigfitResultsFull_RESIDUALS_df,aes(x=mutationLabel_1mer,y=residual_data_minus_model))+
    geom_boxplot(aes(fill=id),alpha=0.8,outlier.shape = NA) + # setting outlier shape to NA so they don't get double plotted
    geom_point(position=position_jitterdodge(jitter.width = .1),aes(group=id,fill=id),size=.5)+
    ylab("residual\n(data proportion - reconstruction proportion)")+
    xlab("")+
    theme_bw()+
    geom_text(data=sigfitResultsFull_RESIDUALS_df[sigfitResultsFull_RESIDUALS_df$outlier=="outlier" & abs(sigfitResultsFull_RESIDUALS_df$residual_data_minus_model)>0.015,],aes(label=paste0(speciesSuperShort,":",gsub("\\.",">",mutationLabelRevComp)),group=id),size=textSize,position=position_dodge(width=1),color="black",vjust=-0.5)+
    theme(legend.title = element_blank(),text=element_text(size=12))+
    ggtitle("Residuals comparison")+
    #ggtitle(paste0("residuals (data proportion - reconstruction proportion)\n",runlabel))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=16))+
    geom_hline(yintercept = 0,color="darkgray",linetype="dashed")+
    theme(legend.position="bottom",legend.title = element_blank()) +
    scale_color_manual(values=modelColors)+
    scale_fill_manual(values=modelColors)
  
  
  
  residualsPlot
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".PERMUTATIONTYPE.BOXPLOT.Facet3merBy1mer.SideBySide.pdf"),residualsPlot,width=5,height=5)
  
  return(residualsPlot)
}






getReconstructionNextToDataTable <- function(sigfitResults){
  reconstructions <- retrieve_pars(sigfitResults,par="reconstructions")$mean
  data <- sigfitResults$data$counts_real
  
  # convert both to proportions:
  # this makes each species' row sum to 1
  reconstructions_proportions = reconstructions/rowSums(reconstructions)
  data_proportions = data/rowSums(data)
  
  reconstructions_proportions_melt <- melt(reconstructions_proportions)
  colnames(reconstructions_proportions_melt) <- c("species","mutationLabel","reconstructionValue")
  # merge:
  # also put data in too:
  data_proportions_melt <- melt(data_proportions)
  colnames(data_proportions_melt) <- c("species","mutationLabel","counts_real")

  combo_reconstructions_data <- merge(reconstructions_proportions_melt,data_proportions_melt,by=c("species","mutationLabel"))
  combo_reconstructions_data <- melt(combo_reconstructions_data)
  
  return(combo_reconstructions_data)
  }
# calculate sum of squared errors SSE = SUM[(model - data)^2]

calculateSSE <- function(sigfitResults){
  reconstructions <- retrieve_pars(sigfitResults,par="reconstructions")$mean
  data <- sigfitResults$data$counts_real
  # convert both to proportions:
  # this makes each species' row sum to 1
  reconstructions_proportions = reconstructions/rowSums(reconstructions)
  data_proportions = data/rowSums(data)
  # SSE: https://www.statology.org/sst-ssr-sse-in-r/
  SqErrorMatrix <- (reconstructions_proportions-data_proportions)^2 # for each cell of matrix,  model - data proportions squared
  # get sum for each species:
  SSE_perSpecies <- rowSums(SqErrorMatrix)
  SSE_perSpecies_melt <- melt(SSE_perSpecies) %>% rownames_to_column(var="species") %>% rename(SSE=value)
  
  return(SSE_perSpecies_melt)
}


########## FIT AGING SIGNATURES TO SPECIES 1MER SPECTRA #############
fit_aging_signatures <- fitSignaturesToData(dataToFitSignaturesTo = spectra_1mer_rescaledByHumanTargets_noCpGTpG_subsetPops,signaturesToFit = agingSignatures,niter=niter,sigfit_model_choice = sigfit_model_choice) 

plot_all(fit_aging_signatures,out_path = paste0(outdir,"/fit_aging_signatures","/")) 
# get cosine similarity:
fit_aging_signatures_COSINE_SIMILARITIES_df <- GetCosineSimilarity_allSamples(fit_aging_signatures)
fit_aging_signatures_COSINE_SIMILARITIES_df$nmodel <- "aging signatures (1-mer-CpG)"
fit_aging_signatures_COSINE_SIMILARITIES_df$nicer_nmodel_name <- "aging signatures (1-mer-CpG)"

# get exposures:
fit_aging_signatures_EXPOSURES_df <- retrieve_pars(fit_aging_signatures,par="exposures")$mean

# plot residuals:
fit_aging_signatures_ResidualsPlot <- makeSigfitResidualsPlot(fit_aging_signatures,runlabel = "aging signatures (1-mer-CpG)",paste0(outdir,"fit_aging_signatures","/"))

# make some tweaks to 1mer residual plot that wouldn't work well for 3mer plot
fit_aging_signatures_ResidualsPlot <- fit_aging_signatures_ResidualsPlot+theme(axis.text.x = element_text(size=16))+ggtitle("")

# save a smaller version of the figure:
ggsave(paste0(outdir,"agingSignatureResiduals.smaller.pdf"),fit_aging_signatures_ResidualsPlot,height=4,width=6)

# facet per spp:
makeSigfitResidualsPlot_perSpp(fit_aging_signatures,runlabel = "aging signatures (1-mer-CpG)",paste0(outdir,"fit_aging_signatures","/"))

# facet per spp with boxplot over all types:
makeSigfitResidualsPlot_perSpp_BoxPlotOverTypes(fit_aging_signatures,runlabel = "aging signatures (1-mer-CpG)",paste0(outdir,"fit_aging_signatures","/"),speciesOrder)

# make pts per type:
makeSigfitResidualsPlot_perSpp_PointsPerType(fit_aging_signatures,runlabel = "aging signatures (1-mer-CpG)",paste0(outdir,"fit_aging_signatures","/"),speciesOrder,sizeOfText =3)

ResidualsBoxPlot_WithOutlierLabels(fit_aging_signatures,runlabel = "aging signatures (1-mer-CpG)",paste0(outdir,"fit_aging_signatures","/"),speciesOrder,textSize  =3)

# sse per species:
SSE_fit_aging_signatures <- calculateSSE(fit_aging_signatures)

SSEPlot_agingsignatures <- ggplot(SSE_fit_aging_signatures,aes(x=species,y=SSE))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))
SSEPlot_agingsignatures
ggsave(paste0(outdir,"SSEPlot_agingsignatures.pdf"),SSEPlot_agingsignatures)


reconstruction_next_to_data_agingSignaturesTable <- getReconstructionNextToDataTable(fit_aging_signatures)

head(reconstruction_next_to_data_agingSignaturesTable)

reconstruction_next_to_data_agingSignaturesPlot <- ggplot(reconstruction_next_to_data_agingSignaturesTable,aes(x=mutationLabel,y=value,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(~species,scales="free_y")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))
reconstruction_next_to_data_agingSignaturesPlot



ggsave(paste0(outdir,"reconstruction_next_to_data_agingSignatures.pdf"),reconstruction_next_to_data_agingSignaturesPlot,width=14,height=7)

# plot just humans:
reconstruction_next_to_data_agingSignaturesPlot_humanonly <- ggplot(reconstruction_next_to_data_agingSignaturesTable[reconstruction_next_to_data_agingSignaturesTable$species=="humans_AFR",],aes(x=mutationLabel,y=value,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(~species,scales="free_y")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))
reconstruction_next_to_data_agingSignaturesPlot_humanonly

ggsave(paste0(outdir,"reconstruction_next_to_data_agingSignatures.humans.only.pdf"),reconstruction_next_to_data_agingSignaturesPlot_humanonly,width=7,height=5)

#### plot with signatures next to reconstructions as well: ####

agingSignatures_melt <- agingSignatures %>%
  rownames_to_column(var="variable") %>%
  melt(.,value.name = "value",variable.name=c("mutationLabel"))

#agingSignatures_melt$label <- "signatures"
#reconstruction_next_to_data_agingSignaturesTable$label <- "data/models"

combo_ReconstructionsAndSignatures_agingModel <- bind_rows(agingSignatures_melt,reconstruction_next_to_data_agingSignaturesTable[reconstruction_next_to_data_agingSignaturesTable$species=="humans_AFR",])


combo_ReconstructionsAndSignatures_agingModel$variable <- factor(combo_ReconstructionsAndSignatures_agingModel$variable,levels=c("counts_real","reconstructionValue","maternal_age","paternal_age","young_parent_13"))

plotSignaturesAndReconstructionsAndEmpiricalAlltogether_agingSignatures <- ggplot(combo_ReconstructionsAndSignatures_agingModel,aes(x=mutationLabel,y=value,fill=variable))+
  geom_col(position="dodge")+
  theme_bw()+
  scale_fill_manual(values=c(counts_real="black",reconstructionValue="darkgoldenrod1",signatureColors_forManuscript[c("maternal_age","paternal_age","young_parent_13")]))+
  theme(axis.text.x = element_text(angle=45,hjust=1))

plotSignaturesAndReconstructionsAndEmpiricalAlltogether_agingSignatures
# plot with signatures next to it as well:
ggsave(paste0(outdir,"plotSignaturesAndReconstructionsAndEmpiricalAlltogether_agingSignatures.pdf"),plotSignaturesAndReconstructionsAndEmpiricalAlltogether_agingSignatures,width=9,height=5)

########## FIT AGING+NOVEL ############

fit_aging_signatures_PlusNovel <- fitSignaturesThenExtractNovel(dataToFitSignaturesTo = spectra_1mer_rescaledByHumanTargets_noCpGTpG_subsetPops,signaturesToFit = agingSignatures,niter=niter,sigfit_model_choice = sigfit_model_choice,num_novel = 1) 

plot_all(fit_aging_signatures_PlusNovel,out_path = paste0(outdir,"/fit_aging_signatures_PlusNovel","/")) 
# get cosine similarity:
fit_aging_signatures_PlusNovel_COSINE_SIMILARITIES_df <- GetCosineSimilarity_allSamples(fit_aging_signatures_PlusNovel)
fit_aging_signatures_PlusNovel_COSINE_SIMILARITIES_df$nmodel <- "aging signatures + novel (1-mer-CpG)"
fit_aging_signatures_PlusNovel_COSINE_SIMILARITIES_df$nicer_nmodel_name <- "aging signatures + novel (1-mer-CpG)"

# get exposures:
fit_aging_signatures_PlusNovel_EXPOSURES_df <- retrieve_pars(fit_aging_signatures_PlusNovel,par="exposures")$mean

# plot residuals:
fit_aging_signatures_PlusNovel_ResidualsPlot <- makeSigfitResidualsPlot(fit_aging_signatures_PlusNovel,runlabel = "aging signatures+novel (1-mer-CpG)",paste0(outdir,"fit_aging_signatures_PlusNovel","/"))

# make some tweaks to 1mer residual plot that wouldn't work well for 3mer plot
fit_aging_signatures_PlusNovel_ResidualsPlot <- fit_aging_signatures_PlusNovel_ResidualsPlot+theme(axis.text.x = element_text(size=16))+ggtitle("")

# save a smaller version of the figure:
ggsave(paste0(outdir,"fit_aging_signatures_PlusNovel/agingSignatureResiduals.smaller.pdf"),fit_aging_signatures_PlusNovel_ResidualsPlot,height=4,width=6)

# facet per spp:
makeSigfitResidualsPlot_perSpp(fit_aging_signatures_PlusNovel,runlabel = "aging signatures+novel (1-mer-CpG)",paste0(outdir,"fit_aging_signatures_PlusNovel","/"))

# make pts per type:
makeSigfitResidualsPlot_perSpp_PointsPerType(fit_aging_signatures_PlusNovel,runlabel = "aging signatures+novel (1-mer-CpG)",paste0(outdir,"fit_aging_signatures_PlusNovel/"),speciesOrder,sizeOfText =3)

## plot novel signature:
fit_aging_signatures_PlusNovel_signatures <- retrieve_pars(fit_aging_signatures_PlusNovel,par = "signatures")$mean
fit_aging_signatures_PlusNovel_signatures # A, B, C are mat, patn, yp13, D is novel signature

novelSignaturePlot <- ggplot(melt(fit_aging_signatures_PlusNovel_signatures["Signature D",]),aes(x=variable,y=value))+
  geom_col()+
  theme_bw()

ggsave(paste0(outdir,"fit_aging_signatures_PlusNovel/novelSiganturePlot.pdf"),novelSignaturePlot)
####### cut out: aging+ sbs5 ############
########## FIT SBS1+5 SIGNATURES TO SPECIES 3MER SPECTRA #############
# this isn't affected by cpg exclusion of 1mers -- still has cpg 3mers present (SBS1 being fit)
fit_SBS1_SBS5_signatures <- fitSignaturesToData(dataToFitSignaturesTo = spectra_3mer_rescaledByHumanTargets_subsetPops,signaturesToFit = cosmic_signatures_v3.2[c("SBS1","SBS5"),],niter=niter,sigfit_model_choice = sigfit_model_choice) 
plot_all(fit_SBS1_SBS5_signatures,out_path = paste0(outdir,"/fit_SBS1_SBS5_signatures","/")) 


# get cosine similarity:
fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df <- GetCosineSimilarity_allSamples(fit_SBS1_SBS5_signatures)
fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df$nmodel <- "SBS1 + SBS5 (3-mer)"
fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df$nicer_nmodel_name <- "SBS1 + SBS5 (3-mer)"

# get exposures
fit_SBS1_SBS5_signatures_EXPOSURES_df <- retrieve_pars(fit_SBS1_SBS5_signatures,par="exposures")$mean

# plot residuals
fit_SBS1_SBS5_signatures_ResidualsPlot <- makeSigfitResidualsPlot(fit_SBS1_SBS5_signatures,runlabel = "SBS1 + SBS5 (3-mer)",paste0(outdir,"fit_SBS1_SBS5_signatures","/"))

fit_SBS1_SBS5_signatures_ResidualsPlot

# facet per spp:
makeSigfitResidualsPlot_perSpp(fit_SBS1_SBS5_signatures,runlabel = "SBS1 + SBS5 (3-mer)",paste0(outdir,"fit_SBS1_SBS5_signatures","/"))

# points per type:
makeSigfitResidualsPlot_perSpp_PointsPerType(fit_SBS1_SBS5_signatures,runlabel = "SBS1 + SBS5 (3-mer)",paste0(outdir,"fit_SBS1_SBS5_signatures","/"),sizeOfText = 1)




########### FIT SBS1+SBS5+NOVEL ############
#*SLOW* because doing 2*niter iterations (20K)
fit_SBS1_SBS5_PlusNovel <- fitSignaturesThenExtractNovel(dataToFitSignaturesTo = spectra_3mer_rescaledByHumanTargets_subsetPops,signaturesToFit = cosmic_signatures_v3.2[c("SBS1","SBS5"),],niter=2*niter,sigfit_model_choice = sigfit_model_choice,num_novel=1) 
# make sure doesn't throw divergent transition errors 
plot_all(fit_SBS1_SBS5_PlusNovel,out_path = paste0(outdir,"/fit_SBS1_SBS5_signatures_PlusNovel","/")) 

# get cosine similarity:
fit_SBS1_SBS5_PlusNovel_signatures_COSINE_SIMILARITIES_df <- GetCosineSimilarity_allSamples(fit_SBS1_SBS5_PlusNovel)
fit_SBS1_SBS5_PlusNovel_signatures_COSINE_SIMILARITIES_df$nmodel <- "SBS1 + SBS5 + novel (3-mer)"
fit_SBS1_SBS5_PlusNovel_signatures_COSINE_SIMILARITIES_df$nicer_nmodel_name <- "SBS1 + SBS5 + novel (3-mer)"


# Warning: There were 4659 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#Warning: Examine the pairs() plot to diagnose sampling problems
# running for enough iterations and upping adapt delta to 0.99 got rid of divergent transition warnings 
########### plot cosine similarities of reconstructions vs cophenetic distance from humans #####
head(phyloDistances)
# file with nice labels:

phyloDistances_fromHuman <- phyloDistances[phyloDistances$Sp1 =="Homo_sapiens" | phyloDistances$Sp2=="Homo_sapiens",]

# species codes to get from latin name to label:
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/analyses/compare_spectrumDist_toPhyloHammingDist/speciesCodes.txt",header = T)
# subset to set I'm using
speciesCodes = speciesCodes[speciesCodes$code %in% speciesToInclude,]
# merge with phylodistances based on spp name (need full species name)
phyloDistances_fromHuman_PlusLabels1 <- merge(phyloDistances_fromHuman,speciesCodes,by.x="Sp1",by.y="species")
phyloDistances_fromHuman_PlusLabels2 <- merge(phyloDistances_fromHuman_PlusLabels1,speciesCodes,by.x="Sp2",by.y="species",suffixes = c(".Sp1",".Sp2"))

phyloDistances_fromHuman_PlusLabels2$NonHumanspeciesLabelForMergingWithSigfit <- phyloDistances_fromHuman_PlusLabels2$code.Sp1

phyloDistances_fromHuman_PlusLabels2[phyloDistances_fromHuman_PlusLabels2$code.Sp1=="humans_AFR",]$NonHumanspeciesLabelForMergingWithSigfit <- phyloDistances_fromHuman_PlusLabels2[phyloDistances_fromHuman_PlusLabels2$code.Sp1=="humans_AFR",]$code.Sp2


#View(phyloDistances_fromHuman_PlusLabels2)

fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS <- merge(fit_aging_signatures_COSINE_SIMILARITIES_df,phyloDistances_fromHuman_PlusLabels2[,c("NonHumanspeciesLabelForMergingWithSigfit","cophenetic_distance")],by.x="species",by.y="NonHumanspeciesLabelForMergingWithSigfit",all.x = T)

# set human - human phylo dist to 0
fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS[fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species=="humans_AFR",]$cophenetic_distance <- 0

plot1_cosineSimilarityVsPhylogeneticDistance_AgingSignatures <- ggplot(fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS,aes(x=cophenetic_distance,y=1-cosine_similarity))+
  geom_point()+
  geom_text(aes(label=species))+
  theme_bw()+
  ggtitle("Aging signatures (1-mer-CpG)")+
  xlab("cophenetic distance from human (based on Upham et al. RAXML tree)\nnote long mouse branch")

plot1_cosineSimilarityVsPhylogeneticDistance_AgingSignatures
ggsave(paste0(outdir,"plot1_cosineSimilarityVsPhylogeneticDistance_AgingSignatures.pdf"),plot1_cosineSimilarityVsPhylogeneticDistance_AgingSignatures)

fit_aging_signatures_COSINE_SIMILARITIES_df$species <- factor(fit_aging_signatures_COSINE_SIMILARITIES_df$species,levels=speciesOrder)

plot2_cosineSimPerSpecies_AgingSignatures <- ggplot(fit_aging_signatures_COSINE_SIMILARITIES_df,aes(x=species,y=cosine_similarity))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("Aging signatures (1-mer-CpG)")+
  xlab("")
plot2_cosineSimPerSpecies_AgingSignatures
ggsave(paste0(outdir,"plot2_cosineSimPerSpecies_AgingSignatures.pdf"),plot2_cosineSimPerSpecies_AgingSignatures)

# plot by ascending distance from human:
fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species <- factor(fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species, levels=fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species[order(fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$cophenetic_distance)], ordered=TRUE)

aging_cosineSimilarity_OrderedByincreasingcopheneticdistance <- ggplot(fit_aging_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS,aes(x=species,y=cosine_similarity))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("model: aging (excluding CpG>TpG)")+
  xlab("species ordered by increasing cophenetic distance from humans")+
  ylab("cosine similarity between\nreconstruction and empirical spectrum")
aging_cosineSimilarity_OrderedByincreasingcopheneticdistance

ggsave(paste0(outdir,"aging_cosineSimilarity_OrderedByincreasingcopheneticdistance.pdf"),aging_cosineSimilarity_OrderedByincreasingcopheneticdistance,height=5,width=8)

# merge with SBS1+5 as well:

fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS <- merge(fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df,phyloDistances_fromHuman_PlusLabels2[,c("NonHumanspeciesLabelForMergingWithSigfit","cophenetic_distance")],by.x="species",by.y="NonHumanspeciesLabelForMergingWithSigfit",all.x = T)

# set human - human phylo dist to 0
fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS[fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species=="humans_AFR",]$cophenetic_distance <- 0
# things to maybe  do
# port over plotting scripts from other script
# maybe write a function to get aitch distance
# need ot bring in phylo distances to human

plot1_cosineSimilarityVsPhylogeneticDistance_SBS1SBS5 <- ggplot(fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS,aes(x=cophenetic_distance,y=1-cosine_similarity))+
  geom_point()+
  geom_text(aes(label=species))+
  theme_bw()+
  ggtitle("SBS1 + SBS5 (3-mer)")+
  xlab("cophenetic distance from human (based on Upham et al. RAXML tree)\nnote long mouse branch")

plot1_cosineSimilarityVsPhylogeneticDistance_SBS1SBS5
ggsave(paste0(outdir,"plot1_cosineSimilarityVsPhylogeneticDistance_SBS1SBS5.pdf"),plot1_cosineSimilarityVsPhylogeneticDistance_SBS1SBS5)

fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species <- factor(fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species,levels=speciesOrder)

plot2_cosineSimPerSpecies_SBS1SBS5 <- ggplot(fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS,aes(x=species,y=cosine_similarity))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("SBS1 + SBS5 (3-mer)")+
  xlab("")
plot2_cosineSimPerSpecies_SBS1SBS5

ggsave(paste0(outdir,"plot2_cosineSimPerSpecies_SBS1SBS5.pdf"),plot2_cosineSimPerSpecies_SBS1SBS5)

######## plot by ascending distance from human: ########
fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species <- factor(fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species, levels=fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species[order(fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$cophenetic_distance)], ordered=TRUE)

SBS1SBS5_cosineSimilarity_OrderedByincreasingcopheneticdistance <- ggplot(fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS,aes(x=species,y=cosine_similarity))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("model: SBS1+SBS5 (3-mer) (including CpG>TpG 3mers)")+
  xlab("species ordered by increasing cophenetic distance from humans")+
  ylab("cosine similarity between\nreconstruction and empirical spectrum")
SBS1SBS5_cosineSimilarity_OrderedByincreasingcopheneticdistance

ggsave(paste0(outdir,"SBS1SBS5_cosineSimilarity_OrderedByincreasingcopheneticdistance.pdf"),SBS1SBS5_cosineSimilarity_OrderedByincreasingcopheneticdistance,height=5,width=8)


### also plot it for the aging+SBS5 (1-mer)
fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS <- merge(fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df,phyloDistances_fromHuman_PlusLabels2[,c("NonHumanspeciesLabelForMergingWithSigfit","cophenetic_distance")],by.x="species",by.y="NonHumanspeciesLabelForMergingWithSigfit",all.x = T)

# set human - human phylo dist to 0
fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS[fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species=="humans_AFR",]$cophenetic_distance <- 0

######## plot by ascending distance from human: ########
fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species <- factor(fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species, levels=fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$species[order(fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS$cophenetic_distance)], ordered=TRUE)

fit_agingMinusCpGPlusSBS5_cosineSimilarity_OrderedByincreasingcopheneticdistance <- ggplot(fit_agingMinusCpGPlusSBS5_signatures_COSINE_SIMILARITIES_df_PLUSPHYLODISTFROMHUMANS,aes(x=species,y=cosine_similarity))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("model: fit_agingMinusCpGPlusSBS5")+
  xlab("species ordered by increasing cophenetic distance from humans")+
  ylab("cosine similarity between\nreconstruction and empirical spectrum")
fit_agingMinusCpGPlusSBS5_cosineSimilarity_OrderedByincreasingcopheneticdistance

ggsave(paste0(outdir,"fit_agingMinusCpGPlusSBS5_cosineSimilarity_OrderedByincreasingcopheneticdistance.pdf"),fit_agingMinusCpGPlusSBS5_cosineSimilarity_OrderedByincreasingcopheneticdistance,height=5,width=8)

##### to do: plot all models cosine sim? for IS? #####

###### >>>>>>> make nice plots for manuscript <<<<<<<<<<<<< ##############

######### nice signature plots ###########
niceplotoutdir=paste0(outdir,"niceplotsforms/")
dir.create(niceplotoutdir)

list_ResultsToPlotTogether_forMansucript=list("SBS1 + SBS5 (3-mer)"=fit_SBS1_SBS5_signatures,"aging signatures (1-mer-CpG)"=fit_aging_signatures,"SBS1 + SBS5 + novel (3-mer)"=fit_SBS1_SBS5_PlusNovel,"aging signatures + novel (1-mer-CpG)"=fit_aging_signatures_PlusNovel)

#modelOrder = c("SBS1 + SBS5 (3-mer)", "aging signatures (1-mer-CpG)","SBS1 + SBS5 + novel (3-mer)","aging signatures + novel (1-mer-CpG)")
modelOrder = c( "aging signatures (1-mer-CpG)","aging signatures + novel (1-mer-CpG)","SBS1 + SBS5 (3-mer)","SBS1 + SBS5 + novel (3-mer)")

# get signatures:
list_ResultsForManuscript_signatures=list("aging signatures (1-mer-CpG)"=data.frame(fit_aging_signatures$data$signatures),"SBS1 + SBS5 (3-mer)"=data.frame(fit_SBS1_SBS5_signatures$data$signatures),"SBS1 + SBS5 + novel (3-mer)"=retrieve_pars(fit_SBS1_SBS5_PlusNovel,par="signatures")$mean,"aging signatures + novel (1-mer-CpG)"=retrieve_pars(fit_aging_signatures_PlusNovel,par="signatures")$mean) ## this may need work 
# for the ones with extracted sigs you need to export signatures 

list_ResultsForManuscript_signatures = lapply(list_ResultsForManuscript_signatures,rownames_to_column,var="signature")

# name the signatures ABC,D based on aging sigs + novel 

list_ResultsForManuscript_signatures$`SBS1 + SBS5 (3-mer)`$SignatureName <- list_ResultsForManuscript_signatures$`SBS1 + SBS5 (3-mer)`$signature

list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`$SignatureName <- list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`$signature

list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`[list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`$signature=="maternal_age",]$SignatureName <- "maternal age"
list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`[list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`$signature=="paternal_age",]$SignatureName <- "paternal age"
list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`[list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`$signature=="young_parent_13",]$SignatureName <- "young parent"

list_ResultsForManuscript_signatures$`aging signatures + novel (1-mer-CpG)`$SignatureName <- c(list_ResultsForManuscript_signatures$`aging signatures (1-mer-CpG)`$SignatureName,"novel (1-mer-CpG)")

list_ResultsForManuscript_signatures$`SBS1 + SBS5 + novel (3-mer)`$SignatureName <- c(list_ResultsForManuscript_signatures$`SBS1 + SBS5 (3-mer)`$signature,"novel (3-mer)")




signatureOrder_forManuscript=c("maternal age","paternal age","young parent","novel (1-mer-CpG)","SBS1","SBS5","novel (3-mer)")

agingPlusNovelAllSigs <- melt(list_ResultsForManuscript_signatures$`aging signatures + novel (1-mer-CpG)`)

agingPlusNovelAllSigs$SignatureName <- factor(agingPlusNovelAllSigs$SignatureName,levels=c("maternal age","paternal age","young parent","novel (1-mer-CpG)"))

AgingAndNovelSignaturesPlot <- ggplot(agingPlusNovelAllSigs,aes(x=gsub("\\.",">",variable),fill=SignatureName,y=value))+
  geom_col(position="dodge")+
  xlab("")+
  ylab("mutation fraction") +
  scale_fill_manual(values=signatureColors_forManuscript_niceNames[c("maternal age","paternal age","young parent","novel (1-mer-CpG)")])+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle=45,hjust=1),legend.title=element_blank())+
  ggtitle("Aging and novel mutation signatures")
ggsave(paste0(niceplotoutdir,"AgingAndNovelSignaturesPlot."))

# Instead of plotting SBS signatures and novel 3mer signatures, use the output of SIGFIT -- is in standard COSMIC format! 
# plot aging signatures:
# agingSignatures_melt <- melt(fit_aging_signatures$data$signatures,varnames = c("Signature","mutation_type"))
# 
# signaturePlot_aging <- ggplot(agingSignatures_melt,aes(x=gsub("\\.",">",mutation_type),y=value,fill=Signature))+
#   geom_col(position="dodge")+
#   theme_bw()+
#   scale_fill_manual(values=signatureColors_forManuscript[c("maternal_age","paternal_age","young_parent_13")])+
#   theme_bw()+
#   ylab("mutation fraction")+
#   xlab("")+
#   ggtitle("Aging Signatures")
# signaturePlot_aging
# ggsave(paste0(niceplotoutdir,"agingSignatures.pdf"),signaturePlot_aging,height=2,width=8.4)
#   # get sbs1 and 5 plots from COSMIC
# 
# # plot SBS1 + SBS5
# SBS1_SBS5_melt <- melt(fit_SBS1_SBS5_signatures$data$signatures,varnames = c("Signature","mutation_type"))
# 
# signaturePlot_SBS1_SBS5 <- ggplot(SBS1_SBS5_melt,aes(x=mutation_type,y=value,fill=Signature))+
#   geom_col(position="dodge")+
#   theme_bw()+
#   scale_fill_manual(values=signatureColors_forManuscript[c("SBS1","SBS5")])+
#   theme_bw()+
#   ylab("mutation fraction")+
#   xlab("")+
#   ggtitle("SBS1 and SBS5 COSMIC Signatures")+
#   theme(axis.text.x = element_text(angle=90,hjust=1,size=4))
# signaturePlot_SBS1_SBS5
# ggsave(paste0(niceplotoutdir,"SBS1_SBS5.pdf"),signaturePlot_SBS1_SBS5,height=2,width=8.4)

############ nice plot of exposures ###########

list_ResultsToPlotTogether_forMansucript_exposures = lapply(list_ResultsToPlotTogether_forMansucript,retrieve_pars,par="exposures")
# get means: 
list_ResultsToPlotTogether_forMansucript_exposures_mean = lapply(list_ResultsToPlotTogether_forMansucript_exposures,"[[","mean")
# make rownames a column: 

list_ResultsToPlotTogether_forMansucript_exposures_mean= lapply(list_ResultsToPlotTogether_forMansucript_exposures_mean,rownames_to_column,var="species")


# need to melt and combine:
list_ResultsToPlotTogether_forMansucript_exposures_mean_melt <- lapply(list_ResultsToPlotTogether_forMansucript_exposures_mean,melt)

# combine:
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt <- bind_rows(list_ResultsToPlotTogether_forMansucript_exposures_mean_melt,.id="id")

df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$species <- factor(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$species,levels=speciesOrder)

df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id <- factor(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id,levels=modelOrder)

# need to give signatures an alt-name 
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$SignatureName <- df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable

# for the aging  + novel model sigs A B and C are mat pat yp , D is novel 1-mer signature

# for the sbs1+5+ novel sigs are sbs1, sbs5, novel 

# get nice signature names:
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$niceSignatureNames <- as.character(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable)

# for aging+novel : sigs A, B and C D are mat/pat/yp age and novel respectively
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="aging signatures + novel (1-mer-CpG)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="Signature A",]$niceSignatureNames <- "maternal age"
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="aging signatures + novel (1-mer-CpG)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="Signature B",]$niceSignatureNames <- "paternal age"
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="aging signatures + novel (1-mer-CpG)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="Signature C",]$niceSignatureNames <- "young parent"
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="aging signatures + novel (1-mer-CpG)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="Signature D",]$niceSignatureNames <- "novel (1-mer-CpG)"

# for sbs1+5+novel sig a = sbs1, b= sbs5, c= novel
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="SBS1 + SBS5 + novel (3-mer)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="Signature A",]$niceSignatureNames <- "SBS1"
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="SBS1 + SBS5 + novel (3-mer)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="Signature B",]$niceSignatureNames <- "SBS5"
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="SBS1 + SBS5 + novel (3-mer)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="Signature C",]$niceSignatureNames <- "novel (3-mer)"

# rename regular aging signatures so they match:
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="aging signatures (1-mer-CpG)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="maternal_age",]$niceSignatureNames <- "maternal age"
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="aging signatures (1-mer-CpG)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="paternal_age",]$niceSignatureNames <- "paternal age"

df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id=="aging signatures (1-mer-CpG)" & df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$variable=="young_parent_13",]$niceSignatureNames <- "young parent"

# separate 1mer and 3mer models:
# add in short names:
df_ResultsToPlotTogether_forMansucript_exposures_mean_melt <- merge(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt,speciesShortNames, by.x="species",by.y = "speciesOriginal",all = T)

df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$speciesSuperShort <- factor(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$speciesSuperShort,levels=speciesOrderShortNames)

exposuresComparisonPlot_forManuscript <- ggplot(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt,aes(x=speciesSuperShort,y=value,fill=niceSignatureNames))+
  geom_col(position="stack")+
  theme_bw()+
  facet_wrap(~id,ncol=1)+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
  scale_fill_manual(values=signatureColors_forManuscript_niceNames[signatureOrder_forManuscript])+
  theme(legend.title=element_blank(),text=element_text(size=14))+
  ylab("exposure fraction")+
  xlab("")+
  ggtitle("Exposures comparison\nNOTE 3mer model still contains CpG.TpG 3mers")+  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16))
exposuresComparisonPlot_forManuscript

ggsave(paste0(niceplotoutdir,"comparisons.exposuresComparisonPlot.ModelsChosenForManuscript.neworder.pdf"),exposuresComparisonPlot_forManuscript
       ,height=8,width=8)

exposuresComparisonPlot_forManuscript_WIDE <- ggplot(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt,aes(x=speciesSuperShort,y=value,fill=niceSignatureNames))+
  geom_col(position="stack")+
  theme_bw()+
  facet_wrap(~id,ncol=2)+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
  scale_fill_manual(values=signatureColors_forManuscript_niceNames[signatureOrder_forManuscript])+
  theme(legend.title=element_blank(),text=element_text(size=14))+
  ylab("exposure fraction")+
  xlab("")+
  ggtitle("Exposures comparison\nNOTE 3mer model still contains CpG.TpG 3mers")+  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16))
exposuresComparisonPlot_forManuscript_WIDE

ggsave(paste0(niceplotoutdir,"comparisons.exposuresComparisonPlot.ModelsChosenForManuscript.neworder.WIDE.pdf"),exposuresComparisonPlot_forManuscript_WIDE
       ,height=8,width=12)




exposuresComparisonPlot_forManuscript_1mer <- ggplot(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id %in% c("aging signatures (1-mer-CpG)","aging signatures + novel (1-mer-CpG)"),],aes(x=speciesSuperShort,y=value,fill=niceSignatureNames))+
  geom_col(position="stack")+
  theme_bw()+
  facet_wrap(~id,ncol=2)+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
  scale_fill_manual(values=signatureColors_forManuscript_niceNames[c( "maternal age" ,     "paternal age"   ,   "young parent"  ,    "novel (1-mer-CpG)")])+
  theme(legend.title=element_blank(),text=element_text(size=14))+
  ylab("exposure fraction")+
  xlab("")+
  theme(strip.background = element_rect(fill="lightblue"))

  #ggtitle("Exposures comparison\nNOTE 3mer model still contains CpG.TpG 3mers")+  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16))
exposuresComparisonPlot_forManuscript_1mer

ggsave(paste0(niceplotoutdir,"comparisons.exposuresComparisonPlot.ModelsChosenForManuscript.neworder.WIDE.1merONLY.pdf"),exposuresComparisonPlot_forManuscript_1mer
       ,height=5,width=12)


exposuresComparisonPlot_forManuscript_3mer <- ggplot(df_ResultsToPlotTogether_forMansucript_exposures_mean_melt[df_ResultsToPlotTogether_forMansucript_exposures_mean_melt$id %in% c("SBS1 + SBS5 (3-mer)","SBS1 + SBS5 + novel (3-mer)"),],aes(x=speciesSuperShort,y=value,fill=niceSignatureNames))+
  geom_col(position="stack")+
  theme_bw()+
  facet_wrap(~id,ncol=2)+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
  scale_fill_manual(values=signatureColors_forManuscript_niceNames[c( "SBS1", "SBS5","novel (3-mer)")])+
  theme(legend.title=element_blank(),text=element_text(size=14))+
  ylab("exposure fraction")+
  xlab("")+
  theme(strip.background = element_rect(fill="lightblue"))
#ggtitle("Exposures comparison\nNOTE 3mer model still contains CpG.TpG 3mers")+  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16))
exposuresComparisonPlot_forManuscript_3mer

ggsave(paste0(niceplotoutdir,"comparisons.exposuresComparisonPlot.ModelsChosenForManuscript.neworder.WIDE.3merONLY.pdf"),exposuresComparisonPlot_forManuscript_3mer
       ,height=5,width=12)
############### nice residuals comparisons ###########
##### add short names: 
residualComparison_aging_agingPlusNovel <- ResidualsBoxPlot_WithOutlierLabels_MultipleModels(list("aging"=fit_aging_signatures,"aging + novel" = fit_aging_signatures_PlusNovel),runlabel = "aging_vs_agingPlusNovel",outdir =niceplotoutdir,textSize = 3,modelColors = c("#1B9E77","#D95F02"),speciesOrderShortNames = speciesOrderShortNames,speciesShortNames = speciesShortNames)
residualComparison_aging_agingPlusNovel

# have box plot be per mut type with species as points: 
residualComparison_aging_agingPlusNovel_PERMUTTYPE <- ResidualsBoxPlotPERMUTATIONTYPE_MultipleModels(list("aging"=fit_aging_signatures,"aging + novel" = fit_aging_signatures_PlusNovel),runlabel = "aging_vs_agingPlusNovel",outdir =niceplotoutdir,textSize = 3,modelColors = c("#1B9E77","#D95F02"),speciesOrderShortNames = speciesOrderShortNames,speciesShortNames = speciesShortNames)


residualComparison_SBS1SBS5_SBS1SBS5PlusNovel <- ResidualsBoxPlot_WithOutlierLabels_MultipleModels(list("SBS1 + SBS5"=fit_SBS1_SBS5_signatures,"SBS1 + SBS5 + novel" = fit_SBS1_SBS5_PlusNovel),runlabel = "SBS1_5_vs_SBSPlusNovel",outdir =niceplotoutdir,textSize = 2,modelColors = c("#1B9E77","#D95F02"),speciesOrderShortNames = speciesOrderShortNames,speciesShortNames = speciesShortNames)

residualComparison_SBS1SBS5_SBS1SBS5PlusNovel
##### plot 3mers per 1mer mutation type
ResidualsBoxPlotPERMUTATIONTYPE_WithOutlierLabels_MultipleModels_SEP3merBy1Mer(list("SBS1 + SBS5"=fit_SBS1_SBS5_signatures,"SBS1 + SBS5 + novel" = fit_SBS1_SBS5_PlusNovel),runlabel = "SBS1_5_vs_SBSPlusNovel",outdir =niceplotoutdir,textSize = 2,modelColors = c("#1B9E77","#D95F02"),speciesOrderShortNames = speciesOrderShortNames,speciesShortNames = speciesShortNames)

############ nice cosine similarity comparisons #########

agingModelsCosineSims <- bind_rows(fit_aging_signatures_COSINE_SIMILARITIES_df,fit_aging_signatures_PlusNovel_COSINE_SIMILARITIES_df)
# merge with short names
agingModelsCosineSims <- merge(agingModelsCosineSims,speciesShortNames,by.x="species",by.y="speciesOriginal",all=T)
agingModelsCosineSims$speciesSuperShort <- factor(agingModelsCosineSims$speciesSuperShort,levels=speciesOrderShortNames)

cosineSimsAgingModelsPlot <- ggplot(agingModelsCosineSims,aes(x=speciesSuperShort,y=cosine_similarity,color=nicer_nmodel_name))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(values= c("#1B9E77","#D95F02"))+
  ylab("cosine similarity between empirical\nspectrum and reconstruction")+
  xlab("")+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.title = element_blank())

cosineSimsAgingModelsPlot
ggsave(paste0(niceplotoutdir,"CosineSims.AgingModels.pdf"),cosineSimsAgingModelsPlot,height=5,width=9)


SBSModelsCosineSims <- bind_rows(fit_SBS1_SBS5_signatures_COSINE_SIMILARITIES_df,fit_SBS1_SBS5_PlusNovel_signatures_COSINE_SIMILARITIES_df)
# merge with short names
SBSModelsCosineSims <- merge(SBSModelsCosineSims,speciesShortNames,by.x="species",by.y="speciesOriginal",all=T)
SBSModelsCosineSims$speciesSuperShort <- factor(SBSModelsCosineSims$speciesSuperShort,levels=speciesOrderShortNames)


cosineSims_SBSModels_Plot <- ggplot(SBSModelsCosineSims,aes(x=speciesSuperShort,y=cosine_similarity,color=nicer_nmodel_name))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(values= c("#1B9E77","#D95F02"))+
  ylab("cosine similarity between empirical\nspectrum and reconstruction")+
  xlab("")+
  theme(axis.text.x = element_text(angle=45,hjust=1),legend.title = element_blank())

cosineSims_SBSModels_Plot
ggsave(paste0(niceplotoutdir,"CosineSims.SBSModels.pdf"),cosineSims_SBSModels_Plot,height=5,width=9)

########### SKIP : functions to carry out mantel test on exposures ##############

# 
# clrTransformExposures_getDistances_forMantelTest <- function(sigfitResultsFull,speciesToInclude,speciesCodes){
#   exposures <- retrieve_pars(sigfitResultsFull,par="exposures")$mean # mean exposures
#   exposures$label <- rownames(exposures)
#   # subset to phylo species: 
#   exposures_phyloSpp <- exposures[exposures$label %in% speciesToInclude,]
#   exposures_phyloSpp_fullSppNames <- merge(exposures_phyloSpp,speciesCodes,by.x="label",by.y="code")
#   rownames(exposures_phyloSpp_fullSppNames) <- exposures_phyloSpp_fullSppNames$species
# 
#   exposures_clr <- clr(select(exposures_phyloSpp_fullSppNames,names(select(exposures,-label)))) # need to switch to young_parent -- will throw an error   # here names of exposures df are the names of the different signatures, subtracting off the label
# 
#   print("clr transform (head):")
#   print(head(exposures_clr))
#   
#   # for Mantel test, keep as a full matrix with full species names:
#   exposures_clr_dist_MatrixForMantel <- as.matrix(dist(exposures_clr,diag=T,upper=T))
#   return(exposures_clr_dist_MatrixForMantel)
#   
# }
# 
# # note for this funciton it will reorder the phylo dist matrix so you don't have to
# # and it will make sure it is in right order 
# # how to get common names in here 
# # this will carry out the square rooting itself: 
# carryOutMantelTest_SQRTDISTANCE <- function(phylo_dist_matrix_NOTYETSQRTED,exposure_dist_matrix_forMantel,npermut,label){
#   # as of 20220902 this is using sqrt phylo distance!! 
#   # reorder phylo dist matrix to match col and row names of exposures:
#   phylo_dist_matrix_NOTYETSQRTED <- as.matrix(phylo_dist_matrix_NOTYETSQRTED)[rownames(exposure_dist_matrix_forMantel),colnames(exposure_dist_matrix_forMantel)]
#   phylo_dist_matrix_sqrt <- sqrt(phylo_dist_matrix_NOTYETSQRTED)
#   if(sum(colnames(phylo_dist_matrix_sqrt)!=colnames(exposure_dist_matrix_forMantel))!=0){
#     print("something is wrong with ordering of colnames")
#     break
#   }
#   if(sum(rownames(phylo_dist_matrix_sqrt)!=rownames(exposure_dist_matrix_forMantel))!=0){
#     print("something is wrong with ordering of rownames")
#     break
#   }
#   
#   # need to reorder the phylo dist to make sure is in same order as distances 
#   mtr = vegan::mantel(xdis=phylo_dist_matrix_sqrt,ydis=exposure_dist_matrix_forMantel,permutations = npermut,method = "pearson")
#   # make a data frame
#   mantel_ExposureToSignatures_df <- data.frame(call="vegan_mantel_SQRTDISTANCE",method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),label=label)
#   
#   return(mantel_ExposureToSignatures_df)
# }
# 
# melt_cophenetic_func_SQRTDISTANCE <- function(cophenetic_dist){
#   cophenetic_dist_sqrt <- sqrt(cophenetic_dist)
#   cophenetic_dist_sqrt_melt <- melt(cophenetic_dist_sqrt)
#   colnames(cophenetic_dist_sqrt_melt) <- c("Sp1","Sp2","sqrt_cophenetic_distance")
#   # get species names # ASSUMES ARE SEPARATED BY "_" eg mus_musculus
#   cophenetic_dist_sqrt_melt$Sp1_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_sqrt_melt$Sp1),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_sqrt_melt$Sp1),"_"),"[",2)))
#   
#   cophenetic_dist_sqrt_melt$Sp2_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_sqrt_melt$Sp2),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_sqrt_melt$Sp2),"_"),"[",2)))
#   
#   # need to deal with reciprocal duplicates
#   # put in alphabetical order so I can get rid of dups:
#   cophenetic_dist_sqrt_melt$comparisonLabel <- paste0(pmin(cophenetic_dist_sqrt_melt$Sp1_species,cophenetic_dist_sqrt_melt$Sp2_species),".",pmax(cophenetic_dist_sqrt_melt$Sp1_species,cophenetic_dist_sqrt_melt$Sp2_species))  # # need to use pmin and pmax to get row by row min max. otherwise gives min max of whoel vector
#   # get rid of self to self comparisons:
#   cophenetic_dist_sqrt_melt_distinct <- cophenetic_dist_sqrt_melt[cophenetic_dist_sqrt_melt$Sp1!=cophenetic_dist_sqrt_melt$Sp2,]
#   # and get rid of reciprocal dups (sp A - B and sp B- A)
#   cophenetic_dist_sqrt_melt_distinct <- cophenetic_dist_sqrt_melt_distinct %>%
#     ungroup() %>% # adding ungroup as a safety measure on 20211201 not sure if ottally necessary though
#     distinct(comparisonLabel,.keep_all=T)  # restrict to just distinct comparison labels. 
#   dim(cophenetic_dist_sqrt_melt_distinct)
#   return(cophenetic_dist_sqrt_melt_distinct)
#   
#   
# }
# 
# meltDistanceCombineWithPhyloToPlot_SQRTDISTANCE <- function(phylo_dist_matrix_NOTYETSQRTED,exposure_dist_matrix_forMantel,speciesCodes){
#   # melt phylo distances and add comparison label
#   phylo_dist_melted_sqrt <- melt_cophenetic_func_SQRTDISTANCE(phylo_dist_matrix_NOTYETSQRTED) #square roots and melts and adds comparison label 
#   exposures_melt <- melt(exposure_dist_matrix_forMantel) # contains symmetrical dups
#   dim(exposures_melt)
#   colnames(exposures_melt) <- c("sp1","sp2","distance")
#   exposures_melt$comparisonLabel <-  paste0(pmin(as.character(exposures_melt$sp1),as.character(exposures_melt$sp1)),".",pmax(as.character(exposures_melt$sp1),as.character(exposures_melt$sp2)))
#   # make unique:
#   exposures_melt_unique <- unique(exposures_melt[,c("comparisonLabel","distance")])
#   dim(exposures_melt_unique)
#   # merge with phylo distances
#   exposures_and_phylo_dist <- merge(exposures_melt_unique,phylo_dist_melted_sqrt,by="comparisonLabel") 
#   return(exposures_and_phylo_dist)
#   
# }
# 
# # this function returns all the different tables as a list (nice!)
# # gives you the distances of expsoures and phylo distances for makign plots
# # as well as a matrix of exposure distances used for mantel test (probably not needed)
# # table of mantel test results 
# # KEY: this iwll output sqrt cophenetic idstance but you want to input the regular dist matrix! keep an eye on this, not very defensively coded. 
# combo_MantelTestPlottingScript_SQRTDISTANCE <- function(sigfitResultsFull,speciesToInclude,speciesCodes,phylo_dist_matrix_NOTYETSQRTED,npermut,label,outdir){
#   # get clr distances based on exposures
#   # get name of results :
#   nameOfInputResults = deparse(substitute(sigfitResultsFull)) # more defensive coding in case you mix up your labels 
#   exposureDistancesForMantelTest <- clrTransformExposures_getDistances_forMantelTest(sigfitResultsFull,speciesToInclude,speciesCodes)
#   # do mantel test:
#   mantelTestResults <- carryOutMantelTest_SQRTDISTANCE(phylo_dist_matrix_NOTYETSQRTED,exposureDistancesForMantelTest,npermut,label=label) # this will convert to sqrt cophenetic distance 
#   write.table(mantelTestResults,paste0(outdir,label,".",nameOfInputResults,".clrDistances.mantelTestResults.sqrtcophentic.txt"),row.names=F,quote=F,sep="\t")
#   # get distances in a format that's easy to plot: 
#   exposureDistancesForPlotting <- meltDistanceCombineWithPhyloToPlot_SQRTDISTANCE(phylo_dist_matrix_NOTYETSQRTED,exposureDistancesForMantelTest,speciesCodes)
#   write.table(exposureDistancesForPlotting,paste0(outdir,label,".",nameOfInputResults,".clrDistances.forplotting.sqrtcophenetic.txt"),row.names=F,quote=F,sep="\t")
#   # make a plot
#   DistancePlot <- ggplot(exposureDistancesForPlotting,aes(x=sqrt_cophenetic_distance,y=distance))+
#     geom_point(size=1,color="dodgerblue")+
#     #geom_text(aes(label=comparisonLabel),size=1)+
#     theme_bw()+
#     ggtitle(paste0(label," ",nameOfInputResults))+
#     geom_text(data=mantelTestResults,aes(x=.5,y=Inf,vjust=2,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test; ",permCount," perm.)")),size=5)+
#     theme(text=element_text(size=12))+
#     xlab("square root of phylogenetic distance (shared branch length)")+
#     ylab("Aitichison distance between\nmutation signature exposures")
#   
#   print(DistancePlot)
#   # save plot:
#   ggsave(paste0(outdir,label,".",nameOfInputResults,".clrDistances.sqrtphylodistance.mantel.pdf"),DistancePlot,height=4,width=6)
#   returnList <- list(all_distancesForPlotting=exposureDistancesForPlotting,mantelResults=mantelTestResults,exposureDistancesForMantelTest=exposureDistancesForMantelTest,plot=DistancePlot)
#   return(returnList) #returns everything as a list so you can replot etc.
#   
# }


######### get raxml matrix for mantel test #######
# process raxml tree here instead of in previous 
raxmlTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")

# want to get rid of the family/order info at end of tip. just keep species : 
raxmlTree_renamedTips <- raxmlTree
raxmlTree_renamedTips$tip.label <- paste0(unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",1)),"_",unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",2)))

speciesToInclude_in_raxmlTree=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Canis_lupus")  # raxml name fmt
# adding in apes 
# restrict to just the species to include
raxmlTree_renamedTips_subset <- keep.tip(raxmlTree_renamedTips,speciesToInclude_in_raxmlTree)

raxml_cophenetic_dist <- cophenetic(raxmlTree_renamedTips_subset)
######### SKIP:  run mantel test on exposures ########
# aging_1merSignaturesExposureMantelComboResults <- combo_MantelTestPlottingScript_SQRTDISTANCE(sigfitResultsFull=fit_aging_signatures,speciesToInclude = speciesToInclude,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "aging (exclue CpG>TpG)",outdir=paste0(outdir,"fit_aging_signatures/"))
# 
# SBS15_3merSignaturesExposureMantelComboResults <- combo_MantelTestPlottingScript_SQRTDISTANCE(sigfitResultsFull=fit_SBS1_SBS5_signatures,speciesToInclude = speciesToInclude,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "SBS1 + SBS5 (3-mer,includes CpG>TpG)",outdir=paste0(outdir,"fit_SBS1_SBS5_signatures/"))
# 
# agingPlus1merSBS5_1merSignaturesExposureMantelComboResults <- combo_MantelTestPlottingScript_SQRTDISTANCE(sigfitResultsFull=fit_agingMinusCpGPlusSBS5_signatures,speciesToInclude = speciesToInclude,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "aging + SBS5 (1-mer,excluding CpG>TpG)",outdir=paste0(outdir,"fit_agingMinusCpGPlusSBS5_signatures/"))

#### functions for mantel test on reconstrutions rather than exposures ########
### adding downsampling: (20230328)
melt_cophenetic_func_SQRTDISTANCE <- function(cophenetic_dist){
    cophenetic_dist_sqrt <- sqrt(cophenetic_dist)
    cophenetic_dist_sqrt_melt <- melt(cophenetic_dist_sqrt)
    colnames(cophenetic_dist_sqrt_melt) <- c("Sp1","Sp2","sqrt_cophenetic_distance")
    # get species names # ASSUMES ARE SEPARATED BY "_" eg mus_musculus
    cophenetic_dist_sqrt_melt$Sp1_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_sqrt_melt$Sp1),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_sqrt_melt$Sp1),"_"),"[",2)))

    cophenetic_dist_sqrt_melt$Sp2_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_sqrt_melt$Sp2),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_sqrt_melt$Sp2),"_"),"[",2)))

    # need to deal with reciprocal duplicates
    # put in alphabetical order so I can get rid of dups:
    cophenetic_dist_sqrt_melt$comparisonLabel <- paste0(pmin(cophenetic_dist_sqrt_melt$Sp1_species,cophenetic_dist_sqrt_melt$Sp2_species),".",pmax(cophenetic_dist_sqrt_melt$Sp1_species,cophenetic_dist_sqrt_melt$Sp2_species))  # # need to use pmin and pmax to get row by row min max. otherwise gives min max of whoel vector
    # get rid of self to self comparisons:
    cophenetic_dist_sqrt_melt_distinct <- cophenetic_dist_sqrt_melt[cophenetic_dist_sqrt_melt$Sp1!=cophenetic_dist_sqrt_melt$Sp2,]
    # and get rid of reciprocal dups (sp A - B and sp B- A)
    cophenetic_dist_sqrt_melt_distinct <- cophenetic_dist_sqrt_melt_distinct %>%
      ungroup() %>% # adding ungroup as a safety measure on 20211201 not sure if ottally necessary though
      distinct(comparisonLabel,.keep_all=T)  # restrict to just distinct comparison labels.
    dim(cophenetic_dist_sqrt_melt_distinct)
    return(cophenetic_dist_sqrt_melt_distinct)


  }
  
clrTransformReconstructions_getDistances_forMantelTest_DOWNSAMPLING <- function(sigfitResultsFull,speciesToInclude,speciesCodes,downsample=TRUE){
  # reconstructions are given as a matrix whereas epxosures come as df 
  reconstructions <- data.frame(retrieve_pars(sigfitResultsFull,par="reconstructions")$mean) # mean reconstructions
  reconstructions$label <- rownames(reconstructions)
  # subset to phylo species: 
  reconstructions_phyloSpp <- reconstructions[reconstructions$label %in% speciesToInclude,]
  reconstructions_phyloSpp_fullSppNames <- merge(reconstructions_phyloSpp,speciesCodes,by.x="label",by.y="code")
  rownames(reconstructions_phyloSpp_fullSppNames) <- reconstructions_phyloSpp_fullSppNames$species
  
  # DOWNSAMPLE:  YOU ARE HERE:::
  # calculate min number of seg sites:
  if(downsample==TRUE){
    print("Downsampling reconstructions to min number of seg sites")
    
    reconstructions_phyloSpp_fullSppNames_melt <- melt(reconstructions_phyloSpp_fullSppNames)
    # note these have been rescaled to human content prior to sigfit
    totalSegSites <- reconstructions_phyloSpp_fullSppNames_melt %>%
    group_by(species,label) %>%
    summarise(totalSegSites=sum(value)) 
  
    subsampleValue=min(totalSegSites$totalSegSites)
  
    subsampleSpecies=totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$label
  
    print(paste0("downsampling to ",subsampleValue," sites (",subsampleSpecies,")"))
    # carry out multinomial sampling:
  
    # get proportions
    reconstructions_phyloSpp_fullSppNames_melt <- reconstructions_phyloSpp_fullSppNames_melt %>%
    group_by(species,label) %>%
    mutate(mutationProportion=value/sum(value))%>%
    ungroup()
    
    # downsample
    reconstructions_phyloSpp_fullSppNames_melt <- reconstructions_phyloSpp_fullSppNames_melt %>%
      group_by(species,label) %>%
      mutate(totalMutations_DOWNSAMPLED=as.numeric(rmultinom(n=1,size=subsampleValue,prob=mutationProportion))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
      ungroup()
    
    # get back into wide
    reconstructions_phyloSpp_fullSppNames_downsamp <- data.frame(pivot_wider(reconstructions_phyloSpp_fullSppNames_melt,names_from=variable,id_cols = c(species,label,common_name,broad_classification),values_from=totalMutations_DOWNSAMPLED))
    
    # need to reset rownames:
    rownames(reconstructions_phyloSpp_fullSppNames_downsamp) <- reconstructions_phyloSpp_fullSppNames_downsamp$species
    
    # set this to be the downsamp version 
    reconstructions_phyloSpp_fullSppNames = reconstructions_phyloSpp_fullSppNames_downsamp # use this downstream
  } else {
    print("NOT downsamling reconstructions prior to CLR transformation")
  }
  

  ## function resumes:

  # here names of reconstructions are the names of the mutation types (1 or 3mer )
  reconstructions_clr <- clr(select(reconstructions_phyloSpp_fullSppNames,names(select(reconstructions,-label)))) # need to switch to young_parent -- will throw an error
  print("clr transform (head):")
  print(head(reconstructions_clr))
  
  # for Mantel test, keep as a full matrix with full species names:
  reconstructions_clr_dist_MatrixForMantel <- as.matrix(dist(reconstructions_clr,diag=T,upper=T))
  return(reconstructions_clr_dist_MatrixForMantel)
  
}

# note for this function it will reorder the phylo dist matrix so you don't have to
# and it will make sure it is in right order 
# how to get common names in here 
# this will carry out the square rooting itself: 
carryOutMantelTest_SQRTDISTANCE_RECONSTRUCTIONS <- function(phylo_dist_matrix_NOTYETSQRTED,reconstructions_clr_dist_MatrixForMantel,npermut,label){
  # as of 20220902 this is using sqrt phylo distance!! 
  # reorder phylo dist matrix to match col and row names of reconstructions:
  phylo_dist_matrix_NOTYETSQRTED <- as.matrix(phylo_dist_matrix_NOTYETSQRTED)[rownames(reconstructions_clr_dist_MatrixForMantel),colnames(reconstructions_clr_dist_MatrixForMantel)]
  phylo_dist_matrix_sqrt <- sqrt(phylo_dist_matrix_NOTYETSQRTED)
  if(sum(colnames(phylo_dist_matrix_sqrt)!=colnames(reconstructions_clr_dist_MatrixForMantel))!=0){
    print("something is wrong with ordering of colnames")
    break
  }
  if(sum(rownames(phylo_dist_matrix_sqrt)!=rownames(reconstructions_clr_dist_MatrixForMantel))!=0){
    print("something is wrong with ordering of rownames")
    break
  }
  
  # need to reorder the phylo dist to make sure is in same order as distances 
  mtr = vegan::mantel(xdis=phylo_dist_matrix_sqrt,ydis=reconstructions_clr_dist_MatrixForMantel,permutations = npermut,method = "pearson")
  # make a data frame
  mantel_ReconstructionBasedonSignatures_df <- data.frame(call="vegan_mantel_SQRTDISTANCE",method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),label=label)
  
  return(mantel_ReconstructionBasedonSignatures_df)
}

### use funciton above for getting phylo dist

meltDistanceCombineWithPhyloToPlot_SQRTDISTANCE_RECONSTRUCTIONS <- function(phylo_dist_matrix_NOTYETSQRTED,reconstructions_clr_dist_MatrixForMantel,speciesCodes){
  # melt phylo distances and add comparison label
  phylo_dist_melted_sqrt <- melt_cophenetic_func_SQRTDISTANCE(phylo_dist_matrix_NOTYETSQRTED) #square roots and melts and adds comparison label 
  reconstructions_melt <- melt(reconstructions_clr_dist_MatrixForMantel) # contains symmetrical dups
  dim(reconstructions_melt)
  colnames(reconstructions_melt) <- c("sp1","sp2","distance")
  reconstructions_melt$comparisonLabel <-  paste0(pmin(as.character(reconstructions_melt$sp1),as.character(reconstructions_melt$sp1)),".",pmax(as.character(reconstructions_melt$sp1),as.character(reconstructions_melt$sp2)))
  # make unique:
  reconstructions_melt_unique <- unique(reconstructions_melt[,c("comparisonLabel","distance")])
  dim(reconstructions_melt_unique)
  # merge with phylo distances
  reconstructions_and_phylo_dist <- merge(reconstructions_melt_unique,phylo_dist_melted_sqrt,by="comparisonLabel") 
  return(reconstructions_and_phylo_dist)
  
}

# this function returns all the different tables as a list (nice!)
# gives you the distances of reconstructions and phylo distances for makign plots
# as well as a matrix of reconstruction distances used for mantel test (probably not needed)
# table of mantel test results 
# KEY: this iwll output sqrt cophenetic idstance but you want to input the regular dist matrix! keep an eye on this, not very defensively coded. 
combo_MantelTestPlottingScript_SQRTDISTANCE_RECONSTRUCTIONS <- function(sigfitResultsFull,speciesToInclude,speciesCodes,phylo_dist_matrix_NOTYETSQRTED,npermut,label,outdir,downsample=TRUE){
  # get clr distances based on reconstructions
  # get name of results :
  nameOfInputResults = deparse(substitute(sigfitResultsFull)) # more defensive coding in case you mix up your labels 
  reconstructionDistancesForMantelTest <- clrTransformReconstructions_getDistances_forMantelTest_DOWNSAMPLING(sigfitResultsFull,speciesToInclude,speciesCodes,downsample = downsample) ## reconstructions! 
  # do mantel test:
  mantelTestResults <- carryOutMantelTest_SQRTDISTANCE_RECONSTRUCTIONS(phylo_dist_matrix_NOTYETSQRTED,reconstructionDistancesForMantelTest,npermut,label=label) # this will convert to sqrt cophenetic distance 
  write.table(mantelTestResults,paste0(outdir,label,".",nameOfInputResults,".clrDistances.mantelTestResults.sqrtcophentic.RECONSTRUCTIONS.downsample",downsample,".txt"),row.names=F,quote=F,sep="\t")
  # get distances in a format that's easy to plot: 
  reconstructionDistancesForPlotting <- meltDistanceCombineWithPhyloToPlot_SQRTDISTANCE_RECONSTRUCTIONS(phylo_dist_matrix_NOTYETSQRTED,reconstructionDistancesForMantelTest,speciesCodes)
  write.table(reconstructionDistancesForPlotting,paste0(outdir,label,".",nameOfInputResults,".clrDistances.forplotting.sqrtcophenetic.RECONSTRUCTIONS.downsample",downsample,".txt"),row.names=F,quote=F,sep="\t")
  # make a plot
  DistancePlot <- ggplot(reconstructionDistancesForPlotting,aes(x=sqrt_cophenetic_distance,y=distance))+
    geom_point(size=1,color="dodgerblue")+
    geom_text(aes(label=comparisonLabel),size=1)+
    theme_bw()+
    ggtitle(paste0("DISTANCES BASED ON RECONSTRUCTIONS",label," ",nameOfInputResults,"\n downsampled =",downsample))+
    geom_text(data=mantelTestResults,aes(x=.5,y=Inf,vjust=2,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test; ",permCount," perm.)")),size=5)+
    theme(text=element_text(size=12))+
    xlab("square root of phylogenetic distance (shared branch length)")+
    ylab("Aitchison distance between\nmutation spectrum reconstructions")
  
  print(DistancePlot)
  # save plot:
  ggsave(paste0(outdir,label,".",nameOfInputResults,".clrDistances.sqrtphylodistance.mantel.RECONSTRUCTIONS.downsampled",downsample,".pdf"),DistancePlot,height=4,width=6)
  returnList <- list(all_distancesForPlotting=reconstructionDistancesForPlotting,mantelResults=mantelTestResults,reconstructionDistancesForMantelTest=reconstructionDistancesForMantelTest,plot=DistancePlot)
  return(returnList) #returns everything as a list so you can replot etc.
  
}
############# run mantel test on reconstructions #############
# apparently you have to re-set seed each time?
set.seed(42)
aging_1merSignaturesRECONSTRUCTIONMantelComboResults <- combo_MantelTestPlottingScript_SQRTDISTANCE_RECONSTRUCTIONS(sigfitResultsFull=fit_aging_signatures,speciesToInclude = speciesToInclude,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "aging (exclude CpG>TpG)",outdir=paste0(outdir,"fit_aging_signatures/"),downsample = TRUE)


set.seed(42)
SBS15_3merSignaturesRECONSTRUCTIONMantelComboResults <- combo_MantelTestPlottingScript_SQRTDISTANCE_RECONSTRUCTIONS(sigfitResultsFull=fit_SBS1_SBS5_signatures,speciesToInclude = speciesToInclude,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "SBS1 + SBS5 (3-mer,includes CpG>TpG)",outdir=paste0(outdir,"fit_SBS1_SBS5_signatures/"),downsample = T)

set.seed(42)
aging_PlusNovel_1merSignaturesRECONSTRUCTIONMantelComboResults <- combo_MantelTestPlottingScript_SQRTDISTANCE_RECONSTRUCTIONS(sigfitResultsFull=fit_aging_signatures_PlusNovel,speciesToInclude = speciesToInclude,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "aging + novel (exclue CpG>TpG)",outdir=paste0(outdir,"fit_aging_signatures_PlusNovel/"),downsample = T)

set.seed(42)
SBS15_PlusNovel_3merSignaturesRECONSTRUCTIONMantelComboResults <- combo_MantelTestPlottingScript_SQRTDISTANCE_RECONSTRUCTIONS(sigfitResultsFull=fit_SBS1_SBS5_PlusNovel,speciesToInclude = speciesToInclude,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "SBS1 + SBS5 + novel (3-mer,includes CpG>TpG)",outdir=paste0(outdir,"fit_SBS1_SBS5_signatures_PlusNovel/"),downsample = T)


####### write out mantel results 
############ save workspace image ###################
save.image(file = paste0(outdir,"RWORKSPACE.",todaysdate,".AtEndOfScript.image.RData"))

