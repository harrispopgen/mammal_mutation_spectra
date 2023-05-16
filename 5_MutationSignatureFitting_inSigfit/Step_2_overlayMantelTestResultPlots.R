########### special plotting code to overlay mantel results from 1/3mers with sigfit recnostruction mantel test results ###############

require(ggplot2)
require(scales)

############ read in mantel test results and clr-trasnformed distance dataframes for plotting ########


# need to write out mantel distances in mantel dir

######### empirical 1mer minus CpG mantel results and distances #######
mantelResults_Empirical1mer_NOCpGTpG <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/compare_spectrumDist_toPhyloHammingDist/20230301.plots.SubsampledIndividualsNotProjection.epsilon.1.sepCpGs.EXCLUDED.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode.9999999.MantelPermutations.CpG.TpGEXCLUDED/mantel_test_results_raxml_SQRTDISTANCE_df.txt",header=T,sep="\t") # for some reason needed quote=T

# aha it's because of the ' in pearson's -- messes up quotes

distances_Empirical1mer_noCpGTpG <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/compare_spectrumDist_toPhyloHammingDist/20230301.plots.SubsampledIndividualsNotProjection.epsilon.1.sepCpGs.EXCLUDED.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode.9999999.MantelPermutations.CpG.TpGEXCLUDED/all_distances_all_conditions_phylo.1merNoCpGTpg.txt",header=T,sep="\t")

distances_Empirical1mer_noCpGTpG_subset <- distances_Empirical1mer_noCpGTpG[distances_Empirical1mer_noCpGTpG$id=="1-mer-CpG spectrum" & distances_Empirical1mer_noCpGTpG$transform_label=="clr" & distances_Empirical1mer_noCpGTpG$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",]

rm(distances_Empirical1mer_noCpGTpG)
######## empirical 3mer mantel results and distances #######

mantelResults_Empirical3mer <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/compare_spectrumDist_toPhyloHammingDist/20221121.plots.SubsampledIndividualsNotProjection.epsilon.1.sepCpGs.no.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode.9999999.MantelPermutations/mantelTest/mantelTest.Results.includesFoldedSpectra.txt",header=T,sep="\t",quote="") # note this includes results from lots of diff spectra including folded spectra, but also has the results we want to plot -- quotes mess it up!!! 

mantelResults_Empirical3mer_subset <- mantelResults_Empirical3mer[mantelResults_Empirical3mer$id=="3-mer spectrum" & mantelResults_Empirical3mer$distance_metric=="raxml_tree_SQRT" & mantelResults_Empirical3mer$call=="vegan_mantel_SQRTPHYLODISTANCE",]

distances_Empirical3mer <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/compare_spectrumDist_toPhyloHammingDist/20221121.plots.SubsampledIndividualsNotProjection.epsilon.1.sepCpGs.no.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode.9999999.MantelPermutations/all_distances_all_conditions_phylo.RESCALED_BY_HUMAN_TARGETS.txt",header=T,sep="\t") # need to subset to

# NOTE you need to sqrt distance! 
# mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS <- restrict to
distances_Empirical3mer_subset  <- distances_Empirical3mer[distances_Empirical3mer$id=="3-mer spectrum" & distances_Empirical3mer$transform_label=="clr" & distances_Empirical3mer$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",] # 78 comparisons (13 choose 2 species)

rm(distances_Empirical3mer) #so you dont accidentally plot with it






######### function to plot ########

PlotEmpiricalOverlayReconstructionMantel <- function(empiricalMantelResults, empiricalDistances,ReconstructionMantelResults,ReconstructionDistances,label,reconstructionPointColor){
  
  
  plot1 <- ggplot(empiricalDistances,aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
    geom_point(aes(color="empirical",shape="empirical"))+
    geom_text(data=empiricalMantelResults,aes(x=0.5,y=Inf,vjust=1.1,label=paste0("[ Empirical data ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)"),color="empirical"),size=5)+
    geom_point(data=ReconstructionDistances,aes(x=sqrt_cophenetic_distance,y=distance,color="reconstructions",shape="reconstructions"))+
    geom_text(data=ReconstructionMantelResults,aes(x=0.5,y=-Inf,vjust=-5,label=paste0(" [ Sigfit reconstructions ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)"),color="reconstructions"),size=5)+
    scale_color_manual(values=c("black",reconstructionPointColor))+
    theme_bw()+
    xlab("sqrt of phylogenetic distance")+
    ylab("Aitchison distance between species' mutation spectra")+
    ggtitle(label)+
    theme(legend.position = "none")
  
  return(plot1)
  
  # note that reconstructions have already converted phylo dist to sqrt, but for empirical I am sqrt'ing it in the above plotting code (manually looks good)
  
}



########## read in sigfit models ############
sigfitMantelDir="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/sigfit/SigfitRuns_notcountsplit_notprojected/20230329.sigfitResults.multinomial.niter.10000.notprojected.CpGTpGExcluded/"


######### aging model ##########
mantelResults_ReconstructionsAging=read.table(paste0(sigfitMantelDir,"/fit_aging_signatures/aging (exclude CpG>TpG).fit_aging_signatures.clrDistances.mantelTestResults.sqrtcophentic.RECONSTRUCTIONS.downsampleTRUE.txt"),header=T,sep="\t",quote="")

distances_ReconstructionsAging=read.table(paste0(sigfitMantelDir,"/fit_aging_signatures/aging (exclude CpG>TpG).fit_aging_signatures.clrDistances.forplotting.sqrtcophenetic.RECONSTRUCTIONS.downsampleTRUE.txt"),header=T,sep="\t")


agingPlot <- PlotEmpiricalOverlayReconstructionMantel(empiricalMantelResults =mantelResults_Empirical1mer_NOCpGTpG,empiricalDistances = distances_Empirical1mer_noCpGTpG_subset,ReconstructionMantelResults =mantelResults_ReconstructionsAging,ReconstructionDistances = distances_ReconstructionsAging,label="Aging model (1-mer-CpG spectrum)",reconstructionPointColor="#1B9E77")

agingPlot


ggsave(paste0(sigfitMantelDir,"OverlayPlot.AgingModel.WithEmpiricalData.pdf"),agingPlot,height=5,width=7)


######### aging + novel model #######
mantelResults_ReconstructionsAging_PlusNovel=read.table(paste0(sigfitMantelDir,"/fit_aging_signatures_PlusNovel/aging + novel (exclue CpG>TpG).fit_aging_signatures_PlusNovel.clrDistances.mantelTestResults.sqrtcophentic.RECONSTRUCTIONS.downsampleTRUE.txt"),header=T,sep="\t",quote="")
distances_ReconstructionsAging_PlusNovel=read.table(paste0(sigfitMantelDir,"/fit_aging_signatures_PlusNovel/aging + novel (exclue CpG>TpG).fit_aging_signatures_PlusNovel.clrDistances.forplotting.sqrtcophenetic.RECONSTRUCTIONS.downsampleTRUE.txt"),header=T,sep="\t")

agingPlusNovelPlot <- PlotEmpiricalOverlayReconstructionMantel(empiricalMantelResults =mantelResults_Empirical1mer_NOCpGTpG,empiricalDistances = distances_Empirical1mer_noCpGTpG_subset,ReconstructionMantelResults =mantelResults_ReconstructionsAging_PlusNovel,ReconstructionDistances = distances_ReconstructionsAging_PlusNovel,label="Aging + Novel model (1-mer-CpG spectrum)",reconstructionPointColor="#D95F02")

agingPlusNovelPlot


ggsave(paste0(sigfitMantelDir,"OverlayPlot.AgingPlusNovelModel.WithEmpiricalData.pdf"),agingPlusNovelPlot,height=5,width=7)



########## SBS 1 + SBS 5 model #########
mantelResults_ReconstructionsSBS1_SBS5=read.table(paste0(sigfitMantelDir,"/fit_SBS1_SBS5_signatures/SBS1 + SBS5 (3-mer,includes CpG>TpG).fit_SBS1_SBS5_signatures.clrDistances.mantelTestResults.sqrtcophentic.RECONSTRUCTIONS.downsampleTRUE.txt"),header=T,sep="\t",quote="")

distances_ReconstructionsSBS1_SBS5 = read.table(paste0(sigfitMantelDir,"/fit_SBS1_SBS5_signatures/SBS1 + SBS5 (3-mer,includes CpG>TpG).fit_SBS1_SBS5_signatures.clrDistances.forplotting.sqrtcophenetic.RECONSTRUCTIONS.downsampleTRUE.txt"),header=T,sep="\t")

SBS1SBS5Plot <- PlotEmpiricalOverlayReconstructionMantel(empiricalMantelResults=mantelResults_Empirical3mer_subset,empiricalDistances = distances_Empirical3mer_subset,ReconstructionMantelResults =mantelResults_ReconstructionsSBS1_SBS5,ReconstructionDistances =distances_ReconstructionsSBS1_SBS5,label="SBS1 + SBS5 model (3-mer spectrum)" ,reconstructionPointColor="#1B9E77" )

SBS1SBS5Plot

ggsave(paste0(sigfitMantelDir,"OverlayPlot.SBS1SBS5Model.WithEmpiricalData.pdf"),SBS1SBS5Plot,height=5,width=7)

######## SBS 1 + SBS5 + novel model ############

mantelResults_ReconstructionsSBS1_SBS5_PlusNovel = read.table(paste0(sigfitMantelDir,"/fit_SBS1_SBS5_signatures_PlusNovel/SBS1 + SBS5 + novel (3-mer,includes CpG>TpG).fit_SBS1_SBS5_PlusNovel.clrDistances.mantelTestResults.sqrtcophentic.RECONSTRUCTIONS.downsampleTRUE.txt"),header=T,sep="\t",quote="")


distances_ReconstructionsSBS1_SBS5_PlusNovel = read.table(paste0(sigfitMantelDir,"/fit_SBS1_SBS5_signatures_PlusNovel/SBS1 + SBS5 + novel (3-mer,includes CpG>TpG).fit_SBS1_SBS5_PlusNovel.clrDistances.forplotting.sqrtcophenetic.RECONSTRUCTIONS.downsampleTRUE.txt"),header=T,sep="\t")


SBS1SBS5PlusNovelPlot <- PlotEmpiricalOverlayReconstructionMantel(empiricalMantelResults=mantelResults_Empirical3mer_subset,empiricalDistances = distances_Empirical3mer_subset,ReconstructionMantelResults =mantelResults_ReconstructionsSBS1_SBS5_PlusNovel,ReconstructionDistances =distances_ReconstructionsSBS1_SBS5_PlusNovel,label="SBS1 + SBS5 + Novel model (3-mer spectrum)" ,reconstructionPointColor="#D95F02" )

SBS1SBS5PlusNovelPlot

ggsave(paste0(sigfitMantelDir,"OverlayPlot.SBS1SBS5PlusNovelModel.WithEmpiricalData.pdf"),SBS1SBS5PlusNovelPlot,height=5,width=7)



########## overlay as one plot ##########

All1merModelsPlot  <- ggplot(distances_Empirical1mer_noCpGTpG_subset,aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color="empirical",shape="empirical"))+
  geom_text(data=mantelResults_Empirical1mer_NOCpGTpG,aes(x=0.5,y=Inf,vjust=1.1,label=paste0("[ Empirical data ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)),"")),color="black",size=5)+
  geom_point(data=distances_ReconstructionsAging,aes(x=sqrt_cophenetic_distance,y=distance,color="aging model",shape="aging model"))+
  geom_text(data=mantelResults_ReconstructionsAging,aes(x=0.5,y=Inf,vjust=2.5,label=paste0(" [ aging model ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)),"")),color="#1B9E77",size=5)+
  geom_point(data=distances_ReconstructionsAging_PlusNovel,aes(x=sqrt_cophenetic_distance,y=distance,color="aging + novel model",shape="aging + novel model"))+
  geom_text(data=mantelResults_ReconstructionsAging_PlusNovel,aes(x=0.5,y=Inf,vjust=3.9,label=paste0("[ aging + novel model ] r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)))),color="#D95F02",size=5)+
  scale_color_manual(name="",values=c("#D95F02","#1B9E77","black"), guide=guide_legend(reverse=T))+
  scale_shape_manual(name="",values=c(15,17,16),guide=guide_legend(reverse=T))+ # for some reason setting name the same fo rthese two scales lets you merge legend
  theme_bw()+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between species' mutation spectra")
All1merModelsPlot

ggsave(paste0(sigfitMantelDir,"OverlayPlot.All1merModelsTogether.pdf"),All1merModelsPlot,height=5,width=9)



All3merModelsPlot  <- ggplot(distances_Empirical3mer_subset,aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color="empirical",shape="empirical"))+
  geom_text(data=mantelResults_Empirical3mer_subset,aes(x=0.5,y=Inf,vjust=1.1,label=paste0("[ Empirical data ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)),"")),color="black",size=5)+
  geom_point(data=distances_ReconstructionsSBS1_SBS5,aes(x=sqrt_cophenetic_distance,y=distance,color="SBS1 + SBS5 model",shape="SBS1 + SBS5 model"))+
  geom_text(data=mantelResults_ReconstructionsSBS1_SBS5,aes(x=0.5,y=Inf,vjust=2.5,label=paste0(" [ SBS1 + SBS5 model ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(round(signif,2)),"")),color="#1B9E77",size=5)+ # using round not sci here bc want to show p vlue as 0.4
  geom_point(data=distances_ReconstructionsSBS1_SBS5_PlusNovel,aes(x=sqrt_cophenetic_distance,y=distance,color="SBS1 + SBS5 + novel model",shape="SBS1 + SBS5 + novel model"))+
  geom_text(data=mantelResults_ReconstructionsSBS1_SBS5_PlusNovel,aes(x=0.5,y=Inf,vjust=3.9,label=paste0("[ SBS1 + SBS5 + novel model ] r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)),"")),color="#D95F02",size=5)+
  scale_color_manual(name="",values=c("black","#1B9E77","#D95F02"),breaks=c("empirical","SBS1 + SBS5 model","SBS1 + SBS5 + novel model"))+
  scale_shape_manual(name="",values=c(16,17,15),breaks=c("empirical","SBS1 + SBS5 model","SBS1 + SBS5 + novel model"))+ # for some reason setting name th # breaks lets you set order
  theme_bw()+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between species' mutation spectra")
All3merModelsPlot

ggsave(paste0(sigfitMantelDir,"OverlayPlot.All3merModelsTogether.pdf"),All3merModelsPlot,height=5,width=9)
