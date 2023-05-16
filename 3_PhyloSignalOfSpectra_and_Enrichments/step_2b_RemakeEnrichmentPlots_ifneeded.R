####### plot enrichments: ########
require(ggplot2)
require(dplyr)
require(viridis)
require(ggrepel)
require(scales)
require(ggpubr)


# file with nice labels:
source("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/colors_and_labels.R") 

######################## If you want to remake any enrichment plots on your computer instead of what was made on sage ############
# SAGE location: plotdir="/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/fishers_exact_test/20220905.plots.projectedCounts.epsilon.1.sepCpGs.no.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode/" # location on sage 
fisherDir="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/compare_spectrumDist_toPhyloHammingDist/20221121.plots.SubsampledIndividualsNotProjection.epsilon.1.sepCpGs.no.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode.9999999.MantelPermutations/fishers_exact_teset_FromSage/" # my computer location
# NEW : not based on projected counts anymore; based on 5-inds per pop
# * on 20220927 I moved this dir over from 20220905.plots.projectedCounts.epsilon.1.sepCpGs.no.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode; it's based on undownsampled counts so is the same from run to run.
# name of dir should match where the files came from for the enrichment (but if you didn't change anything upstream of them being generated it's fine not to rerun this if you're just changing downstream stuff )

############## READ IN FISHERS EXACT TEST RESULTS  #################
all7merSpectra_FishersExactTestResults <- read.table(paste0(fisherDir,"FISHERS.EXACT.TEST.RESULTS.7merSpectrum.notrescaled_notdownsampled_notregularized.txt"),header=T)

all5merSpectra_FishersExactTestResults <- read.table(paste0(fisherDir,"FISHERS.EXACT.TEST.RESULTS.5merSpectrum.notrescaled_notdownsampled_notregularized.txt"),header=T)

all3merSpectra_FishersExactTestResults <- read.table(paste0(fisherDir,"FISHERS.EXACT.TEST.RESULTS.3merSpectrum.notrescaled_notdownsampled_notregularized.txt"),header=T)
# 3mer results excluding cpgs:

all3merSpectra_FishersExactTestResults_excludingCpG <- read.table(paste0(fisherDir,"FISHERS.EXACT.TEST.RESULTS.3merSpectrum.EXCLUDINGCpGs.notrescaled_notdownsampled_notregularized.txt"),header=T)


######### determine what enrichment level enrichment 3mer CpGs>TpGs have and use as a threshold for finding interesting enrichments ###########
# 3mer CpGs:
CpG_TpG_3merTypes=c("ACG.ATG","CCG.CTG","GCG.GTG","TCG.TTG")

#all3merSpectra_FishersExactTestResults[all3merSpectra_FishersExactTestResults$mutation_label %in% CpG_TpG_3merTypes,] %>% View()

all3merSpectra_FishersExactTestResults_GetCpGEnrichment <- data.frame(all3merSpectra_FishersExactTestResults)

all3merSpectra_FishersExactTestResults_GetCpGEnrichment$isCpG_TpG3mer <- "no"

all3merSpectra_FishersExactTestResults_GetCpGEnrichment[all3merSpectra_FishersExactTestResults_GetCpGEnrichment$mutation_label %in% CpG_TpG_3merTypes, ]$isCpG_TpG3mer <- "yes"

# sum up
CpGCountsPerSpp <- all3merSpectra_FishersExactTestResults_GetCpGEnrichment %>%
  group_by(species,population,label) %>%
  filter(isCpG_TpG3mer=="yes") %>%
  summarise(total_mutations_CpG=sum(as.numeric(total_mutations)),total_target_count_CpG=sum(as.numeric(total_target_count)),count_div_by_targets_CpG=total_mutations_CpG/total_target_count_CpG) %>%
  ungroup()

# need to add in per species C>T rates:
CTCoutnsDivByTargetsPerSpp <- all3merSpectra_FishersExactTestResults_GetCpGEnrichment %>%
  group_by(species,population,label) %>%
  filter(mutation_1mer=="C.T") %>% 
  select(species,population,label,count_div_by_targets_1mer) %>%
  unique() %>%
  ungroup()

# merge them:
CpG_Rates <- merge(CpGCountsPerSpp,CTCoutnsDivByTargetsPerSpp,by=c("species","population","label"))

CpG_Rates$CpGEnrichmentOverCT <- CpG_Rates$count_div_by_targets_CpG/CpG_Rates$count_div_by_targets_1mer

write.table(CpG_Rates,paste0(fisherDir,"CpG.EnrichmentRates.PerSpecies.txt"),row.names = F,quote=F,sep="\t")
################## add cpg rate info #################
# so for each species let's find out which of their kmers exceed their respective cpg threshold
all3merSpectra_FishersExactTestResults_PlusCpGRateInfo <- merge(all3merSpectra_FishersExactTestResults,select(CpG_Rates,c("population","species","label","CpGEnrichmentOverCT")),by=c("species","population","label"))

all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$greaterOrEqualToSpeciesCpGRate <- "no"
all3merSpectra_FishersExactTestResults_PlusCpGRateInfo[all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$kmerEnrichmentOverCentralBP>=all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$CpGEnrichmentOverCT,]$greaterOrEqualToSpeciesCpGRate <- 'yes'

# find any kmers that aren't cpg>tpg but have greater enrichment than cpgs:

pvalue_threshold_3mers = 0.05/dim(all3merSpectra_FishersExactTestResults_PlusCpGRateInfo)[1] # number o mut types * number of species

# label if is cpg?tpg mutation: 

all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$isCpG_TpG <- "no"


all3merSpectra_FishersExactTestResults_PlusCpGRateInfo[grepl("CG",all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$ancestral3mer) & all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$mutation_1mer =="C.T", ]$isCpG_TpG <- "yes"


greaterThanCpG_3mers <- all3merSpectra_FishersExactTestResults_PlusCpGRateInfo[all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$greaterOrEqualToSpeciesCpGRate=="yes" & all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$isCpG_TpG=="no" & all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$p_value_fisher <pvalue_threshold_3mers,]

greaterThanCpG_3mers 
### ^^^ there are no 3mers that are higher than CpGs. ###


###### 5mers #######


all5merSpectra_FishersExactTestResults_PlusCpGRateInfo <- merge(all5merSpectra_FishersExactTestResults,select(CpG_Rates,c("population","species","label","CpGEnrichmentOverCT")),by=c("species","population","label"))
all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$greaterOrEqualToSpeciesCpGRate <- "no"

pvalue_threshold_5mers = 0.05/dim(all5merSpectra_FishersExactTestResults_PlusCpGRateInfo)[1]

all5merSpectra_FishersExactTestResults_PlusCpGRateInfo[all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$kmerEnrichmentOverCentralBP>=all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$CpGEnrichmentOverCT,]$greaterOrEqualToSpeciesCpGRate <- 'yes'

all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$ancestral3mer <- substr(all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$ancestral5mer,2,4)

all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$isCpG_TpG <- "no"


all5merSpectra_FishersExactTestResults_PlusCpGRateInfo[grepl("CG",all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$ancestral3mer) & all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$mutation_1mer =="C.T", ]$isCpG_TpG <- "yes"




greaterThanCpG_5mers <- all5merSpectra_FishersExactTestResults_PlusCpGRateInfo[all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$greaterOrEqualToSpeciesCpGRate=="yes" &  all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$isCpG_TpG=="no" & all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$p_value_fisher < pvalue_threshold_5mers,]

write.table(greaterThanCpG_5mers,paste0(fisherDir,"5mersWithEnrichmentsGreaterThanCpGThatAreSignificant.txt"),row.names = F,quote=F,sep="\t")

# p-value significance threshold 
# this code : (grepl("CG",all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$ancestral3mer) & all5merSpectra_FishersExactTestResults_PlusCpGRateInfo$mutation_1mer=="C.T") excludes CpG>TpG mutations but allows other CpG mutations to stay in 

# restrict to species to include and pvalue threshold 
### get 7mer results
all7merSpectra_FishersExactTestResults_PlusCpGRateInfo <- merge(all7merSpectra_FishersExactTestResults,select(CpG_Rates,c("population","species","label","CpGEnrichmentOverCT")),by=c("species","population","label"))

all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$greaterOrEqualToSpeciesCpGRate <- "no"
all7merSpectra_FishersExactTestResults_PlusCpGRateInfo[all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$kmerEnrichmentOverCentralBP>=all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$CpGEnrichmentOverCT,]$greaterOrEqualToSpeciesCpGRate <- 'yes'

all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$ancestral3mer <- substr(all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$ancestral7mer,3,5)


all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$isCpG_TpG <- "no"


all7merSpectra_FishersExactTestResults_PlusCpGRateInfo[grepl("CG",all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$ancestral3mer) & all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$mutation_1mer =="C.T", ]$isCpG_TpG <- "yes"



pvalue_threshold_7mers = 0.05/dim(all7merSpectra_FishersExactTestResults_PlusCpGRateInfo)[1]

greaterThanCpG_7mers <- all7merSpectra_FishersExactTestResults_PlusCpGRateInfo[all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$greaterOrEqualToSpeciesCpGRate=="yes" & all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$isCpG_TpG=="no" & all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$p_value_fisher < pvalue_threshold_7mers,]

length(unique(greaterThanCpG_7mers$mutation_label))

write.table(greaterThanCpG_7mers,paste0(fisherDir,"7mersWithEnrichmentsGreaterThanCpGThatAreSignificant.txt"),row.names = F,quote=F,sep="\t")

# get some stats
# greaterThanCpG_7mers %>%
#   group_by(species,population,label,mutation_1mer) %>%
#   summarise(totalPerCentralBP=n() ) %>%
#   ungroup() %>%
#   View()
# 
# # how many per species
# greaterThanCpG_7mers %>%
#   group_by(species,population,label) %>%
#   summarise(totalPerSpecies=n() ) %>%
#   ungroup() %>%
#   View()

# do all have elevated TTTAAA?
greaterThanCpG_7mers[greaterThanCpG_7mers$mutation_label=="TTTAAAA.TTTTAAA",]
# could also just focus on species that are in our plots:


write.table(greaterThanCpG_7mers,paste0(fisherDir,"7mersWithEnrichmentsGreaterThanCpGThatAreSignificant.txt"),row.names = F,quote=F,sep="\t")

#### figure out better way to summarise . don't have to put in that much. 
# maybe figure out which kmers are consistently cool/hot


######### want to make U shaped plots like Kelleys where x axis is OR and y axis is -log10 pvalue #########
# color by central 1mer and facet over species
# relabel p=0 to be <1e-322
############## u shaped plot for 3mers ##########
speciesToIncludeInPlot=c("humans_EUR","humans_EAS","humans_AMR","Mus_musculus_Mmd","Mus_spretus_Ms","brown_bear_ABC","polar_bear_PB","vaquita","fin_whale_GOC","humans_AFR","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","Pongo_abelii","Pongo_pygmaeus","wolves")


# order for plot labels (includes somes pecies that aren't in species to include list )
#labelorder <- c("Mus_spretus_Ms","Mus_musculus_Mmd","Mus_musculus_Mmc","Mus_musculus_Mmm","Pongo_abelii","Pongo_pygmaeus","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","humans_AFR","humans_EUR","humans_SAS","humans_EAS","humans_AMR","wolves","brown_bear_ABC","brown_bear_EUR","polar_bear_PB","fin_whale_GOC","fin_whale_ENP","vaquita") # update this if you change species
# add nice labels: (which was sourced from colors and labels file)
all3merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels <- merge(all3merSpectra_FishersExactTestResults_PlusCpGRateInfo,niceLabels,by="label",all.x = T,all.y=F)

niceLabelorder <- c("mouse (Ms)" ,"mouse (Mmd)","mouse (Mmc)"  , "mouse (Mmm)"  ,   "orangutan (Sumatran)","orangutan (Bornean)","gorilla", "bonobo"   ,"chimpanzee" ,"human (AFR)"  , "human (EUR)"  , "human (SAS)" , "human (EAS)"    ,"human (AMR)"  ,"wolf"  ,"brown bear (ABC)","brown bear (EUR)" , "polar bear" ,"fin whale (GOC)",   "fin whale (ENP)" ,"vaquita")

all3merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$niceLabel <-factor(all3merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$niceLabel,levels=niceLabelorder)


# making bigger dots for 3mers
plotEnrichments_UShapedPlot_3mer <- ggplot(all3merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels,aes(x=kmerEnrichmentOverCentralBP,y=-log10(p_value_fisher),color=mutation_1mer,shape=isCpG_TpG))+
  geom_point(size=4)+
  facet_wrap(~niceLabel)+
  scale_x_log10()+
  theme_bw()+
  scale_color_brewer(palette = "Set2",name="central 1-mer")+
  geom_vline(aes(xintercept = CpGEnrichmentOverCT),color="red",linetype="dashed")+ # species specific
  geom_hline(yintercept = pvalue_threshold_3mers,linetype="dashed")+ #dataset specific+
  #geom_text_repel(data=all3merSpectra_FishersExactTestResults_PlusCpGRateInfo[all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$p_value_fisher<pvalue_threshold_3mers & all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$kmerEnrichmentOverCentralBP>all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$CpGEnrichmentOverCT & all3merSpectra_FishersExactTestResults_PlusCpGRateInfo$isCpG_TpG=="no",],aes(label=mutation_label),size=3,color="black") + # label ones that are sig and have greater than CpG enrichments
  #ggtitle("dashed line is significance threshold; red line is CpG>T/C>T enrichment level (species-specific)\nlabelling non CpG>TpG sites")+
  scale_shape_manual(values=c(16,8),name="is CpG>TpG?")+
  theme(strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  xlab("enrichment over central 1-mer mutation rate (log10 scale)")

plotEnrichments_UShapedPlot_3mer
ggsave(paste0(fisherDir,"3merEnrichment.NewPlotUshaped.allSppAllPops.pdf"),plotEnrichments_UShapedPlot_3mer,height=12,width=14)
ggsave(paste0(fisherDir,"3merEnrichment.NewPlotUshaped.allSppAllPops.png"),plotEnrichments_UShapedPlot_3mer,height=12,width=14,dpi = 300) # also save as png


############## u shaped plot for 5mers ##########

all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels <- merge(all5merSpectra_FishersExactTestResults_PlusCpGRateInfo,niceLabels,by="label",all.x = T,all.y=F)


all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$niceLabel <-factor(all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$niceLabel,levels=niceLabelorder)

plotEnrichments_UShapedPlot_5mer <- ggplot(all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels,aes(x=kmerEnrichmentOverCentralBP,y=-log10(p_value_fisher),color=mutation_1mer,shape=isCpG_TpG))+
  geom_point(size=2,alpha=0.6)+
  facet_wrap(~niceLabel)+
  theme_bw()+
  scale_x_log10()+
  #scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2",name="central 1-mer")+
  geom_vline(aes(xintercept = CpGEnrichmentOverCT),color="red",linetype="dashed")+ # species specific
  geom_hline(yintercept = pvalue_threshold_5mers,linetype="dashed")+ #dataset specific+
  geom_text_repel(data=all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels[all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$p_value_fisher<pvalue_threshold_5mers & all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$kmerEnrichmentOverCentralBP>all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$CpGEnrichmentOverCT & all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$isCpG_TpG=="no",],aes(label=mutation_label,color=mutation_1mer),size=3,min.segment.length = 0,direction="y") + # label ones that are sig and have greater than CpG enrichments
  #ggtitle("dashed line is significance threshold; red line is CpG>T/C>T enrichment level (species-specific)\nlabelling non CpG>TpG sites")+
  scale_shape_manual(values=c(16,8),name="is CpG>TpG?")+
  theme(strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  xlab("enrichment over central 1-mer mutation rate (log10 scale)")

plotEnrichments_UShapedPlot_5mer
ggsave(paste0(fisherDir,"5merEnrichment.NewPlotUshaped.allSppAllPops.pdf"),plotEnrichments_UShapedPlot_5mer,height=12,width=14)
ggsave(paste0(fisherDir,"5merEnrichment.NewPlotUshaped.allSppAllPops.png"),plotEnrichments_UShapedPlot_5mer,height=12,width=14,dpi = 300) # also save as png


############## u shaped plot for 7mers ##########


all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels <- merge(all7merSpectra_FishersExactTestResults_PlusCpGRateInfo,niceLabels,by="label",all.x = T,all.y=F)


all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$niceLabel <-factor(all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$niceLabel,levels=niceLabelorder)


#all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$label <- factor(all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$label,levels=labelorder)

plotEnrichments_UShapedPlot_7mer <- ggplot(all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels,aes(x=kmerEnrichmentOverCentralBP,y=-log10(p_value_fisher),color=mutation_1mer,shape=isCpG_TpG))+
  geom_point(size=2,alpha=0.6)+
  facet_wrap(~niceLabel)+
  scale_x_log10()+
  theme_bw()+
  scale_color_brewer(palette = "Set2",name="central 1-mer")+
  geom_vline(aes(xintercept = CpGEnrichmentOverCT),color="red",linetype="dashed")+ # species specific
  geom_hline(yintercept = pvalue_threshold_7mers,linetype="dashed")+ #dataset specific+
  geom_text_repel(data=all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels[all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$p_value_fisher<pvalue_threshold_7mers & all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$kmerEnrichmentOverCentralBP>all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$CpGEnrichmentOverCT & all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$isCpG_TpG=="no",],aes(label=mutation_label),size=2,min.segment.length = 0,direction="y",max.overlaps = 15) + # label ones that are sig and have greater than CpG enrichments
  #ggtitle("dashed line is significance threshold; red line is CpG>T/C>T enrichment level (species-specific)\nlabelling non CpG>TpG sites")+
  scale_shape_manual(values=c(16,8),name="is CpG>TpG?")+
  theme(strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  xlab("enrichment over central 1-mer mutation rate (log10 scale)")

plotEnrichments_UShapedPlot_7mer
ggsave(paste0(fisherDir,"7merEnrichment.NewPlotUshaped.allSppAllPops.pdf"),plotEnrichments_UShapedPlot_7mer,height=12,width=14)
ggsave(paste0(fisherDir,"7merEnrichment.NewPlotUshaped.allSppAllPops.png"),plotEnrichments_UShapedPlot_7mer,height=12,width=14,dpi = 300) # also save as png


############## make main text figure with 5/7mers from a small subset of species #######
mainTextSpecies=c("vaquita","humans_AFR","fin_whale_GOC","Mus_musculus_Mmd") # examples
# want to plot 5 and 7mer together 

all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$id <- "7-mer spectrum"
all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$id <- "5-mer spectrum"

# will grid arrange them. 

## want to make inf p-values say ">300"

volcanoPlotForMainText_5mer <- ggplot(all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels[all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$label %in% mainTextSpecies,],aes(x=kmerEnrichmentOverCentralBP,y=-log10(p_value_fisher),color=mutation_1mer,shape=isCpG_TpG))+
  geom_point(size=2)+
  facet_grid(id~niceLabel,switch="y")+
  scale_x_log10()+
  theme_bw()+
  scale_color_brewer(palette = "Set2",name="central 1-mer")+
  geom_vline(aes(xintercept = CpGEnrichmentOverCT),color="red",linetype="dashed")+ # species specific
  geom_hline(yintercept = pvalue_threshold_5mers,linetype="dashed")+ #dataset specific+
  geom_text_repel(data=all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels[all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$label %in% mainTextSpecies & all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$p_value_fisher<pvalue_threshold_5mers & all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$kmerEnrichmentOverCentralBP>all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$CpGEnrichmentOverCT & all5merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$isCpG_TpG=="no",],aes(label=mutation_label),size=4,min.segment.length = 0) + # label ones that are sig and have greater than CpG enrichments
  #ggtitle("dashed line is significance threshold; red line is CpG>T/C>T enrichment level (species-specific)\nlabelling non CpG>TpG sites")+
  scale_shape_manual(values=c(16,8),name="is CpG>TpG?")+
  theme(strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  xlab("enrichment over central 1-mer mutation rate (log10 scale)")

volcanoPlotForMainText_5mer
ggsave(paste0(fisherDir,"mainTextVolcanoGrid.5meronly.pdf"),volcanoPlotForMainText_5mer,height=5,width=10)
ggsave(paste0(fisherDir,"mainTextVolcanoGrid.5meronly.png"),volcanoPlotForMainText_5mer,height=5,width=10,dpi = 300) # also save as png



volcanoPlotForMainText_7mer <- ggplot(all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels[all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$label %in% mainTextSpecies,],aes(x=kmerEnrichmentOverCentralBP,y=-log10(p_value_fisher),color=mutation_1mer,shape=isCpG_TpG,))+
  geom_point(size=2)+
  facet_grid(id~niceLabel,switch="y")+
  scale_x_log10()+
  theme_bw()+
  scale_color_brewer(palette = "Set2",name="central 1-mer")+
  geom_vline(aes(xintercept = CpGEnrichmentOverCT),color="red",linetype="dashed")+ # species specific
  geom_hline(yintercept = pvalue_threshold_7mers,linetype="dashed")+ #dataset specific+
  geom_text_repel(data=all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels[all7merSpectra_FishersExactTestResults_PlusCpGRateInfo$label %in% mainTextSpecies & all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$p_value_fisher<pvalue_threshold_7mers & all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$kmerEnrichmentOverCentralBP>all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$CpGEnrichmentOverCT & all7merSpectra_FishersExactTestResults_PlusCpGRateInfo_plusLabels$isCpG_TpG=="no",],aes(label=mutation_label),size=3,min.segment.length = 0,direction="y") + # label ones that are sig and have greater than CpG enrichments
  #ggtitle("dashed line is significance threshold; red line is CpG>T/C>T enrichment level (species-specific)\nlabelling non CpG>TpG sites")+
  scale_shape_manual(values=c(16,8),name="is CpG>TpG?")+
  theme(strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  xlab("enrichment over central 1-mer mutation rate (log10 scale)")

volcanoPlotForMainText_7mer
ggsave(paste0(fisherDir,"mainTextVolcanoGrid.7meronly.pdf"),volcanoPlotForMainText_7mer,height=5,width=10)
ggsave(paste0(fisherDir,"mainTextVolcanoGrid.7meronly.png"),volcanoPlotForMainText_7mer,height=5,width=10,dpi = 300) # also save as png

mainTextVolcanoGrid <- ggarrange(volcanoPlotForMainText_5mer,volcanoPlotForMainText_7mer,common.legend = T,nrow = 2,legend = "right")
# then grid arrange them into one 
mainTextVolcanoGrid
ggsave(paste0(fisherDir,"mainTextVolcanoGrid.pdf"),mainTextVolcanoGrid,height=8,width=11)
ggsave(paste0(fisherDir,"mainTextVolcanoGrid.png"),mainTextVolcanoGrid,height=8,width=11,dpi = 300) # also save as png


############## analysis to figure out which kmers are cool/hot across all species #########
###### 3mers ############
all3merSpectra_FishersExactTestResults_SIGNIFICANTONLY <- all3merSpectra_FishersExactTestResults[all3merSpectra_FishersExactTestResults$p_value_fisher < pvalue_threshold_3mers,] 



significantlyGreaterThan1InAllSpeciesAndPops_3mer <- all3merSpectra_FishersExactTestResults_SIGNIFICANTONLY %>%
  group_by(mutation_label) %>%
  filter(kmerEnrichmentOverCentralBP>1) %>%
  summarise(totalGT1=n()) %>%
  filter(totalGT1==21) %>% # 21 species/pops
  ungroup()

significantlyGreaterThan1InAllSpeciesAndPops_3mer$ancestral3mer <- substr(significantlyGreaterThan1InAllSpeciesAndPops_3mer$mutation_label,1,3)

significantlyGreaterThan1InAllSpeciesAndPops_3mer$central1mer <- paste0(substr(significantlyGreaterThan1InAllSpeciesAndPops_3mer$mutation_label,2,2),".",substr(significantlyGreaterThan1InAllSpeciesAndPops_3mer$mutation_label,6,6))

significantlyGreaterThan1InAllSpeciesAndPops_3mer$isCpG_TpG <- "no"

significantlyGreaterThan1InAllSpeciesAndPops_3mer[grepl("CG",significantlyGreaterThan1InAllSpeciesAndPops_3mer$ancestral3mer) & significantlyGreaterThan1InAllSpeciesAndPops_3mer$central1mer =="C.T", ]$isCpG_TpG <- "yes"

dim(significantlyGreaterThan1InAllSpeciesAndPops_3mer)

sum(significantlyGreaterThan1InAllSpeciesAndPops_3mer$isCpG_TpG=="yes")

write.table(all3merSpectra_FishersExactTestResults_SIGNIFICANTONLY[all3merSpectra_FishersExactTestResults_SIGNIFICANTONLY$mutation_label %in% significantlyGreaterThan1InAllSpeciesAndPops_3mer$mutation_label,],paste0(fisherDir,"3mers.significantlyGT1.inallspeciesandpops.txt"),row.names=F,quote=F,sep="\t")


significantlyLessThan1InAllSpeciesAndPops_3mer <- all3merSpectra_FishersExactTestResults_SIGNIFICANTONLY %>%
  group_by(mutation_label) %>%
  filter(kmerEnrichmentOverCentralBP<1) %>%
  summarise(totalLT1=n()) %>%
  filter(totalLT1==21) %>%
  ungroup()

significantlyLessThan1InAllSpeciesAndPops_3mer$ancestral3mer <- substr(significantlyLessThan1InAllSpeciesAndPops_3mer$mutation_label,1,3)

significantlyLessThan1InAllSpeciesAndPops_3mer$central1mer <- paste0(substr(significantlyLessThan1InAllSpeciesAndPops_3mer$mutation_label,2,2),".",substr(significantlyLessThan1InAllSpeciesAndPops_3mer$mutation_label,6,6))


write.table(all3merSpectra_FishersExactTestResults_SIGNIFICANTONLY[all3merSpectra_FishersExactTestResults_SIGNIFICANTONLY$mutation_label %in% significantlyLessThan1InAllSpeciesAndPops_3mer$mutation_label,],paste0(fisherDir,"3mers.significantlyLT1.inallspeciesandpops.txt"),row.names=F,quote=F,sep="\t")


View(significantlyLessThan1InAllSpeciesAndPops_3mer)

#### what about switching 3mers?
significantlyGreaterThanInSomeAndLessThanSomeInAllSpeciesAndPops_3mer <- all3merSpectra_FishersExactTestResults_SIGNIFICANTONLY %>%
  group_by(mutation_label) %>%
  filter(kmerEnrichmentOverCentralBP<1) %>%
  ungroup()


######### 5mers #######
# should be significant
all5merSpectra_FishersExactTestResults_SIGNIFICANTONLY <- all5merSpectra_FishersExactTestResults[all5merSpectra_FishersExactTestResults$p_value_fisher < pvalue_threshold_5mers,] 
# this contains >1 and <1 than 1. 

significantlyGreaterThan1InAllSpeciesAndPops_5mer <- all5merSpectra_FishersExactTestResults_SIGNIFICANTONLY %>%
  group_by(mutation_label) %>%
  filter(kmerEnrichmentOverCentralBP>1) %>%
  summarise(totalGT1=n()) %>%
  filter(totalGT1==21) %>%
  ungroup

significantlyGreaterThan1InAllSpeciesAndPops_5mer$ancestral3mer <- substr(significantlyGreaterThan1InAllSpeciesAndPops_5mer$mutation_label,2,4)

significantlyGreaterThan1InAllSpeciesAndPops_5mer$central1mer <- paste0(substr(significantlyGreaterThan1InAllSpeciesAndPops_5mer$mutation_label,3,3),".",substr(significantlyGreaterThan1InAllSpeciesAndPops_5mer$mutation_label,9,9))

significantlyGreaterThan1InAllSpeciesAndPops_5mer$isCpG_TpG <- "no"

significantlyGreaterThan1InAllSpeciesAndPops_5mer[grepl("CG",significantlyGreaterThan1InAllSpeciesAndPops_5mer$ancestral3mer) & significantlyGreaterThan1InAllSpeciesAndPops_5mer$central1mer =="C.T", ]$isCpG_TpG <- "yes"
dim(significantlyGreaterThan1InAllSpeciesAndPops_5mer)

sum(significantlyGreaterThan1InAllSpeciesAndPops_5mer$isCpG_TpG=="yes")

write.table(all5merSpectra_FishersExactTestResults_SIGNIFICANTONLY[all5merSpectra_FishersExactTestResults_SIGNIFICANTONLY$mutation_label %in% significantlyGreaterThan1InAllSpeciesAndPops_5mer$mutation_label,],paste0(fisherDir,"5mers.significantlyGT1.inallspeciesandpops.txt"),row.names=F,quote=F,sep="\t")

######## now get ones that are consistently lower than 1 across species ###########
significantlyLessThan1InAllSpeciesAndPops_5mer <- all5merSpectra_FishersExactTestResults_SIGNIFICANTONLY %>%
  group_by(mutation_label) %>%
  filter(kmerEnrichmentOverCentralBP<1) %>%
  summarise(totalLT1=n()) %>%
  filter(totalLT1==21) %>%
  ungroup

significantlyLessThan1InAllSpeciesAndPops_5mer$ancestral3mer <- substr(significantlyLessThan1InAllSpeciesAndPops_5mer$mutation_label,2,4)

significantlyLessThan1InAllSpeciesAndPops_5mer$central1mer <- paste0(substr(significantlyLessThan1InAllSpeciesAndPops_5mer$mutation_label,3,3),".",substr(significantlyLessThan1InAllSpeciesAndPops_5mer$mutation_label,9,9))

# many are likely C>T because we find that cpgs drie the rest of c.ts down. 

sum(significantlyLessThan1InAllSpeciesAndPops_5mer$central1mer=="C.T")

write.table(all5merSpectra_FishersExactTestResults_SIGNIFICANTONLY[all5merSpectra_FishersExactTestResults_SIGNIFICANTONLY$mutation_label %in% significantlyLessThan1InAllSpeciesAndPops_5mer$mutation_label,],paste0(fisherDir,"5mers.significantlyLT1.inallspeciesandpops.5mer.txt"),row.names=F,quote=F,sep="\t")



#################### plotting functions #######################
# code to run fisher's test itself is in previous script to run on sage 
speciesToIncludeInPlot=c("humans_EUR","humans_EAS","humans_AMR","Mus_musculus_Mmd","Mus_spretus_Ms","brown_bear_ABC","polar_bear_PB","vaquita","fin_whale_GOC","humans_AFR","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","Pongo_abelii","Pongo_pygmaeus","wolves")

plotEnrichmentsAndPValues <- function(resultsOfFishersExactTest,speciesToIncludeInPlot,fisherDir){
  
  spectrumName=unique(resultsOfFishersExactTest$spectrumName) # this includes whether CPGs were included or excluded
  # order for plot labels (includes somes pecies that aren't in species to include list )
  labelorder <- c("Mus_spretus_Ms","Mus_musculus_Mmd","Mus_musculus_Mmc","Mus_musculus_Mmm","Pongo_abelii","Pongo_pygmaeus","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","humans_AFR","humans_EUR","humans_SAS","humans_EAS","humans_AMR","wolves","brown_bear_ABC","brown_bear_EUR","polar_bear_PB","fin_whale_GOC","fin_whale_ENP","vaquita") # update this if you change species
  labelorder <- labelorder[labelorder %in% unique(resultsOfFishersExactTest$label)] # subset to just the species I want
   
  # order species labels:
  resultsOfFishersExactTest$label <- factor(resultsOfFishersExactTest$label,levels=labelorder)
  
  
  significanceCutoff = 0.05/dim(resultsOfFishersExactTest)[1] # 
  
  ###### plot straight enrichment odds ratios ########
  enrichmentPlot1 <- ggplot(resultsOfFishersExactTest[resultsOfFishersExactTest$label %in% speciesToIncludeInPlot,],aes(x=mutation_label,y=kmerEnrichmentOverCentralBP,color=log10(total_target_count)))+
    geom_point()+
    facet_grid(mutation_1mer~label)+
    theme_bw()+theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
    scale_color_viridis()+
    theme(axis.text.x=element_blank())+
    #geom_text_repel(data=resultsOfFishersExactTest[resultsOfFishersExactTest$kmerEnrichmentOverCentralBP_DOWNSAMPLED>enrichmentThreshold & resultsOfFishersExactTest$label %in% speciesToIncludeInPlot,],aes(label=mutation_label),size=3)+
    geom_hline(yintercept = 1,color="grey")+
    xlab("kmer")+
    ggtitle(paste0(spectrumName,"Enrichment of kmer relative rate over 1-mer relative rate;\nnot rescaled, downsampled, or regularized"))+
    theme(legend.position = "none")
  ggsave(paste0(fisherDir,"enrichmentplot.",spectrumName,".notRescaledDownsampledOrRegularized.pdf"),enrichmentPlot1,height=12,width=18)
  
  

    ##### plot p-values from fisher's exact test ############
    enrichmentPlot2 <- ggplot(resultsOfFishersExactTest[resultsOfFishersExactTest$label %in% speciesToIncludeInPlot,],aes(x=mutation_label,y=-log10(p_value_fisher),color=log2(kmerEnrichmentOverCentralBP)))+
      geom_point()+
      facet_grid(mutation_1mer~label)+
      theme_bw()+theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
      #scale_color_viridis()+
      scale_color_gradient2(midpoint=0,low="blue",high="red",mid="gray")+
      theme(axis.text.x=element_blank())+
      geom_text_repel(data=resultsOfFishersExactTest[resultsOfFishersExactTest$p_value_fisher<significanceCutoff & resultsOfFishersExactTest$label %in% speciesToIncludeInPlot,],aes(label=mutation_label),size=3)+
      geom_hline(yintercept = -log10(significanceCutoff),color="black")+
      xlab("kmer")+
      ggtitle(paste0("p-values from Fisher's exact test; Enrichment of kmer relative rate over 1-mer relative rate;\nnot rescaled, downsampled or regularized;\nsignificance cutoff = ", scientific(significanceCutoff,2)))+
      theme(legend.position = "none")
    ggsave(paste0(fisherDir,"enrichmentplot.",spectrumName,".PVALUESFROMFISHERSEXACTTEST.notDownsampledOrRegularized.pdf"),enrichmentPlot2,height=12,width=18)

  # not returning anything, this is just a plotting function 
}


############# >>>> plot 3mers <<<<< ##############


plotEnrichmentsAndPValues(all3merSpectra_FishersExactTestResults,fisherDir=fisherDir,speciesToIncludeInPlot = speciesToIncludeInPlot)

####### >>>> plot 3mers with CPG 3mers excluded  <<<<< #############

plotEnrichmentsAndPValues(all3merSpectra_FishersExactTestResults_excludingCpG,fisherDir=fisherDir,speciesToIncludeInPlot = speciesToIncludeInPlot)


######## >>>>>> plots exploring TCC pulse ORs <<<<<<<<<<<< ######
tccPlot1 <- ggplot(all3merSpectra_FishersExactTestResults[all3merSpectra_FishersExactTestResults$mutation_1mer=="C.T" & all3merSpectra_FishersExactTestResults$species=="humans",],aes(x=mutation_label,y=kmerEnrichmentOverCentralBP,fill=label,color=label))+
  theme_bw()+
  geom_hline(yintercept = 1)+
  geom_col(position="dodge")+
  theme(axis.text.x = element_text(angle=90))
ggsave(paste0(fisherDir,"kmerEnrichment.C.T.HumansOnly.OddsRatio.pdf"),tccPlot1)

# without cpgs:
tccPlot2 <- ggplot(all3merSpectra_FishersExactTestResults_excludingCpG[all3merSpectra_FishersExactTestResults_excludingCpG$mutation_1mer=="C.T" &all3merSpectra_FishersExactTestResults_excludingCpG$species=="humans",],aes(x=mutation_label,y=kmerEnrichmentOverCentralBP,fill=label,color=label))+
  theme_bw()+
  geom_hline(yintercept = 1)+
  geom_col(position="dodge")+
  theme(axis.text.x = element_text(angle=90))
ggsave(paste0(fisherDir,"kmerEnrichment.C.T.HumansOnly.OddsRatio.ExcludingCpGs.pdf"),tccPlot2)

# try as a heatmap
tccPlot2b <- ggplot(all3merSpectra_FishersExactTestResults_excludingCpG[all3merSpectra_FishersExactTestResults_excludingCpG$mutation_1mer=="C.T" & all3merSpectra_FishersExactTestResults_excludingCpG$species=="humans",],aes(y=mutation_label,x=label,fill=kmerEnrichmentOverCentralBP))+
  geom_tile()+
  scale_fill_gradient2(midpoint=1,low="blue",high="red",mid="white")+
  geom_text(aes(label=round(kmerEnrichmentOverCentralBP,2)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))+
  ggtitle("odds ratio of 3mer rate divided by 1mer rate (excluding CpGs)")
ggsave(paste0(fisherDir,"kmerEnrichment.C.T.HumansOnly.OddsRatio.ExcludingCpGs.heatmap.pdf"),tccPlot2b)

tccPlot2c <- ggplot(all3merSpectra_FishersExactTestResults_excludingCpG[all3merSpectra_FishersExactTestResults_excludingCpG$mutation_1mer=="C.T",],aes(y=mutation_label,x=label,fill=kmerEnrichmentOverCentralBP))+
  geom_tile()+
  scale_fill_gradient2(midpoint=1,low="blue",high="red",mid="white")+
  geom_text(aes(label=round(kmerEnrichmentOverCentralBP,2)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))+
  ggtitle("odds ratio of 3mer rate divided by 1mer rate (excluding CpGs)")
ggsave(paste0(fisherDir,"kmerEnrichment.C.T.AllSpecies.OddsRatio.ExcludingCpGs.heatmap.pdf"),tccPlot2c,height=8,width=16)

############## >>>> plot 5/7mers <<<<<<<<<<< ##############

###### 5mer

plotEnrichmentsAndPValues(all5merSpectra_FishersExactTestResults,fisherDir=fisherDir,speciesToIncludeInPlot = speciesToIncludeInPlot)

###### 7mer: (VERY SLOW)

# then plot results:
plotEnrichmentsAndPValues(all7merSpectra_FishersExactTestResults,fisherDir=fisherDir,speciesToIncludeInPlot = speciesToIncludeInPlot)

