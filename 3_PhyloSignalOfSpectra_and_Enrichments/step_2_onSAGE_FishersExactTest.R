####### plot enrichments: ########
require(ggplot2)
require(dplyr)
require(viridis)
require(ggrepel)
require(scales)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (wd).n", call.=FALSE)
} 



plotdir=paste0(args[1],"/") # get from command line and add / in case that wasn't supplied
# name of dir should match where the files came from for the enrichment (but if you didn't change anything upstream of them being generated it's fine not to rerun this if you're just changing downstream stuff )
fisherDir=paste0(plotdir,"fishers_exact_test/")
dir.create(fisherDir,showWarnings = F)
############## READ IN SUMMED UP SPECTRA WITH TARGET INFO FOR 3,5,7mers #################
all7merSpectra_plusTargets <- read.table(paste0(plotdir,"all7merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),header=T)

all5merSpectra_plusTargets <- read.table(paste0(plotdir,"all5merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),header=T)

all3merSpectra_plusTargets <- read.table(paste0(plotdir,"all3merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),header=T)
# don't need 1mers

#################### set up functions #######################
fishersexacttest_nondownsampledData <- function(spectradf,fisherDir,spectrumName,excludeCpGs="no"){

  if(excludeCpGs=="yes"){
    print("excluding CpGs")

    spectradf <- spectradf[!grepl(pattern = "CG",spectradf$ancestral3mer),]
    spectrumName <- paste0(spectrumName,".EXCLUDINGCpGs")
  }
  print("am running enrichments tests on non-scaled by human targets, non-downsampled, non-regularized counts")
  #### calculate things based on downsampled counts but also based on non-downsampled counts
  spectradf$count_div_by_targets <- spectradf$total_mutations/spectradf$total_target_count # note should be divided by this species NOT by human target count ; not downsampled, regularized or rescaled

  
  # get 1mer counts/rates
  SUMMARYOF1mercounts <- spectradf %>%
    group_by(species,population,label,mutation_1mer) %>%
    summarise(total_mutations_1mer = sum(total_mutations),total_target_count_1mer=sum(total_target_count)) %>% # summing up these targets work because you've faceted by 1mer type 
    ungroup()


  SUMMARYOF1mercounts$count_div_by_targets_1mer <- SUMMARYOF1mercounts$total_mutations_1mer / SUMMARYOF1mercounts$total_target_count_1mer  #should be div by each species target count NOT by human target count!

  # now want to merge it back
  mergedf <- merge(spectradf,SUMMARYOF1mercounts,by=c("species","population","label","mutation_1mer"))


  # this takes ratio of a given kmer mutation type rate over its central 1mer rate
  mergedf$kmerEnrichmentOverCentralBP <- mergedf$count_div_by_targets / mergedf$count_div_by_targets_1mer

  ##### TO DO: Fisher's exact test p-values #########
  # null hypo: ratio of TTTTAAAA>T/A>T is 1 (etc.)
  # two-sided alternate hypothesis that it's either enriched or depleted
  
  print('carrying out fishers exact test')
  mergedf <- mergedf %>%
      rowwise() %>%
      mutate(p_value_fisher = fisher.test(as.table(rbind(c(round(total_mutations,digits=0),round(total_mutations_1mer,digits=0)), c(total_target_count, total_target_count_1mer))),alternative ="two.sided")$p.value)
# NOTE: fisher's text requires integers (will work otherwise but gives a lot of warnings)

  mergedf$spectrumName <- spectrumName
  
  write.table(mergedf,paste0(fisherDir,"FISHERS.EXACT.TEST.RESULTS.",spectrumName,".notrescaled_notdownsampled_notregularized.txt"),row.names=F,quote=F,sep="\t")
  
  return(mergedf)
}


#tt=fishersexacttest_nondownsampledData(head(all3merSpectra_plusTargets),plotdir,"test3merSpectrum")
                                         
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

################# CARRY OUT FISHERS EXACT TEST ####################

############# >>>> 3mers <<<<< ##############
# test: 
print("starting 3mers, including CpGs")
# plotEnrichmentsAndPValues(tt,speciesToIncludeInPlot,plotdir)
speciesToIncludeInPlot=c("humans_EUR","humans_EAS","humans_AMR","Mus_musculus_Mmd","Mus_spretus_Ms","brown_bear_ABC","polar_bear_PB","vaquita","fin_whale_GOC","humans_AFR","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","Pongo_abelii","Pongo_pygmaeus","wolves")

all3merSpectra_FishersExactTestResults <- fishersexacttest_nondownsampledData(all3merSpectra_plusTargets,fisherDir=fisherDir,spectrumName = "3merSpectrum",excludeCpGs="no")

# then plot results:
plotEnrichmentsAndPValues(all3merSpectra_FishersExactTestResults,fisherDir=fisherDir,speciesToIncludeInPlot = speciesToIncludeInPlot)

####### >>>> 3mers with CPG 3mers excluded  <<<<< #############
print("starting 3mers, excluding CpGs")

all3merSpectra_FishersExactTestResults_excludingCpG <- fishersexacttest_nondownsampledData(all3merSpectra_plusTargets,fisherDir=fisherDir,spectrumName = "3merSpectrum",excludeCpGs="yes")

# make plots
plotEnrichmentsAndPValues(all3merSpectra_FishersExactTestResults_excludingCpG,fisherDir=fisherDir,speciesToIncludeInPlot = speciesToIncludeInPlot)


######## >>>>>> exploring TCC pulse ORs <<<<<<<<<<<< ######
print("exploring TCC pulse related Odds ratios")
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

############## >>>> Fishers test on 5/7mers <<<<<<<<<<< ##############

###### 5mer
print('starting 5mers')
all5merSpectra_FishersExactTestResults <- fishersexacttest_nondownsampledData(all5merSpectra_plusTargets,fisherDir=fisherDir,spectrumName = "5merSpectrum",excludeCpGs="no")

# then plot results:
plotEnrichmentsAndPValues(all5merSpectra_FishersExactTestResults,fisherDir=fisherDir,speciesToIncludeInPlot = speciesToIncludeInPlot)

###### 7mer: (VERY SLOW)
print('starting 7mers')

all7merSpectra_FishersExactTestResults <- fishersexacttest_nondownsampledData(all7merSpectra_plusTargets,fisherDir=fisherDir,spectrumName = "7merSpectrum",excludeCpGs="no")

# then plot results:
plotEnrichmentsAndPValues(all7merSpectra_FishersExactTestResults,fisherDir=fisherDir,speciesToIncludeInPlot = speciesToIncludeInPlot)

