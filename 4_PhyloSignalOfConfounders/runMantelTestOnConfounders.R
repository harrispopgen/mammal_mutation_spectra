require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(compositions)  # for clr
require(broom) # for tidy
require(ggrepel)
require(ggfortify) # for pca autoplot
require(gridExtra)
require(sjPlot) # for grids of plots
require(ggpmisc) # for nicer plotting of linear models
require(ape)
require(ggpubr)
set.seed(42)
outdir="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/confoundersMantelTest/"
mantelTestPermutationCount = 99999
#speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/analyses/compare_spectrumDist_toPhyloHammingDist/speciesCodes.txt",header = T)
#head(speciesCodes) # for getting phylo distances 

########### species to include ######################
#speciesList=c('humans', 'Mus_musculus','Mus_spretus' ,'polar_bear','brown_bear', 'fin_whale', 'vaquita', 'Gorilla_gorilla', 'Pan_troglodytes', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus','wolves') # overarching species list 

speciesToInclude_in_raxmlTree=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Canis_lupus")  # raxml name fmt
# taking Ms out: "Mus_spretus_Ms",

######  get phylo distances from upham et al Raxml tree:  #########
# process raxml tree here instead of in previous 
raxmlTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")

# want to get rid of the family/order info at end of tip. just keep species : 
raxmlTree_renamedTips <- raxmlTree
raxmlTree_renamedTips$tip.label <- paste0(unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",1)),"_",unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",2)))

# restrict to just the species to include
raxmlTree_renamedTips_subset <- keep.tip(raxmlTree_renamedTips,speciesToInclude_in_raxmlTree)

raxml_cophenetic_dist <- cophenetic(raxmlTree_renamedTips_subset)
# save this as an object
#saveRDS(raxml_cophenetic_dist,file=paste0(plotdir,"raxml_cophenetic_dist.dist"))
## ^^ this is what you need in matrix format for the mantel test
# and then make a tidy list for plotting correlations:
melt_cophenetic_func <- function(cophenetic_dist){
  cophenetic_dist_melt <- melt(cophenetic_dist)
  colnames(cophenetic_dist_melt) <- c("Sp1","Sp2","cophenetic_distance")
  # get species names # ASSUMES ARE SEPARATED BY "_" eg mus_musculus
  cophenetic_dist_melt$Sp1_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp1),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp1),"_"),"[",2)))
  
  cophenetic_dist_melt$Sp2_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp2),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp2),"_"),"[",2)))
  
  # need to deal with reciprocal duplicates
  # put in alphabetical order so I can get rid of dups:
  cophenetic_dist_melt$comparisonLabel <- paste0(pmin(cophenetic_dist_melt$Sp1_species,cophenetic_dist_melt$Sp2_species),".",pmax(cophenetic_dist_melt$Sp1_species,cophenetic_dist_melt$Sp2_species))  # # need to use pmin and pmax to get row by row min max. otherwise gives min max of whoel vector
  # get rid of self to self comparisons:
  cophenetic_dist_melt_distinct <- cophenetic_dist_melt[cophenetic_dist_melt$Sp1!=cophenetic_dist_melt$Sp2,]
  # and get rid of reciprocal dups (sp A - B and sp B- A)
  cophenetic_dist_melt_distinct <- cophenetic_dist_melt_distinct %>%
    ungroup() %>% # adding ungroup as a safety measure on 20211201 not sure if ottally necessary though
    distinct(comparisonLabel,.keep_all=T)  # restrict to just distinct comparison labels. 
  dim(cophenetic_dist_melt_distinct)
  return(cophenetic_dist_melt_distinct)
  
  
}


phyloDistances = melt_cophenetic_func(raxml_cophenetic_dist)

################### read in confounders ############
#confounders <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/information/ALT_FinWhaleENP_confoundersToPlotPhyloSignal.GenomeN50.etc.updatedthetavals_updatedCoverage_alt.txt",sep="\t",header=T) # 20221129 updated values of theta and coverage to be based on the 5 inds that were picked at random per pop
confounders <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/information/confoundersToPlotPhyloSignal.GenomeN50.etc.updatedthetavals.updatedCoverage.txt",sep="\t",header=T)
# this is an alt file that subs in fin whale ENP information for GOC since I am using ENP in this analysis.
# 20230516 -- going back to using GOC since GOC was used in all analyses, not ENP.

head(confounders)

#confoundersIWantToTest=c("wattersons_theta","avg_coverage_per_individual") # just redoing coverage and theta which were recalc'd based on non projected data
confoundersIWantToTest=c("scaffold_N50","contig_N50","avg_coverage_per_individual","Rspan_d","AFR_d","wattersons_theta") # maybe I could add a 0/1 thing for whether they were processed/sequenced together somehow? How would that work with mantel test though?

## remember to use sqrt of phylo distance
# you can input your full confounder dataframe and then give it the particular confounder you want to use (note must be in my specific format though)
confounders_comboFunction_mantelTest_SQRTDISTANCE <- function(confounder_full_dataframe, nameOfConfounderYouWant,phylo_dist,speciesToInclude,mantelTestPermutationCount){
  # restrict just to species you want:
  confounders_restrictedSpecies <- confounder_full_dataframe[confounder_full_dataframe$full_name_to_match_RAXML %in% speciesToInclude,] # in case you want to exclude some species from the analyses. for now it should be all 13
  # then restrict to the particular nameOfConfounderYouWant you want : 
  confounders_wide = pivot_wider(select(confounders_restrictedSpecies,c("full_name_to_match_RAXML",nameOfConfounderYouWant)),values_from = all_of(nameOfConfounderYouWant),names_from = "full_name_to_match_RAXML" )
  # get pairwise euclidean distances which just reduce down to abs(x-y) when it's just two numbers (did manual checks)
  # want to keep the diag and upper tri of matrix for mantel test even though redundant
  confounder_pairwise_distances = dist(t(confounders_wide),diag=T,upper=T) # want full matrix for mantel test  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(confounder_pairwise_distances)),colnames(as.matrix(confounder_pairwise_distances))]
  #head(phylo_dist_reordered)
  # TAKE SQUARE ROOT of phylo distance
  phylo_dist_reordered_sqrt <- sqrt(phylo_dist_reordered)
  # run mantel test
  mtr = vegan::mantel(xdis=phylo_dist_reordered_sqrt,ydis=confounder_pairwise_distances,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel_SQRTPHYLODISTANCE",confounder=nameOfConfounderYouWant,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}

plotMantelResults <- function(confounder_full_dataframe, nameOfConfounderYouWant,phylo_dist_TIDIED_notmatrix,speciesToInclude,mantelTestResult_df){
  # restrict just to species you want:
  confounders_restrictedSpecies <- confounder_full_dataframe[confounder_full_dataframe$full_name_to_match_RAXML %in% speciesToInclude,] # in case you want to exclude some species from the analyses. for now it should be all 13
  # then restrict to the particular nameOfConfounderYouWant you want : 
  confounders_wide = pivot_wider(select(confounders_restrictedSpecies,c("full_name_to_match_RAXML",nameOfConfounderYouWant)),values_from = all_of(nameOfConfounderYouWant),names_from = "full_name_to_match_RAXML" )
  # get pairwise euclidean distances which just reduce down to abs(x-y) when it's just two numbers (did manual checks)
  # want to keep the diag and upper tri of matrix for mantel test even though redundant
  # for plotting don't want to keep upper tri or diag and want to tidy it 
  confounder_pairwise_distances_forPlotting = tidy(dist(t(confounders_wide),diag=F,upper=F)) # want full matrix for mantel test  
  
  # make alphabetical order comparison labels so I can merge with phylo dist df:
  # put in alphabetical order so I can get rid of dups:
  confounder_pairwise_distances_forPlotting$comparisonLabel <- paste0(pmin(as.character(confounder_pairwise_distances_forPlotting$item1),as.character(confounder_pairwise_distances_forPlotting$item2)),".",pmax(as.character(confounder_pairwise_distances_forPlotting$item1),as.character(confounder_pairwise_distances_forPlotting$item2)))  # # need to use pmin and pmax to get row by row min max. otherwise gives min max of
# want to do this for every confounder and save plots
  # want to merge this with tidied phylo dist not the matrix  
  merge_confounderDist_plusPhyloDist <- merge(confounder_pairwise_distances_forPlotting,phylo_dist_TIDIED_notmatrix,by="comparisonLabel",all = T)
  # check that dimensions didn't change 
  if(dim(merge_confounderDist_plusPhyloDist)[1]!=dim(confounder_pairwise_distances_forPlotting)[1])
  {stop("something went wrong with merge")}
  
  # plot SQUARE ROOT(phylo dist) and confounder dist:
  mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_confounder <- ggplot(merge_confounderDist_plusPhyloDist,aes(x=sqrt(cophenetic_distance),y=distance))+
    geom_point(size=1,color="dodgerblue")+
    geom_text(data=mantelTestResult_df,aes(x=0.3,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=4)+
    ggtitle(paste0(nameOfConfounderYouWant,"\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantelTestResult_df$permCount)," permutations-- SQRT of phylo distance"))+
    theme_bw()+
    theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
    xlab("square root of phylogenetic distance (shared branch length)\nbranch lengths represent expected substitutions per site")+
    ylab(paste0(nameOfConfounderYouWant," abs. diff."))
  # making strips not gray
  
  return(mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_confounder)
  
}

plotMantelResults_AddingLabels <- function(confounder_full_dataframe, nameOfConfounderYouWant,phylo_dist_TIDIED_notmatrix,speciesToInclude,mantelTestResult_df){
  # restrict just to species you want:
  confounders_restrictedSpecies <- confounder_full_dataframe[confounder_full_dataframe$full_name_to_match_RAXML %in% speciesToInclude,] # in case you want to exclude some species from the analyses. for now it should be all 13
  # then restrict to the particular nameOfConfounderYouWant you want : 
  confounders_wide = pivot_wider(select(confounders_restrictedSpecies,c("full_name_to_match_RAXML",nameOfConfounderYouWant)),values_from = all_of(nameOfConfounderYouWant),names_from = "full_name_to_match_RAXML" )
  # get pairwise euclidean distances which just reduce down to abs(x-y) when it's just two numbers (did manual checks)
  # want to keep the diag and upper tri of matrix for mantel test even though redundant
  # for plotting don't want to keep upper tri or diag and want to tidy it 
  confounder_pairwise_distances_forPlotting = tidy(dist(t(confounders_wide),diag=F,upper=F)) # want full matrix for mantel test  
  
  # make alphabetical order comparison labels so I can merge with phylo dist df:
  # put in alphabetical order so I can get rid of dups:
  confounder_pairwise_distances_forPlotting$comparisonLabel <- paste0(pmin(as.character(confounder_pairwise_distances_forPlotting$item1),as.character(confounder_pairwise_distances_forPlotting$item2)),".",pmax(as.character(confounder_pairwise_distances_forPlotting$item1),as.character(confounder_pairwise_distances_forPlotting$item2)))  # # need to use pmin and pmax to get row by row min max. otherwise gives min max of
  # want to do this for every confounder and save plots
  # want to merge this with tidied phylo dist not the matrix  
  merge_confounderDist_plusPhyloDist <- merge(confounder_pairwise_distances_forPlotting,phylo_dist_TIDIED_notmatrix,by="comparisonLabel",all = T)
  # check that dimensions didn't change 
  if(dim(merge_confounderDist_plusPhyloDist)[1]!=dim(confounder_pairwise_distances_forPlotting)[1])
  {stop("something went wrong with merge")}
  
  # plot SQUARE ROOT(phylo dist) and confounder dist:
  mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_confounder_label <- ggplot(merge_confounderDist_plusPhyloDist,aes(x=sqrt(cophenetic_distance),y=distance))+
    geom_point(size=1,color="dodgerblue")+
    geom_text_repel(aes(label=comparisonLabel),size=1.3,max.overlaps=100)+
    geom_text(data=mantelTestResult_df,aes(x=0.3,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=4)+
    ggtitle(paste0(nameOfConfounderYouWant,"\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantelTestResult_df$permCount)," permutations-- SQRT of phylo distance"))+
    theme_bw()+
    theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
    xlab("square root of phylogenetic distance (shared branch length)\nbranch lengths represent expected substitutions per site")+
    ylab(paste0(nameOfConfounderYouWant," abs. diff."))
  # making strips not gray
  
  return(mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_confounder_label)
  
}

###### do for every confounder:
for(confounder in confoundersIWantToTest){
  # first run the mantel test
  # note here you use the matrix of raxml distances
  mantelResultsPerConfounder <- 
    confounders_comboFunction_mantelTest_SQRTDISTANCE(confounders,confounder,raxml_cophenetic_dist,speciesToInclude_in_raxmlTree,mantelTestPermutationCount)
 
   # write out table
  write.table(mantelResultsPerConfounder,paste0(outdir,confounder,".mantelTestResults.txt"),row.names = F,quote=F,sep="\t")
  
  # then plot using results of mantel test and tidied phylo distances
  mantelPlotPerConfounder <- plotMantelResults(confounders,confounder,phyloDistances,speciesToInclude_in_raxmlTree,mantelResultsPerConfounder)

  # then save the plot
  ggsave(paste0(outdir,confounder,".mantelTestResults.CorrelationPlot.pdf"),mantelPlotPerConfounder)
  
  # also make a version with labels:
  mantelPlotPerConfounder_label <- plotMantelResults_AddingLabels(confounders,confounder,phyloDistances,speciesToInclude_in_raxmlTree,mantelResultsPerConfounder)
  
  # then save the plot
  ggsave(paste0(outdir,confounder,".mantelTestResults.CorrelationPlot.labels.pdf"),mantelPlotPerConfounder_label)
  
}



