############# gather all spectra and calculate rescaled mutation rates  ###########
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
require(spgs) # for rev comp. be VERY CAREFUL not to use its rev comp on a vector of 1mers
require(viridis)
require(scales)
require(lsa) # for cosine similarity function 

set.seed(42) # adding seed
########### note: dplyr is a bit dangerous when dividing by sum(X) -- if the df has ever been grouped by things in the past, it maintains that grouping. This can manifest as sneaky bugs for fraction of segregating sites if you don't ungroup() after grouping.

epsilon=1 # amount to add to all mutation types to deal with the 0-entries
# a key paper to cite for dealing with regularization of compositional data: https://www.sciencedirect.com/science/article/pii/S0169743921000162#bib8
separateCpGs="yes" # yes or no
todaysdate=format(Sys.Date(),"%Y%m%d")
maskLabel="maskALL" # if species get different mask labels this will be a pain.
mantelTestPermutationCount=9999999 # trying with more (10M now)
# will be quite slow
flagOfAnalysis=paste0(todaysdate,".plots.SubsampledIndividualsNotProjection.epsilon.",epsilon,".sepCpGs.",separateCpGs,".masked",maskLabel,".RESCALEDTOHUMANTARGETS.SimplifiedCode.",as.character(mantelTestPermutationCount),".MantelPermutations.1merCpGOnly") # a label that will be in title of plots and outdir 
plotdir=paste0("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/compare_spectrumDist_toPhyloHammingDist/",flagOfAnalysis,"/")
dir.create(plotdir,showWarnings = F)

# file with nice labels:
source("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/colors_and_labels.R") 

######### CHOOSE VARIABLE FOR MANTEL TEST #############
variableToUseInMantelTest="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" # for now but can change

trblshootdir=paste0(plotdir,"troubleshooting/")
dir.create(trblshootdir,showWarnings = F)

############## all spectra based on 5 subset individuals #############
kmersize="7"
wd=paste0("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/allspecies_mutyper_results_unified_SubsettingIndividuals/",kmersize,"_mer/")

######### TO DO: read in spectra here

allSpectra <- read.table(paste0(wd,"/concatenated_pop_level_spectra_allSpecies_basedOnSubsetOfIndividuals/masked_maskALL/including_ucla_wolves_comparison/allSpectra.PopulationLevel.BasedOn5individualsPerPop.notprojected.usethis.txt"),sep="\t",header=T) # note this contains UCLA wolves but am not using for generating main results; can change down in species to include list
###### this is NOT projected, but is instead pop-level spectra based on drawing 5 individuals at random per population. population spectra are randomized within pops (but not between pops) and fixed derived sites within a population are removed.
head(allSpectra)
dim(allSpectra)
############### species codes (for going form latin name to common ) ##################
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/analyses/compare_spectrumDist_toPhyloHammingDist/speciesCodes.txt",header = T)
head(speciesCodes) # for getting phylo distances 

########### species to include ######################
speciesList=c('humans', 'Mus_musculus','Mus_spretus' ,'polar_bear','brown_bear', 'fin_whale', 'vaquita', 'Gorilla_gorilla', 'Pan_troglodytes', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus','wolves') # overarching species list 


speciesToInclude_in_raxmlTree=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Canis_lupus")  # raxml name fmt
# taking Ms out: "Mus_spretus_Ms",

speciesToInclude_phylo=c("Mus_musculus_Mmd","Mus_spretus_Ms","brown_bear_ABC","polar_bear_PB","vaquita","fin_whale_GOC","humans_AFR","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","Pongo_abelii","Pongo_pygmaeus","wolves") # human_AFR# switching to GOC fin whale for now due to ksfs weirdness
# note am excluding ucla wolves here



############ read in confounders dataframe ##############
confounders <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/information/confoundersToPlotPhyloSignal.GenomeN50.etc.updatedthetavals.updatedCoverage.txt",sep="\t",header=T) # updated values of theta and coverage to be based on the 5 inds that were picked at random per pop

head(confounders)

confoundersIWantToTest=c("scaffold_N50","contig_N50","avg_coverage_per_individual","Rspan_d","AFR_d","wattersons_theta") # maybe I could add a 0/1 thing for whether they were processed/sequenced together somehow? How would that work with mantel test though?

######  get phylo distances from upham et al Raxml tree:  #########
# used to get these distances in a separate script. now doing it here. 

# process raxml tree here instead of in previous 
raxmlTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")

# want to get rid of the family/order info at end of tip. just keep species : 
raxmlTree_renamedTips <- raxmlTree
raxmlTree_renamedTips$tip.label <- paste0(unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",1)),"_",unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",2)))

# restrict to just the species to include
raxmlTree_renamedTips_subset <- keep.tip(raxmlTree_renamedTips,speciesToInclude_in_raxmlTree)

raxml_cophenetic_dist <- cophenetic(raxmlTree_renamedTips_subset)
# save this as an object
saveRDS(raxml_cophenetic_dist,file=paste0(plotdir,"raxml_cophenetic_dist.dist"))

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


######### read in time-calibrated tree cophenetic distances -- this is from UpHam et al. may want to use TimeTree Instead in future!  ########3

timeTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/timeTreePhyloTree/listOfSpecies.ForTimeTree.pluswolves.20220204.nwk") # updated to include wolf 
# get cophenetic distances:
timeTree_cophenetic_dist <- cophenetic(timeTree)
timeDistances = melt_cophenetic_func(timeTree_cophenetic_dist)

# now do same processing:

##### 1mers are dealt with in the if statement below this (dealing with CpGs) 
# get ancestral kmers:
allSpectra$ancestral7mer <- substr(allSpectra$mutation_type,1,7)
allSpectra$ancestral5mer <- substr(allSpectra$mutation_type,2,6)
allSpectra$ancestral3mer <- substr(allSpectra$mutation_type,3,5)

# get central mutation types:
allSpectra$mutation_7mer <- allSpectra$mutation_type
allSpectra$mutation_5mer <- paste0(substr(allSpectra$mutation_type,2,6),".",substr(allSpectra$mutation_type,10,14))
allSpectra$mutation_3mer <- paste0(substr(allSpectra$mutation_type,3,5),".",substr(allSpectra$mutation_type,11,13))

# want to get ancestral 1mers for all data frames so I can plot things with central mutattion type separated later:

#### SEPARATE OUT CPGS  -- OPTIONAL ############
####### prior to 20221208 this would separate all kinds of CpG>N muations (CpG.ApG etc) but to match with sigfit we want to just separate CpG.TpG out and leave the rest as part of the overall C>A C>G 1mer categories.
# but need to change target sizes below to accomodate this
if(separateCpGs=="yes"){
  #### label CpGs and specify them as ancestral target as well
  print("separating out CpG sites!")
  # make a noncpg label first:
  allSpectra$mutation_1mer <- paste0(substr(allSpectra$mutation_type,4,4),".",substr(allSpectra$mutation_type,12,12)) # if you don't want to separate out CpGs then this is just the regular center bp
  
  # label CpG.TpG mutations (ignore other kinds of CpG>NpG mutations ; they'll stay in their respective categories)
  allSpectra[grepl("CG",allSpectra$ancestral3mer) & allSpectra$mutation_1mer=="C.T",]$mutation_1mer <- "CpG.TpG"
  
  # changing how I label ancestral 1mers for different mutation types. since we aren't separating out non C>T CpG mutations I want them to have the full "C" target size (CnonCpG+C_cpg) whereas noncpg C>T mtuations should just have nonCpG_Cs as their targets and cpg C>Ts should have CpGs as target size (NOTE this is different from how I've been doing it! -- redo the CpG-separated distance script to incorporate this as well)
  # get regular ancestral 1mer info and then update it:
  
  allSpectra$ancestral1mer <- substr(allSpectra$mutation_type,4,4) # if you don't want to separate out CpGs then this is just the regular center bp
  
  allSpectra[allSpectra$mutation_1mer=="C.T",]$ancestral1mer <- "C_nonCpG" # non-cpg C>Ts
  allSpectra[allSpectra$mutation_1mer=="CpG.TpG",]$ancestral1mer <- "C_CpG"
  # note that other C>N mutations will have "C" as ancestral 1mer because all kinds of Cs are valid (non cpg + cpg)
  
  
  
} else if(separateCpGs=="no"){
  print("note that CpGs are not separated out!")
  allSpectra$mutation_1mer <- paste0(substr(allSpectra$mutation_type,4,4),".",substr(allSpectra$mutation_type,12,12)) # if you don't want to separate out CpGs then this is just the regular center bp
  # still get ancestral1mer (but don't specify if CpG or not)
  allSpectra$ancestral1mer <- substr(allSpectra$mutation_type,4,4) # if you don't want to separate out CpGs then this is just the regular center bp
  
  
  
} else {
  print("invalid choice for separateCpGs")
  break
}



# get spectrum summed up per type
# all7merSpectraOnly <- allSpectra %>%
#   # adding groupbing by central 1mer but that doesn't actually change the sums since 7mers already group by those 
#   # I just want to keep it in the dataframe 
#   group_by(species,population,label,ancestral7mer,mutation_7mer,ancestral1mer,mutation_1mer) %>%
#   summarise(total_mutations = sum(totalSites)) %>%
#   ungroup() # this will be identical to original spectrum (7mer), just doing for consistency of formatting
# 
# all5merSpectraOnly <- allSpectra %>%
#   group_by(species,population,label,ancestral5mer,mutation_5mer,ancestral1mer,mutation_1mer) %>%
#   summarise(total_mutations = sum(totalSites)) %>%
#   ungroup()

# all3merSpectraOnly <- allSpectra %>%
#   group_by(species,population,label,ancestral3mer,mutation_3mer,ancestral1mer,mutation_1mer) %>%
#   summarise(total_mutations = sum(totalSites)) %>%
#   ungroup()

all1merSpectraOnly <- allSpectra %>%
  group_by(species,population,label,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>%
  ungroup()
# if you don't want to partition over CpGs you can use mutation_1mer_CpGNotLabeled instead
# okay NOTE HERE: "C" counts are nonCpG and then C_CpG are CpG ancestral targets


########## read in targets ##############
#all7merTargets=data.frame()
#all5merTargets=data.frame()
#all3merTargets = data.frame()
all1merTargets = data.frame()

for(species in speciesList){
  # some 'samples' may contain multipops like mice and bears
  targetdir=paste0(wd,"allspecies_summed_up_over_intervals_forTransfer/",species,"/mutyper_results_masked_",maskLabel,"/mutyper_target_files/")
  targetsfilename = paste0(targetdir,species,".summedup.mutyper.targets.SeeLogForFilters.",maskLabel,".",kmersize,"mer.txt")
  
  targets=read.table(targetsfilename,header=T) 
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  targets$species <- species
  # get central 5mer, 3mer and 1mer:
  targets <- targets %>%
    mutate(target_5mer=substr(target_7mer,2,6),target_3mer=substr(target_7mer,3,5),target_1mer= substr(target_7mer,4,4)) # get substrings 
  
  ######### relabel target_1mer with CpGs if you are separating out CpGs ###########
  # 20220112 -- note that grepping central CpG is fine because there are no CGX ancestral 3mers because they would be revcomped to YCG so you're good here!  
  # but we want C>A and C>G mutatoins to have non CpG and CpG sites as their targets since we don't separate those into CpG or not anymore
  if(separateCpGs=="yes"){
    targets$CpGLabel <- ""
    targets[grepl("CG",targets$target_3mer),]$CpGLabel <- "_CpG"
    targets[targets$target_1mer=="C" & !grepl("CG",targets$target_3mer),]$CpGLabel <- "_nonCpG" # label the nonCpG C's as well.
    # so now nothing is labeled as plain "C"
    targets$target_1mer_CpGsNotLabeled <- targets$target_1mer # update with this 
    targets$target_1mer <- paste0(targets$target_1mer_CpGsNotLabeled,targets$CpGLabel)
    
  } # don't need an else statement, just keep going with targets target$target_1mer as is.
  
  
  
  # TARGETS per 3mer type #
  #targets_7mer <- targets %>% group_by(species,target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup() # original targets file, just doing to be consistent but it's identical to original targets
  #targets_3mer <- targets %>% group_by(species,target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup()
  #targets_5mer <- targets %>% group_by(species,target_5mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup()
  targets_1mer <- targets %>% group_by(species,target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup()

  # make all species targets:
  
  
  # combine  all targets:
  #all7merTargets = bind_rows(all7merTargets,targets_7mer)
  #all5merTargets = bind_rows(all5merTargets,targets_5mer)
  #all3merTargets = bind_rows(all3merTargets,targets_3mer)
  all1merTargets = bind_rows(all1merTargets,targets_1mer)
  
  
}

############ get target proportions (new as of 20220727) ##############
all1merTargets <- all1merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()


# all3merTargets <- all3merTargets %>%
#   group_by(species) %>%
#   mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
#   ungroup()
# 
# all5merTargets <- all5merTargets %>%
#   group_by(species) %>%
#   mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
#   ungroup()
# 
# 
# all7merTargets <- all7merTargets %>%
#   group_by(species) %>%
#   mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
#   ungroup()
# plot targets:

if(separateCpGs=="yes"){
  print("adding in an extra row to targets that is total 'C's : the sum of CpGs + nonCpGs - use this for target for C>A and C>G mutations\nnote this was done after calculating proportions so that it wasn't included when calcuating proportions; non CpG and CpG proportions were added together to get total C proportion ")
  
  all1merTargets_CTotals_summingCpGandNonCpG <- all1merTargets %>%
    group_by(species) %>%
    filter(target_1mer %in% c("C_CpG","C_nonCpG")) %>%
    summarise(total_target_count=sum(total_target_count),target_proportion=sum(target_proportion)) %>%
    mutate(target_1mer="C") %>%
    ungroup()
  
  all1merTargets_plusCTotals <- bind_rows(all1merTargets_CTotals_summingCpGandNonCpG,all1merTargets)
  
  all1merTargets_plusCTotals <- arrange(all1merTargets_plusCTotals,species) # sorts by species so that when you look at it, you see all the types together.
  
  head(all1merTargets_plusCTotals)
  
} else if(separateCpGs=="no") {
  all1merTargets_plusCTotals=all1merTargets # dont need to change anything
  head(all1merTargets_plusCTotals)
}


########### fill in 0-entry types #############
# all7merSpectraOnly_filledin <- all7merSpectraOnly %>%
#   ungroup() %>%
#   complete(nesting(mutation_7mer,ancestral7mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations=0))%>%
#   mutate(mutation_label=mutation_7mer)
# 
# if(dim(all7merSpectraOnly_filledin)[1]==dim(all7merSpectraOnly)[1]){
#   print("filling in didn't work!")
# }
# all5merSpectraOnly_filledin <- all5merSpectraOnly %>%
#   ungroup() %>% # need to pre-ungroup ; got grouped up above
#   complete(nesting(mutation_5mer,ancestral5mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations=0))%>%
#   mutate(mutation_label=mutation_5mer) # adding mutation label for use in functions
# 
# all3merSpectraOnly_filledin <- all3merSpectraOnly %>%
#   ungroup() %>%
#   complete(nesting(mutation_3mer,ancestral3mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations=0))%>%
#   mutate(mutation_label=mutation_3mer)

all1merSpectraOnly_filledin <- all1merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_1mer,ancestral1mer),nesting(species,population,label),fill=list(total_mutations=0)) %>%
  mutate(mutation_label=mutation_1mer)
# fill in should only be needed for 7mer (won't change dim of others)

# so when I was merging before all species were combined there was a bug where missing mtuation types wouldn't get targets filled in which led to them being NA downstream 
# okay now fill in missing mutation types : as 0s PRIOR TO MERGING WITH TARGETS
#  note that targets are the same within species (mouse1 and mouse2 have same targets; bear 1 and bear2 have targets)


############## merge with targets ###############

# all7merSpectra_plusTargets <- merge(all7merSpectraOnly_filledin, all7merTargets,by.x=c("species","ancestral7mer"),by.y=c("species","target_7mer"))
# 
# if(dim(all7merSpectra_plusTargets)[1] != dim(all7merSpectraOnly_filledin)[1]){
#   print("something went wrong with merge!")
#   
# }

# write out for enrichments
#write.table(all7merSpectra_plusTargets,paste0(plotdir,"all7merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")


#all5merSpectra_plusTargets <- merge(all5merSpectraOnly_filledin, all5merTargets,by.x=c("species","ancestral5mer"),by.y=c("species","target_5mer"))

#write.table(all5merSpectra_plusTargets,paste0(plotdir,"all5merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")

#all3merSpectra_plusTargets <- merge(all3merSpectraOnly_filledin, all3merTargets,by.x=c("species","ancestral3mer"),by.y=c("species","target_3mer"))

#write.table(all3merSpectra_plusTargets,paste0(plotdir,"all3merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")


all1merSpectra_plusTargets <- merge(all1merSpectraOnly_filledin, all1merTargets_plusCTotals,by.x=c("species","ancestral1mer"),by.y=c("species","target_1mer"))
# don't need to write out the 1mers because they aren't used as numerator in enrichments 
# you are ehre
# NOW CpGs are the target size for CpG.TpG mtuations, but C+CpG are the target size for C>A and C>G
# you can see:
# CHECK:
# all1merSpectra_plusTargets[all1merSpectra_plusTargets$label=="brown_bear_ABC" & all1merSpectra_plusTargets$mutation_1mer %in% c("C.G","C.A","C.T","CpG.TpG"),] # 15031288 CpGs sites (target size for CpG.TpG + 254275125 non CpGs sites (target size for C>T)  = 269306413 which is the target size for C>A and C>G (works!)

############### >>>> RESCALE TO HUMAN TARGETS <<<<< #################
# rescaling all counts to be scaled by human targets total_mutations_of_typeX_in_SpeciesA* (target_proportion_HUMAN /target_proportion_speciesA)

humanTargetProportions_1mer <- all1merTargets_plusCTotals[all1merTargets_plusCTotals$species=="humans",]
# note that these proportions add up to more than one because both C and CPG and nonCPG C are included
# but that works out okay because you are going to merge based on the ancestral 1mer
# so C>A mtuations will be multiplied up by the frac of C in human / frac of C in species A where C frac is CpG+nonCPG Cs, whereas CpG.TPG muts will be scaled by cpg frac in human/cpg frac in species A. so that all is okay.
#humanTargetProportions_3mer <- all3merTargets[all3merTargets$species=="humans",]
#humanTargetProportions_5mer <- all5merTargets[all5merTargets$species=="humans",]
#humanTargetProportions_7mer <- all7merTargets[all7merTargets$species=="humans",]

# note these don't yet have the human counts; need to run processSpectra_AdjustForHumanTargets on them (is done in combo distance script at bottom fo the script )
# going to keep these names 
all1merSpectra_plusHumanInfoNotYetRescaled <- merge(all1merSpectra_plusTargets,humanTargetProportions_1mer[,c("target_1mer","total_target_count","target_proportion")],by.x=c("ancestral1mer"),by.y=c("target_1mer"),suffixes=c("_thisSpecies","_HUMAN"))

# all3merSpectra_plusHumanInfoNotYetRescaled <- merge(all3merSpectra_plusTargets,humanTargetProportions_3mer[,c("target_3mer","total_target_count","target_proportion")],by.x=c("ancestral3mer"),by.y=c("target_3mer"),suffixes=c("_thisSpecies","_HUMAN"))
# 
# all5merSpectra_plusHumanInfoNotYetRescaled <- merge(all5merSpectra_plusTargets,humanTargetProportions_5mer[,c("target_5mer","total_target_count","target_proportion")],by.x=c("ancestral5mer"),by.y=c("target_5mer"),suffixes=c("_thisSpecies","_HUMAN"))
# 
# all7merSpectra_plusHumanInfoNotYetRescaled <- merge(all7merSpectra_plusTargets,humanTargetProportions_7mer[,c("target_7mer","total_target_count","target_proportion")],by.x=c("ancestral7mer"),by.y=c("target_7mer"),suffixes=c("_thisSpecies","_HUMAN"))

# okay you want to get spectrum and the mutation rate and want to add 1 to all counts to deal with 0 entries (assuming every mutation type would happen if you sampled enough -- though note this will have a bigger i)
print("rescaling to human targets")

# NOTE! need to be dividing by HUMAN target sizes for mutation rates!!! 
#### scale counts to human targets: 

all1merSpectra_HumanRescaled <- all1merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

# all3merSpectra_HumanRescaled <- all3merSpectra_plusHumanInfoNotYetRescaled %>%
#   group_by(species,population,label) %>%
#   mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
#   ungroup() %>%
#   rename(total_mutations_ORIGINAL=total_mutations)
# 
# all5merSpectra_HumanRescaled <- all5merSpectra_plusHumanInfoNotYetRescaled %>%
#   group_by(species,population,label) %>%
#   mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
#   ungroup() %>%
#   rename(total_mutations_ORIGINAL=total_mutations)
# 
# all7merSpectra_HumanRescaled <- all7merSpectra_plusHumanInfoNotYetRescaled %>%
#   group_by(species,population,label) %>%
#   mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
#   ungroup() %>%
#   rename(total_mutations_ORIGINAL=total_mutations)

########### >>>> DOWNSAMPLING <<<<< #######################
print("downsampling!")

# do separately for every spectrum
# need to add plotdir to write files to 
downsampleFunc_RescaledByHumanTargets <- function(spectradf,plotdir){
  
  spectrumName <- deparse(substitute(spectradf)) # gets the name of the object that was put in
###### trying experiment 
  totalSegSites <- spectradf %>%
   group_by(species,label) %>%
    summarise(totalSegSites=sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    ungroup()
  
  subsampleValue=min(totalSegSites$totalSegSites)
  
  subsampleSpecies=totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species

  print(paste0("downsampling to ",subsampleValue," sites (",subsampleSpecies,") [ counts already rescaled by human targets ]"))
  
  
  multinomlogfile <- file(paste0(plotdir,"multinomialLogFile.",spectrumName,".Downsampling.RESCALED_BY_HUMAN_TARGETS.txt"))
  
  writeLines(paste0("downsampling to ",min(totalSegSites$totalSegSites)," sites (",totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species," [after rescaling by human targets ])"),multinomlogfile)
  
  close(multinomlogfile)
  
  # first need to get frac segregating sites:
  spectradf <- spectradf %>% 
    group_by(species,population,label) %>%
    mutate(fractionOfSegregatingSites_RESCALED_BY_HUMAN_TARGETS = total_mutations_RESCALED_BY_HUMAN_TARGETS / sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    mutate(total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS=as.numeric(rmultinom(n=1,size=subsampleValue,prob=fractionOfSegregatingSites_RESCALED_BY_HUMAN_TARGETS))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
    ungroup() %>% # ungroup just in case here 
    
    return(spectradf)
  }

##### ok!! These have now been downsampled !!!!!!!!! ###########
all1merSpectra <- downsampleFunc_RescaledByHumanTargets(all1merSpectra_HumanRescaled,plotdir)
# all3merSpectra <- downsampleFunc_RescaledByHumanTargets(all3merSpectra_HumanRescaled,plotdir)
# all5merSpectra <- downsampleFunc_RescaledByHumanTargets(all5merSpectra_HumanRescaled,plotdir)
# all7merSpectra <- downsampleFunc_RescaledByHumanTargets(all7merSpectra_HumanRescaled,plotdir)


# write these out:
write.table(all1merSpectra,paste0(plotdir,"all1merSpectra.plusCpGs.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
# write.table(all3merSpectra,paste0(plotdir,"all3merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
# write.table(all5merSpectra,paste0(plotdir,"all5merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
# write.table(all7merSpectra,paste0(plotdir,"all7merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)


####### define some functions ###########

######## PROCESS SPECTRA NO LONGER DOES DOWNSAMPLING (!!)
# HAPPENS BEFORE HAND 
######### 20220829: no longer want to normalize mutation rates (!!) - since they're human targets and downsampled don't need to normalize them. and that might introduce biases. 
processSpectra <- function(rescaled_and_downsampled_spectradf,epsilon){ 
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +epsilon (1) to every count to deal with 0 entries ; taking out adding epsilon to target counts since I fixed the empty target bug
  # no longer is doing downsampling (doing it beforehand)
  # and 
 # don't need to get frac seg sites here, doing down below # but am renaming so that it doesn't get accidentally used 
  # not keeping original values because going to downsample based on human counts
  # don't take away original columns bc want to be able to compare original stuff
  # maybe rename them though
  # now am dividing by human target counts to get rates (same for all species)
  # figure out kmer size to act as a label for writing out:

  
  spectradf <- rescaled_and_downsampled_spectradf %>%
    mutate(total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=(total_mutations_RESCALED_BY_HUMAN_TARGETS+epsilon),
           mutation_rate_RESCALED_BY_HUMAN_TARGETS=(total_mutations_RESCALED_BY_HUMAN_TARGETS/total_target_count_HUMAN),
           mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=(total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS/total_target_count_HUMAN))  # %>%
    
    spectradf <- spectradf %>%
    # now am adding epsilon and doing all transforms to the downsampled counts 
    mutate(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS = (total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS+epsilon),
           mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS = (total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS/total_target_count_HUMAN),
           mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS/total_target_count_HUMAN))
    
    # getting rid of normalization
   # 20220420 fixing frac seg sites plus epsilon bug -- need to use plus epsilon version 
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +epsilon (1) to every count to deal with 0 entries ; taking out adding epsilon to target counts since I fixed the empty target bug
  # do want to also calculate frac seg sites plus epsilon for both the downsampled and non downsampled verisons: 
    spectradf <- spectradf %>%
      group_by(species,population,label) %>%
      mutate(fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=(total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS/sum(total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS)),fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS/sum(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS))) %>% 
      ungroup()
    # now that things are rescaled to human targets, I am thinking that frac seg sites may not have the same issues it used to have of baking in genome composition. 
  return(spectradf)

}
# test function:
#test = processSpectra(all7merSpectra,epsilon) # ,plotdir,writeOut=F
#test3mers = processSpectra(all3merSpectra,epsilon) # ,plotdir,writeOut=F
test1mers = processSpectra(all1merSpectra,epsilon) # ,plotdir,writeOut=F

#sum(test[test$label=="vaquita",]$total_mutations_ORIGINAL)
#sum(test[test$label=="vaquita",]$total_mutations_RESCALED_BY_HUMAN_TARGETS)
# vaquita went to 139377 instead of 129668.9 Went down when resclaed by human target sizes. 


# need to pivot it wider? 
pivotSpectra_perVariable <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("label"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}

# this is just for mantel test:
pivotSpectra_perVariable_KEEPFULLSPECIESNAMES <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("species.sppCodes"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}

# test function
testvar="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" # no longer using the normalized version
testpivoted = pivotSpectra_perVariable(test1mers,testvar)
head(testpivoted)
######## CHANGE THIS IF CHANGE TEST KMER: 


  # want to speed up with only a subset of species.
# ilr/clr calcs require no 0 entries
# can also specify no transform 
clr_or_ilr_orNoTransform_calculations <- function(pivotedspectradf_plusEpsilonVersionOfVariables,ilr_or_clr_or_none){
  # remove character vectors:
  tableWithoutIDVars <- select(pivotedspectradf_plusEpsilonVersionOfVariables, -c("variable","mutation_label"))
  # AHA so if you run clr just on a matrix it goes ROW-WISE not column wise! annoying! so need to transpose
  #clr_matrix <- clr(t(tableWithoutIDVars)) # if you do it this way it matches when I do it just on on un-tranposed column by itself. if you do it without transposing then you get weird row-wise results that are WRONG. so be cautious here. 
  # so alternative without transposing is to use sapply -- this is a lot nicer I think
  if(ilr_or_clr_or_none=="clr"){
    print("carrying out clr transformation")
    df <- data.frame(sapply(tableWithoutIDVars,clr)) # this way you don't have to transpose and you get a much nicer shaped output
    return(df)
  } else if(ilr_or_clr_or_none=="ilr"){
    print("carrying out ilr transformation")
    
    df <- data.frame(sapply(tableWithoutIDVars,ilr)) # this way you don't have to transpose and you get a much nicer shaped output
    return(df)
    
  } else if(ilr_or_clr_or_none=="none") {
    print("carrying out no transformation")
    df <- tableWithoutIDVars # just return the table without id vars and no transformation
    return(df)
  } else {
    print("invalid option")
    break
  }
  
}


euclideanDistance <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- tidy(dist(t(pivotedspectradf_noIDVars))) # transpose for dist otherwise it goes along rows and is wrong
  colnames(distances) <- c("item1","item2","spectrum_distance")
  return(distances)
}

# match with species codes for raxml tree distances
addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf <- function(distancesdf,speciesCodesdf){
  merge1 <- merge(distancesdf,speciesCodesdf,by.x="item1",by.y="code",all.x=TRUE,all.y=FALSE) # only want comparisons which have phylo distances (?)
  merge2 <- merge(merge1,speciesCodesdf,by.x="item2",by.y="code",all.x=TRUE,all.y=FALSE,suffixes=c(".item1",".item2"))
  print("if there are NAs in species name, it's because I don't have it in raxml tree")
  return(merge2)
}

# add in phylo distances
addPhyloDistances_excludeNas <- function(dfWithSpeciesCodesAdded,phylodistancesdf) {
  # remove NAs (no phylo distance present)
  dfWithSpeciesCodesAdded <- na.omit(dfWithSpeciesCodesAdded)
  # add alphabetical comparison label:
  dfWithSpeciesCodesAdded$comparisonLabel <- paste0(pmin(dfWithSpeciesCodesAdded$species.item1,dfWithSpeciesCodesAdded$species.item2),".",pmax(dfWithSpeciesCodesAdded$species.item1,dfWithSpeciesCodesAdded$species.item2))
  # add a common name comparison label (note order may be diff from latin name because alphabetical)
  dfWithSpeciesCodesAdded$comparisonLabel_common_alphabetical <-  paste0(pmin(dfWithSpeciesCodesAdded$common_name.item1,dfWithSpeciesCodesAdded$common_name.item2),".",pmax(dfWithSpeciesCodesAdded$common_name.item1,dfWithSpeciesCodesAdded$common_name.item2))
  
  dfWithSpeciesCodesAdded$comparisonLabel_broad_alphabetical <-  paste0(pmin(dfWithSpeciesCodesAdded$broad_classification.item1,dfWithSpeciesCodesAdded$broad_classification.item2),".",pmax(dfWithSpeciesCodesAdded$broad_classification.item1,dfWithSpeciesCodesAdded$broad_classification.item2))
  
  merge1 <- merge(dfWithSpeciesCodesAdded,phylodistancesdf,by="comparisonLabel")
  return(merge1)
}


# keep it as distance mat
euclideanDistance_keepMatrix <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- dist(t(pivotedspectradf_noIDVars)) # transpose for dist otherwise it goes along rows and is wrong
  #colnames(distances) <- c("item1","item2","spectrum_distance")
  return(distances)
}


euclideanDistance_keepMatrix_UPPERANDLOWER <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- dist(t(pivotedspectradf_noIDVars),diag=T,upper=T) # transpose for dist otherwise it goes along rows and is wrong
  #colnames(distances) <- c("item1","item2","spectrum_distance")
  return(distances)
}


###### Get ILR and CLR distances ########
# this is slow so restrict just to species you want
##### this subsets just to the species you want above ^^ 
############# clr mutation rate ##############

####### combine functions ##########
# input: spectrum, variable (frac seg sites or mut rate), label (1mer , 5mer etc, and choice of ilr or clr)
combo_function_phylo <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,speciesCodes,phyloDistances,epsilon){
  distance_dataframe <- 
    processSpectra(spectrumdf,epsilon) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
    euclideanDistance(.) %>%
    addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
    addPhyloDistances_excludeNas(.,phyloDistances) %>%
    select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
    mutate(variable=variable,transform_label=ilr_or_clr_or_none)
  
  
  return(distance_dataframe)
  
  
}



############## run combo functions ######################
# do some loops
# all Distances
# 20220617 making labels nicer
listOfDFs = list(`1-mer+CpG spectrum`=all1merSpectra)
# lapply applies combo_function_phylo to each of the above spectra ^^ 
# if you set the names of a list, they get preserved in the output of lapply! cool!


listOfVars=c("mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS") # no longer using the normalized versions ; what about nondownsampled frac seg sites? adding in total mutations as well
listOfTransforms=c("clr") # clr only here=

########### run combo function on ALL data frames with all conditions; do for phylo and hamming ############
all_distances_all_conditions_phylo = data.frame()
all_distances_all_conditions_time = data.frame()
for(variable in listOfVars){
  for(transform in listOfTransforms){
    # get for all dfs: 
    ########## phylo distances ###############
    list_of_distance_dfs_phylo <- lapply(listOfDFs,combo_function_phylo, variable =variable, ilr_or_clr_or_none = transform,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances,epsilon=epsilon)
    # this will return a list of named results. How to deal with that? 
    # can bind them together and add some ids:
    distance_dfs_phylo <- bind_rows(list_of_distance_dfs_phylo, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_distances_all_conditions_phylo <- bind_rows(all_distances_all_conditions_phylo,distance_dfs_phylo)
    rm(distance_dfs_phylo,list_of_distance_dfs_phylo) # to keep memory open
    
 
    ########### time distances (from time calibrated phylogeny )#########
    list_of_distance_dfs_time <- lapply(listOfDFs,combo_function_phylo, variable =variable, ilr_or_clr_or_none = transform,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = timeDistances,epsilon=epsilon) # just put in time distances instead of raw phylo distances 
    # this will return a list of named results. How to deal with that? 
    # can bind them together and add some ids:
    distance_dfs_time <- bind_rows(list_of_distance_dfs_time, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_distances_all_conditions_time <- bind_rows(all_distances_all_conditions_time,distance_dfs_time)
    rm(distance_dfs_time,list_of_distance_dfs_time) # to keep memory open
    
  }
}





# need to keep in labels 
# note now as of 20220729 this doesn't re-downsample the spectra (good!)
combo_function_JustGetTransformedSpectra_notDistances <- function(spectrumdf,variable,clr_or_none,speciesToInclude,epsilon){
  procSpec <- processSpectra(spectrumdf,epsilon) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) # just get the cols for the species you want

  # do clr transform
  if(clr_or_none=="clr"){
    procSpec_transform <- data.frame(sapply(select(procSpec,-c(variable,mutation_label)),clr)) # this way you don't have to transpose and you get a much nicer shaped output
    procSpec_transform$mutation_label <- procSpec$mutation_label
    procSpec <- procSpec_transform 
  } else {
    
  }
  # want mutation labels to stay in here 
  return(procSpec)
}

all_processed_spectra_to_write_out = data.frame()

for(variable in listOfVars){
  for(transform in c("clr")){
    print("writing out processed spectra")
    list_of_processed_spectra_notdistances <- lapply(listOfDFs,combo_function_JustGetTransformedSpectra_notDistances, variable =variable, clr_or_none = transform,speciesToInclude=unique(all1merSpectra$label),epsilon=epsilon) # unique(all1merSpectra$label) is just to get all the species names. isn't 1mer specific 
    all_spectra <- bind_rows(list_of_processed_spectra_notdistances, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_spectra$variable <- variable
    all_spectra$transform_id <- transform 
    all_processed_spectra_to_write_out <- bind_rows(all_processed_spectra_to_write_out,all_spectra)
    rm(all_spectra,list_of_processed_spectra_notdistances) # to keep memory open
    
  }}

# then write this out.
write.table(all_processed_spectra_to_write_out,paste0(plotdir,'AllProcessedSpectra.AllConditions.AllVariables.AllSpecies.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.PROCESSED.txt'),quote=F,sep="\t",row.names = F)



################### MANTEL TEST ########################
# need dist matrix and phylo dist matrix of phylo distances
# may want to sq euclidean distances (bc expected to sale linearly with brownian motion apparently)


############# RUN MANTEL TEST ##################
manteldir=paste0(plotdir,"mantelTest/")
dir.create(manteldir,showWarnings = F)
comboFunction_mantelTest <- function(spectrum, phylo_dist,variable,speciesToInclude,mantelTestPermutationCount,epsilon,speciesCodes){
  spectrum_dist <- spectrum %>% 
    processSpectra(.,epsilon) %>% 
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(phylo_dist_reordered)

  mtr = vegan::mantel(xdis=phylo_dist_reordered,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}
# according to Hardy and Pavoine you get most power from euclidean dist vs sqrt(patristic)
# harmon and glor do euclidean squared vs patristic. either way you want one term to be squared or rooted because that is expectation for linearity under brownian motion
comboFunction_mantelTest_SQRTDISTANCE <- function(spectrum, phylo_dist,variable,speciesToInclude,mantelTestPermutationCount,epsilon,speciesCodes){
  spectrum_dist <- spectrum %>% 
    processSpectra(.,epsilon) %>% 
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(phylo_dist_reordered)
  phylo_dist_reordered_sqrt <- sqrt(phylo_dist_reordered)
  mtr = vegan::mantel(xdis=phylo_dist_reordered_sqrt,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel_SQRTPHYLODISTANCE",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}
# switching here to using a variable to set which variable to use 
mantel_test_results_raxml_list = lapply(listOfDFs,comboFunction_mantelTest,phylo_dist= raxml_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
# bind into a df:
mantel_test_results_raxml_df = bind_rows(mantel_test_results_raxml_list,.id="id")
mantel_test_results_raxml_df$distance_metric <- "raxml_tree"

mantel_test_results_raxml_SQRTDISTANCE_list =  lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE,phylo_dist= raxml_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
mantel_test_results_raxml_SQRTDISTANCE_df = bind_rows(mantel_test_results_raxml_SQRTDISTANCE_list,.id="id")
mantel_test_results_raxml_SQRTDISTANCE_df$distance_metric <- "raxml_tree_SQRT"

# you are here
mantel_test_results_time_list = lapply(listOfDFs,comboFunction_mantelTest,phylo_dist= timeTree_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
mantel_test_results_time_df = bind_rows(mantel_test_results_time_list,.id="id")
mantel_test_results_time_df$distance_metric <- "time_tree"

mantel_test_results_time_SQRTDISTANCE_list = lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE,phylo_dist= timeTree_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
mantel_test_results_time_SQRTDISTANCE_df = bind_rows(mantel_test_results_time_SQRTDISTANCE_list,.id="id")
mantel_test_results_time_SQRTDISTANCE_df$distance_metric <- "time_tree"






##################### plot correlations with mantel coefficients #############

######## making this look nice for manuscript :
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df,aes(x=0.5,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
  xlab("square root of phylogenetic distance (shared branch length)\nbranch lengths represent expected substitutions per site")+
  ylab("Aitchison distance between mutation spectra")
 # making strips not gray

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE,height=8,width=11)


###### plot just 1-mer spectrum for manuscript: 
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_Just1merSpectrumForMs <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$id=="1-mer+CpG spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id=="1-mer+CpG spectrum",],aes(x=0.4,y=Inf,vjust=1.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=6)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=16))+
  xlab("square root of phylogenetic distance (shared branch length)\nbranch lengths represent expected substitutions per site")+
  ylab("Aitchison distance between mutation spectra")

# making strips not gray

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_Just1merSpectrumForMs
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.1merOnly.ForManucript.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_Just1merSpectrumForMs,height=6,width=8)



######## plot mantel results with labels : #######
# add in nice labels:
# nice labels are loaded in from colors_and_labels file
all_distances_all_conditions_phylo <- merge(all_distances_all_conditions_phylo,niceLabels,by.x="item1",by.y="label")
all_distances_all_conditions_phylo <- merge(all_distances_all_conditions_phylo,niceLabels,by.x="item2",by.y="label",suffixes=c(".item1",".item2"))

all_distances_all_conditions_phylo$niceLabelForComparison <- paste0(all_distances_all_conditions_phylo$shortSpeciesLabel.item1,"-",all_distances_all_conditions_phylo$shortSpeciesLabel.item2)

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(color="gray")+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df,aes(x=0.35,y=Inf,vjust=1.3,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  geom_text(aes(label=niceLabelForComparison),size=2)+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=16))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between species' mutation spectra")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb,height=8,width=11)

####### plot jsut 1mer with labels for my manually labeling
mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr" &all_distances_all_conditions_phylo$id=="1-mer+CpG spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id =="1-mer=CpG spectrum",],aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=2)

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.Just1merToHelpMakeManualLabels.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels,height=8,width=11)


# compare with and without sqrt patterns:

########### mantel plot: time tree ###################
mantelTestCorrelationsPlot_time <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (",permCount," permutations)")))+
  geom_text(data=mantel_test_results_time_df,aes(x=150,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
  theme_bw()
mantelTestCorrelationsPlot_time
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.pdf"),mantelTestCorrelationsPlot_time,height=8,width=11)

########### timetree sqrt distance mantel plot #############
mantelTestCorrelationsPlot_time_SQRTDISTANCE <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (",permCount," permutations)")))+
  geom_text(data=mantel_test_results_time_SQRTDISTANCE_df,aes(x=5,y=Inf,vjust=1.15,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_time_SQRTDISTANCE_df$permCount)," permutations"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
  xlab("square root of phylogenetic distance (shared branch length)\nbranch length represent millions of years bp")+
  ylab("Aitchison distance between mutation spectra")
mantelTestCorrelationsPlot_time_SQRTDISTANCE
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.SQRTDISTANCE.pdf"),mantelTestCorrelationsPlot_time_SQRTDISTANCE,height=8,width=11)

######### >>> Mantel test comparing confounders to spectrum distance without correction for phylo relationships <<< #########
# NOTE: in a separate script, I also do mantel test to determine phylogenetic signal of confoudners (correlation with sqrt phylo distance)
# spectrum, phylo_dist,variable,speciesToInclude,mantelTestPermutationCount,epsilon
# note: not taking sqrt of distance, that's just for brownian motion of phylo distance. here am just correlating straight distances
# don't need species codes

confoundermantelTestPermutationCount=99999 # 100K instead of 1M  or 10M (phylomantel is quite slow -- it doesn't make a difference (I checked) unless any are hitting the pvalue ceiling)
print("note that confounder mantel test run with fewer perms")

confounderManteldir=paste0(plotdir,"/confoundersMantelTest",as.character(confoundermantelTestPermutationCount),"Permutations/")
dir.create(confounderManteldir,showWarnings = F)

########## FUNCTION TO RUN MANTEL TEST ON CONFOUNDERS: 
runMantelTest_onConfounders_vs_SpectrumDistances <- function(confounder_full_dataframe, nameOfConfounderYouWant,spectrum,speciesToInclude,variable,epsilon,confoundermantelTestPermutationCount){
  # restrict just to species you want:
  confounders_restrictedSpecies <- confounder_full_dataframe[confounder_full_dataframe$label  %in% speciesToInclude,] # in case you want to exclude some species from the analyses. for now it should be all 13
  # note that when you're doing this to match the phylogenetic tree you use the full_name_to_match_RAXML
  # but here you can use "label"! because both spectrum and confounders have same labels -- nice!
  # then restrict to the particular nameOfConfounderYouWant you want : 
  confounders_wide = pivot_wider(select(confounders_restrictedSpecies,c("label",nameOfConfounderYouWant)),values_from = all_of(nameOfConfounderYouWant),names_from = "label" )
  # get pairwise euclidean distances which just reduce down to abs(x-y) when it's just two numbers (did manual checks)
  # want to keep the diag and upper tri of matrix for mantel test even though redundant
  confounder_pairwise_distances = dist(t(confounders_wide),diag=T,upper=T) # want full matrix for mantel test  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  
  ##### get spectrum distances #######
  spectrum_dist <- spectrum %>% 
    processSpectra(.,epsilon) %>% 
    # merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% # dont need species codes because both confoudners and spectra have labels! only need when dealing with phylo dist stuff.
    filter(label %in% speciesToInclude) %>% # here filter on SPECIES from sppCodes not on label ! 
    pivotSpectra_perVariable(variable) %>% # don't need to keep full names, labels are fine
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  
  spectrum_dist_reordered <- as.matrix(spectrum_dist)[rownames(as.matrix(confounder_pairwise_distances)),colnames(as.matrix(confounder_pairwise_distances))] # get in same order ! 
  #head(phylo_dist_reordered)
  # run mantel test
  mtr = vegan::mantel(xdis=spectrum_dist_reordered,ydis=confounder_pairwise_distances,permutations = confoundermantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel_confounders",confounder=nameOfConfounderYouWant,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  
  return(mantelTestResult_df)
}

############# phylogenetically informed mantel test function ###################
# https://www.rdocumentation.org/packages/evolqg/versions/0.2-9/topics/PhyloMantel
# evolqg has : tree

# phylogenetic tree. Tip labels must match names in input matrices
# matrix.1 and matrix.2 with names that match phylo tree ; oh okay so tree isn't input as a matrix; k determines influence of phylo (1= strong influence, larger values converge to traditional mantel test )

runPHYLOMantelTest_onConfounders_vs_SpectrumDistances_CorrectsForPhyloSignal <- function(confounder_full_dataframe, nameOfConfounderYouWant,spectrum,speciesToInclude,variable,epsilon,confoundermantelTestPermutationCount,raxmlTree_renamedTips_subset){
  # restrict just to species you want:
  confounders_restrictedSpecies <- confounder_full_dataframe[confounder_full_dataframe$label  %in% speciesToInclude,] # note am using latin names here! 
  # note that when you're doing this to match the phylogenetic tree you use the full_name_to_match_RAXML
  # but here you can use "label"! because both spectrum and confounders have same labels -- nice!
  # then restrict to the particular nameOfConfounderYouWant you want : 
  # get names from "full_name_to_match_RAXML" instad of "label" because you're matching it to raxml tree! 
  confounders_wide = pivot_wider(select(confounders_restrictedSpecies,c("full_name_to_match_RAXML",nameOfConfounderYouWant)),values_from = all_of(nameOfConfounderYouWant),names_from = "full_name_to_match_RAXML" )
  # get pairwise euclidean distances which just reduce down to abs(x-y) when it's just two numbers (did manual checks)
  # want to keep the diag and upper tri of matrix for mantel test even though redundant
  confounder_pairwise_distances = dist(t(confounders_wide),diag=T,upper=T) # want full matrix for mantel test  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  
  ##### get spectrum distances #######
  spectrum_dist <- spectrum %>% 
    processSpectra(.,epsilon) %>% 
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% # getting full species names in 
    filter(label %in% speciesToInclude) %>% #  
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>% #KEEPING FULL SPECIES NAMES 
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  
  spectrum_dist_reordered <- as.matrix(spectrum_dist)[rownames(as.matrix(confounder_pairwise_distances)),colnames(as.matrix(confounder_pairwise_distances))] # get in same order ! 
  
  if(sum(colnames(as.matrix(confounder_pairwise_distances))!=colnames(as.matrix(spectrum_dist_reordered)))!=0){
    stop("something is wrong with your matrices column/row names!")
  }
  
  if(sum(!raxmlTree_renamedTips_subset$tip.label %in% colnames(as.matrix(confounder_pairwise_distances)))!=0){
    stop("tree contains species that aren't in matrices!")
  }
  
  if(sum(!colnames(as.matrix(confounder_pairwise_distances)) %in% raxmlTree_renamedTips_subset$tip.label)!=0){
    stop("matrices contain species that aren't in tree!")
  }
  # run mantel test
  #mtr = vegan::mantel(xdis=spectrum_dist_reordered,ydis=confounder_pairwise_distances,permutations = confoundermantelTestPermutationCount,method = "pearson")
  phylo_mtr= evolqg::PhyloMantel(tree=raxmlTree_renamedTips_subset,matrix.1 = as.matrix(spectrum_dist_reordered),matrix.2=as.matrix(confounder_pairwise_distances),permutations = confoundermantelTestPermutationCount,k = 1) # comparison func default is cor which has pearson as default
  # make a data frame
  phylomantelTestResult_df <- data.frame(call="evolqg_PhyloMantel",confounder=nameOfConfounderYouWant,statistic=as.numeric(as.list(phylo_mtr)$rho),permCount=as.numeric(confoundermantelTestPermutationCount),signif=((((as.numeric(as.list(phylo_mtr)$Probability))*as.numeric(confoundermantelTestPermutationCount))+1)/(as.numeric(confoundermantelTestPermutationCount)+1)))
  
  # seems to be called rho arbitrarily but is pearsons r. 
  # fixed how p values are calc'd from phylo mantel to match 
  ## AHA! this is an issue! phylomantel calculates probability as num with value > empirical/permutations
  # whereas vegan mantel calcs it as 1+num with value >empirical/permutations+1 . That extra 1 in the numerator makes a big difference, leading to phylomantel having more sig pvalues. So we need to rescale it
  # since phylomantel prob (f) = x / permuation_count, then x (num having value >emp) = f * perm_count.
  # So then can recalc p-value as (f*perm count)+ 1 / (perm count + 1)
  # this has big impacts at sig p values and little to no impact for non sig p values 
  
  return(phylomantelTestResult_df)
}


getConfounderAndSpectrumDistancesForPlotting <- function(confounder_full_dataframe, nameOfConfounderYouWant,spectrum,speciesToInclude,variable,epsilon){
  
  confounders_restrictedSpecies <- confounder_full_dataframe[confounder_full_dataframe$label  %in% speciesToInclude,] # in case you want to exclude some species from the analyses. for now it should be all 13
  # note that when you're doing this to match the phylogenetic tree you use the full_name_to_match_RAXML
  # but here you can use "label"! because both spectrum and confounders have same labels -- nice!
  # then restrict to the particular nameOfConfounderYouWant you want : 
  confounders_wide = pivot_wider(select(confounders_restrictedSpecies,c("label",nameOfConfounderYouWant)),values_from = all_of(nameOfConfounderYouWant),names_from = "label" )
  
  confounder_pairwise_distances_forPlotting = tidy(dist(t(confounders_wide),diag=F,upper=F)) # want full matrix for mantel test  
  colnames(confounder_pairwise_distances_forPlotting) <- c("item1","item2","confounder_distance")
  ##### get spectrum distances in tidy format for potting: 
  spectrum_dist_forPlotting <- spectrum %>% 
    processSpectra(.,epsilon) %>% 
    filter(label %in% speciesToInclude) %>% # here filter on SPECIES from sppCodes not on label ! 
    pivotSpectra_perVariable(variable) %>% # don't need to keep full names, labels are fine
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance(.)  
  
  ### add alphabetical labels (will be unique because you re-got distances using diag=f and upper=f: 
  confounder_pairwise_distances_forPlotting$comparisonLabel <- paste0(pmin(as.character(confounder_pairwise_distances_forPlotting$item1),as.character(confounder_pairwise_distances_forPlotting$item2)),".",pmax(as.character(confounder_pairwise_distances_forPlotting$item1),as.character(confounder_pairwise_distances_forPlotting$item2)))  # # need to use pmin and pmax to get row by row min max. otherwise gives min max of
  
  spectrum_dist_forPlotting$comparisonLabel <- paste0(pmin(as.character(spectrum_dist_forPlotting$item1),as.character(spectrum_dist_forPlotting$item2)),".",pmax(as.character(spectrum_dist_forPlotting$item1),as.character(spectrum_dist_forPlotting$item2)))
  # 20220905 I checked whether there are duplicates which would be bad, and there not.  
  
  # merge:
  merge_confounderDist_plusSpectrumDist <- merge(confounder_pairwise_distances_forPlotting,spectrum_dist_forPlotting,by="comparisonLabel",all = T)
  
  
  if(dim(merge_confounderDist_plusSpectrumDist)[1]!=dim(confounder_pairwise_distances_forPlotting)[1]){stop("something went wrong with merge")}
  
  return(merge_confounderDist_plusSpectrumDist)
  
}



### CAUTION: Make sure parameters are the same between the two functions !! 
all_confounder_spectrum_mantel_results = data.frame()
for(confounder in confoundersIWantToTest){
  
  mantel_test_results_comparingConfoundersToSpectra_withoutphylo_LIST <- lapply(listOfDFs,runMantelTest_onConfounders_vs_SpectrumDistances,confounder_full_dataframe=confounders,speciesToInclude=speciesToInclude_phylo,nameOfConfounderYouWant=confounder,confoundermantelTestPermutationCount=confoundermantelTestPermutationCount,epsilon=epsilon,variable = variableToUseInMantelTest) ### IF YOU CHANGE ANY OF THESE PARAMETERS BE SURE TO CHANGE THEM IN FUNCTION BELOW AS WELL!!!! 
  
  # turn into dataframe with spectrum in id column 
  mantel_test_results_comparingConfoundersToSpectra_withoutphylo_df = bind_rows(mantel_test_results_comparingConfoundersToSpectra_withoutphylo_LIST,.id="id") # bind over spectra
  
  # add to all-confounders dataframe:
  all_confounder_spectrum_mantel_results <- bind_rows(all_confounder_spectrum_mantel_results,mantel_test_results_comparingConfoundersToSpectra_withoutphylo_df) # has confounder as a column already # going to write this out at the end
  
  # now run PHYLO CORRECTED mantel test:
  print("starting PHYLO corrected mantel test of confounders")
  PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_LIST <- lapply(listOfDFs,runPHYLOMantelTest_onConfounders_vs_SpectrumDistances_CorrectsForPhyloSignal,confounder_full_dataframe=confounders,speciesToInclude=speciesToInclude_phylo,nameOfConfounderYouWant=confounder,confoundermantelTestPermutationCount=confoundermantelTestPermutationCount,epsilon=epsilon,variable = variableToUseInMantelTest,raxmlTree_renamedTips_subset=raxmlTree_renamedTips_subset)
  
  # make it a dataframe:
  PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_df <- bind_rows(PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_LIST,.id="id") # bind over spectra
  
  all_confounder_spectrum_mantel_results <- bind_rows(all_confounder_spectrum_mantel_results,PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_df) # has confounder as a column already # going to write this out at the end
  
  
  ######### MAKE *SURE* you're using all the same settings below as in the function above, otherwise plot and mantel values will be silently discordant (not ideal defensive coding, but works okay)
  distances_spectra_and_confounders_forPlottingWithMantelResults_LIST <- lapply(listOfDFs,getConfounderAndSpectrumDistancesForPlotting,confounder_full_dataframe=confounders,speciesToInclude=speciesToInclude_phylo,nameOfConfounderYouWant=confounder,epsilon=epsilon,variable = variableToUseInMantelTest) 
  
  # turn into dataframe with spectrum in id column 
  distances_spectra_and_confounders_forPlottingWithMantelResults_df = bind_rows(distances_spectra_and_confounders_forPlottingWithMantelResults_LIST,.id="id") # bind over spectra
  
  # make a plot without labels:
  # note: don't use sqrt of either distance! that's just a phylo distance thing
  mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_nolabels <- ggplot(distances_spectra_and_confounders_forPlottingWithMantelResults_df,aes(y=spectrum_distance,x=confounder_distance))+ # note need aes_string here so you can call confounder distance
    geom_point(size=1,color="dodgerblue")+
    # am assigning text position as middle value of x-axis confounder value so that is what that median( ) is doing there
    geom_text(data=mantel_test_results_comparingConfoundersToSpectra_withoutphylo_df,aes(x=median(distances_spectra_and_confounders_forPlottingWithMantelResults_df$confounder_distance),y=Inf,hjust=0.2,vjust=1.1,label=paste0("Mantel test (no phylo correction): ",confounder," Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)))),size=4)+
    geom_text(data=PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_df,aes(x=median(distances_spectra_and_confounders_forPlottingWithMantelResults_df$confounder_distance),y=Inf,hjust=0.2,vjust=2.5,label=paste0("Mantel test (phylo correction): ",confounder," Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)))),size=4)+
    
    ggtitle(paste0("Corrrelation between spectrum distance and ",confounder,"\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_comparingConfoundersToSpectra_withoutphylo_df$permCount)," permutations\n With and without correction for phylogenetic signal"))+
    theme_bw()+
    theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
    xlab(paste0(confounder," abs. diff."))+
    ylab("Aitchison distance")+
    facet_wrap(~id,scales="free_y")
  
  ggsave(paste0(confounderManteldir,confounder,".mantelTestResults.allspectra.pdf"),mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_nolabels,height=9,width=9)
  
  # make a plot with labels:
  mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_withlabels <- mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_nolabels+
    geom_text_repel(aes(label=comparisonLabel),size=1.3,max.overlaps=100)
  
  ggsave(paste0(confounderManteldir,confounder,".mantelTestResults.allspectra.labels.pdf"),mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_withlabels,height=9,width=9)
  
  
}

# write out:
write.table(all_confounder_spectrum_mantel_results,paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.forAnSITable.txt"),row.names=F,quote=F,sep="\t")

### plot confounders PLUS phylo signal mantel results: #####
all_confounder_spectrum_mantel_results_PLUSphylodist <- bind_rows(all_confounder_spectrum_mantel_results,mantel_test_results_raxml_SQRTDISTANCE_df)
all_confounder_spectrum_mantel_results_PLUSphylodist$confounder <- as.character(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder) # need to be characters so I can add a new one below

all_confounder_spectrum_mantel_results_PLUSphylodist[is.na(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder),]$confounder <- "sqrt_phylogenetic_distance"
all_confounder_spectrum_mantel_results_PLUSphylodist$confounder <- factor(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder,levels=c("sqrt_phylogenetic_distance","avg_coverage_per_individual","scaffold_N50","contig_N50","wattersons_theta"  ,"Rspan_d","AFR_d"))

# nice labels
all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting <- ""
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$call=="evolqg_PhyloMantel",]$niceLabelForPlotting <- "phyloMantel"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$call=="vegan_mantel_confounders",]$niceLabelForPlotting <- "uncorrected Mantel"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$call=="vegan_mantel_SQRTPHYLODISTANCE",]$niceLabelForPlotting <- "phylogenetic signal"

all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting <- factor(all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting,levels=c("phylogenetic signal","uncorrected Mantel","phyloMantel"))

### make nice confounder labels:
all_confounder_spectrum_mantel_results_PLUSphylodist$niceConfounderLabel <- ""
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$confounder=="scaffold_N50",]$niceConfounderLabel <- "scaffold N50"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$confounder=="contig_N50",]$niceConfounderLabel <- "contig N50"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$confounder=="avg_coverage_per_individual",]$niceConfounderLabel <- "coverage"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$confounder=="Rspan_d",]$niceConfounderLabel <- "reproductive span"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$confounder=="AFR_d",]$niceConfounderLabel <- "age at reproduction"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$confounder=="wattersons_theta",]$niceConfounderLabel <- "genetic diversity"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$confounder=="sqrt_phylogenetic_distance",]$niceConfounderLabel <- "phylogeny"

# order them:
all_confounder_spectrum_mantel_results_PLUSphylodist$niceConfounderLabel <- factor(all_confounder_spectrum_mantel_results_PLUSphylodist$niceConfounderLabel, levels=c("phylogeny","coverage","scaffold N50","contig N50","genetic diversity","reproductive span","age at reproduction"))

allConfounderMantelResultsPlot_plusPhyloPValues <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist),aes(x=niceConfounderLabel,y=-log10(signif),fill=niceLabelForPlotting))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  theme(legend.title = element_blank())+
  #geom_text(position = position_dodge(width = 1),aes(label=paste0("r =\n",round(statistic,2))),size=3)+
  theme(strip.background = element_rect(fill="lightblue"))+
  facet_grid(~id,switch="y")+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=10,"Paired"))))+
  ggtitle(paste0("Correlation of confounders with spectrum distance\n",confoundermantelTestPermutationCount," permutations"))
allConfounderMantelResultsPlot_plusPhyloPValues
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.phyloMantel.pdf"),allConfounderMantelResultsPlot_plusPhyloPValues,height=7, width=12)
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.phyloMantel.png"),allConfounderMantelResultsPlot_plusPhyloPValues,height=7, width=12,dpi = 300)

#### separate phylomantel test results into 1/3mer and 5/7mer #######

############### nice confounder plots for manuscript main text #############
# make a version JUST with non-phylo Mantel version and color by confounder! 
allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist)[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel",],aes(x=niceConfounderLabel,y=-log10(signif),fill=niceConfounderLabel))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
 #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  theme(legend.title = element_blank())+
  #geom_text(position = position_dodge(width = 1),aes(label=paste0("r =\n",round(statistic,2))),size=3)+
  theme(plot.margin = margin(13,10,10,10),strip.background = element_rect(fill="lightblue"))+
  facet_grid(~id,switch="y")+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=6,"Dark2"))))+
  ggtitle(paste0("Correlation of confounders with spectrum distance\n",confoundermantelTestPermutationCount," permutations"))
allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.NOPHYLOMANTEL.pdf"),allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL,height=6, width=12)


####### plot R values ###########

# annotate p values;
all_confounder_spectrum_mantel_results_PLUSphylodist$signif_star <- ""
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<0.05,]$signif_star <- "*"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<0.01,]$signif_star <- "**"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<0.001,]$signif_star <- "***"

all_confounder_spectrum_mantel_results_PLUSphylodist$signif_annotation <- "NS"
# annotate those below threshold as significant:
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),]$signif_annotation <- paste0(scientific(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),]$signif,1),"\n",all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),]$signif_star)

######### plot r values: ###########
# get ylim and ymax limits for plot (+20% to max so annotations show up)

ylimmin_rValues_1mer_nonsep=min(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel" & all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("1-mer+CpG spectrum"),]$statistic)

ylimmax_rValues_1mer_nonsep=max(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel" & all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("1-mer+CpG spectrum"),]$statistic)
# add 20% so that annotations will show up
ylimmax_rValues_1mer_nonsep=ylimmax_rValues_1mer_nonsep+.2*ylimmax_rValues_1mer_nonsep

allConfounderMantelResultsPlot_rValues_NOPHYLOMANTEL_1meronly <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel" & all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("1-mer+CpG spectrum"),]),aes(x=niceConfounderLabel,y=statistic,fill=niceConfounderLabel))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14),strip.background = element_rect(fill="lightblue"))+
  xlab("")+
  ylab("Pearson's r")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=20))+
  facet_wrap(~id,ncol=1)+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=6,"Dark2"))))+
  geom_text(aes(label = signif_annotation, x = niceConfounderLabel, y = statistic), position = position_dodge(width = 0.8), size=3,vjust=-0.05) +
  scale_y_continuous(limits=c(min(0,ylimmin_rValues_1mer_nonsep),ylimmax_rValues_1mer_nonsep))

allConfounderMantelResultsPlot_rValues_NOPHYLOMANTEL_1meronly

ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.rvaluesOnly.PLUSPHYLODISTFORCOMPARISON.NOPHYLOMANTEL.1merplusCpGonly.rvalues.pdf"),allConfounderMantelResultsPlot_rValues_NOPHYLOMANTEL_1meronly,height=4, width=6)



# allConfounderMantelResultsPlot_plusPhyloPValues_RVALUES <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist),aes(x=confounder,y=statistic,fill=niceLabelForPlotting))+
#   geom_col(position="dodge")+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
#   xlab("")+
#   ylab("Pearson's r")+
#   facet_grid(~id,switch="y")+
#   scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=10,"Paired"))))+
#   theme(legend.title = element_blank(),strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
#   ggtitle(paste0("Pearson's r values from Mantel test with ",confoundermantelTestPermutationCount," permutations"))
# 
# allConfounderMantelResultsPlot_plusPhyloPValues_RVALUES
# ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.RValues.PLUSPHYLODISTFORCOMPARISON.pdf"),allConfounderMantelResultsPlot_plusPhyloPValues_RVALUES,height=6, width=12)


############# Calculate pairwise cosine similarities between non-CLR transformed mutation spectra COUNTS (not proportions(?) ###############
# this is new
# uses sigfit's cosine_sim function to keep things similar to my sigfit cosine similarity calculations

# use the human-rescaled and downsampled spectra (don't need 'process' spectra or to do ilr or clr transform (whole point is to not do those so I compare how cosine simlarity performs))

head(all1merSpectra)

cosineDISSimilarity_keepMatrix_UPPERANDLOWER <- function(pivotedspectradf_noIDVars){
    # make it a matrix with each row as a species (transposed and all numeric)
    # variable in question:
    # distances are row-wise nto column wise so need to transpose 
    cosineSimilarity <- lsa::cosine(as.matrix(pivotedspectradf_noIDVars)) #  # transpose for dist otherwise it goes along rows and is wrong
    cosineDISSimilarity <- 1-cosineSimilarity
    # want cosine DISsimiliarity (1-cosine similarity) so that it scales up with phylo distance the same way the other distance measures do
    # this results in cosine disimilarity matrix
    return(cosineDISSimilarity)
  }

cosineDISSimilarity_tidiedIntoTable <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  cosineSimilarity <- lsa::cosine(as.matrix(pivotedspectradf_noIDVars)) #  # transpose for dist otherwise it goes along rows and is wrong
  cosineDISSimilarity <- 1-cosineSimilarity
  cosineDISSimilarity_tidy <- tidy(as.dist(cosineDISSimilarity)) # makes it a distance matrix (doesn't change values, I checked and then tidies it into a dataframe using broom)
  # want cosine DISsimiliarity (1-cosine similarity) so that it scales up with phylo distance the same way the other distance measures do
  # this results in cosine disimilarity matrix
  # rename as item 1 item2 cosine_dis
  colnames(cosineDISSimilarity_tidy) <- c("item1","item2","cosine_dissimilarity")
  return(cosineDISSimilarity_tidy)
}


comboFunction_mantelTest_SQRTDISTANCE_USINGCOSINEDISSIMILARITY <- function(spectrum, phylo_dist,variable,speciesToInclude,mantelTestPermutationCount,epsilon,speciesCodes){
  spectrum_dist <- spectrum %>% 
    processSpectra(.,epsilon) %>% 
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,"none") %>% # fixing as NO transformation; this step just removes variable columns
    cosineDISSimilarity_keepMatrix_UPPERANDLOWER(.)   # calculates cosine DISsimilarity (1-cosine similarity) as a matrix 
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(phylo_dist_reordered)
  phylo_dist_reordered_sqrt <- sqrt(phylo_dist_reordered)
  mtr = vegan::mantel(xdis=phylo_dist_reordered_sqrt,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel_SQRTPHYLODISTANCE_COSINEDISIMILARITY",variable="cosine_disimilarity",method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}


mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_list =  lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE_USINGCOSINEDISSIMILARITY,phylo_dist= raxml_cophenetic_dist,variable="total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS",speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
# you are here ! 

mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_list
mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_df = bind_rows(mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_list,.id="id")
mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_df$distance_metric <- "raxml_tree_SQRT"


# get cosine similarities to plot:
getCosineDisimilaritiesForPlotting <- function(spectrum, variable,speciesToInclude,epsilon,speciesCodes){
  
  cosine_table <- spectrum %>% 
    processSpectra(.,epsilon) %>% 
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,"none") %>% # fixing as NO transformation; this step just removes variable columns
    cosineDISSimilarity_tidiedIntoTable(.)
    
  
  
  return(cosine_table)
  
}

cosineDisimilarityTables <- lapply(listOfDFs,getCosineDisimilaritiesForPlotting,variable="total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS",speciesToInclude=speciesToInclude_phylo,epsilon=epsilon,speciesCodes=speciesCodes)

# collapse into one dataframe:

cosineDisimilarityTables_id <- bind_rows(cosineDisimilarityTables,.id="id")

# merge with phylo distances:
cosineDisimilarityTables_id$comparisonLabel <- paste0(pmin(as.character(cosineDisimilarityTables_id$item1),as.character(cosineDisimilarityTables_id$item2)),".",pmax(as.character(cosineDisimilarityTables_id$item1),as.character(cosineDisimilarityTables_id$item2)))  # alphabetical

# check that all are in phylo comp labels:
sum(!cosineDisimilarityTables_id$comparisonLabel %in% phyloDistances$comparisonLabel) # cool they are all in there.

cosineDisimilarityTables_mergedWithPhyloDistance <- merge(cosineDisimilarityTables_id,phyloDistances,by="comparisonLabel") # dim stays the same! 

# get cosine similarity as well:
cosineDisimilarityTables_mergedWithPhyloDistance$cosine_similarity <- 1-cosineDisimilarityTables_mergedWithPhyloDistance$cosine_dissimilarity

# now can plot:
cosine_dissimilarity_plot_SQRTPHYLODIST <- ggplot(cosineDisimilarityTables_mergedWithPhyloDistance,aes(x=sqrt(cophenetic_distance),y=cosine_dissimilarity))+
  geom_point(size=1)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_df,aes(x=0.5,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
  xlab("square root of phylogenetic distance (shared branch length)\nbranch lengths represent expected substitutions per site")+
  ylab("Cosine dissimilarty between mutation spectra")

ggsave(paste0(manteldir,"mantelTest.COSINESIMILARITY.SQRTDISTANCE.pdf"),cosine_dissimilarity_plot_SQRTPHYLODIST,height=8,width=11)

# plot as violin plot:
head(cosineDisimilarityTables_mergedWithPhyloDistance)


# all types:
cosineSimilarityViolinPlot1 <- ggplot(cosineDisimilarityTables_mergedWithPhyloDistance,aes(x=id,y=cosine_similarity))+
  geom_violin(fill="lightblue",alpha=0.5)+
  theme_bw()+
  ylab("Pairwise cosine similarity between\nspecies' empirical mutation spectra")+
  theme(text=element_text(size=14),axis.text.x = element_text(angle=45,hjust=1))+
  xlab("")+
  scale_y_continuous(breaks=seq(0.75,1,by=0.05))
cosineSimilarityViolinPlot1
ggsave(paste0(plotdir,"COSINESIMILARITY.allTypes.pdf"),cosineSimilarityViolinPlot1,height=4,width=8)



############ save workspace image ###################
save.image(file = paste0(plotdir,"RWORKSPACE.",todaysdate,".AtEndOfScript.image.RData"))

