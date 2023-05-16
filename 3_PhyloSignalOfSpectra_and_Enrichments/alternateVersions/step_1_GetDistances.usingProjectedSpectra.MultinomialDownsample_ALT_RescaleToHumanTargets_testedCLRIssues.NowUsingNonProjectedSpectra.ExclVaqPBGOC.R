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
#separateCpGs="no" # yes or no -- 20221208 note I have removed the code that separates cpgs from this script; it's in a separate script now. re-ran the start of the script to make sure it still works; it does.
todaysdate=format(Sys.Date(),"%Y%m%d")
maskLabel="maskALL" # if species get different mask labels this will be a pain.
mantelTestPermutationCount=9999999 # trying with more (10M now)
# will be quite slow
flagOfAnalysis=paste0(todaysdate,".plots.SubsampledIndividualsNotProjection.epsilon.",epsilon,".sepCpGs.no.masked",maskLabel,".RESCALEDTOHUMANTARGETS.SimplifiedCode.",as.character(mantelTestPermutationCount),".MantelPermutations.ExclVaqPBFwGOC") # a label that will be in title of plots and outdir 
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
speciesList=c('humans', 'Mus_musculus','Mus_spretus' ,'brown_bear', 'fin_whale', 'Gorilla_gorilla', 'Pan_troglodytes', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus','wolves') # overarching species list  # taking out vaq and pb


speciesToInclude_in_raxmlTree=c("Balaenoptera_physalus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Canis_lupus")  # raxml name fmt
# # taking out pb and vaq

speciesToInclude_phylo=c("Mus_musculus_Mmd","Mus_spretus_Ms","brown_bear_ABC","fin_whale_ENP","humans_AFR","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","Pongo_abelii","Pongo_pygmaeus","wolves") # human_AFR# switching to ENP fin whale for more diversity; excluding ucla wolves.



############ read in confounders dataframe ##############
confounders <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/information/ALT_FinWhaleENP_confoundersToPlotPhyloSignal.GenomeN50.etc.updatedthetavals_updatedCoverage_alt.txt",sep="\t",header=T) # updated values of theta and coverage to be based on teh 5 inds that were picked at random per pop
# this is an alt file that subs in fin whale ENP information for GOC since I am using ENP in this analysis.
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
  print("note that CpGs are not separated out!")
  allSpectra$mutation_1mer <- paste0(substr(allSpectra$mutation_type,4,4),".",substr(allSpectra$mutation_type,12,12)) # if you don't want to separate out CpGs then this is just the regular center bp
  # still get ancestral1mer (but don't specify if CpG or not)
  allSpectra$ancestral1mer <- substr(allSpectra$mutation_type,4,4) # if you don't want to separate out CpGs then this is just the regular center bp




# get spectrum summed up per type
all7merSpectraOnly <- allSpectra %>%
  # adding groupbing by central 1mer but that doesn't actually change the sums since 7mers already group by those 
  # I just want to keep it in the dataframe 
  group_by(species,population,label,ancestral7mer,mutation_7mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>%
  ungroup() # this will be identical to original spectrum (7mer), just doing for consistency of formatting

all5merSpectraOnly <- allSpectra %>%
  group_by(species,population,label,ancestral5mer,mutation_5mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>%
  ungroup()

all3merSpectraOnly <- allSpectra %>%
  group_by(species,population,label,ancestral3mer,mutation_3mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>%
  ungroup()

all1merSpectraOnly <- allSpectra %>%
  group_by(species,population,label,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>%
  ungroup()
# if you don't want to partition over CpGs you can use mutation_1mer_CpGNotLabeled instead
# okay NOTE HERE: "C" counts are nonCpG and then C_CpG are CpG ancestral targets


########## read in targets ##############
all7merTargets=data.frame()
all5merTargets=data.frame()
all3merTargets = data.frame()
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
  

  # TARGETS per 3mer type #
  targets_7mer <- targets %>% group_by(species,target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup() # original targets file, just doing to be consistent but it's identical to original targets
  targets_3mer <- targets %>% group_by(species,target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup()
  targets_5mer <- targets %>% group_by(species,target_5mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup()
  targets_1mer <- targets %>% group_by(species,target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup()

  # make all species targets:
  
  
  # combine  all targets:
  all7merTargets = bind_rows(all7merTargets,targets_7mer)
  all5merTargets = bind_rows(all5merTargets,targets_5mer)
  all3merTargets = bind_rows(all3merTargets,targets_3mer)
  all1merTargets = bind_rows(all1merTargets,targets_1mer)
  
  
}

############ get target proportions (new as of 20220727) ##############
all1merTargets <- all1merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()


all3merTargets <- all3merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()

all5merTargets <- all5merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()


all7merTargets <- all7merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()
# plot targets:
target1merPlot <- ggplot(all1merTargets,aes(y=species,x=total_target_count,fill=target_1mer))+
  geom_col(position="dodge")
target1merPlot
ggsave(paste0(plotdir,"targetPlot.1mers.png"),target1merPlot,height=9,width=5)

# plot proportions:
target1merPlot_frac <-  ggplot(all1merTargets,aes(y=species,x=target_proportion,fill=target_1mer))+geom_col(position='dodge')+facet_wrap(~target_1mer,scales="free")
target1merPlot_frac
ggsave(paste0(plotdir,"targetPlot.1mers.frac.png"),target1merPlot_frac,height=9,width=15)

########### fill in 0-entry types #############
all7merSpectraOnly_filledin <- all7merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_7mer,ancestral7mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations=0))%>%
  mutate(mutation_label=mutation_7mer)

if(dim(all7merSpectraOnly_filledin)[1]==dim(all7merSpectraOnly)[1]){
  print("filling in didn't work!")
}
all5merSpectraOnly_filledin <- all5merSpectraOnly %>%
  ungroup() %>% # need to pre-ungroup ; got grouped up above
  complete(nesting(mutation_5mer,ancestral5mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations=0))%>%
  mutate(mutation_label=mutation_5mer) # adding mutation label for use in functions

all3merSpectraOnly_filledin <- all3merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_3mer,ancestral3mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations=0))%>%
  mutate(mutation_label=mutation_3mer)

all1merSpectraOnly_filledin <- all1merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_1mer,ancestral1mer),nesting(species,population,label),fill=list(total_mutations=0)) %>%
  mutate(mutation_label=mutation_1mer)
# fill in should only be needed for 7mer (won't change dim of others)

# so when I was merging before all species were combined there was a bug where missing mtuation types wouldn't get targets filled in which led to them being NA downstream 
# okay now fill in missing mutation types : as 0s PRIOR TO MERGING WITH TARGETS
#  note that targets are the same within species (mouse1 and mouse2 have same targets; bear 1 and bear2 have targets)


############## merge with targets ###############

all7merSpectra_plusTargets <- merge(all7merSpectraOnly_filledin, all7merTargets,by.x=c("species","ancestral7mer"),by.y=c("species","target_7mer"))

if(dim(all7merSpectra_plusTargets)[1] != dim(all7merSpectraOnly_filledin)[1]){
  print("something went wrong with merge!")
  
}

# write out for enrichments
write.table(all7merSpectra_plusTargets,paste0(plotdir,"all7merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")


all5merSpectra_plusTargets <- merge(all5merSpectraOnly_filledin, all5merTargets,by.x=c("species","ancestral5mer"),by.y=c("species","target_5mer"))

write.table(all5merSpectra_plusTargets,paste0(plotdir,"all5merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")

all3merSpectra_plusTargets <- merge(all3merSpectraOnly_filledin, all3merTargets,by.x=c("species","ancestral3mer"),by.y=c("species","target_3mer"))

write.table(all3merSpectra_plusTargets,paste0(plotdir,"all3merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")


all1merSpectra_plusTargets <- merge(all1merSpectraOnly_filledin, all1merTargets,by.x=c("species","ancestral1mer"),by.y=c("species","target_1mer"))
# don't need to write out the 1mers because they aren't used as numerator in enrichments 
# you are ehre
########### WRITE THESE OUT FOR TRANSfERRING TO SAGE to carry out fishers exact test #################

############### >>>> RESCALE TO HUMAN TARGETS <<<<< #################
# rescaling all counts to be scaled by human targets total_mutations_of_typeX_in_SpeciesA* (target_proportion_HUMAN /target_proportion_speciesA)

humanTargetProportions_1mer <- all1merTargets[all1merTargets$species=="humans",]
humanTargetProportions_3mer <- all3merTargets[all3merTargets$species=="humans",]
humanTargetProportions_5mer <- all5merTargets[all5merTargets$species=="humans",]
humanTargetProportions_7mer <- all7merTargets[all7merTargets$species=="humans",]

# note these don't yet have the human counts; need to run processSpectra_AdjustForHumanTargets on them (is done in combo distance script at bottom fo the script )
# going to keep these names 
all1merSpectra_plusHumanInfoNotYetRescaled <- merge(all1merSpectra_plusTargets,humanTargetProportions_1mer[,c("target_1mer","total_target_count","target_proportion")],by.x=c("ancestral1mer"),by.y=c("target_1mer"),suffixes=c("_thisSpecies","_HUMAN"))

all3merSpectra_plusHumanInfoNotYetRescaled <- merge(all3merSpectra_plusTargets,humanTargetProportions_3mer[,c("target_3mer","total_target_count","target_proportion")],by.x=c("ancestral3mer"),by.y=c("target_3mer"),suffixes=c("_thisSpecies","_HUMAN"))

all5merSpectra_plusHumanInfoNotYetRescaled <- merge(all5merSpectra_plusTargets,humanTargetProportions_5mer[,c("target_5mer","total_target_count","target_proportion")],by.x=c("ancestral5mer"),by.y=c("target_5mer"),suffixes=c("_thisSpecies","_HUMAN"))

all7merSpectra_plusHumanInfoNotYetRescaled <- merge(all7merSpectra_plusTargets,humanTargetProportions_7mer[,c("target_7mer","total_target_count","target_proportion")],by.x=c("ancestral7mer"),by.y=c("target_7mer"),suffixes=c("_thisSpecies","_HUMAN"))

# okay you want to get spectrum and the mutation rate and want to add 1 to all counts to deal with 0 entries (assuming every mutation type would happen if you sampled enough -- though note this will have a bigger i)
print("rescaling to human targets")

# NOTE! need to be dividing by HUMAN target sizes for mutation rates!!! 
#### scale counts to human targets: 

all1merSpectra_HumanRescaled <- all1merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

all3merSpectra_HumanRescaled <- all3merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

all5merSpectra_HumanRescaled <- all5merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

all7merSpectra_HumanRescaled <- all7merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

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
all3merSpectra <- downsampleFunc_RescaledByHumanTargets(all3merSpectra_HumanRescaled,plotdir)
all5merSpectra <- downsampleFunc_RescaledByHumanTargets(all5merSpectra_HumanRescaled,plotdir)
all7merSpectra <- downsampleFunc_RescaledByHumanTargets(all7merSpectra_HumanRescaled,plotdir)


# write these out:
write.table(all1merSpectra,paste0(plotdir,"all1merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
write.table(all3merSpectra,paste0(plotdir,"all3merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
write.table(all5merSpectra,paste0(plotdir,"all5merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
write.table(all7merSpectra,paste0(plotdir,"all7merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)

############ also make a set of folded spectra ##################
### as of 20220728 this needs the human rescaled version to be present
# also getting this to work on downsampled version 
foldSpectrum <- function(allXmerSpectra){
  # get ancestral and derived states by splitting X.Y
  allXmerSpectra$ancestral <- unlist(lapply(strsplit(allXmerSpectra$mutation_label,"\\."),"[",1))
  allXmerSpectra$derived <- unlist(lapply(strsplit(allXmerSpectra$mutation_label,"\\."),"[",2))
  # already have central bp as "mutation_1mer"
  # set up the matches (opposite mutation direction so C-->D match is D -->C )
  # but if this will result in a non A or C central bp then you need to rev comp the sequences
  allXmerSpectra$match <- ""
  # get matches for C>A mutations which don't have issues with rev comping 
  allXmerSpectra[allXmerSpectra$mutation_1mer %in% c("A.C","C.A"),]$match <- paste0(allXmerSpectra[allXmerSpectra$mutation_1mer %in% c("A.C","C.A"),]$derived,".",allXmerSpectra[allXmerSpectra$mutation_1mer %in% c("A.C","C.A"),]$ancestral)
  
  # get matches for all other mutation types:
  # THIS DOES NOT WORK FOR 1mers!!!!!! gives you totally the wrong thing because it rev comps the whole vector of 1mers . so need do manually -- so annoying 
  if(allXmerSpectra[1,]$mutation_1mer==allXmerSpectra[1,]$mutation_label){
    print("is this a 1mer spectrum? note that spgs can't accurately rev comp these types so I'm doing it manually in this function. check results carefully!")
    
    # manually set matches: (A.C and C.A are already done above correclty)
    allXmerSpectra[allXmerSpectra$mutation_label=="A.G",]$match <- "C.T" # match is G>A but then need to rev comp so it's C>T
    allXmerSpectra[allXmerSpectra$mutation_label=="C.T",]$match <- "A.G" # match is T>C but then need to rev comp so it's A>G
    allXmerSpectra[allXmerSpectra$mutation_label=="A.T",]$match <- "A.T" # match is T>A but then need to rev comp so it's A>T again! these are already totally collapsed
    allXmerSpectra[allXmerSpectra$mutation_label=="C.G",]$match <- "C.G" # match is G>C but then need to rev comp so it's C>G again! these are already totally collapsed
    
    
    
  } else if(nchar(allXmerSpectra[1,]$mutation_label)>3){ # if mut label (X.Y) is >3 then it's not a 1mer spectrum
    print("this isn't a 1-mer spectrum, correct? then it's okay to use spgs.")
    allXmerSpectra[!(allXmerSpectra$mutation_1mer %in% c("A.C","C.A")),]$match <- paste0(reverseComplement(allXmerSpectra[!(allXmerSpectra$mutation_1mer %in% c("A.C","C.A")),]$derived,case = "upper"),".",reverseComplement(allXmerSpectra[!(allXmerSpectra$mutation_1mer %in% c("A.C","C.A")),]$ancestral,case="upper"))
    
  }
  # okay if dealing with 1mers have to be really careful with spgs because it will treat 1mers as one lnog sequence and reverse the whole thing
  
  # okay then want to match them up with a label that will the the same for the matches
  allXmerSpectra$foldingGroup <- paste0(pmin(allXmerSpectra$mutation_label,allXmerSpectra$match),"_plus_",pmax(allXmerSpectra$mutation_label,allXmerSpectra$match))
  
  # if match == mutation label then there will be only 1 folding group (see notes below)
  
  print(head(allXmerSpectra))
  # for 7mers:
  # there are 128 mutation types where reverse is totally symmetric eg TTTCAAA > TTTGAAA : allXmerSpectra[allXmerSpectra$match==allXmerSpectra$mutation_label,]
  # so those will end up with just one folding group even after filling in because  they are already totally collapsed with their reverse. each species has 128 of them. folding group looks like "CACCGTG.CACGGTG_plus_CACCGTG.CACGGTG" where it's just the same mutation type repeated. I think it's okay to leave it like that rather than changing the label. should still work out. 
  # NOTE that I'm NOT double adding these types together, that's just a little quirk of the foldingGroup naming process. 
  # note that you end up with 24,576 original types but 128 are already collapsed to matches so you end up with (24576-128 that don't have a match/2)  + 128 (the ones that don't match). this works out to 12352 and that works! ok cool. so it's not exactly 24,576/2 because of those 128 collapsed types that are symmetrical.
  # okay need to make this work with the other types too
  # for 5mers there are 32 types that don't match so (1536-32)/2 + 32 = 784 folding groups (checked)
  # for 3mers there are 8 types that don't match so (96-8)/2 + 8 = 52 (checked)
  # for 1mers there are 2 types that don't match (A>T and C>G) so (6-2)/2 + 2 = 4 groups (good)
  
  # okay now want to sum over the groups. do I also want to sum the targets? I think yes?
  ######## THIS IS MODIFIED TO WORK ON HUMAN COUNTS INSTEAD
  # also getting original counts and downsampled counts as well
  allXmerSpectra_folded <- allXmerSpectra %>%
    group_by(species,label,population,foldingGroup) %>%
    summarise(total_mutations_ORIGINAL=sum(total_mutations_ORIGINAL),total_mutations_RESCALED_BY_HUMAN_TARGETS= sum(total_mutations_RESCALED_BY_HUMAN_TARGETS),total_target_count_HUMAN=sum(total_target_count_HUMAN),total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS=sum(total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS)) %>% 
    ungroup() %>%
    rename(mutation_label = foldingGroup)  #keep names the same for downstream 
  # need to rename foldingGroup to mutation_label so it works with things downstream 
  # I think this works for targets because you want to know the target space of each of the matches in the fold group so if you have X.Y that has X target count and Y.X with Y target count you'd add up the X+Y targets for the full target space (?) where these mutations could ahve occurred? Do we need to scale it in any way?
  
  # and need to rename 
  # let's check that total mutations have remained the same
  totalSites_folded <- allXmerSpectra_folded %>%
    group_by(species,population,label) %>%
    summarise(totalSites_RESCALED_BY_HUMAN_TARGETS=sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    ungroup()
  
  totalSites_unfolded <- allXmerSpectra %>%
    group_by(species,population,label) %>%
    summarise(totalSites_RESCALED_BY_HUMAN_TARGETS=sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    ungroup()
  
  totalSites_compareFoldUnfold <- merge(totalSites_folded,totalSites_unfolded,by=c("species","population","label"),all=T,suffixes = c(".folded",".unfolded"))
  
  if(!isTRUE(all.equal(totalSites_compareFoldUnfold$totalSites_RESCALED_BY_HUMAN_TARGETS.folded,totalSites_compareFoldUnfold$totalSites_RESCALED_BY_HUMAN_TARGETS.unfolded))){
    print("something's wrong -- different number of sites at end") # need to use all.equal to deal with floating point issues
    
  }
  
  return(allXmerSpectra_folded) # returns folded spectrum ! 
}

####### fold the 4 spectra ############
listOfDFs_forFolding = list(`1-mer spectrum`=all1merSpectra,`3-mer spectrum`=all3merSpectra,`5-mer spectrum`=all5merSpectra,`7-mer spectrum` = all7merSpectra)

listOfFoldedSpectra <- lapply(listOfDFs_forFolding,foldSpectrum) 

# checked BB ABC A>C+ C>A manually, added up correctly!
# write these out:
write.table(listOfFoldedSpectra$`1-mer spectrum`,paste0(plotdir,"all1merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.FOLDED.txt"),quote=F,sep="\t",row.names = F)
write.table(listOfFoldedSpectra$`3-mer spectrum`,paste0(plotdir,"all3merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.FOLDED.RESCALED_BY_HUMAN_TARGETS.txt"),quote=F,sep="\t",row.names = F)
write.table(listOfFoldedSpectra$`5-mer spectrum`,paste0(plotdir,"all5merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.FOLDED.RESCALED_BY_HUMAN_TARGETS.txt"),quote=F,sep="\t",row.names = F)
write.table(listOfFoldedSpectra$`7-mer spectrum`,paste0(plotdir,"all7merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.FOLDED.txt"),quote=F,sep="\t",row.names = F)



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
test = processSpectra(all7merSpectra,epsilon) # ,plotdir,writeOut=F
test3mers = processSpectra(all3merSpectra,epsilon) # ,plotdir,writeOut=F
test1mers = processSpectra(all1merSpectra,epsilon) # ,plotdir,writeOut=F

sum(test[test$label=="vaquita",]$total_mutations_ORIGINAL)
sum(test[test$label=="vaquita",]$total_mutations_RESCALED_BY_HUMAN_TARGETS)
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
testpivoted = pivotSpectra_perVariable(test,testvar)
head(testpivoted)
######## CHANGE THIS IF CHANGE TEST KMER: 
testpivoted$centralBP <- paste0(substr(testpivoted$mutation_label,4,4),".",substr(testpivoted$mutation_label,12,12))
ggplot(testpivoted,aes(x=wolves,y=Mus_musculus_Mmd,color="wolf-mmd"))+
  geom_point()+
  geom_abline()+
  geom_point(data=testpivoted,aes(x=wolves,y=Mus_spretus_Ms,color="wolf-ms"))+
  #geom_text(data=testpivoted,aes(label=mutation_label))+
  #geom_point(data=testpivoted,aes(x=wolves,y=fin_whale_ENP,color="wolf-fwENP"))+
  #geom_point(data=testpivoted,aes(x=wolves,y=fin_whale_GOC,color="wolf-gocENP"))+
  geom_point(data=testpivoted,aes(x=wolves,y=brown_bear_ABC,color="wolf-BB"))+
  facet_wrap(~centralBP,scales="free")


bearplot1 <- ggplot(testpivoted,aes(x=polar_bear_PB,y=brown_bear_ABC,color="PB-BB"))+
  geom_point()+
  #geom_point(aes(x=polar_bear_PB,y=Mus_musculus_Mmd,color="PB-Mmd"))+
  geom_abline()+
  facet_wrap(~centralBP,scales="free")+
  ggtitle(testvar)
bearplot1
ggsave(paste0(trblshootdir,"comparing.bears.png"),bearplot1)
  
bearplot2 <- ggplot(testpivoted,aes(x=polar_bear_PB,y=brown_bear_ABC,color=centralBP))+
  geom_point()+
  geom_abline()+
  ggtitle(testvar)
bearplot2
ggsave(paste0(trblshootdir,"comparing.bears2.png"),bearplot2)

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
listOfDFs = list(`1-mer spectrum`=all1merSpectra,`3-mer spectrum`=all3merSpectra,`5-mer spectrum`=all5merSpectra,`7-mer spectrum` = all7merSpectra)
# lapply applies combo_function_phylo to each of the above spectra ^^ 
# if you set the names of a list, they get preserved in the output of lapply! cool!


listOfVars=c("mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS") # no longer using the normalized versions ; what about nondownsampled frac seg sites? adding in total mutations as well
listOfTransforms=c("none","clr") # skipping ilr for now (exhausts vector memory for 7mers, maybe just show 1/3mer results.)

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


############# Get distances for FOLDED spectra ################### 
all_distances_all_conditions_phylo_FOLDED = data.frame()
all_distances_all_conditions_time_FOLDED = data.frame()
for(variable in listOfVars){
  for(transform in listOfTransforms){
    # get for all dfs: 
    ########## phylo distances ###############
    list_of_distance_dfs_phylo <- lapply(listOfFoldedSpectra,combo_function_phylo, variable =variable, ilr_or_clr_or_none = transform,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances,epsilon=epsilon)
    # this will return a list of named results. How to deal with that? 
    # can bind them together and add some ids:
    distance_dfs_phylo <- bind_rows(list_of_distance_dfs_phylo, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_distances_all_conditions_phylo_FOLDED <- bind_rows(all_distances_all_conditions_phylo_FOLDED,distance_dfs_phylo)
    rm(distance_dfs_phylo,list_of_distance_dfs_phylo) # to keep memory open
    

    
    ########### time distances (from time calibrated phylogeny )#########
    list_of_distance_dfs_time <- lapply(listOfFoldedSpectra,combo_function_phylo, variable =variable, ilr_or_clr_or_none = transform,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = timeDistances,epsilon=epsilon) # just put in time distances instead of raw phylo distances 
    # this will return a list of named results. How to deal with that? 
    # can bind them together and add some ids:
    distance_dfs_time <- bind_rows(list_of_distance_dfs_time, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_distances_all_conditions_time_FOLDED <- bind_rows(all_distances_all_conditions_time_FOLDED,distance_dfs_time)
    rm(distance_dfs_time,list_of_distance_dfs_time) # to keep memory open
    
  }
}

write.table(all_distances_all_conditions_phylo,paste0(plotdir,"all_distances_all_conditions_phylo.RESCALED_BY_HUMAN_TARGETS.txt"),row.names=F,quote=F,sep="\t")



write.table(all_distances_all_conditions_time,paste0(plotdir,"all_distances_all_conditions_timeCalibrated.RESCALED_BY_HUMAN_TARGETS.txt"),row.names=F,quote=F,sep="\t")


##### write out folded versions as well #####
write.table(all_distances_all_conditions_phylo_FOLDED,paste0(plotdir,"all_distances_all_conditions_phylo.FOLDEDSPECTRUM.RESCALED_BY_HUMAN_TARGETS.txt"),row.names=F,quote=F,sep="\t")



write.table(all_distances_all_conditions_time_FOLDED,paste0(plotdir,"all_distances_all_conditions_timeCalibrated.FOLDEDSPECTRUM.RESCALED_BY_HUMAN_TARGETS.txt"),row.names=F,quote=F,sep="\t")


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
  for(transform in c("clr","none")){
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

######### plot CLR values and geomeans -- a troublehsooting section (not necessary) ########
# NOT plotting normalized version 
all_processed_spectra_to_write_out_melt <- melt(all_processed_spectra_to_write_out,variable.name = "label")
ggplot(all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_processed_spectra_to_write_out_melt$label %in% speciesToInclude_phylo & all_processed_spectra_to_write_out_melt$id=="7-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="clr",],aes(x=mutation_label,y=value,color=label))+
  theme(axis.text.x=element_blank())+geom_point()+ggtitle("not using normalized mut rate anymore")
# fixed bug in this plot that was plotting multiple variables ! 
# no bands in human scaled version! that is promising...

# plot a subset
muttypesall=unique(all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_processed_spectra_to_write_out_melt$label %in% speciesToInclude_phylo & all_processed_spectra_to_write_out_melt$id=="7-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="clr",]$mutation_label)
randomsubsetofmutations=sample(muttypesall,100)



### PLOT SUBSET
ggplot(all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$mutation_label %in% randomsubsetofmutations & all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_processed_spectra_to_write_out_melt$label %in% speciesToInclude_phylo & all_processed_spectra_to_write_out_melt$id=="7-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="clr",],aes(x=mutation_label,y=value,fill=label,color=label))+
  theme(axis.text.x=element_blank())+geom_point(size=2)+ggtitle("7mers - random 100 subset")

# want to plot just vaquita and human and just plot the mutations that were 0 in both *originally* (so should be the same rate)
missingFromVaqAndHuman_7mers <- test[test$label %in% c("humans_AFR","vaquita") & test$total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS==0,]$mutation_label[duplicated(test[test$label %in% c("humans_AFR","vaquita") & test$total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS==0,]$mutation_label)]

# just plot those clr values of types that are missing in vaq and human so should be *identical* rates pre-clr
ggplot(all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$mutation_label %in% missingFromVaqAndHuman_7mers & all_processed_spectra_to_write_out_melt$label %in% c("vaquita","humans_AFR") & all_processed_spectra_to_write_out_melt$id=="7-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="clr" & all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",],aes(x=mutation_label,y=value,color=label))+
  geom_point(aes(shape=label))+
  scale_shape_manual(values=c(1,3))+
  theme(axis.text.x=element_blank())+
  ggtitle("7mers that were not observed in vaquita or human after downsampling\nAfter regularization (+1) should have identical mutation rates.\nBut due to differences in geometric means between vaq and human they are slightly different")

#plot NON CLR rates (UNTRANSFORMED)
ggplot(all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$mutation_label %in% missingFromVaqAndHuman_7mers & all_processed_spectra_to_write_out_melt$label %in% c("vaquita","humans_AFR") & all_processed_spectra_to_write_out_melt$id=="7-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="none" & all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",],aes(x=mutation_label,y=value,color=label))+
  geom_point(aes(shape=label))+
  scale_shape_manual(values=c(1,3))+
  theme(axis.text.x=element_blank())+
  ggtitle("UNTRANSFORMED 7mers that were not observed in vaquita or human after downsampling\nAfter regularization (+1) should have identical mutation rates.\nBut due to differences in geometric means between vaq and human they are slightly different")


## plot the CLR transformed rates that AREN"T 0 in both vaq/human
ggplot(all_processed_spectra_to_write_out_melt[!(all_processed_spectra_to_write_out_melt$mutation_label %in% missingFromVaqAndHuman_7mers) & all_processed_spectra_to_write_out_melt$label %in% c("vaquita","humans_AFR") & all_processed_spectra_to_write_out_melt$id=="7-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="clr" & all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",],aes(x=mutation_label,y=value,color=label))+
  geom_point(aes(shape=label))+
  scale_shape_manual(values=c(1,3))+
  theme(axis.text.x=element_blank())+
  ggtitle("7mers that were not observed in vaquita or human after downsampling\nAfter regularization (+1) should have identical mutation rates.\nBut due to differences in geometric means between vaq and human they are slightly different")



# try with 3mers
ggplot(all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_processed_spectra_to_write_out_melt$label %in% speciesToInclude_phylo & all_processed_spectra_to_write_out_melt$id=="3-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="clr",],aes(x=mutation_label,y=value,color=label))+
  theme(axis.text.x=element_blank())+geom_point()
# plot geometric means per species for 7mers:
# regularized
test$label <- factor(test$label,levels=c("wolves","brown_bear_EUR","brown_bear_ABC","polar_bear_PB","fin_whale_GOC","fin_whale_ENP","vaquita","humans_EUR"  ,"humans_EAS","humans_SAS","humans_AMR","humans_AFR" , "Pan_paniscus", "Pan_troglodytes" ,"Gorilla_gorilla"  , "Pongo_pygmaeus"  , "Pongo_abelii","Mus_musculus_Mmd","Mus_musculus_Mmm","Mus_musculus_Mmc","Mus_spretus_Ms" ))
test %>% group_by(label) %>% summarise(geometricmeanmutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=exp(mean(log(mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS)))) %>% 
  ggplot(.,aes(x=label,y=geometricmeanmutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS))+
  geom_col()+
  theme(text=element_text(size=14),axis.text.x=element_text(angle=90))+
  ggtitle("geometric means -- scaled to human targets prior to downsampling; regularized usign +1 to all types.")


test3mers$label <- factor(test3mers$label,levels=c("wolves","brown_bear_EUR","brown_bear_ABC","polar_bear_PB","fin_whale_GOC","fin_whale_ENP","vaquita","humans_EUR"  ,"humans_EAS","humans_SAS","humans_AMR","humans_AFR" , "Pan_paniscus", "Pan_troglodytes" ,"Gorilla_gorilla"  , "Pongo_pygmaeus"  , "Pongo_abelii","Mus_musculus_Mmd","Mus_musculus_Mmm","Mus_musculus_Mmc","Mus_spretus_Ms" ))
test3mers %>% group_by(label) %>% summarise(geometricmeanmutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=exp(mean(log(mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS)))) %>% 
  ggplot(.,aes(x=label,y=geometricmeanmutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS))+
  geom_col()+
  theme(text=element_text(size=14),axis.text.x=element_text(angle=90))+
  ggtitle("geometric means -- based on 3mers; scaled to human targets prior to downsampling; regularized usign +1 to all types.")



test1mers$label <- factor(test1mers$label,levels=c("wolves","brown_bear_EUR","brown_bear_ABC","polar_bear_PB","fin_whale_GOC","fin_whale_ENP","vaquita","humans_EUR"  ,"humans_EAS","humans_SAS","humans_AMR","humans_AFR" , "Pan_paniscus", "Pan_troglodytes" ,"Gorilla_gorilla"  , "Pongo_pygmaeus"  , "Pongo_abelii","Mus_musculus_Mmd","Mus_musculus_Mmm","Mus_musculus_Mmc","Mus_spretus_Ms" ))
test1mers %>% group_by(label) %>% summarise(geometricmeanmutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=exp(mean(log(mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS)))) %>% 
  ggplot(.,aes(x=label,y=geometricmeanmutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS))+
  geom_col()+
  theme(text=element_text(size=14),axis.text.x=element_text(angle=90))+
  ggtitle("geometric means -- based on 1mers; scaled to human targets prior to downsampling; regularized using +1 to all types.")


# note: geomean(x/y) = geomean(x)/geomean(y) apparenlty
# so let's see if that's what could be driving differences?
head(test)
geomeanFunction <- function(x){
  geomean=exp(mean(log(x)))
}

geomeans7mers <- test %>%
  group_by(label) %>%
  summarise(geoMeanTotalMutationsDownsampledPlusEpsilonHumanTargets=geomeanFunction(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS),geoMeanTargetSizeHumanScaled=geomeanFunction(total_target_count_HUMAN),geoMeanMutationRatedownsampledPlusEpsilsonHumanTargets=geomeanFunction(mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS),arithmeticMeanTotalMutationsDownsampledPlusEpsilonHumanTargets=mean(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS)) %>%
  ungroup()

geomeans1mers <- test1mers  %>%  
  group_by(label) %>%
  summarise(geoMeanTotalMutationsDownsampledPlusEpsilonHumanTargets=geomeanFunction(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS),geoMeanTargetSizeHumanScaled=geomeanFunction(total_target_count_HUMAN),geoMeanMutationRatedownsampledPlusEpsilsonHumanTargets=geomeanFunction(mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS),arithmeticMeanTotalMutationsDownsampledPlusEpsilonHumanTargets=mean(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS)) %>%
  ungroup()

geomeans3mers <- test3mers  %>%  
  group_by(label) %>%
  summarise(geoMeanTotalMutationsDownsampledPlusEpsilonHumanTargets=geomeanFunction(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS),geoMeanTargetSizeHumanScaled=geomeanFunction(total_target_count_HUMAN),geoMeanMutationRatedownsampledPlusEpsilsonHumanTargets=geomeanFunction(mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS),arithmeticMeanTotalMutationsDownsampledPlusEpsilonHumanTargets=mean(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS)) %>%
  ungroup()

speciesToTest=c("wolves","Mus_musculus_Mmd","vaquita","humans_AFR")
ggplot(test1mers[test1mers$label %in% speciesToTest,],aes(x=mutation_label,y=total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS,fill=label))+
  geom_col(position='dodge')+
  geom_hline(data=geomeans1mers[geomeans1mers$label %in% speciesToTest,],aes(yintercept=geoMeanTotalMutationsDownsampledPlusEpsilonHumanTargets,color=label),size=2)+
  geom_hline(data=geomeans1mers[geomeans1mers$label %in% speciesToTest,],aes(yintercept=arithmeticMeanTotalMutationsDownsampledPlusEpsilonHumanTargets,color=label),linetype="dashed",size=1)+
  
  ggtitle("Counts (downsampled after rescaling to human targets)\nsubset of species for illustration\ngeometric means shown as solid lines\narithmetic mean shown as dashed lines (all are equal)")


# plot that with 7mers but just a couple species
speciesToTest=c("vaquita","humans_AFR")

ggplot(test[test$label %in% speciesToTest,],aes(x=mutation_label,y=total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS,fill=label))+
  geom_col(position='dodge')+
  theme(axis.text.x = element_blank())+
  geom_label_repel(data=test[test$total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS>100 & test$label %in% speciesToTest,],aes(label=mutation_label))+
  geom_point(data=test[test$total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS>100 & test$label %in% speciesToTest,],aes(x=mutation_label,y=total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS,color=label))+
  ggtitle("plot showing unevenness in vaquita spectrum -- mutation counts")+
  geom_hline(data=geomeans7mers[geomeans7mers$label %in% speciesToTest,],aes(yintercept=geoMeanTotalMutationsDownsampledPlusEpsilonHumanTargets,color=label),size=1)+
  scale_y_log10()+
  ggtitle("Counts (downsampled after rescaling to human targets)\nsubset of species for illustration\ngeometric means shown as solid lines\narithmetic mean shown as dashed lines (all are equal)")


# try for gorilla and human
speciesToTest=c("Gorilla_gorilla","humans_AFR")

ggplot(test3mers[test3mers$label %in% speciesToTest,],aes(x=mutation_label,y=total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS,fill=label))+
  geom_col(position='dodge')+
  theme(axis.text.x = element_blank())+
  #geom_label_repel(data=test3mers[test3mers$total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS>1000 & test3mers$label %in% speciesToTest,],aes(label=mutation_label))+
  #geom_point(data=test3mers[test3mers$total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS>100 & test3mers$label %in% speciesToTest,],aes(x=mutation_label,y=total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS,color=label))+
  ggtitle("plot showing unevenness in vaquita spectrum -- mutation counts")+
  geom_hline(data=geomeans3mers[geomeans3mers$label %in% speciesToTest,],aes(yintercept=geoMeanTotalMutationsDownsampledPlusEpsilonHumanTargets,color=label),size=1)+
  #scale_y_log10()+
  ggtitle("Counts (downsampled after rescaling to human targets)\nsubset of species for illustration\ngeometric means shown as solid lines\narithmetic mean shown as dashed lines (all are equal)")

############# done with that troubleshooting section (could comment out) ############
########## CASE STUDIES want to see which mutation types are contributing most to diffs between GOR-HUMAN 3mers pre and post CLR #########
speciesToTest_humanGorilla = c("humans_AFR","Gorilla_gorilla")
# these are clr transformed 3mers of humans/gorilla 
casestudydf_clr <- all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_processed_spectra_to_write_out_melt$label %in% speciesToTest_humanGorilla & all_processed_spectra_to_write_out_melt$id=="3-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="clr",]

casestudydf_untransformed <- all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_processed_spectra_to_write_out_melt$label %in% speciesToTest_humanGorilla & all_processed_spectra_to_write_out_melt$id=="3-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="none",]

casestudydf_clr_wide <- pivot_wider(casestudydf_clr,id_cols = c(id,transform_id,mutation_label),names_from = label,values_from = value)


casestudydf_untransformed_wide <- pivot_wider(casestudydf_untransformed,id_cols = c(id,transform_id,mutation_label),names_from = label,values_from = value)

casestudydf_clr_wide$squaredDiff <- (casestudydf_clr_wide$humans_AFR-casestudydf_clr_wide$Gorilla_gorilla)^2
#sqrt(sum(casestudydf_clr_wide$squaredDiff)) # matches distance in dist plot, cool.

casestudydf_untransformed_wide$squaredDiff <- (casestudydf_untransformed_wide$humans_AFR-casestudydf_untransformed_wide$Gorilla_gorilla)^2
#sqrt(sum(casestudydf_untransformed_wide$squaredDiff)) # matches untransformed distance in plot, cool.  

casestudydf_clr_wide$squaredDiffFracOfTotal <- casestudydf_clr_wide$squaredDiff/sum(casestudydf_clr_wide$squaredDiff)
casestudydf_untransformed_wide$squaredDiffFracOfTotal <- casestudydf_untransformed_wide$squaredDiff/sum(casestudydf_untransformed_wide$squaredDiff)

casestudydf_combo <- bind_rows(casestudydf_clr_wide,casestudydf_untransformed_wide)

casestudydf_combo$niceDistanceLabel <- ""
casestudydf_combo[casestudydf_combo$transform_id=="clr",]$niceDistanceLabel <- "Aitchison Distance"
casestudydf_combo[casestudydf_combo$transform_id=="none",]$niceDistanceLabel <- "Euclidean Distance"

casestudyplot1 <- ggplot(casestudydf_combo,aes(x=mutation_label,y=squaredDiffFracOfTotal,fill=transform_id,color=transform_id))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),strip.text.x=element_text(size=14),legend.position = "none")+
  ylab("fraction of squared distance\n[squared difference/sum(squared differences)]")+
  #ggtitle('Contributions to distance between human and gorilla -- with and without CLR')+
  facet_wrap(~niceDistanceLabel,scales="free_y",ncol=1)
casestudyplot1
ggsave(paste0(plotdir,"casestudyplot1.humangorilladistancesper3mer.pdf"),casestudyplot1,height=7,width=16)
# let's see what mutation types contribute most to that distance

###### add cosine similarity to case study 1:
#### want to plot individual components that add to cosine distance ####
# human/gorilla casestudy
# https://en.wikipedia.org/wiki/Cosine_similarity 
# cosine similarity is calculated as ( sum(Ai*Bi) over all i ) / ( sqrt(sum(Ai^2) over all i )+sqrt(sum(Bi^2) over all i) )
# so could get components for each type as Ai*Bi / ( sqrt(sum(Ai^2) over all i )+sqrt(sum(Bi^2) over all i) 
# and then do 1 - that to get distance rather than similarity
# use allXmerspectra: 
getCosineSimilarityComponents_CaseStudy <- function(spectra,epsilon,speciesALabel,speciesBLabel){
  
  # get your two species:
  subsetSpectraForCaseStudy <- spectra %>%
    processSpectra(.,epsilon) %>%
    pivotSpectra_perVariable(.,"total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS") %>%
    select("mutation_label",speciesALabel,speciesBLabel)
  
  colnames(subsetSpectraForCaseStudy) <- c("mutation_label","speciesA","speciesB")
  
  # get the sqaure root of the sum of squared values for each
  spA_sqrt_sum_of_squares = sqrt(sum(subsetSpectraForCaseStudy$speciesA^2))
  spB_sqrt_sum_of_squares = sqrt(sum(subsetSpectraForCaseStudy$speciesB^2))
  
  product_spA_sqrt_sum_of_squares_spB_sqrt_sum_of_squares = spA_sqrt_sum_of_squares*spB_sqrt_sum_of_squares
  
  subsetSpectraForCaseStudy <- subsetSpectraForCaseStudy %>%
    mutate(Ai_x_Bi=speciesA*speciesB)
  
  subsetSpectraForCaseStudy$contributionToCosineSIMILARITY <- subsetSpectraForCaseStudy$Ai_x_Bi / product_spA_sqrt_sum_of_squares_spB_sqrt_sum_of_squares
  
  # note total cosine similarity is sum(subsetSpectraForCaseStudy$contributionToCosineSIMILARITY) # 0.9872027 *checked manually that it matches the lsa value
  
  # note that contributin to dissimilarity is 1 - (similarity)
  # need to distributte that 1 across all the mutation types (e.g 1 - (x+y+z) = (1/3-x + 1/3-y + 1/3-z))
  nmutationtypes=dim(subsetSpectraForCaseStudy)[1] # total mutation types
  
  subsetSpectraForCaseStudy <- subsetSpectraForCaseStudy %>%
    mutate(contributionToCosineDISSIMILARITY=(1/nmutationtypes)-contributionToCosineSIMILARITY) # this works!
  
  # sum(subsetSpectraForCaseStudy$contributionToCosineDISSIMILARITY) = 0.01279735 which matches the lsa value
  
  return(subsetSpectraForCaseStudy)
}

humanGorilla3merCosineDISTANCEContributions <- getCosineSimilarityComponents_CaseStudy(all3merSpectra,epsilon,"humans_AFR","Gorilla_gorilla")

humanGorilla1merCosineDISTANCEContributions <- getCosineSimilarityComponents_CaseStudy(all1merSpectra,epsilon,"humans_AFR","Gorilla_gorilla")


humanGorilla3merCosineDISTANCEContributions$distanceLabel <- "cosine similarity"
humanGorilla1merCosineDISTANCEContributions$distanceLabel <- "cosine similarity"

casestudyplot1_cosineSim <- ggplot(humanGorilla3merCosineDISTANCEContributions,aes(x=mutation_label,y=contributionToCosineSIMILARITY))+
  geom_col(fill="orange")+
  theme_bw()+
  facet_wrap(~distanceLabel)+
  theme(axis.text.x=element_text(angle=45,hjust=1),strip.text.x=element_text(size=14))+
  ylab("contribution to cosine similarity")#+
  #ggtitle('Contributions to cosine similarity between human and gorilla')
casestudyplot1_cosineSim
ggsave(paste0(plotdir,"casestudyplot1.humangorilladistancesper3mer.COSINESIMILARITY.pdf"),casestudyplot1_cosineSim,height=4,width=16)
# let's see what mutation types contribute most to that distance

# 1mer plot:
casestudyplot1_cosineSim_1mer <- ggplot(humanGorilla1merCosineDISTANCEContributions,aes(x=mutation_label,y=contributionToCosineSIMILARITY))+
  geom_col(fill="orange")+
  theme_bw()+
  facet_wrap(~distanceLabel)+
  theme(axis.text.x=element_text(angle=45,hjust=1),strip.text.x=element_text(size=14))+
  ylab("contribution to cosine similarity")#+
#ggtitle('Contributions to cosine similarity between human and gorilla')
casestudyplot1_cosineSim_1mer
ggsave(paste0(plotdir,"casestudyplot1.humangorilladistancesper1mer.COSINESIMILARITY.pdf"),casestudyplot1_cosineSim_1mer,height=4,width=8)


# show that most abundant types contribute most:
casestudyplot1_cosineSimAbudance <- ggplot(humanGorilla3merCosineDISTANCEContributions,aes(x=speciesA,y=contributionToCosineSIMILARITY))+
  geom_point()+
  theme_bw()+
  xlab("3-mer abundance in human (downsampled count)")
casestudyplot1_cosineSimAbudance
ggsave(paste0(plotdir,"casestudyplot1.humangorilladistancesper1mer.COSINESIMILARITYABUDNANCE.pdf"),casestudyplot1_cosineSimAbudance,height=4,width=8)

############## case study 2: vaquita and human

# now without CLR are those same ones contributing the most?
speciesToTest_VaquitaHuman <- c("vaquita","humans_AFR")
casestudydf2_clr <- all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_processed_spectra_to_write_out_melt$label %in% speciesToTest_VaquitaHuman & all_processed_spectra_to_write_out_melt$id=="7-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="clr",]

casestudydf2_untransformed <- all_processed_spectra_to_write_out_melt[all_processed_spectra_to_write_out_melt$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_processed_spectra_to_write_out_melt$label %in% speciesToTest_VaquitaHuman & all_processed_spectra_to_write_out_melt$id=="7-mer spectrum" & all_processed_spectra_to_write_out_melt$transform_id=="none",]


casestudydf2_clr_wide <- pivot_wider(casestudydf2_clr,id_cols = c(id,transform_id,mutation_label),names_from = label,values_from = value)

casestudydf2_untransformed_wide <- pivot_wider(casestudydf2_untransformed,id_cols = c(id,transform_id,mutation_label),names_from = label,values_from = value)

casestudydf2_clr_wide$squaredDiff <- (casestudydf2_clr_wide$humans_AFR-casestudydf2_clr_wide$vaquita)^2
sqrt(sum(casestudydf2_clr_wide$squaredDiff)) # matches distance in dist plot, cool.

casestudydf2_clr_wide$squaredDiffFracOfTotal <- casestudydf2_clr_wide$squaredDiff/sum(casestudydf2_clr_wide$squaredDiff)

casestudydf2_untransformed_wide$squaredDiff <- (casestudydf2_untransformed_wide$humans_AFR-casestudydf2_untransformed_wide$vaquita)^2
sqrt(sum(casestudydf2_untransformed_wide$squaredDiff)) # matches distance in dist plot, cool.

casestudydf2_untransformed_wide$squaredDiffFracOfTotal <- casestudydf2_untransformed_wide$squaredDiff/sum(casestudydf2_untransformed_wide$squaredDiff)



casestudydf2_clr_wide$colorLabel <- "not missing in both species"
casestudydf2_clr_wide[casestudydf2_clr_wide$mutation_label %in% missingFromVaqAndHuman_7mers,]$colorLabel <- "missing in both species"


casestudydf2_untransformed_wide$colorLabel <- "not missing in both species"
casestudydf2_untransformed_wide[casestudydf2_untransformed_wide$mutation_label %in% missingFromVaqAndHuman_7mers,]$colorLabel <- "missing in both species"

casestudydf2_combo <- bind_rows(casestudydf2_clr_wide,casestudydf2_untransformed_wide)

casestudy2plot1 <- ggplot(casestudydf2_combo,aes(x=mutation_label,y=squaredDiffFracOfTotal,fill=colorLabel,color=colorLabel))+
  geom_point()+
  theme(axis.text.x=element_blank())+
  ylab("fraction of squared distance\n[squared difference/sum(squared differences)]")+
  ggtitle('Contributions to distance between human and vaquita 7mer spectra')+
  facet_wrap(~transform_id,scales="free_y",ncol=1)
casestudy2plot1
ggsave(paste0(plotdir,"casestudy2plot1.humanvaquitadistancesper7mer.pdf"),casestudy2plot1,height=7,width=16)

#casestudydf2_combo[casestudydf2_combo$squaredDiffFracOfTotal>0.01,]
# let's see what mutation types contribute most to that distance
################ returning to main script: phylo plots #################
##### be careful that you aren't falsely combining things 
# adding sqrt cophenetic now 
phylo_plot1 <- ggplot(all_distances_all_conditions_phylo,aes(x=sqrt(cophenetic_distance),y=spectrum_distance,color=variable))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=1.8)+
  ggtitle(paste0("phylo distances\n",flagOfAnalysis))

phylo_plot1
ggsave(paste0(plotdir,"phylo_plot1_transforms.allconditions.sqrt.pdf"),height=12,width=24,phylo_plot1)

# plot each individual statistic:
### plot individual stats:
for(variable in listOfVars){
  phylo_plot1b <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variable,],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
    facet_wrap(~transform_label~id,scales="free",ncol=4)+
    geom_point()+
    theme_bw()+
    geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=2)+
    ggtitle(paste0(variable,"\n",flagOfAnalysis))
  ggsave(paste0(plotdir,"phylo_plot1b_transforms.",variable,".sqrt.pdf"),height=8,width=16,phylo_plot1b)
  
  # also without labels:
  phylo_plot1b <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variable,],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
    facet_wrap(~transform_label~id,scales="free",ncol=4)+
    geom_point()+
    theme_bw()+
    #geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=2)+
    ggtitle(paste0(variable,"\n",flagOfAnalysis))
  ggsave(paste0(plotdir,"phylo_plot1b_transforms.",variable,".sqrt.noabels.pdf"),height=8,width=16,phylo_plot1b)
  
  
}

####### plot both frac seg sites and mut rate with and without downsamp ######

pairsOfVars=list(mutrate=c("mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS"),fracSegSites=c("fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS"),comparingDownsampFrac_with_MutRate=c("mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS"),comparingNonDownsampFrac_with_MutRate=c("mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS"),comparingDownsampFrac_with_MutRat_and_Counts=c("mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS"),comparingNonDownsampFrac_with_MutRat_and_Counts=c("mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS"))
for(index in c("mutrate","fracSegSites","comparingDownsampFrac_with_MutRate","comparingNonDownsampFrac_with_MutRate","comparingDownsampFrac_with_MutRat_and_Counts","comparingNonDownsampFrac_with_MutRat_and_Counts")){
  variablepair=pairsOfVars[index]
  phylo_plot1c <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable %in% unlist(variablepair),],aes(x=sqrt(cophenetic_distance),y=spectrum_distance,color=variable))+
    facet_wrap(~transform_label~id,scales="free",ncol=4)+
    geom_point()+
    theme_bw()+
    geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=.8)+
    ggtitle(paste0(index,"\n",flagOfAnalysis))+
    geom_smooth(method="lm",se=F)
  ggsave(paste0(plotdir,"phylo_plot1c_transforms.comparingPairsOfVariables.",index,".sqrt.png"),height=10,width=24,phylo_plot1c)
}

### Make a plot for the manuscript that is just CLR 3-mers showing all 3 variables in one row (PCA will go above it) ####
# to show that it doesn't change:
subsetToShowVariableDoesntMatter <- all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$id=="3-mer spectrum" & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$variable %in% c("mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS"),]

# making simple labels for the SI plot 
subsetToShowVariableDoesntMatter$niceVariableLabel <- ""
subsetToShowVariableDoesntMatter[subsetToShowVariableDoesntMatter$variable=="total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",]$niceVariableLabel <- "mutation counts"
subsetToShowVariableDoesntMatter[subsetToShowVariableDoesntMatter$variable=="fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",]$niceVariableLabel <- "mutation proportions"
subsetToShowVariableDoesntMatter[subsetToShowVariableDoesntMatter$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",]$niceVariableLabel <- "mutation 'rates' (count/target size)"

subsetToShowVariableDoesntMatter$niceVariableLabel <- factor(subsetToShowVariableDoesntMatter$niceVariableLabel,levels=c("mutation counts","mutation proportions","mutation 'rates' (count/target size)"))

variablesDontMatterIfHumanRescaledPlot <- ggplot(subsetToShowVariableDoesntMatter,aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+  
  geom_point()+
  theme_bw()+
  facet_wrap(~niceVariableLabel,nrow=1)+
  theme(strip.background = element_rect(fill="lightblue"))+
  ggtitle("Comparing spectrum distances for 3-mer spectrum based on different ways of scaling the data")+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance")
variablesDontMatterIfHumanRescaledPlot                                                                                    
ggsave(paste0(plotdir,"Comparing3merDistancesWithDifferentVariables.ForSI.pdf"),variablesDontMatterIfHumanRescaledPlot,height=3,width=8)


############# time plots: time calibrated (should make mice an outlier ) ############
time_plot1 <- ggplot(all_distances_all_conditions_time,aes(x=cophenetic_distance,y=spectrum_distance,color=variable))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  xlab("cophenetic distance (time-calibrated tree)")+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=0.8)+
  ggtitle(paste0("time distances\n",flagOfAnalysis))

time_plot1
ggsave(paste0(plotdir,"time_plot1_transforms.allconditions.png"),height=6,width=16,time_plot1)

# plot each individual statistic:
time_plot1b <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",],aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  xlab("cophenetic distance (time-calibrated tree)")+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=.8,max.overlaps = 10)+
  ggtitle(paste0("time distances\nmutation_rate_plusEpsilon\n",flagOfAnalysis))+
  geom_smooth(method="lm")

time_plot1b
ggsave(paste0(plotdir,"time_plot1b_transforms.mutationrate.png"),height=6,width=16,time_plot1b)

# color by comparison type:
time_plot1c <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",],aes(x=cophenetic_distance,y=spectrum_distance,color=comparisonLabel_broad_alphabetical))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  xlab("cophenetic distance (time-calibrated tree)")+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=.8)+
  ggtitle(paste0("time distances\nmutation_rate_plusEpsilon\n",flagOfAnalysis))# +
  #geom_smooth(method="lm")

time_plot1c
ggsave(paste0(plotdir,"time_plot1c_transforms.mutationrate.ColoredByType.png"),height=6,width=16,time_plot1c)



# got rid of pca -- if you need it it's in older script versions. 




############ try to make nj tree based on spectrum distances instead of genetic dsitance #######
# from ape, need distance matrix 
require(ape)
#ape::nj()


combo_function_njTree <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,epsilon){
  distance_dataframe <- 
    processSpectra(spectrumdf,epsilon) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
    euclideanDistance_keepMatrix(.)
  
  return(distance_dataframe)
  
  
}
distanceMatsFromSpectrumDistancesForMakingTrees_mutRate <- lapply(listOfDFs,combo_function_njTree,variable="mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",speciesToInclude=speciesToInclude_phylo,epsilon=epsilon)
# I think something is wrong here. Why are Mmd and Ms so different in the distance matrix. 

distanceMatsFromSpectrumDistancesForMakingTrees_mutRate_multinomDownsamp <- lapply(listOfDFs,combo_function_njTree,variable="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",speciesToInclude=c(speciesToInclude_phylo),epsilon=epsilon)


trees_based_on_spectra <- lapply(distanceMatsFromSpectrumDistancesForMakingTrees_mutRate,ape::nj)
trees_based_on_spectra_multinomDownsamp <- lapply(distanceMatsFromSpectrumDistancesForMakingTrees_mutRate_multinomDownsamp,ape::nj)

treesBasedOnSpectraPlotDir=paste0(plotdir,"/treesBasedOnSpectrumDistance/")
dir.create(treesBasedOnSpectraPlotDir,showWarnings = F)
# make tree plots:
for(spectrum in names(trees_based_on_spectra)){
  png(paste0(treesBasedOnSpectraPlotDir,"nj.tree.basedon.distances.clr.mutationrate.",spectrum,".png"),height=8,width=8,units="in",res=300)
  plot.phylo(trees_based_on_spectra[[spectrum]],type="unrooted")
  title(paste0("neighbor joining tree based on clr distance (mut rate)\n",spectrum))
  dev.off()
  
}

for(spectrum in names(trees_based_on_spectra_multinomDownsamp)){
  png(paste0(treesBasedOnSpectraPlotDir,"nj.tree.basedon.distances.clr.mutationrate.multinomDownsampled.",spectrum,".png"),height=8,width=8,units="in",res=300)
  plot.phylo(trees_based_on_spectra_multinomDownsamp[[spectrum]],type="unrooted",)
  title(paste0("neighbor joining tree based on clr distance ( mut rate with multinom downsampling)\n",spectrum))
  dev.off()
  
}

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


mantel_test_results_time_list = lapply(listOfDFs,comboFunction_mantelTest,phylo_dist= timeTree_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
mantel_test_results_time_df = bind_rows(mantel_test_results_time_list,.id="id")
mantel_test_results_time_df$distance_metric <- "time_tree"

mantel_test_results_time_SQRTDISTANCE_list = lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE,phylo_dist= timeTree_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
mantel_test_results_time_SQRTDISTANCE_df = bind_rows(mantel_test_results_time_SQRTDISTANCE_list,.id="id")
mantel_test_results_time_SQRTDISTANCE_df$distance_metric <- "time_tree"


####### run on folded spectrum!!!!!!!!!! #####
# just doing it with sqrt of phylo time (skipping time-time and hamming distances etc.)
mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_list =  lapply(listOfFoldedSpectra,comboFunction_mantelTest_SQRTDISTANCE,phylo_dist= raxml_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_df = bind_rows(mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_list,.id="id")
mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_df$distance_metric <- "raxml_tree_SQRT_FOLDEDSPECTRUM"

mantel_results_all <- bind_rows(mantel_test_results_raxml_SQRTDISTANCE_df,mantel_test_results_time_SQRTDISTANCE_df,mantel_test_results_raxml_df,mantel_test_results_time_df,mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_df)
write.table(mantel_results_all,paste0(manteldir,"mantelTest.Results.includesFoldedSpectra.txt"),row.names=F,quote=F,sep="\t")
# did this work?



# you are here
############ plot mantel r values:  ############
# just raxml: 
mantel_r_plot_raxml <- 
  ggplot(mantel_test_results_raxml_df,aes(x=id,y=statistic,fill=distance_metric))+
  geom_col(position="dodge")+
  #facet_wrap(distance_label~transform_label,ncol=2)+
  theme_bw() +
  scale_fill_manual(values=c("#1F78B4"))+
  xlab("")+
  theme(text=element_text(size=14),legend.title = element_blank())+
  ggtitle("Mantel Test: r statistic (Pearson)")+
  ylab("correlation")
mantel_r_plot_raxml
ggsave(paste0(manteldir,"mantelTest.r.stat.PerSpectrum.raxmlTreeOnly.pdf"),mantel_r_plot_raxml,height=4,width=7)
# want to add mantel test information to plots: 

######## raxml and time :
mantel_r_plot_raxml_both <-   ggplot(mantel_results_all,aes(x=id,y=statistic,fill=distance_metric))+
  geom_col(position="dodge")+
  #facet_wrap(distance_label~transform_label,ncol=2)+
  theme_bw() +
  scale_fill_brewer(palette = "Paired",direction = -1)+
  xlab("")+
  theme(text=element_text(size=14),legend.title = element_blank())+
  ggtitle("Mantel Test: r statistic (Pearson)")+
  ylab("correlation")
mantel_r_plot_raxml_both
ggsave(paste0(manteldir,"mantelTest.r.stat.PerSpectrum.raxml.and.time.pdf"),mantel_r_plot_raxml_both,height=4,width=7)
# want to add mantel test information to plots: 




##################### plot correlations with mantel coefficients #############
####### PLOTTING BOTh SQRT AND NON SQRT DISTANCE ala Hardy and Pavoine. SQRT cophenetic is best!
# mantelTestCorrelationsPlot_phylo <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance))+
#   geom_point()+
#   facet_wrap(~id,scales="free")+
#   geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
#   ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
#   theme_bw()
# 
# mantelTestCorrelationsPlot_phylo
# ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.pdf"),mantelTestCorrelationsPlot_phylo,height=8,width=11)

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
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_Just1merSpectrumForMs <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$id=="1-mer spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id=="1-mer spectrum",],aes(x=0.4,y=Inf,vjust=1.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=6)+
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

###### mantel results with labels: plot 1/3 and 5/7mer separately: ######
mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_13merOnly <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$id %in% c("1-mer spectrum","3-mer spectrum") & all_distances_all_conditions_phylo$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(color="gray")+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id %in% c("1-mer spectrum","3-mer spectrum"),],aes(x=0.45,y=Inf,vjust=1.3,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  geom_text(aes(label=niceLabelForComparison),size=2)+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=16))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between\nspecies' mutation spectra")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_13merOnly
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.13merOnly.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_13merOnly,height=5,width=11)
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.13merOnly.png"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_13merOnly,height=5,width=11,dpi=300)


mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_57merOnly <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$id %in% c("5-mer spectrum","7-mer spectrum") & all_distances_all_conditions_phylo$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(color="gray")+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=0.45,y=Inf,vjust=1.3,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  geom_text(aes(label=niceLabelForComparison),size=2)+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=16))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between\nspecies' mutation spectra")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_57merOnly
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.57merOnly.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_57merOnly,height=5,width=11)
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.57merOnly.png"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_57merOnly,height=5,width=11,dpi=300)


####### plot jsut 1mer with labels for my manually labeling
mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr" &all_distances_all_conditions_phylo$id=="1-mer spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id =="1-mer spectrum",],aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=2)

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.Just1merToHelpMakeManualLabels.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels,height=8,width=11)


# compare with and without sqrt patterns:

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color="sqrt"))+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance,color="nonsqrt"))+
  facet_wrap(~id,scales="free")+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  xlab("sqrt or non-sqrt cophenetic distance")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.COMPARISON.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc,height=8,width=11)

########### plot mantel test based on folded spectra ##########
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA <- ggplot(all_distances_all_conditions_phylo_FOLDED[all_distances_all_conditions_phylo_FOLDED$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo_FOLDED$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_df,aes(x=0.5,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance, FOLDED spectra"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between folded mutation spectra")
# making strips not gray

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.FOLDEDSPECTRA.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA,height=8,width=11)


########## folded mantel separated by 1/3mer 5/7mer ##########
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA_13mer <- ggplot(all_distances_all_conditions_phylo_FOLDED[all_distances_all_conditions_phylo_FOLDED$id %in% c("1-mer spectrum","3-mer spectrum") & all_distances_all_conditions_phylo_FOLDED$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo_FOLDED$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1,alpha=0.85)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_df[mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_df$id %in% c("1-mer spectrum","3-mer spectrum"), ],aes(x=0.5,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance, FOLDED spectra"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between folded mutation spectra")
# making strips not gray




mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA_13mer
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.FOLDEDSPECTRA.13mer.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA_13mer,height=4,width=11)
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.FOLDEDSPECTRA.13mer.png"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA_13mer,height=4,width=11,dpi=300)


mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA_57mer <- ggplot(all_distances_all_conditions_phylo_FOLDED[all_distances_all_conditions_phylo_FOLDED$id %in% c("5-mer spectrum","7-mer spectrum") & all_distances_all_conditions_phylo_FOLDED$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo_FOLDED$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1,alpha=0.85)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_df[mantel_test_results_raxml_SQRTDISTANCE_FOLDEDSPECTRA_df$id %in% c("5-mer spectrum","7-mer spectrum"), ],aes(x=0.5,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance, FOLDED spectra"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between folded mutation spectra")
# making strips not gray




mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA_57mer
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.FOLDEDSPECTRA.57mer.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA_57mer,height=4,width=11)
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.FOLDEDSPECTRA.57mer.png"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_FOLDEDSPECTRA_57mer,height=4,width=11,dpi=300)


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

########### timetree sqrt distance mantel plot: separate 1/3 and 5/7mers #############
mantelTestCorrelationsPlot_time_SQRTDISTANCE_13mer <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$id %in% c("1-mer spectrum","3-mer spectrum") & all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (",permCount," permutations)")))+
  geom_text(data=mantel_test_results_time_SQRTDISTANCE_df[mantel_test_results_time_SQRTDISTANCE_df$id %in% c("1-mer spectrum","3-mer spectrum"),],aes(x=5,y=Inf,vjust=1.15,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_time_SQRTDISTANCE_df$permCount)," permutations"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
  xlab("square root of phylogenetic distance (shared branch length)\nbranch length represent millions of years bp")+
  ylab("Aitchison distance between mutation spectra")
mantelTestCorrelationsPlot_time_SQRTDISTANCE_13mer
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.SQRTDISTANCE.13mer.pdf"),mantelTestCorrelationsPlot_time_SQRTDISTANCE_13mer,height=5,width=11)
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.SQRTDISTANCE.13mer.png"),mantelTestCorrelationsPlot_time_SQRTDISTANCE_13mer,height=5,width=11,dpi=300)


mantelTestCorrelationsPlot_time_SQRTDISTANCE_57mer <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$id %in% c("5-mer spectrum","7-mer spectrum") & all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (",permCount," permutations)")))+
  geom_text(data=mantel_test_results_time_SQRTDISTANCE_df[mantel_test_results_time_SQRTDISTANCE_df$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=5,y=Inf,vjust=1.15,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_time_SQRTDISTANCE_df$permCount)," permutations"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
  xlab("square root of phylogenetic distance (shared branch length)\nbranch length represent millions of years bp")+
  ylab("Aitchison distance between mutation spectra")
mantelTestCorrelationsPlot_time_SQRTDISTANCE_57mer
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.SQRTDISTANCE.57mer.pdf"),mantelTestCorrelationsPlot_time_SQRTDISTANCE_57mer,height=5,width=11)
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.SQRTDISTANCE.57mer.png"),mantelTestCorrelationsPlot_time_SQRTDISTANCE_57mer,height=5,width=11,dpi=300)

########### mantel test with mouse labels ##########
mice=c("Mus_musculus_Mmd","Mus_spretus_Ms","Mus_musculus_Mmm","Mus_musculus_Mmc") # note phylo nly has Mmd and Ms
all_distances_all_conditions_time$mouseLabel <- "no mouse"
all_distances_all_conditions_time[all_distances_all_conditions_time$item1 %in% mice,]$mouseLabel <- "mouse"
all_distances_all_conditions_time[all_distances_all_conditions_time$item2 %in% mice,]$mouseLabel <- "mouse"

mantelTestCorrelationsPlot_time_mouselabel <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point(aes(color=mouseLabel,shape=mouseLabel))+
  facet_wrap(~id,scales="free")+
  #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (",permCount," permutations)")))+
  geom_text(data=mantel_test_results_time_df,aes(x=150,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
  theme_bw()+
  scale_color_manual(values=c("blue","black"))+
  scale_shape_manual(values=c(5,16))+
  theme(legend.position = "none")
mantelTestCorrelationsPlot_time_mouselabel 
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.mouseLabel.pdf"),mantelTestCorrelationsPlot_time_mouselabel,height=8,width=11)


mantelTestCorrelationsPlot_time_mouselabel_SQRTDISTANCE <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color=mouseLabel,shape=mouseLabel))+
  facet_wrap(~id,scales="free")+
  ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
  theme_bw()+
  scale_color_manual(values=c("blue","black"))+
  scale_shape_manual(values=c(5,16))+
  geom_text(data=mantel_test_results_time_SQRTDISTANCE_df,aes(x=12,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  theme(legend.position = "none")
mantelTestCorrelationsPlot_time_mouselabel_SQRTDISTANCE
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.mouseLabel.SQRTDISTANCE.pdf"),mantelTestCorrelationsPlot_time_mouselabel_SQRTDISTANCE,height=8,width=11)


####### make a heatmap ##########
# use the downsampled metric
pivot1merForHeatmap <- all1merSpectra %>%
  processSpectra(.,epsilon ) %>%
  pivotSpectra_perVariable(variableToUseInMantelTest)
pivot3merForHeatmap <- all3merSpectra %>%
  processSpectra(.,epsilon ) %>%
  pivotSpectra_perVariable(variableToUseInMantelTest)
pivot5merForHeatmap <- all5merSpectra %>%
  processSpectra(.,epsilon ) %>%
  pivotSpectra_perVariable(variableToUseInMantelTest)
pivot7merForHeatmap <- all7merSpectra %>%
  processSpectra(.,epsilon ) %>%
  pivotSpectra_perVariable(variableToUseInMantelTest)


#### 1mer heatmap:
pivot1merForHeatmap$wolves_over_humans <- pivot1merForHeatmap$wolves/pivot1merForHeatmap$humans_AFR
pivot1merForHeatmap$mmd_over_humans <- pivot1merForHeatmap$Mus_musculus_Mmd / pivot1merForHeatmap$humans_AFR
pivot1merForHeatmap$ms_over_humans <- pivot1merForHeatmap$Mus_spretus_Ms / pivot1merForHeatmap$humans_AFR
pivot1merForHeatmap$wolves_over_vaquita <- pivot1merForHeatmap$wolves / pivot1merForHeatmap$vaquita
pivot1merForHeatmap$wolves_over_mmd <- pivot1merForHeatmap$wolves / pivot1merForHeatmap$Mus_musculus_Mmd
pivot1merForHeatmap$wolves_over_ms <- pivot1merForHeatmap$wolves / pivot1merForHeatmap$Mus_spretus_Ms
pivot1merForHeatmap$vaquita_over_ms <- pivot1merForHeatmap$vaquita / pivot1merForHeatmap$Mus_spretus_Ms
pivot1merForHeatmap$vaquita_over_humans <- pivot1merForHeatmap$vaquita / pivot1merForHeatmap$humans_AFR

ggplot(melt(pivot1merForHeatmap[,c("mutation_label","vaquita_over_humans","wolves_over_humans","mmd_over_humans","ms_over_humans","wolves_over_vaquita","wolves_over_mmd","wolves_over_ms","vaquita_over_ms")]),aes(x=mutation_label,y=variable,fill=value))+
  geom_tile()+
  scale_fill_gradient2(high="red",mid="white",low="blue",midpoint = 1)+
  geom_text(aes(label=round(value,2)))
  
##### make a heatmap comparison of a few species
# ape/human
pivot3merForHeatmap$Pan_troglodytes_over_humans_AFR <- pivot3merForHeatmap$Pan_troglodytes / pivot3merForHeatmap$humans_AFR
# ape / mouse 
pivot3merForHeatmap$Pan_troglodytes_over_Mus_musculus_Mmd <- pivot3merForHeatmap$Pan_troglodytes / pivot3merForHeatmap$Mus_musculus_Mmd
# ape / whale
pivot3merForHeatmap$Pan_troglodytes_over_vaquita <- pivot3merForHeatmap$Pan_troglodytes / pivot3merForHeatmap$vaquita
#human/ mouse
pivot3merForHeatmap$humans_AFR_over_Mus_musculus_Mmd<- pivot3merForHeatmap$humans_AFR / pivot3merForHeatmap$Mus_musculus_Mmd
# vaquita/mouse
pivot3merForHeatmap$vaquita_over_Mus_musculus_Mmd<- pivot3merForHeatmap$vaquita / pivot3merForHeatmap$Mus_musculus_Mmd
# polar bear / mouse 
pivot3merForHeatmap$polar_bear_PB_over_Mus_musculus_Mmd<- pivot3merForHeatmap$polar_bear_PB / pivot3merForHeatmap$Mus_musculus_Mmd
# polar bear / bb
pivot3merForHeatmap$polar_bear_PB_over_brown_bear_ABC <- pivot3merForHeatmap$polar_bear_PB / pivot3merForHeatmap$brown_bear_ABC

# bears and wolves
pivot3merForHeatmap$wolves_over_brown_bear_ABC <- pivot3merForHeatmap$wolves / pivot3merForHeatmap$brown_bear_ABC

# try other mouse
pivot3merForHeatmap$Mus_musculus_Mmd_over_Mus_spretus_Ms <- pivot3merForHeatmap$Mus_musculus_Mmd / pivot3merForHeatmap$Mus_spretus_Ms


# try wolves:
pivot3merForHeatmap$wolves_over_miceMs <- pivot3merForHeatmap$wolves / pivot3merForHeatmap$Mus_spretus_Ms
pivot3merForHeatmap$Mus_musculus_Mmd_over_wolves <- pivot3merForHeatmap$Mus_musculus_Mmd / pivot3merForHeatmap$wolves
pivot3merForHeatmap$Mus_musculus_Mmm_over_wolves <- pivot3merForHeatmap$Mus_musculus_Mmm / pivot3merForHeatmap$wolves

# wolves human
pivot3merForHeatmap$wolves_over_humans <- pivot3merForHeatmap$wolves / pivot3merForHeatmap$humans_AFR
pivot3merForHeatmap$Mmd_over_humans <- pivot3merForHeatmap$Mus_musculus_Mmd / pivot3merForHeatmap$humans_AFR
pivot3merForHeatmap$brown_bear_ABC_over_humans <- pivot3merForHeatmap$brown_bear_ABC / pivot3merForHeatmap$humans_AFR

#comparisons=c("apes_Pan_troglodytes_over_humans_AFR","apes_Pan_troglodytes_over_Mus_musculus_Mmd","apes_Pan_troglodytes_over_vaquita","apes_humans_AFR_over_Mus_musculus_Mmd","vaquita_over_Mus_musculus_Mmd","polar_bear_PB_over_Mus_musculus_Mmd","polar_bear_PB_over_brown_bear_ABC","Mus_spretus_Ms_over_Mus_musculus_Mmd","Mus_spretus_Ms_over_wolves","Mus_musculus_Mmd_over_wolves")
comparisons = c("Mus_musculus_Mmd_over_Mus_spretus_Ms","wolves_over_miceMs","Mmd_over_humans","wolves_over_humans","brown_bear_ABC_over_humans")
# make heatmap:
pivot3merForHeatmap$centralBP <- paste0(substr(pivot3merForHeatmap$mutation_label,2,2),".",substr(pivot3merForHeatmap$mutation_label,6,6))

pivot3merForHeatmap$five_prime_flanking_base = substr(pivot3merForHeatmap$mutation_label,1,1)
pivot3merForHeatmap$three_prime_flanking_base = substr(pivot3merForHeatmap$mutation_label,3,3)

pivot3merForHeatmap_melt_subset <- melt(pivot3merForHeatmap[,c("mutation_label","centralBP","five_prime_flanking_base","three_prime_flanking_base",comparisons)])

pivot3merForHeatmap_melt_subset$nicerCentralBPLabel <- gsub("\\.",">",pivot3merForHeatmap_melt_subset$centralBP,)

heatmap_lotsOfComparisons <- ggplot(pivot3merForHeatmap_melt_subset,aes(x=three_prime_flanking_base,y=five_prime_flanking_base,fill=value))+
  geom_tile()+
  facet_grid(~nicerCentralBPLabel~variable,switch="y")+
  scale_fill_gradient2(high="red",mid="white",low="blue",midpoint = 1)+
  theme_bw()+
  geom_text(aes(label=mutation_label),size=1)
heatmap_lotsOfComparisons
ggsave(paste0(plotdir,"heatmap.LotsOfComparisons.pdf"),heatmap_lotsOfComparisons,height=10,width=16)

######## SEPARATE DISTANCES BY CENTRAL BP #####################
# as of 20220729 changed order of operations 
combo_function_separateByMutationType <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,speciesCodes,phyloDistances,centralMutationType,epsilon){
  distance_dataframe <- spectrumdf %>%
    filter(mutation_1mer==centralMutationType) %>% # filter first!!! 
    processSpectra(.,epsilon) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
    euclideanDistance(.) %>%
    addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
    addPhyloDistances_excludeNas(.,phyloDistances) %>%
    select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
    mutate(variable=variable,transform_label=ilr_or_clr_or_none,centralMutationType=centralMutationType)
  
  
  return(distance_dataframe)
  
  
}
######## mantel test faceted by central bp: ###############
combo_function_separateByMutationType_MantelTest <- function(spectrum, phylo_dist,variable,speciesToInclude,centralMutationType,mantelTestPermutationCount,epsilon,speciesCodes){
  spectrum_dist <- spectrum %>% 
    filter(mutation_1mer==centralMutationType) %>% # filter first~
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
  print(paste0("carrying out mantel test with",mantelTestPermutationCount," permutations"))
  
  mtr = vegan::mantel(xdis=phylo_dist_reordered,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),centralMutationType=centralMutationType)
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}

combo_function_separateByMutationType_MantelTest_SQRTDISTANCE <- function(spectrum, phylo_dist,variable,speciesToInclude,centralMutationType,mantelTestPermutationCount,epsilon,speciesCodes){
  spectrum_dist <- spectrum %>% 
    filter(mutation_1mer==centralMutationType) %>% # filter first!
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
  print(paste0("carrying out mantel test with",mantelTestPermutationCount," permutations"))
  phylo_dist_reordered_SQRTDISTANCE <- sqrt(phylo_dist_reordered)
  mtr = vegan::mantel(xdis=phylo_dist_reordered_SQRTDISTANCE,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel_SQRTDISTANCE",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),centralMutationType=centralMutationType)
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}


#centralMutationTypes=c("C.T","C.A","C.G","A.C","A.T","A.G")
# for now haven't been separating out CpGs. 
centralMutationTypes=unique(all1merSpectra$mutation_1mer) # forn ow doesn't sep CpGs out

##### separate 3mer, 5mer and 7mer distances by central bp: ##########
all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP = data.frame()
all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP = data.frame()
all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST = data.frame()
for(centralMutationType in centralMutationTypes){
  
  list_of_distance_dfs_phylo_SEPBYCENTRALBP <- lapply(listOfDFs[c("3-mer spectrum","5-mer spectrum","7-mer spectrum")]
,combo_function_separateByMutationType, variable =variableToUseInMantelTest, ilr_or_clr_or_none = "clr",speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances,centralMutationType=centralMutationType,epsilon=epsilon)
  
  distance_dfs_phylo_sepByCentralBP <- bind_rows(list_of_distance_dfs_phylo_SEPBYCENTRALBP, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
  
  all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP <- bind_rows(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,distance_dfs_phylo_sepByCentralBP)
  
  rm(list_of_distance_dfs_phylo_SEPBYCENTRALBP)
  
  # carry out mantel test:
  list_of_mantel_dfs_phylo_SEPBYCENTRALBP <- lapply(listOfDFs[c("3-mer spectrum","5-mer spectrum","7-mer spectrum")],combo_function_separateByMutationType_MantelTest, variable =variableToUseInMantelTest, speciesToInclude = speciesToInclude_phylo,phylo_dist = raxml_cophenetic_dist,centralMutationType=centralMutationType,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
  
  mantel_dfs_phylo_sepByCentralBP <- bind_rows(list_of_mantel_dfs_phylo_SEPBYCENTRALBP, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
  
  all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP <- bind_rows(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,mantel_dfs_phylo_sepByCentralBP)
  
  rm(list_of_mantel_dfs_phylo_SEPBYCENTRALBP)
  
  # carry out mantel test with sqrt phylo distance
  # NOTE: you don't need to supply sqrt distance here! the funciton L combo_function_separateByMutationType_MantelTest_SQRTDISTANCE will sqrt it for you! (it would double sqrt it if you pre-sqrt it, so don't do that!)
  list_of_mantel_dfs_phylo_SEPBYCENTRALBP_SQRTDIST <- lapply(listOfDFs[c("3-mer spectrum","5-mer spectrum","7-mer spectrum")],combo_function_separateByMutationType_MantelTest_SQRTDISTANCE, variable =variableToUseInMantelTest, speciesToInclude = speciesToInclude_phylo,phylo_dist = raxml_cophenetic_dist,centralMutationType=centralMutationType,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
  
  mantel_dfs_phylo_sepByCentralBP_SQRSTDIST <- bind_rows(list_of_mantel_dfs_phylo_SEPBYCENTRALBP_SQRTDIST, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
  
  all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST <- bind_rows(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST,mantel_dfs_phylo_sepByCentralBP_SQRSTDIST)
  
  rm(list_of_mantel_dfs_phylo_SEPBYCENTRALBP_SQRTDIST)
  
}

write.table(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST,paste0(manteldir,"mantelResults.facetedByCentralBP.",variableToUseInMantelTest,".includessqrtdist.txt"),row.names = F,quote=F,sep="\t")
write.table(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,paste0(manteldir,"mantelResults.facetedByCentralBP",variableToUseInMantelTest,".notsqrtdist.txt"),row.names = F,quote=F,sep="\t")


allDistances_facetedByCentralBP_plot <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  facet_grid(id~centralMutationType,scales="free_y")+
  theme_bw()+
  ggtitle(paste0(variableToUseInMantelTest," -- clr"))+
  geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=2)

allDistances_facetedByCentralBP_plot

ggsave(paste0(manteldir,"FacetedByCentralBP.",variableToUseInMantelTest,".clr.png"),allDistances_facetedByCentralBP_plot,height=8,width=16)

allDistances_facetedByCentralBP_plot_SQRTDIST <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_grid(id~centralMutationType,scales="free_y")+
  theme_bw()+
  ggtitle(paste0(variableToUseInMantelTest," -- clr; SQRT phylo dist"))+
  geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=2)


allDistances_facetedByCentralBP_plot_SQRTDIST

ggsave(paste0(manteldir,"FacetedByCentralBP.",variableToUseInMantelTest,".clr.SQRTDIST.png"),allDistances_facetedByCentralBP_plot_SQRTDIST,height=8,width=16)


########## Add the 3mer 5mer and 7mer all-sites plots to this plot ! 
# can I combine? 


allSitesSpectra_notFaceted_3mer_5mer_7mer_clr_variableToUseInMantelTest_NOTSEP <- all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$id %in% c("3-mer spectrum","5-mer spectrum","7-mer spectrum") & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$variable==variableToUseInMantelTest,]

allSitesSpectra_notFaceted_3mer_5mer_7mer_clr_variableToUseInMantelTest_NOTSEP$centralMutationType <- "all mutation types"

# let's combine htese
all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES <- bind_rows(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,allSitesSpectra_notFaceted_3mer_5mer_7mer_clr_variableToUseInMantelTest_NOTSEP)


# combine mantel results for sep by central bp and all sites: 
all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults <- bind_rows(mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id!="1-mer spectrum",], all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST)
# make sure these are both for sqrt distance
all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults[is.na(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType),]$centralMutationType <- "all mutation types"
########### Plot sep by central bp *WITH* the all sites versions so they are comparable: ######
# set order"

# replace . with > 
all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$centralMutationType <- gsub("\\.",">",all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$centralMutationType)

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType <- gsub("\\.",">",all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType)
# order facets:

all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$centralMutationType <- factor(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$centralMutationType,levels=c("all mutation types","A>C","A>G","A>T","C>A","C>G","C>T"))

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType <- factor(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType,levels=c("all mutation types","A>C","A>G","A>T","C>A","C>G","C>T"))

######## make and combine plots showing all sites and faceted by central bp ################
# exclude all-sites:
facetedByCentralBpPlot_niceForMs <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES[all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$centralMutationType!="all mutation types",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_grid(id~centralMutationType,scales="free_y",switch="y")+
  theme_bw()+
  geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults[all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType!="all mutation types",],aes(x=0.5,y=Inf,label=paste0("r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)))),size=4,vjust=1.1)+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=14))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance")

facetedByCentralBpPlot_niceForMs


JUSTPLOTTINGALLSITES_357_tobecombinedplot <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES[all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$centralMutationType=="all mutation types",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_grid(id~centralMutationType,scales="free_y",switch = "y")+ # switch puts on left side of plot
  theme_bw()+
  geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults[all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType=="all mutation types",],aes(x=0.5,y=Inf,label=paste0("r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)))),size=5,vjust=1.5)+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=14))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance")

JUSTPLOTTINGALLSITES_357_tobecombinedplot



all_sites_plus_facet_plot <- ggarrange(JUSTPLOTTINGALLSITES_357_tobecombinedplot,facetedByCentralBpPlot_niceForMs,nrow=1,widths=c(0.25,0.9))
all_sites_plus_facet_plot

ggsave(paste0(manteldir,"AllSites.Plus.FacetedByCentralBP.",variableToUseInMantelTest,".clr.SQRTDIST.pdf"),all_sites_plus_facet_plot,height=7,width=16)


############ separate 1/3mer from 5/7mer plots for ms ##################

# want to add an "all mutation types" flag so that they match
# just 1/3mers: 
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_13merONLY <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$id %in% c("1-mer spectrum" ,"3-mer spectrum"),],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id %in% c("1-mer spectrum" ,"3-mer spectrum"),],aes(x=0.5,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=14))+
  xlab("square root of phylogenetic distance (shared branch length)\nbranch lengths represent expected substitutions per site")+
  ylab("Aitchison distance between mutation spectra")
# making strips not gray

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_13merONLY
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.1-3merONLY.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_13merONLY,height=5,width=11)

# faceted by central bp: 3mer only:
facetedByCentralBpPlot_niceForMs_3merONLY <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES[all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$centralMutationType!="all mutation types" & all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$id=="3-mer spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1)+
  facet_grid(id~centralMutationType,scales="free_y",switch="y")+
  theme_bw()+
  geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults[all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType!="all mutation types" & all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$id=="3-mer spectrum",],aes(x=0.5,y=Inf,label=paste0("r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)))),size=4,vjust=1.1)+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=14))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance")

facetedByCentralBpPlot_niceForMs_3merONLY
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.FACETEDBYCENTRALBP.3merONLY.pdf"),facetedByCentralBpPlot_niceForMs_3merONLY,height=3,width=11)

###### arrange 1/3mer plot for main text ########

all_sites_plus_facet_plot_13merONLY <- ggarrange(mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_13merONLY,facetedByCentralBpPlot_niceForMs_3merONLY,nrow=2,heights = c(0.6,0.4))
all_sites_plus_facet_plot_13merONLY

ggsave(paste0(manteldir,"AllSites.Plus.FacetedByCentralBP.",variableToUseInMantelTest,".clr.SQRTDIST.1-3merONLY.pdf"),all_sites_plus_facet_plot_13merONLY,height=8,width=12)



# just 5/7mers:
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_57merONLY <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$id %in% c("5-mer spectrum" ,"7-mer spectrum"),],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1)+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id %in% c("5-mer spectrum" ,"7-mer spectrum"),],aes(x=0.5,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  #ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=14))+
  xlab("square root of phylogenetic distance (shared branch length)\nbranch lengths represent expected substitutions per site")+
  ylab("Aitchison distance between mutation spectra")
# making strips not gray

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_57merONLY
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.5-7merONLY.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_57merONLY,height=5,width=11)

# faceted by central bp: 5-7mer only:
facetedByCentralBpPlot_niceForMs_57merONLY <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES[all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$centralMutationType!="all mutation types" & all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PLUSALLSITES$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(size=1)+
  facet_grid(id~centralMutationType,scales="free_y",switch="y")+
  theme_bw()+
  geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults[all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$centralMutationType!="all mutation types" & all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST_PlusAllSitesResults$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=0.5,y=Inf,label=paste0("r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2)))),size=4,vjust=1.1)+
  theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=14))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance")

facetedByCentralBpPlot_niceForMs_57merONLY
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.FACETEDBYCENTRALBP.5-7merONLY.pdf"),facetedByCentralBpPlot_niceForMs_57merONLY,height=6,width=11)

###### arrange 5/7mer plot for main text ########
all_sites_plus_facet_plot_57merONLY <- ggarrange(mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_57merONLY,facetedByCentralBpPlot_niceForMs_57merONLY,nrow=2,heights = c(0.43,0.42))
all_sites_plus_facet_plot_57merONLY

ggsave(paste0(manteldir,"AllSites.Plus.FacetedByCentralBP.",variableToUseInMantelTest,".clr.SQRTDIST.5-7merONLY.pdf"),all_sites_plus_facet_plot_57merONLY,height=11,width=12)



###### plot highlighting the wolf comparisons ############
all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP$wolf_mouse_label <- "non"
all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP[all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP$comparisonLabel_broad_alphabetical =="canid.mouse",]$wolf_mouse_label <- "wolf-mouse"

allDistances_facetedByCentralBP_plot_highlightwolfmouse <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point(aes(color=wolf_mouse_label))+
  facet_grid(id~centralMutationType,scales="free_y")+
  theme_bw()+
  ggtitle(paste0(variableToUseInMantelTest,"-- clr"))+
  geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=2)


allDistances_facetedByCentralBP_plot_highlightwolfmouse

ggsave(paste0(manteldir,"FacetedByCentralBP.",variableToUseInMantelTest,".clr.COLORWOLF_MOUSE.png"),allDistances_facetedByCentralBP_plot_highlightwolfmouse,height=8,width=16)

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP$negLog10Pvalue <- -log10(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP$signif)

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST$negLog10Pvalue <- -log10(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST$signif)

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PlotPearsonR <- ggplot(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,aes(x=centralMutationType,y=statistic,fill=id))+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle(paste0("correlation between spectrum distance",variableToUseInMantelTest,"; clr) and cophenetic distance "))+
  ylab("Pearson's r")
all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PlotPearsonR

ggsave(paste0(manteldir,"FacetedByCentralBP.",variableToUseInMantelTest,".clr.mantelTestPearsonValues.png"),all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PlotPearsonR,height=8,width=16)

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP$negLog10Pvalue <- -log10(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP$signif)

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PlotPearsonR_SQRTDIST <- ggplot(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST,aes(x=centralMutationType,y=statistic,fill=id))+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle(paste0("correlation between spectrum distance (",variableToUseInMantelTest,"; clr) and cophenetic distance "))+
  ylab("Pearson's r")
all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PlotPearsonR_SQRTDIST

ggsave(paste0(manteldir,"FacetedByCentralBP.",variableToUseInMantelTest,".clr.mantelTestPearsonValues.SQRTDIST.png"),all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_PlotPearsonR_SQRTDIST,height=8,width=16)



all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_Plotlog10pvalue <- ggplot(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP,aes(x=centralMutationType,y=negLog10Pvalue,fill=id))+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle(paste0("correlation between spectrum distance (",variableToUseInMantelTest,"; clr) and cophenetic distance\ndashed line is 0.05/6 to correct for multiple testing (6 mutation types)"))+
  ylab("-log10 p-value")+
  geom_hline(yintercept = -log10(0.05/6),linetype="dashed")

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_Plotlog10pvalue

ggsave(paste0(manteldir,"FacetedByCentralBP.",variableToUseInMantelTest,".clr.mantelTestPValues.png"),all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_Plotlog10pvalue,height=8,width=16)


all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_Plotlog10pvalue_SQRTDIST <- ggplot(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST,aes(x=centralMutationType,y=negLog10Pvalue,fill=id))+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle(paste0("correlation between spectrum distance (",variableToUseInMantelTest,"; clr) and cophenetic distance\ndashed line is 0.05/6 to correct for multiple testing (6 mutation types)"))+
  ylab("-log10 p-value")+
  geom_hline(yintercept = -log10(0.05/6),linetype="dashed")

all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_Plotlog10pvalue_SQRTDIST

ggsave(paste0(manteldir,"FacetedByCentralBP.",variableToUseInMantelTest,".clr.mantelTestPValues.SQRTDIST.png"),all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_Plotlog10pvalue_SQRTDIST,height=8,width=16)


# ######### BGC categories ##########
# BGC_WS = c("A.C","A.G")
# BGC_conserved = c("A.T","C.G")
# BGC_SW = c("C.A","C.T")
# BGC_Categories = list(BGC_WS=BGC_WS,BGC_conserved=BGC_conserved,BGC_SW=BGC_SW)
# 
# combo_function_separateByBGCCategory <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,speciesCodes,phyloDistances,BGC_Categories,BGC_Category,epsilon){
#   BGC_category_mutationTypes=BGC_Categories[[BGC_Category]]
#   # 20220729 changing order of operations! need to filter out types first 
#   distance_dataframe <- spectrumdf %>%
#     filter(mutation_1mer %in% BGC_category_mutationTypes) %>% # changed order of ops
#     processSpectra(.,epsilon) %>%
#     pivotSpectra_perVariable(.,variable)  %>%
#     select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
#     clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
#     euclideanDistance(.) %>%
#     addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
#     addPhyloDistances_excludeNas(.,phyloDistances) %>%
#     select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
#     mutate(variable=variable,transform_label=ilr_or_clr_or_none,BGC_Category=BGC_Category)
#   
#   
#   return(distance_dataframe)
#   
#   
# }
# 
# ######## mantel test faceted by central bp: ###############
# ### NEED TO CHANGE ORDER OF OPERATIONS !!! 
# # as of 20220729 changed order of operations so fitering comes first !
# combo_function_separateByBGCCategory_MantelTest <- function(spectrum, phylo_dist,variable,speciesToInclude,BGC_Categories,BGC_Category,mantelTestPermutationCount,epsilon,speciesCodes){
#   BGC_category_mutationTypes=BGC_Categories[[BGC_Category]]
#   
#   spectrum_dist <- spectrum %>% 
#     filter(mutation_1mer %in% BGC_category_mutationTypes) %>% ## filter first!!! 
#     processSpectra(.,epsilon) %>% 
#     merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
#     filter(label %in% speciesToInclude) %>%
#     pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
#     clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
#     euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
#   # the matrices need to be ordered in the same way 
#   # reorder the time matrix in the same way as the spectrum matrix:
#   phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
#   head(phylo_dist_reordered)
#   print(paste0("carrying out mantel test with",mantelTestPermutationCount," permutations"))
#   
#   mtr = vegan::mantel(xdis=phylo_dist_reordered,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
#   # make a data frame
#   mantelTestResult_df <- data.frame(call="vegan_mantel",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),centralMutationType=centralMutationType,BGC_Category=BGC_Category)
#   # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
#   
#   return(mantelTestResult_df)
#   
# }
# 
# combo_function_separateByBGCCategory_MantelTest_SQRTDISTANCE <- function(spectrum, phylo_dist,variable,speciesToInclude,BGC_Categories,BGC_Category,mantelTestPermutationCount,epsilon,speciesCodes){
#   BGC_category_mutationTypes=BGC_Categories[[BGC_Category]]
#   
#   spectrum_dist <- spectrum %>% 
#     filter(mutation_1mer %in% BGC_category_mutationTypes) %>% # filter first!!
#     processSpectra(.,epsilon) %>% 
#     merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
#     filter(label %in% speciesToInclude) %>%
#     pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
#     clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
#     euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
#   # the matrices need to be ordered in the same way 
#   # reorder the time matrix in the same way as the spectrum matrix:
#   phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
#   head(phylo_dist_reordered)
#   phylo_dist_reordered_SQRT <- sqrt(phylo_dist_reordered)
#   print(paste0("carrying out mantel test with",mantelTestPermutationCount," permutations"))
#   
#   mtr = vegan::mantel(xdis=phylo_dist_reordered_SQRT,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
#   # make a data frame
#   mantelTestResult_df <- data.frame(call="vegan_mantel_SQRTDIST",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),centralMutationType=centralMutationType,BGC_Category=BGC_Category)
#   # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
#   
#   return(mantelTestResult_df)
#   
# }
# 
# ##### separate 3mer, 5mer and 7mer distances by bgc category: ##########
# all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC = data.frame()
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC = data.frame()
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_SQRTDIST = data.frame()
# for(BGC_Category in c("BGC_WS","BGC_conserved","BGC_SW")){
#   print(BGC_Category)  
#   list_of_distance_dfs_phylo_SEPBYBGC<- lapply(listOfDFs[c("3-mer spectrum","5-mer spectrum","7-mer spectrum")],combo_function_separateByBGCCategory, variable =variableToUseInMantelTest, ilr_or_clr_or_none = "clr",speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances,BGC_Categories,BGC_Category,epsilon=epsilon)
#   
#   distance_dfs_phylo_sepByBGC <- bind_rows(list_of_distance_dfs_phylo_SEPBYBGC, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
#   
#   all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC <- bind_rows(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC,distance_dfs_phylo_sepByBGC)
#   
#   rm(list_of_distance_dfs_phylo_SEPBYBGC)
#   
#   # carry out mantel test:
#   list_of_mantel_dfs_phylo_SEPBYBGC <- lapply(listOfDFs[c("3-mer spectrum","5-mer spectrum","7-mer spectrum")],combo_function_separateByBGCCategory_MantelTest, variable =variableToUseInMantelTest, speciesToInclude = speciesToInclude_phylo,phylo_dist = raxml_cophenetic_dist,BGC_Categories=BGC_Categories,BGC_Category=BGC_Category,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
#   
#   mantel_dfs_phylo_sepByBGC <- bind_rows(list_of_mantel_dfs_phylo_SEPBYBGC, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
#   
#   all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC <- bind_rows(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC,mantel_dfs_phylo_sepByBGC)
#   
#   rm(list_of_mantel_dfs_phylo_SEPBYBGC)
#   
#   # mantel test with sqrt(phylo):
#   # note don't give this funciton sqrt of distances -- it will do it interally 
#   list_of_mantel_dfs_phylo_SEPBYBGC_SQRTDIST <- lapply(listOfDFs[c("3-mer spectrum","5-mer spectrum","7-mer spectrum")],combo_function_separateByBGCCategory_MantelTest_SQRTDISTANCE, variable =variableToUseInMantelTest, speciesToInclude = speciesToInclude_phylo,phylo_dist = raxml_cophenetic_dist,BGC_Categories=BGC_Categories,BGC_Category=BGC_Category,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
#   
#   mantel_dfs_phylo_sepByBGC_SQRTDIST <- bind_rows(list_of_mantel_dfs_phylo_SEPBYBGC_SQRTDIST, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
#   
#   all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_SQRTDIST <- bind_rows(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_SQRTDIST,mantel_dfs_phylo_sepByBGC_SQRTDIST)
#   
#   rm(list_of_mantel_dfs_phylo_SEPBYBGC_SQRTDIST)
#   
#   
#   
# }
# 
# write.table(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC,paste0(manteldir,"mantelresults.sepbybgc.notsqrtdistance.txt"),row.names=F,quote=F,sep="\t")
# write.table(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_SQRTDIST,paste0(manteldir,"mantelresults.sepbybgc.sqrtdistance.txt"),row.names=F,quote=F,sep="\t")
# 
# allDistances_facetedByBGC <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC,aes(x=cophenetic_distance,y=spectrum_distance))+
#   geom_point()+
#   facet_grid(id~BGC_Category,scales="free_y")+
#   theme_bw()+
#   ggtitle(paste0(variableToUseInMantelTest," -- clr"))+
#   geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=2)
# 
# allDistances_facetedByBGC
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.png"),allDistances_facetedByBGC,height=8,width=12)
# 
# 
# # note that you plot the same distance df and just take sqrt of cophenetic for plotting; but need to use mantel results from sqrt dist (different dataframe)
# allDistances_facetedByBGC_SQRTDIST <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC,aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
#   geom_point()+
#   facet_grid(id~BGC_Category,scales="free_y",switch="y")+
#   theme_bw()+
#   #ggtitle("mutation_rate_multinom_downsampled_normalized_plusEpsilon -- clr; SQRT dist")+
#   geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC,aes(x=0.45,y=Inf,vjust=1.3,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=4)+
#   theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
#   xlab("sqrt of phylogenetic distance")+
#   ylab("Aitchison distance between species' mutation spectra")
# 
# 
# allDistances_facetedByBGC_SQRTDIST # making nicer 
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.sqrtdist.pdf"),allDistances_facetedByBGC_SQRTDIST,height=8,width=12)
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.sqrtdist.png"),allDistances_facetedByBGC_SQRTDIST,height=8,width=12,dpi=300)
# 
# ############# mantel test faceted by bgc: separate 1/3mer 5/7mer ########
# # there is no 1mer because these are faceted at higher dimensional level 
# allDistances_facetedByBGC_SQRTDIST_3mer <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC[all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC$id=="3-mer spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
#   geom_point()+
#   facet_grid(id~BGC_Category,scales="free_y",switch="y")+
#   theme_bw()+
#   #ggtitle("mutation_rate_multinom_downsampled_normalized_plusEpsilon -- clr; SQRT dist")+
#   geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC[all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC$id=="3-mer spectrum",],aes(x=0.45,y=Inf,vjust=1.3,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=4)+
#   theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
#   xlab("sqrt of phylogenetic distance")+
#   ylab("Aitchison distance between species' mutation spectra")
# 
# 
# allDistances_facetedByBGC_SQRTDIST_3mer
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.sqrtdist.3meronly.pdf"),allDistances_facetedByBGC_SQRTDIST_3mer,height=4,width=12)
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.sqrtdist.3meronly.png"),allDistances_facetedByBGC_SQRTDIST_3mer,height=4,width=12,dpi=300)
# 
# 
# allDistances_facetedByBGC_SQRTDIST_57mer <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC[all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
#   geom_point()+
#   facet_grid(id~BGC_Category,scales="free_y",switch="y")+
#   theme_bw()+
#   #ggtitle("mutation_rate_multinom_downsampled_normalized_plusEpsilon -- clr; SQRT dist")+
#   geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC[all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=0.45,y=Inf,vjust=1.3,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=4)+
#   theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
#   xlab("sqrt of phylogenetic distance")+
#   ylab("Aitchison distance between species' mutation spectra")
# 
# 
# allDistances_facetedByBGC_SQRTDIST_57mer
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.sqrtdist.57meronly.pdf"),allDistances_facetedByBGC_SQRTDIST_57mer,height=7,width=12)
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.sqrtdist.57meronly.png"),allDistances_facetedByBGC_SQRTDIST_57mer,height=7,width=12,dpi=300)
# 
# ######## plot faceted by bgc, labeling wolf-mouse ########
# all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC$wolf_mouse_label <- "non"
# all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC[all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC$comparisonLabel_broad_alphabetical =="canid.mouse",]$wolf_mouse_label <- "wolf-mouse"
# 
# 
# allDistances_facetedByBGC_labelwolfmouse <- ggplot(all_distances_phylo_variableToUseInMantelTest_clr_SEPBYBGC,aes(x=cophenetic_distance,y=spectrum_distance))+
#   geom_point(aes(color=wolf_mouse_label))+
#   facet_grid(id~BGC_Category,scales="free_y")+
#   theme_bw()+
#   ggtitle(paste0(variableToUseInMantelTest," -- clr"))+
#   geom_text(data=all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=2)
# 
# allDistances_facetedByBGC_labelwolfmouse
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.LABELWOLFMOUSE.png"),allDistances_facetedByBGC_labelwolfmouse,height=8,width=12)
# 
# 
# 
# 
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC$negLog10Pvalue <- -log10(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC$signif)
# 
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_SQRTDIST$negLog10Pvalue <- -log10(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_SQRTDIST$signif)
# 
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_PlotPearsonR <- ggplot(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC,aes(x=BGC_Category,y=statistic,fill=id))+
#   geom_col(position="dodge")+
#   theme_bw()+
#   ggtitle(paste0("correlation between spectrum distance (",variableToUseInMantelTest,"; clr) and cophenetic distance "))+
#   ylab("Pearson's r")
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_PlotPearsonR
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.mantelTestPearsonValues.png"),all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_PlotPearsonR,height=5,width=8)
# 
# 
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_PlotPearsonR_SQRTDIST <- ggplot(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_SQRTDIST,aes(x=BGC_Category,y=statistic,fill=id))+
#   geom_col(position="dodge")+
#   theme_bw()+
#   ggtitle(paste0("correlation between spectrum distance (",variableToUseInMantelTest,"; clr) and cophenetic distance "))+
#   ylab("Pearson's r")
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_PlotPearsonR_SQRTDIST
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.mantelTestPearsonValues.sqrtdist.png"),all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_PlotPearsonR_SQRTDIST,height=5,width=8)
# 
# 
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_Plotlog10pvalue <- ggplot(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC,aes(x=BGC_Category,y=negLog10Pvalue,fill=id))+
#   geom_col(position="dodge")+
#   theme_bw()+
#   ggtitle(paste0("correlation between spectrum distance (",variableToUseInMantelTest,"; clr) and cophenetic distance\ndashed line is 0.05/9 to correct for multiple testing (3 mutation types * 3 types of 3mer)"))+
#   ylab("-log10 p-value")+
#   geom_hline(yintercept = -log10(0.05/(3*3)),linetype="dashed")
# 
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_Plotlog10pvalue
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.mantelTestPValues.png"),all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_Plotlog10pvalue,height=5,width=8)
# 
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_Plotlog10pvalue_SQRTDIST <- ggplot(all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_SQRTDIST,aes(x=BGC_Category,y=negLog10Pvalue,fill=id))+
#   geom_col(position="dodge")+
#   theme_bw()+
#   ggtitle(paste0("correlation between spectrum distance (",variableToUseInMantelTest,"; clr) and cophenetic distance\ndashed line is 0.05/9 to correct for multiple testing (3 mutation types * 3 types of 3mer)"))+
#   ylab("-log10 p-value")+
#   geom_hline(yintercept = -log10(0.05/(3*3)),linetype="dashed")
# 
# all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_Plotlog10pvalue_SQRTDIST
# 
# ggsave(paste0(manteldir,"FacetedByBGC.",variableToUseInMantelTest,".clr.mantelTestPValues.sqrtdist.png"),all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYBGC_Plotlog10pvalue_SQRTDIST,height=5,width=8)

##################### generating fake control spectra #####################
controldatadir=paste0(plotdir,"controlDataMultinomialDownsampledWithinCentralKmer/")
dir.create(controldatadir,showWarnings = F)
# goal: create randomized (fake) control spectra which apportion total mutations in a 3mer mutation cateogry across all 5mer mutation types proportionally by target size
# then do mantel test 
# if it fits equally well as 5mer mantel test (real) then 5mers aren't adding anything beyond 3mers. If however the the control spectrum has less phylo signal then 5mers are adding something.
# see if control does as well as real -- in which case 5mers don't add much beyond 3mers!
# will also do for 7mers 
# want to do this on non-downsampled human-scaled spectrum but then process it (downsample,regularize) as before using the 'fake' counts
# order of operations: we want to do the mulitnomial sampling AFTER human rescaling so that genome content doens't get baked into the resampling at all. 
# order of operations: start with HUMAN RESCALED spectra: all5merSpectra_HumanRescaled
# then reshuffle the counts 
# then downsample (!)
# so going to use all5merSpectra_HumanRescaled and all7merSpectra_HumanRescaled which are human-rescaled but not yet downsampled
################## fake 5mer spectrum ##############################
makeFakeControl5merSpectrumFunction_humanRescaled <- function(spectradf_5mer){
  ##### spectra should have target info and human info in it #####
  # want to use all5merSpectra_plusHumanInfo #
  if(nchar(spectradf_5mer[1,]$mutation_label)!=11){
    stop("this isn't a 5mer spectrum!")
    
  }
  # get central 3mer 
  spectradf_5mer$central3mer <- paste0(substr(spectradf_5mer$mutation_5mer,2,4),".",substr(spectradf_5mer$mutation_5mer,8,10))
  # group by central 3mer and get target fractions of 5mers within that 3mer type
  # UPDATE: doing this with HUMAN TARGETS because want to shuffle HUMAN-RESCALED counts!!
  spectradf_5mer <- spectradf_5mer %>% 
    group_by(species,population,label,central3mer) %>%
    mutate(targetFractionOf5merTypeWithinCentral3mer_HUMAN = total_target_count_HUMAN/sum(total_target_count_HUMAN),total_mutations_InOverallCentral3merCategory_AfterHumanRescaling=sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    ungroup()
  
  # also need to get total number of sites in 
  spectradf_5mer <- spectradf_5mer %>% 
    group_by(species,population,label,central3mer) %>% # must group by central 3mer!!
    mutate(FAKE_total_mutations_MultinomDrawBasedOn5merTargetFracWithin3merCategory_AfterHumanRescaling=as.numeric(rmultinom(n=1,size=total_mutations_InOverallCentral3merCategory_AfterHumanRescaling,prob=targetFractionOf5merTypeWithinCentral3mer_HUMAN))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
    ungroup()# ungroup just in case here 
  
  ### NOTE  : here you are not yet using downsampled counts ; that happens in the next function! am using full counts prior to downsampling for each species ! (then downsample outside of this function) <- important! 
  # subset to just the parts that were due to multinomial sampling: 
  # am aggressively labeling as "FAKE" so it doesn't ever go into any other analyses 
  FAKE_spectradf_5mer_FAKE <- select(spectradf_5mer,-c(central3mer,targetFractionOf5merTypeWithinCentral3mer_HUMAN,total_mutations_InOverallCentral3merCategory_AfterHumanRescaling,total_mutations_RESCALED_BY_HUMAN_TARGETS,total_mutations_ORIGINAL)) %>% # get rid of these columns that are new or being replaced
    rename(total_mutations_RESCALED_BY_HUMAN_TARGETS=FAKE_total_mutations_MultinomDrawBasedOn5merTargetFracWithin3merCategory_AfterHumanRescaling) # renaming so that it matches the format of my other spectra
  # this worked!
  # okay so this apportioned total 3mers proportionally across 5mers based on 5mer target size, with mulitnomial samplign to add some randomness. 
  print("returning 5mer spectrum that has been human-rescaled already with counts that are a multinomial draw within each 3mer based on human-rescaled target fraction of each 5mer type -- these are NOT empirical counts! note that fake counts are renamed to be 'total_mutations_RESCALED_BY_HUMAN_TARGETS' to work with downstream functions")
  return(FAKE_spectradf_5mer_FAKE) # return the control spectrum
  
}


FAKE_all5merSpectra_humanRescaled <- makeFakeControl5merSpectrumFunction_humanRescaled(all5merSpectra_HumanRescaled) # aggressivley labelled it as FAKE so that it doesnt' end up in any empirical analyses 
# checking that it worked:
# kept central3mer in for testing (then went back to taking it out)
# checked that values are the same:
#sum(FAKE_all5merSpectra_humanRescaled[FAKE_all5merSpectra_humanRescaled$label=="brown_bear_ABC"& FAKE_all5merSpectra_humanRescaled$central3mer=="TAA.TTA",]$total_mutations_RESCALED_BY_HUMAN_TARGETS) # 17039
#all5merSpectra_HumanRescaled_forchecking <- all5merSpectra_HumanRescaled

#all5merSpectra_HumanRescaled_forchecking$central3mer <- paste0(substr(all5merSpectra_HumanRescaled_forchecking$mutation_5mer,2,4),".",substr(all5merSpectra_HumanRescaled_forchecking$mutation_5mer,8,10)) # 17039.13 

# now downsample the control spectrum:
FAKE_all5merSpectra <- downsampleFunc_RescaledByHumanTargets(FAKE_all5merSpectra_humanRescaled,controldatadir)

# this is now downsampled!
write.table(FAKE_all5merSpectra,paste0(controldatadir,"FAKE_all5merSpectra_rescaledToHumanTargetsBeforeShuffling.Downsampled.txt"),quote=F,sep="\t",row.names = F)

FAKE_all5merSpectra_distances_phylo <- combo_function_phylo(FAKE_all5merSpectra,variable =variableToUseInMantelTest, ilr_or_clr_or_none = "clr",speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances,epsilon=epsilon)


write.table(FAKE_all5merSpectra_distances_phylo,paste0(controldatadir,"FAKE_all5merSpectra.Distances.txt"),quote=F,sep="\t",row.names = F)

# run mantel test on control spectrum:
FAKE_all5merSpectra_mantelTestResults <- comboFunction_mantelTest_SQRTDISTANCE(FAKE_all5merSpectra,phylo_dist= raxml_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)

write.table(FAKE_all5merSpectra_mantelTestResults,paste0(controldatadir,"FAKE_all5merSpectra.MantelTest.SqrtPhylo.txt"),quote=F,sep="\t",row.names = F)

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKE5merSpectrum <- ggplot(FAKE_all5merSpectra_distances_phylo[FAKE_all5merSpectra_distances_phylo$variable==variableToUseInMantelTest & FAKE_all5merSpectra_distances_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color="randomized control (multinomial sampled w/in 3-mer)"))+ # fake distances based on multinomial sample within 3mer
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest& all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$id=="5-mer spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance,color="empirical"))+ # real distances
  geom_text(data=FAKE_all5merSpectra_mantelTestResults,aes(x=0.5,y=Inf,vjust=3.5,label=paste0("[ RANDOMIZED CONTROL DATA ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)"),color="randomized control (multinomial sampled w/in 3-mer)"),size=4)+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id=="5-mer spectrum",],aes(x=0.5,y=Inf,vjust=1.1,label=paste0("[ EMPIRICAL DATA ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)"),color="empirical"),size=4)+ # empirical mantel results 
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance\n[Based on randomized 5mer counts multinomial sampled based on 5mer target size within a 3mer (after human rescaling)]"))+
  theme_bw()+
  theme(legend.title = element_blank())+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between species' mutation spectra")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKE5merSpectrum

ggsave(paste0(controldatadir,"manteltest.sqrtPhylo.exploded5merAnalyses.ComparingFakeEmpirical.FakeGeneratedByMultinomSamplOverTargetSizesWithin3mer.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKE5merSpectrum,height=5,width=9)


################## fake 7mer spectrum ##############################

################ repeat exploded analysis for 7mers:
### different function:
makeFakeControl7merSpectrumFunction_humanRescaled <- function(spectradf_7mer){
  ##### spectra should have target info and human info in it #####
  # want to use all7merSpectra_plusHumanInfo #
  if(nchar(spectradf_7mer[1,]$mutation_label)!=15){
    stop("this isn't a 7mer spectrum!")
    
  }
  # get central 5mer  (not 3mer)
  spectradf_7mer$central5mer <- paste0(substr(spectradf_7mer$mutation_7mer,2,6),".",substr(spectradf_7mer$mutation_7mer,10,14))
  spectradf_7mer <- spectradf_7mer %>% 
    group_by(species,population,label,central5mer) %>%
    mutate(targetFractionOf7merTypeWithinCentral5mer_HUMAN = total_target_count_HUMAN/sum(total_target_count_HUMAN),total_mutations_InOverallCentral5merCategory_AfterHumanRescaling=sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    ungroup()
  
  # also need to get total number of sites in 
  spectradf_7mer <- spectradf_7mer %>% 
    group_by(species,population,label,central5mer) %>% # must group by central 5mer!!
    mutate(FAKE_total_mutations_MultinomDrawBasedOn7merTargetFracWithin5merCategory_AfterHumanRescaling=as.numeric(rmultinom(n=1,size=total_mutations_InOverallCentral5merCategory_AfterHumanRescaling,prob=targetFractionOf7merTypeWithinCentral5mer_HUMAN))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
    ungroup()# ungroup just in case here 
  
  
  # subset to just the parts that were due to multinomial sampling: 
  # am aggressively labeling as "FAKE" so it doesn't ever go into any other analyses 
  FAKE_spectradf_7mer_FAKE <- select(spectradf_7mer,-c(central5mer,targetFractionOf7merTypeWithinCentral5mer_HUMAN,total_mutations_InOverallCentral5merCategory_AfterHumanRescaling,total_mutations_RESCALED_BY_HUMAN_TARGETS,total_mutations_ORIGINAL)) %>% # get rid of these columns that are new or being replaced
    rename(total_mutations_RESCALED_BY_HUMAN_TARGETS=FAKE_total_mutations_MultinomDrawBasedOn7merTargetFracWithin5merCategory_AfterHumanRescaling) # renaming so that it matches the format of my other spectra
  # this worked!
  # okay so this apportioned total 5mers proportionally across 7mers based on 7mer target size, with mulitnomial samplign to add some randomness. 
  print("returning 7mer spectrum that has been human-rescaled already with counts that are a multinomial draw within each 5mer based on human-rescaled target fraction of each 7mer type -- these are NOT empirical counts! note that fake counts are renamed to be 'total_mutations_RESCALED_BY_HUMAN_TARGETS' to work with downstream functions")
  return(FAKE_spectradf_7mer_FAKE) # return the control spectrum
  
}

FAKE_all7merSpectra_humanRescaled <- makeFakeControl7merSpectrumFunction_humanRescaled(all7merSpectra_HumanRescaled) # 
# now downsample the control spectrum:
FAKE_all7merSpectra <- downsampleFunc_RescaledByHumanTargets(FAKE_all7merSpectra_humanRescaled,controldatadir)

# this is now downsampled!
write.table(FAKE_all7merSpectra,paste0(controldatadir,"FAKE_all7merSpectra_rescaledToHumanTargetsBeforeShuffling.Downsampled.txt"),quote=F,sep="\t",row.names = F)

FAKE_all7merSpectra_distances_phylo <- combo_function_phylo(FAKE_all7merSpectra,variable =variableToUseInMantelTest, ilr_or_clr_or_none = "clr",speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances,epsilon=epsilon)


write.table(FAKE_all7merSpectra_distances_phylo,paste0(controldatadir,"FAKE_all7merSpectra.Distances.txt"),quote=F,sep="\t",row.names = F)

# run mantel test on control spectrum:
FAKE_all7merSpectra_mantelTestResults <- comboFunction_mantelTest_SQRTDISTANCE(FAKE_all7merSpectra,phylo_dist= raxml_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)

write.table(FAKE_all7merSpectra_mantelTestResults,paste0(controldatadir,"FAKE_all7merSpectra.MantelTest.SqrtPhylo.txt"),quote=F,sep="\t",row.names = F)

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKE7merSpectrum <- ggplot(FAKE_all7merSpectra_distances_phylo[FAKE_all7merSpectra_distances_phylo$variable==variableToUseInMantelTest & FAKE_all7merSpectra_distances_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color="randomized control (multinomial sampled w/in 5-mer)"))+ # fake distances based on multinomial sample within 5mer
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest& all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$id=="7-mer spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance,color="empirical"))+ # real distances
  geom_text(data=FAKE_all7merSpectra_mantelTestResults,aes(x=0.5,y=Inf,vjust=3.5,label=paste0("[ RANDOMIZED CONTROL DATA ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)"),color="randomized control (multinomial sampled w/in 5-mer)"),size=4)+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id=="7-mer spectrum",],aes(x=0.5,y=Inf,vjust=1.1,label=paste0("[ EMPIRICAL DATA ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)"),color="empirical"),size=4)+ # empirical mantel results 
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance\n[Based on randomized 7mer counts multinomial sampled based on 7mer target size within a 5mer (after human rescaling)]"))+
  theme_bw()+
  theme(legend.title = element_blank())+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between species' mutation spectra")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKE7merSpectrum

ggsave(paste0(controldatadir,"manteltest.sqrtPhylo.exploded7merAnalyses.ComparingFakeEmpirical.FakeGeneratedByMultinomSamplOverTargetSizesWithin5mer.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKE7merSpectrum,height=5,width=9)

######## plot exploded kmer plots side by side for manuscript ###############

FAKE_all5merSpectra_distances_phylo$id <- "5-mer spectrum"
FAKE_all7merSpectra_distances_phylo$id <- "7-mer spectrum"

FAKE_all5merSpectra_mantelTestResults$id <- "5-mer spectrum"
FAKE_all7merSpectra_mantelTestResults$id <- "7-mer spectrum"


FAKE_distances_5mer_7mer_combo <- bind_rows(FAKE_all5merSpectra_distances_phylo,FAKE_all7merSpectra_distances_phylo)

# combo
FAKE_all5mer_and7mer_Spectra_mantelTestResults_combo <- bind_rows(FAKE_all5merSpectra_mantelTestResults,FAKE_all7merSpectra_mantelTestResults)

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKECONTROLSpectrum_COMBOPLOT57mers <- ggplot(FAKE_distances_5mer_7mer_combo[FAKE_distances_5mer_7mer_combo$variable==variableToUseInMantelTest & FAKE_distances_5mer_7mer_combo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color="randomized control (multinomial sampled w/in 3-mer)"),size=1)+ # fake distances based on multinomial sample within 3mer
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest& all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=sqrt(cophenetic_distance),y=spectrum_distance,color="empirical"),size=1)+ # real distances
  facet_wrap(~id,scales="free_y")+
  geom_text(data=FAKE_all5mer_and7mer_Spectra_mantelTestResults_combo,aes(x=0.5,y=Inf,vjust=3.5,label=paste0("[ CONTROL DATA ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)"),color="randomized control (multinomial sampled w/in 3-mer)"),size=4)+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=0.5,y=Inf,vjust=1.1,label=paste0("[ EMPIRICAL DATA ] Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)"),color="empirical"),size=4)+ # empirical mantel results 
  #ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance\n[Based on randomized 5mer counts multinomial sampled based on 5mer target size within a 3mer (after human rescaling)]"))+
  theme_bw()+
  theme(legend.title = element_blank(),legend.position = "none",text=element_text(size=14),strip.background = element_rect(fill="lightblue"))+
  xlab("sqrt of phylogenetic distance")+
  ylab("Aitchison distance between species' mutation spectra")+
  scale_color_manual(values=c("black","red"))

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKECONTROLSpectrum_COMBOPLOT57mers

ggsave(paste0(controldatadir,"manteltest.sqrtPhylo.exploded5mer7merAnalyses.SideBySide.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKECONTROLSpectrum_COMBOPLOT57mers,height=5,width=11)
ggsave(paste0(controldatadir,"manteltest.sqrtPhylo.exploded5mer7merAnalyses.SideBySide.png"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_BasedOnFAKECONTROLSpectrum_COMBOPLOT57mers,height=5,width=11,dpi=300)


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
  print(confounder)
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
allConfounderMantelResultsPlot_plusPhyloPValues_13mer <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("1-mer spectrum","3-mer spectrum"),]),aes(x=niceConfounderLabel,y=-log10(signif),fill=niceLabelForPlotting))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  theme(legend.title = element_blank())+
  #geom_text(position = position_dodge(width = 1),aes(label=paste0("r =\n",round(statistic,2))),size=3)+
  theme(strip.background = element_rect(fill="lightblue"))+
  facet_grid(~id,switch="y")+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=10,"Paired"))))
  #ggtitle(paste0("Correlation of confounders with spectrum distance\n",confoundermantelTestPermutationCount," permutations"))
allConfounderMantelResultsPlot_plusPhyloPValues_13mer
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.phyloMantel.13mer.pdf"),allConfounderMantelResultsPlot_plusPhyloPValues_13mer,height=4, width=7)
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.phyloMantel.13mer.png"),allConfounderMantelResultsPlot_plusPhyloPValues_13mer,height=4, width=7,dpi=300)


allConfounderMantelResultsPlot_plusPhyloPValues_57mer <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("5-mer spectrum","7-mer spectrum"),]),aes(x=niceConfounderLabel,y=-log10(signif),fill=niceLabelForPlotting))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  theme(legend.title = element_blank())+
  #geom_text(position = position_dodge(width = 1),aes(label=paste0("r =\n",round(statistic,2))),size=3)+
  theme(strip.background = element_rect(fill="lightblue"))+
  facet_grid(~id,switch="y")+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=10,"Paired"))))
#ggtitle(paste0("Correlation of confounders with spectrum distance\n",confoundermantelTestPermutationCount," permutations"))
allConfounderMantelResultsPlot_plusPhyloPValues_57mer
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.phyloMantel.57mer.pdf"),allConfounderMantelResultsPlot_plusPhyloPValues_57mer,height=4, width=7)
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.phyloMantel.57mer.png"),allConfounderMantelResultsPlot_plusPhyloPValues_57mer,height=4, width=7,dpi=300)

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

# plot 1/3mer and 5/7mer separately:
allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL_13meronly <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist)[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel" & all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("1-mer spectrum","3-mer spectrum"),],aes(x=niceConfounderLabel,y=-log10(signif),fill=niceConfounderLabel))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14),strip.background = element_rect(fill="lightblue"))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  theme(legend.title = element_blank())+
  #geom_text(position = position_dodge(width = 1),aes(label=paste0("r =\n",round(statistic,2))),size=3)+
  theme(text=element_text(size=20))+
  facet_wrap(~id,ncol=1)+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=6,"Dark2"))))
  #ggtitle(paste0("Correlation of confounders with spectrum distance\n",confoundermantelTestPermutationCount," permutations"))
allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL_13meronly
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.NOPHYLOMANTEL.13meronly.pdf"),allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL_13meronly,height=6, width=6)


allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL_57meronly <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist)[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel" & all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=niceConfounderLabel,y=-log10(signif),fill=niceConfounderLabel))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14),strip.background = element_rect(fill="lightblue"))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  theme(legend.title = element_blank())+
  #geom_text(position = position_dodge(width = 1),aes(label=paste0("r =\n",round(statistic,2))),size=3)+
  theme(text=element_text(size=20))+
  facet_wrap(~id,ncol=1)+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=6,"Dark2"))))
allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL_57meronly
ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.NOPHYLOMANTEL.57meronly.pdf"),allConfounderMantelResultsPlot_plusPhyloPValues_NOPHYLOMANTEL_57meronly,height=6, width=6)

####### plot r values (5/7mer only) ######

# annotate p values;
all_confounder_spectrum_mantel_results_PLUSphylodist$signif_star <- ""
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<0.05,]$signif_star <- "*"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<0.01,]$signif_star <- "**"
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<0.001,]$signif_star <- "***"

all_confounder_spectrum_mantel_results_PLUSphylodist$signif_annotation <- "NS"
# annotate those below threshold as significant:
all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),]$signif_annotation <- paste0(scientific(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),]$signif,1),"\n",all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_PLUSphylodist$confounder))),]$signif_star)


ylimmin_rValues_57mer_nonsep=min(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel" & all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("5-mer spectrum","7-mer spectrum"),]$statistic)

ylimmax_rValues_57mer_nonsep=max(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel" & all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("5-mer spectrum","7-mer spectrum"),]$statistic)
# add 20% so that annotations will show up
ylimmax_rValues_57mer_nonsep=ylimmax_rValues_57mer_nonsep+.2*ylimmax_rValues_57mer_nonsep

allConfounderMantelResultsPlot_rValues_NOPHYLOMANTEL_57meronly <- ggplot(data.frame(all_confounder_spectrum_mantel_results_PLUSphylodist[all_confounder_spectrum_mantel_results_PLUSphylodist$niceLabelForPlotting!="phyloMantel" & all_confounder_spectrum_mantel_results_PLUSphylodist$id %in% c("5-mer spectrum","7-mer spectrum"),]),aes(x=niceConfounderLabel,y=statistic,fill=niceConfounderLabel))+
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
  scale_y_continuous(limits=c(min(0,ylimmin_rValues_57mer_nonsep),ylimmax_rValues_57mer_nonsep))

allConfounderMantelResultsPlot_rValues_NOPHYLOMANTEL_57meronly

ggsave(paste0(confounderManteldir,"allMantelResults.allConfounders.AllSpectra.rvaluesOnly.PLUSPHYLODISTFORCOMPARISON.NOPHYLOMANTEL.57meronly.rvalues.pdf"),allConfounderMantelResultsPlot_rValues_NOPHYLOMANTEL_57meronly,height=6, width=6)

###### also plot r values:
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

############### Mantel test on confounders, separate spectrum by central BP ###########
###### this takes too long to do with 10M permutations (and that's overkill ) -- redo with fewer perms

# you are here
# key is to filter the spectrum before processing it (same order of operations as above)
# this is the function to run mantel test on 3/5/7mers separated by central 1mer 
confounderManteldir_sepbymutationtype=paste0(confounderManteldir,"sepByMutationType/")
dir.create(confounderManteldir_sepbymutationtype,showWarnings = F)

runMantelTest_onConfounders_vs_SpectrumDistances_separateByMutationType <- function(confounder_full_dataframe, nameOfConfounderYouWant,spectrum,speciesToInclude,variable,epsilon,confoundermantelTestPermutationCount,centralMutationType){
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
    filter(mutation_1mer==centralMutationType) %>% # filter first!!! 
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
  mantelTestResult_df <- data.frame(call="vegan_mantel_confounders",confounder=nameOfConfounderYouWant,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),centralMutationType=centralMutationType) # adding central mut type as a column
  
  return(mantelTestResult_df)
}


#### phylo-corrected mantel, separated by central bp #######

runPHYLOMantelTest_onConfounders_vs_SpectrumDistances_CorrectsForPhyloSignal_separateByMutationType <- function(confounder_full_dataframe, nameOfConfounderYouWant,spectrum,speciesToInclude,variable,epsilon,confoundermantelTestPermutationCount,raxmlTree_renamedTips_subset,centralMutationType){
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
    filter(mutation_1mer==centralMutationType) %>% # filter first!!! 
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
  phylomantelTestResult_df <- data.frame(call="evolqg_PhyloMantel_SEPBYCENTRALBP",confounder=nameOfConfounderYouWant,statistic=as.numeric(as.list(phylo_mtr)$rho),permCount=as.numeric(confoundermantelTestPermutationCount),signif=((((as.numeric(as.list(phylo_mtr)$Probability))*as.numeric(confoundermantelTestPermutationCount))+1)/(as.numeric(confoundermantelTestPermutationCount)+1)),centralMutationType=centralMutationType) # seems to be called rho arbitrarily (don't think it's spearmans rho, )
  
  return(phylomantelTestResult_df)
}

#### function to get distances for plotting, sep by central bp
getConfounderAndSpectrumDistancesForPlotting_separateByMutationType <- function(confounder_full_dataframe, nameOfConfounderYouWant,spectrum,speciesToInclude,variable,epsilon,centralMutationType){
  
  confounders_restrictedSpecies <- confounder_full_dataframe[confounder_full_dataframe$label  %in% speciesToInclude,] # in case you want to exclude some species from the analyses. for now it should be all 13
  # note that when you're doing this to match the phylogenetic tree you use the full_name_to_match_RAXML
  # but here you can use "label"! because both spectrum and confounders have same labels -- nice!
  # then restrict to the particular nameOfConfounderYouWant you want :
  confounders_wide = pivot_wider(select(confounders_restrictedSpecies,c("label",nameOfConfounderYouWant)),values_from = all_of(nameOfConfounderYouWant),names_from = "label" )
  
  confounder_pairwise_distances_forPlotting = tidy(dist(t(confounders_wide),diag=F,upper=F)) # want full matrix for mantel test
  colnames(confounder_pairwise_distances_forPlotting) <- c("item1","item2","confounder_distance")
  ##### get spectrum distances in tidy format for potting:
  spectrum_dist_forPlotting <- spectrum %>%
    filter(mutation_1mer==centralMutationType) %>% # filter first!!!
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
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE = data.frame()
for(confounder in confoundersIWantToTest){
  for(centralMutationType in centralMutationTypes){
    # dont' run on 1mer spectrum
    mantel_test_results_comparingConfoundersToSpectra_withoutphylo_SEPBYMUTYTPE_LIST <- lapply(listOfDFs[c("3-mer spectrum" ,"5-mer spectrum" ,"7-mer spectrum")],runMantelTest_onConfounders_vs_SpectrumDistances_separateByMutationType,confounder_full_dataframe=confounders,speciesToInclude=speciesToInclude_phylo,nameOfConfounderYouWant=confounder,confoundermantelTestPermutationCount=confoundermantelTestPermutationCount,epsilon=epsilon,variable = variableToUseInMantelTest,centralMutationType=centralMutationType) ### IF YOU CHANGE ANY OF THESE PARAMETERS BE SURE TO CHANGE THEM IN FUNCTION BELOW AS WELL!!!! 
    
    # turn into dataframe with spectrum in id column 
    mantel_test_results_comparingConfoundersToSpectra_withoutphylo_SEPBYMUTYTPE_df = bind_rows(mantel_test_results_comparingConfoundersToSpectra_withoutphylo_SEPBYMUTYTPE_LIST,.id="id") # bind over spectra
    
    # add to all-confounders dataframe:
    all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE <- bind_rows(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE,mantel_test_results_comparingConfoundersToSpectra_withoutphylo_SEPBYMUTYTPE_df) # has confounder as a column already # going to write this out at the end
    
    
    ###### run phylo corrected mantel test:
    PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_SEPBYMUTYTPE_LIST <- lapply(listOfDFs[c("3-mer spectrum" ,"5-mer spectrum" ,"7-mer spectrum")],runPHYLOMantelTest_onConfounders_vs_SpectrumDistances_CorrectsForPhyloSignal_separateByMutationType,confounder_full_dataframe=confounders,speciesToInclude=speciesToInclude_phylo,nameOfConfounderYouWant=confounder,confoundermantelTestPermutationCount=confoundermantelTestPermutationCount,epsilon=epsilon,variable = variableToUseInMantelTest,centralMutationType=centralMutationType,raxmlTree_renamedTips_subset=raxmlTree_renamedTips_subset) ### IF YOU CHANGE ANY OF THESE PARAMETERS BE SURE TO CHANGE THEM IN FUNCTION BELOW AS WELL!!!! 
    
    # turn into dataframe with spectrum in id column 
    PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_SEPBYMUTYTPE_df = bind_rows(PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_SEPBYMUTYTPE_LIST,.id="id") # bind over spectra
    
    # add to all-confounders dataframe:
    all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE <- bind_rows(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE,PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_SEPBYMUTYTPE_df) # has confounder as a column already # going to write this out at the end
    
    ######### MAKE *SURE* you're using all the same settings below as in the function above, otherwise plot and mantel values will be silently discordant (not ideal defensive coding, but works okay)
    distances_spectra_and_confounders_forPlottingWithMantelResults_separateByMutationType_LIST <- lapply(listOfDFs[c("3-mer spectrum" ,"5-mer spectrum" ,"7-mer spectrum")],getConfounderAndSpectrumDistancesForPlotting_separateByMutationType,confounder_full_dataframe=confounders,speciesToInclude=speciesToInclude_phylo,nameOfConfounderYouWant=confounder,epsilon=epsilon,variable = variableToUseInMantelTest,centralMutationType=centralMutationType) 
    
    # turn into dataframe with spectrum in id column 
    distances_spectra_and_confounders_forPlottingWithMantelResults_separateByMutationType_df = bind_rows(distances_spectra_and_confounders_forPlottingWithMantelResults_separateByMutationType_LIST,.id="id") # bind over spectra
    
    # make a plot without labels:
    # note: don't use sqrt of either distance! that's just a phylo distance thing
    mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_SEPBYMUTTYPE_nolabels <- ggplot(distances_spectra_and_confounders_forPlottingWithMantelResults_separateByMutationType_df,aes(y=spectrum_distance,x=confounder_distance))+ # note need aes_string here so you can call confounder distance
      geom_point(size=1,color="dodgerblue")+
      # am assigning text position as middle value of x-axis confounder value so that is what that median( ) is doing there
      geom_text(data=mantel_test_results_comparingConfoundersToSpectra_withoutphylo_SEPBYMUTYTPE_df,aes(x=median(distances_spectra_and_confounders_forPlottingWithMantelResults_separateByMutationType_df$confounder_distance),y=Inf,hjust=0.2,vjust=1.1,label=paste0(confounder," Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=4)+
      geom_text(data=PHYLOmantel_test_results_comparingConfoundersToSpectra_withPhyloCorrection_SEPBYMUTYTPE_df,aes(x=median(distances_spectra_and_confounders_forPlottingWithMantelResults_separateByMutationType_df$confounder_distance),y=Inf,hjust=0.2,vjust=2.3,label=paste0(confounder," Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (phyloMantel Test)")),size=4)+
      ggtitle(paste0(centralMutationType,"\nCorrrelation between spectrum distance and ",confounder,"\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_comparingConfoundersToSpectra_withoutphylo_SEPBYMUTYTPE_df$permCount)," permutations\n with and without correcting for phylogenetic signal"))+
      theme_bw()+
      theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
      xlab(paste0(confounder," abs. diff."))+
      ylab("Aitchison distance")+
      facet_wrap(~id,scales="free_y")
    
    ggsave(paste0(confounderManteldir_sepbymutationtype,confounder,".mantelTestResults.allspectra.",centralMutationType,".pdf"),mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_SEPBYMUTTYPE_nolabels,height=5,width=12)
    
    # make a plot with labels:
    mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_SEPBYMUTTYPE_withlabels <- mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_SEPBYMUTTYPE_nolabels+
      geom_text_repel(aes(label=comparisonLabel),size=1.3,max.overlaps=100)
    
    ggsave(paste0(confounderManteldir_sepbymutationtype,confounder,".mantelTestResults.allspectra.labels.",centralMutationType,".pdf"),mantelTestCorrelationsPlot_spectrum_confounder_notusingsqrtonpurpose_SEPBYMUTTYPE_withlabels,height=5,width=9)
    
    
  }
}

######### summarize p-values:
# write out:
write.table(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE,paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.ALLSEPBYCENTRALMUTTYPE.txt"),row.names=F,quote=F,sep="\t")

# add in phylo dist results and plot together
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist <- bind_rows(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE,all_mantelresults_phylo_variableToUseInMantelTest_clr_SEPBYCENTRALBP_SQRTDIST)
# make confounder a character so I can add another one
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder <- as.character(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder)

all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[is.na(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder),]$confounder <- "sqrt_phylogenetic_distance"
# reorganize order of confounders for plot:
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder <- factor(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder,levels=c("sqrt_phylogenetic_distance","avg_coverage_per_individual","scaffold_N50","contig_N50","wattersons_theta"  ,"Rspan_d","AFR_d"))


# nice labels
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceLabelForPlotting <- ""
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$call=="evolqg_PhyloMantel_SEPBYCENTRALBP",]$niceLabelForPlotting <- "phyloMantel"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$call=="vegan_mantel_confounders",]$niceLabelForPlotting <- "uncorrected Mantel"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$call=="vegan_mantel_SQRTDISTANCE",]$niceLabelForPlotting <- "phylogenetic signal"

all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceLabelForPlotting <- factor(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceLabelForPlotting,levels=c("phylogenetic signal","uncorrected Mantel","phyloMantel"))

######### nice confounder plot sep by bp for manuscript #############
### make nice confounder labels:
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceConfounderLabel <- ""
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder=="scaffold_N50",]$niceConfounderLabel <- "scaffold N50"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder=="contig_N50",]$niceConfounderLabel <- "contig N50"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder=="avg_coverage_per_individual",]$niceConfounderLabel <- "coverage"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder=="Rspan_d",]$niceConfounderLabel <- "reproductive span"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder=="AFR_d",]$niceConfounderLabel <- "age at reproduction"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder=="wattersons_theta",]$niceConfounderLabel <- "genetic diversity"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder=="sqrt_phylogenetic_distance",]$niceConfounderLabel <- "phylogeny"


# order them:
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceConfounderLabel <- factor(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceConfounderLabel, levels=c("phylogeny","coverage","scaffold N50","contig N50","genetic diversity","reproductive span","age at reproduction"))

######### annotate p-values as to whether they fall above or below alpha cutoff after correction for multiple testing ##########
# annotate p-values as above or below sig line
# cutoff: 0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6

# * = <0.05, ** < 0.01, *** < 0.001
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif_star <- ""
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif<0.05,]$signif_star <- "*"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif<0.01,]$signif_star <- "**"
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif<0.001,]$signif_star <- "***"

all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif_annotation <- "NS"

#all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),]$signif_annotation <- "*"

# or annotate with pvalue for significant ones::
all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),]$signif_annotation <- paste0(scientific(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),]$signif,1),"\n",all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$signif<(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),]$signif_star)

# plot summary of mantel results:
allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues <- ggplot(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist),aes(x=niceConfounderLabel,y=-log10(signif),fill=niceLabelForPlotting))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  facet_grid(~id~centralMutationType,switch="y")+
  theme(legend.title = element_blank(),strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  ggtitle(paste0("-log10(p-values) from Mantel test with ",confoundermantelTestPermutationCount," permutations\nblack sig line is 0.05/total confounders/6 mutation types (tests per spectrum)"))+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=10,"Paired"))))+
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))
allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues
ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.pdf"),height=10, width=18)

######### separate phylomantel results to 3mer and 5/7mer #######
# there's no 1mer faceted by central bp 
allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_3mer <- ggplot(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist)[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$id=="3-mer spectrum",],aes(x=niceConfounderLabel,y=-log10(signif),fill=niceLabelForPlotting))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  facet_grid(~id~centralMutationType,switch="y")+
  theme(legend.title = element_blank(),strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  #ggtitle(paste0("-log10(p-values) from Mantel test with ",confoundermantelTestPermutationCount," permutations\nblack sig line is 0.05/total confounders/6 mutation types (tests per spectrum)"))+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=10,"Paired"))))+
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_3mer

ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.3mer.pdf"),allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_3mer,height=5, width=18)
ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.3mer.png"),allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_3mer,height=5, width=18,dpi=300)


allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_57mer <- ggplot(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist)[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=niceConfounderLabel,y=-log10(signif),fill=niceLabelForPlotting))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  facet_grid(~id~centralMutationType,switch="y")+
  theme(legend.title = element_blank(),strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  #ggtitle(paste0("-log10(p-values) from Mantel test with ",confoundermantelTestPermutationCount," permutations\nblack sig line is 0.05/total confounders/6 mutation types (tests per spectrum)"))+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=10,"Paired"))))+
  theme(plot.margin = margin(1,1,1.5,1.2, "cm"))

allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_57mer

ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.57mer.pdf"),allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_57mer,height=8, width=18)
ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.57mer.png"),allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_57mer,height=8, width=18,dpi=300)

# make a version of this plot that is 3mer only with uncorrected mantel -- save the bigger plot for the SI 
allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_3merOnly_noPhyloMantel_forms <- ggplot(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist)[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$id=="3-mer spectrum" & all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceLabelForPlotting!="phyloMantel",],aes(x=niceConfounderLabel,y=-log10(signif),fill=niceConfounderLabel))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  facet_grid(id~centralMutationType,switch="y")+
  theme(legend.title = element_blank(),strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  #ggtitle(paste0("-log10(p-values) from Mantel test with ",confoundermantelTestPermutationCount," permutations\nblack sig line is 0.05/total confounders/6 mutation types (tests per spectrum)"))+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=8,"Dark2"))))+
  theme(legend.position="none")
allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_3merOnly_noPhyloMantel_forms

ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.3meronly.nophylomantel.forms.pdf"),allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_3merOnly_noPhyloMantel_forms,height=4, width=12)


# make a version of this plot that is 5/7mer only with uncorrected mantel -- save the bigger plot for the SI 
allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_57merOnly_noPhyloMantel_forms <- ggplot(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist)[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$id %in% c("7-mer spectrum","5-mer spectrum") & all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceLabelForPlotting!="phyloMantel",],aes(x=niceConfounderLabel,y=-log10(signif),fill=niceConfounderLabel))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  geom_hline(yintercept = -log10(0.05/length(unique(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$confounder))/6),linetype="dashed")+ # dividing by total confounder tests per spectrum not by all tests (since 3mer and 5mer results arent ind. tests)
  #geom_hline(yintercept = -log10(1/(confoundermantelTestPermutationCount+1)),linetype="dashed",color="red")+ # max possible pvalue
  xlab("")+
  ylab("-log10(p-value)")+
  facet_grid(id~centralMutationType,switch="y")+
  theme(legend.title = element_blank(),strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  #ggtitle(paste0("-log10(p-values) from Mantel test with ",confoundermantelTestPermutationCount," permutations\nblack sig line is 0.05/total confounders/6 mutation types (tests per spectrum)"))+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=8,"Dark2"))))+
  theme(legend.position="none")
allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_57merOnly_noPhyloMantel_forms

ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.57meronly.nophylomantel.forms.pdf"),allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_57merOnly_noPhyloMantel_forms,height=5, width=11)

##########  plot R values faceted per central bp (5/7 only) :##########
ylimmin_rValues_57mer=min(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$id %in% c("7-mer spectrum","5-mer spectrum") & all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceLabelForPlotting!="phyloMantel",])$statistic)

ylimmax_rValues_57mer=max(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$id %in% c("7-mer spectrum","5-mer spectrum") & all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceLabelForPlotting!="phyloMantel",])$statistic) 
# and add 10%:
ylimmax_rValues_57mer=ylimmax_rValues_57mer+.1*ylimmax_rValues_57mer


allConfounderMantelResultsPlot_SEPBYCENTRALBP_rValues_57merOnly_forms <- ggplot(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist[all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$id %in% c("7-mer spectrum","5-mer spectrum") & all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist$niceLabelForPlotting!="phyloMantel",]),aes(x=niceConfounderLabel,y=statistic,fill=niceConfounderLabel))+
  geom_col(position="dodge")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
  xlab("")+
  ylab("Pearson's r")+
  facet_grid(id~centralMutationType,switch="y")+
  theme(legend.title = element_blank(),strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
  scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=8,"Dark2"))))+
  theme(legend.position="none")+
  geom_text(aes(label = signif_annotation, x = niceConfounderLabel, y = statistic), position = position_dodge(width = 0.8), size=2,vjust = -0.1)+
  scale_y_continuous(limits=c(min(0,ylimmin_rValues_57mer),ylimmax_rValues_57mer),breaks=seq(round(ylimmin_rValues_57mer,1),round(ylimmax_rValues_57mer,1),0.2))


allConfounderMantelResultsPlot_SEPBYCENTRALBP_rValues_57merOnly_forms

ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.rvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.57meronly.nophylomantel.forms.rvals.pdf"),allConfounderMantelResultsPlot_SEPBYCENTRALBP_rValues_57merOnly_forms,height=6, width=11)



# R values: 
# note this includes both methods! r vals are the same between phylo and regular mantel
# allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_rValues <- ggplot(data.frame(all_confounder_spectrum_mantel_results_SEPARATED_BY_MUTATION_TYPE_PLUSphylodist),aes(x=confounder,y=statistic,fill=niceLabelForPlotting))+
#   geom_col(position="dodge")+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=14))+
#   xlab("")+
#   ylab("Pearson's r")+
#   facet_grid(~id~centralMutationType,switch="y")+
#   scale_fill_manual(values=c("black",c(RColorBrewer::brewer.pal(n=10,"Paired"))))+
#   theme(legend.title = element_blank(),strip.background = element_rect(fill="lightblue"),text=element_text(size=14))+
#   ggtitle(paste0("Pearson's r values from Mantel test with ",confoundermantelTestPermutationCount," permutations"))
# 
# allConfounderMantelResultsPlot_SEPBYCENTRALBP_plusPhyloPValues_rValues
# ggsave(paste0(confounderManteldir_sepbymutationtype,"allMantelResults.allConfounders.AllSpectra.pvaluesOnly.PLUSPHYLODISTFORCOMPARISON.SEPARATEDBYBP.RVALUES.pdf"),height=10, width=16)

############# Calculate pairwise cosine similarities between non-CLR transformed mutation spectra COUNTS (not proportions(?) ###############
# this is new
# uses sigfit's cosine_sim function to keep things similar to my sigfit cosine similarity calculations

# use the human-rescaled and downsampled spectra (don't need 'process' spectra or to do ilr or clr transform (whole point is to not do those so I compare how cosine simlarity performs))

# head(all1merSpectra)
# 
# cosineDISSimilarity_keepMatrix_UPPERANDLOWER <- function(pivotedspectradf_noIDVars){
#     # make it a matrix with each row as a species (transposed and all numeric)
#     # variable in question:
#     # distances are row-wise nto column wise so need to transpose 
#     cosineSimilarity <- lsa::cosine(as.matrix(pivotedspectradf_noIDVars)) #  # transpose for dist otherwise it goes along rows and is wrong
#     cosineDISSimilarity <- 1-cosineSimilarity
#     # want cosine DISsimiliarity (1-cosine similarity) so that it scales up with phylo distance the same way the other distance measures do
#     # this results in cosine disimilarity matrix
#     return(cosineDISSimilarity)
#   }
# 
# cosineDISSimilarity_tidiedIntoTable <- function(pivotedspectradf_noIDVars){
#   # make it a matrix with each row as a species (transposed and all numeric)
#   # variable in question:
#   # distances are row-wise nto column wise so need to transpose 
#   cosineSimilarity <- lsa::cosine(as.matrix(pivotedspectradf_noIDVars)) #  # transpose for dist otherwise it goes along rows and is wrong
#   cosineDISSimilarity <- 1-cosineSimilarity
#   cosineDISSimilarity_tidy <- tidy(as.dist(cosineDISSimilarity)) # makes it a distance matrix (doesn't change values, I checked and then tidies it into a dataframe using broom)
#   # want cosine DISsimiliarity (1-cosine similarity) so that it scales up with phylo distance the same way the other distance measures do
#   # this results in cosine disimilarity matrix
#   # rename as item 1 item2 cosine_dis
#   colnames(cosineDISSimilarity_tidy) <- c("item1","item2","cosine_dissimilarity")
#   return(cosineDISSimilarity_tidy)
# }
# 
# 
# comboFunction_mantelTest_SQRTDISTANCE_USINGCOSINEDISSIMILARITY <- function(spectrum, phylo_dist,variable,speciesToInclude,mantelTestPermutationCount,epsilon,speciesCodes){
#   spectrum_dist <- spectrum %>% 
#     processSpectra(.,epsilon) %>% 
#     merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
#     filter(label %in% speciesToInclude) %>%
#     pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
#     clr_or_ilr_orNoTransform_calculations(.,"none") %>% # fixing as NO transformation; this step just removes variable columns
#     cosineDISSimilarity_keepMatrix_UPPERANDLOWER(.)   # calculates cosine DISsimilarity (1-cosine similarity) as a matrix 
#   # the matrices need to be ordered in the same way 
#   # reorder the time matrix in the same way as the spectrum matrix:
#   phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
#   head(phylo_dist_reordered)
#   phylo_dist_reordered_sqrt <- sqrt(phylo_dist_reordered)
#   mtr = vegan::mantel(xdis=phylo_dist_reordered_sqrt,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
#   # make a data frame
#   mantelTestResult_df <- data.frame(call="vegan_mantel_SQRTPHYLODISTANCE_COSINEDISIMILARITY",variable="cosine_disimilarity",method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
#   # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
#   
#   return(mantelTestResult_df)
#   
# }
# 
# 
# mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_list =  lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE_USINGCOSINEDISSIMILARITY,phylo_dist= raxml_cophenetic_dist,variable="total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS",speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
# # you are here ! 
# 
# mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_list
# mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_df = bind_rows(mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_list,.id="id")
# mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_df$distance_metric <- "raxml_tree_SQRT"
# 
# 
# # get cosine similarities to plot:
# getCosineDisimilaritiesForPlotting <- function(spectrum, variable,speciesToInclude,epsilon,speciesCodes){
#   
#   cosine_table <- spectrum %>% 
#     processSpectra(.,epsilon) %>% 
#     merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
#     filter(label %in% speciesToInclude) %>%
#     pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
#     clr_or_ilr_orNoTransform_calculations(.,"none") %>% # fixing as NO transformation; this step just removes variable columns
#     cosineDISSimilarity_tidiedIntoTable(.)
#     
#   
#   
#   return(cosine_table)
#   
# }
# 
# cosineDisimilarityTables <- lapply(listOfDFs,getCosineDisimilaritiesForPlotting,variable="total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS",speciesToInclude=speciesToInclude_phylo,epsilon=epsilon,speciesCodes=speciesCodes)
# 
# # collapse into one dataframe:
# 
# cosineDisimilarityTables_id <- bind_rows(cosineDisimilarityTables,.id="id")
# 
# # merge with phylo distances:
# cosineDisimilarityTables_id$comparisonLabel <- paste0(pmin(as.character(cosineDisimilarityTables_id$item1),as.character(cosineDisimilarityTables_id$item2)),".",pmax(as.character(cosineDisimilarityTables_id$item1),as.character(cosineDisimilarityTables_id$item2)))  # alphabetical
# 
# # check that all are in phylo comp labels:
# sum(!cosineDisimilarityTables_id$comparisonLabel %in% phyloDistances$comparisonLabel) # cool they are all in there.
# 
# cosineDisimilarityTables_mergedWithPhyloDistance <- merge(cosineDisimilarityTables_id,phyloDistances,by="comparisonLabel") # dim stays the same! 
# 
# # get cosine similarity as well:
# cosineDisimilarityTables_mergedWithPhyloDistance$cosine_similarity <- 1-cosineDisimilarityTables_mergedWithPhyloDistance$cosine_dissimilarity
# 
# # now can plot:
# cosine_dissimilarity_plot_SQRTPHYLODIST <- ggplot(cosineDisimilarityTables_mergedWithPhyloDistance,aes(x=sqrt(cophenetic_distance),y=cosine_dissimilarity))+
#   geom_point(size=1)+
#   facet_wrap(~id,scales="free")+
#   geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_df,aes(x=0.5,y=Inf,vjust=1.1,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
#   ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_COSINEDISIMILARITY_df$permCount)," permutations-- SQRT of phylo distance"))+
#   theme_bw()+
#   theme(strip.background = element_rect(fill="lightblue"),strip.text=element_text(size=16),text=element_text(size=12))+
#   xlab("square root of phylogenetic distance (shared branch length)\nbranch lengths represent expected substitutions per site")+
#   ylab("Cosine dissimilarty between mutation spectra")
# 
# ggsave(paste0(manteldir,"mantelTest.COSINESIMILARITY.SQRTDISTANCE.pdf"),cosine_dissimilarity_plot_SQRTPHYLODIST,height=8,width=11)
# 
# # plot as violin plot:
# head(cosineDisimilarityTables_mergedWithPhyloDistance)
# 
# 
# # all types:
# cosineSimilarityViolinPlot1 <- ggplot(cosineDisimilarityTables_mergedWithPhyloDistance,aes(x=id,y=cosine_similarity))+
#   geom_violin(fill="lightblue",alpha=0.5)+
#   theme_bw()+
#   ylab("Pairwise cosine similarity between\nspecies' empirical mutation spectra")+
#   theme(text=element_text(size=14),axis.text.x = element_text(angle=45,hjust=1))+
#   xlab("")+
#   scale_y_continuous(breaks=seq(0.75,1,by=0.05))
# cosineSimilarityViolinPlot1
# ggsave(paste0(plotdir,"COSINESIMILARITY.allTypes.pdf"),cosineSimilarityViolinPlot1,height=4,width=8)
# 
# # then maybe just 1/3mers:
# cosineSimilarityViolinPlot2 <- ggplot(cosineDisimilarityTables_mergedWithPhyloDistance[cosineDisimilarityTables_mergedWithPhyloDistance$id %in% c("1-mer spectrum","3-mer spectrum"),],aes(x=id,y=cosine_similarity))+
#   geom_violin(fill="lightblue",alpha=0.5,draw_quantiles = c(0.5))+
#   theme_bw()+
#   ylab("Pairwise cosine similarity between\nspecies' empirical mutation spectra")+
#   theme(text=element_text(size=14))+
#   xlab("")+
#   scale_y_continuous(breaks=seq(0.92,1,by=0.01))
# cosineSimilarityViolinPlot2
# ggsave(paste0(plotdir,"COSINESIMILARITY.1-3merOnly.pdf"),cosineSimilarityViolinPlot2,height=4,width=5)
# 
# # and 5/7mers
# cosineSimilarityViolinPlot3 <- ggplot(cosineDisimilarityTables_mergedWithPhyloDistance[cosineDisimilarityTables_mergedWithPhyloDistance$id %in% c("5-mer spectrum","7-mer spectrum"),],aes(x=id,y=cosine_similarity))+
#   geom_violin(fill="lightblue",alpha=0.5,draw_quantiles = c(0.5))+
#   theme_bw()+
#   ylab("Pairwise cosine similarity between\nspecies' empirical mutation spectra")+
#   theme(text=element_text(size=14))+
#   xlab("")+
#   scale_y_continuous(breaks=seq(0.75,1,by=0.05))
# cosineSimilarityViolinPlot3
# ggsave(paste0(plotdir,"COSINESIMILARITY.5-7merOnly.pdf"),cosineSimilarityViolinPlot3,height=4,width=5)

############ save workspace image ###################
save.image(file = paste0(plotdir,"RWORKSPACE.",todaysdate,".AtEndOfScript.image.RData"))

