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
flagOfAnalysis=paste0(todaysdate,".plots.SubsampledIndividualsNotProjection.epsilon.",epsilon,".sepCpGs.no.masked",maskLabel,".RESCALEDTOHUMANTARGETS.SimplifiedCode.",as.character(mantelTestPermutationCount),".MantelPermutations.ILR1-3mer") # a label that will be in title of plots and outdir 
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
#write.table(all7merSpectra_plusTargets,paste0(plotdir,"all7merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")


all5merSpectra_plusTargets <- merge(all5merSpectraOnly_filledin, all5merTargets,by.x=c("species","ancestral5mer"),by.y=c("species","target_5mer"))

#write.table(all5merSpectra_plusTargets,paste0(plotdir,"all5merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")

all3merSpectra_plusTargets <- merge(all3merSpectraOnly_filledin, all3merTargets,by.x=c("species","ancestral3mer"),by.y=c("species","target_3mer"))

#write.table(all3merSpectra_plusTargets,paste0(plotdir,"all3merSpectra_plusTargets.ToTransferToSageForEnrichments.txt"),row.names=F,quote=F,sep="\t")


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
#all5merSpectra <- downsampleFunc_RescaledByHumanTargets(all5merSpectra_HumanRescaled,plotdir)
#all7merSpectra <- downsampleFunc_RescaledByHumanTargets(all7merSpectra_HumanRescaled,plotdir)


# write these out:
write.table(all1merSpectra,paste0(plotdir,"all1merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
write.table(all3merSpectra,paste0(plotdir,"all3merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
#write.table(all5merSpectra,paste0(plotdir,"all5merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)
#write.table(all7merSpectra,paste0(plotdir,"all7merSpectra.CountsAndTargets.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.txt"),quote=F,sep="\t",row.names = F)





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
listOfDFs = list(`1-mer spectrum`=all1merSpectra,`3-mer spectrum`=all3merSpectra) # 1/3 only
# lapply applies combo_function_phylo to each of the above spectra ^^ 
# if you set the names of a list, they get preserved in the output of lapply! cool!


listOfVars=c("mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS") # no longer using the normalized versions ; what about nondownsampled frac seg sites? adding in total mutations as well
listOfTransforms=c("ilr") # doing ilr here. 

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



write.table(all_distances_all_conditions_phylo,paste0(plotdir,"all_distances_all_conditions_phylo.RESCALED_BY_HUMAN_TARGETS.txt"),row.names=F,quote=F,sep="\t")



write.table(all_distances_all_conditions_time,paste0(plotdir,"all_distances_all_conditions_timeCalibrated.RESCALED_BY_HUMAN_TARGETS.txt"),row.names=F,quote=F,sep="\t")


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
  for(transform in c("ilr")){
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

# cut here


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
    clr_or_ilr_orNoTransform_calculations(.,"ilr") %>% # fixing as ILR here (this script only!) -- could set as an option in the function tho
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
    clr_or_ilr_orNoTransform_calculations(.,"ilr") %>% # fixing as clr -- could set as an option in the function tho
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
#mantel_test_results_raxml_list = lapply(listOfDFs,comboFunction_mantelTest,phylo_dist= raxml_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
# bind into a df:
#mantel_test_results_raxml_df = bind_rows(mantel_test_results_raxml_list,.id="id")
#mantel_test_results_raxml_df$distance_metric <- "raxml_tree"

mantel_test_results_raxml_SQRTDISTANCE_list =  lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE,phylo_dist= raxml_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
mantel_test_results_raxml_SQRTDISTANCE_df = bind_rows(mantel_test_results_raxml_SQRTDISTANCE_list,.id="id")
mantel_test_results_raxml_SQRTDISTANCE_df$distance_metric <- "raxml_tree_SQRT"


#mantel_test_results_time_list = lapply(listOfDFs,comboFunction_mantelTest,phylo_dist= timeTree_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
#mantel_test_results_time_df = bind_rows(mantel_test_results_time_list,.id="id")
#mantel_test_results_time_df$distance_metric <- "time_tree"

mantel_test_results_time_SQRTDISTANCE_list = lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE,phylo_dist= timeTree_cophenetic_dist,variable=variableToUseInMantelTest,speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount,epsilon=epsilon,speciesCodes=speciesCodes)
mantel_test_results_time_SQRTDISTANCE_df = bind_rows(mantel_test_results_time_SQRTDISTANCE_list,.id="id")
mantel_test_results_time_SQRTDISTANCE_df$distance_metric <- "time_tree"



mantel_results_all <- bind_rows(mantel_test_results_raxml_SQRTDISTANCE_df,mantel_test_results_time_SQRTDISTANCE_df)
write.table(mantel_results_all,paste0(manteldir,"mantelTest.Results.nononswrtandnofolded.txt"),row.names=F,quote=F,sep="\t")
# did this work?

############ plot mantel r values:  ############
# just raxml: 
# mantel_r_plot_raxml <- 
#   ggplot(mantel_test_results_raxml_df,aes(x=id,y=statistic,fill=distance_metric))+
#   geom_col(position="dodge")+
#   #facet_wrap(distance_label~transform_label,ncol=2)+
#   theme_bw() +
#   scale_fill_manual(values=c("#1F78B4"))+
#   xlab("")+
#   theme(text=element_text(size=14),legend.title = element_blank())+
#   ggtitle("Mantel Test: r statistic (Pearson)")+
#   ylab("correlation")
# mantel_r_plot_raxml
# ggsave(paste0(manteldir,"mantelTest.r.stat.PerSpectrum.raxmlTreeOnly.pdf"),mantel_r_plot_raxml,height=4,width=7)
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


######## making this look nice for manuscript :
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="ilr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
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
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_Just1merSpectrumForMs <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="ilr" & all_distances_all_conditions_phylo$id=="1-mer spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
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

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_distances_all_conditions_phylo$transform_label=="ilr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
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
mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_13merOnly <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$id %in% c("1-mer spectrum","3-mer spectrum") & all_distances_all_conditions_phylo$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_distances_all_conditions_phylo$transform_label=="ilr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
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





####### plot jsut 1mer with labels for my manually labeling
mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="ilr" &all_distances_all_conditions_phylo$id=="1-mer spectrum",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df[mantel_test_results_raxml_SQRTDISTANCE_df$id =="1-mer spectrum",],aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=2)

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.Just1merToHelpMakeManualLabels.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb_Just1merWithLabels,height=8,width=11)


# compare with and without sqrt patterns:

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="ilr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color="sqrt"))+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="ilr",],aes(x=cophenetic_distance,y=spectrum_distance,color="nonsqrt"))+
  facet_wrap(~id,scales="free")+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  xlab("sqrt or non-sqrt cophenetic distance")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.COMPARISON.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc,height=8,width=11)




########### mantel plot: time tree ###################
# mantelTestCorrelationsPlot_time <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="ilr",],aes(x=cophenetic_distance,y=spectrum_distance))+
#   geom_point()+
#   facet_wrap(~id,scales="free")+
#   #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (",permCount," permutations)")))+
#   geom_text(data=mantel_test_results_time_df,aes(x=150,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",signif(statistic,2),"\np-value: ",as.character(scientific(signif,2))," (Mantel Test)")),size=5)+
#   ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
#   theme_bw()
# mantelTestCorrelationsPlot_time
# ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.pdf"),mantelTestCorrelationsPlot_time,height=8,width=11)

########### timetree sqrt distance mantel plot #############
mantelTestCorrelationsPlot_time_SQRTDISTANCE <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="ilr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
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
mantelTestCorrelationsPlot_time_SQRTDISTANCE_13mer <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$id %in% c("1-mer spectrum","3-mer spectrum") & all_distances_all_conditions_time$variable==variableToUseInMantelTest & all_distances_all_conditions_time$transform_label=="ilr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
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



############ separate 1/3mer from 5/7mer plots for ms ##################

# want to add an "all mutation types" flag so that they match
# just 1/3mers: 
mantelTestCorrelationsPlot_phylo_SQRTDISTANCE_13merONLY <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variableToUseInMantelTest & all_distances_all_conditions_phylo$transform_label=="ilr" & all_distances_all_conditions_phylo$id %in% c("1-mer spectrum" ,"3-mer spectrum"),],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
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


############ save workspace image ###################
save.image(file = paste0(plotdir,"RWORKSPACE.",todaysdate,".AtEndOfScript.image.RData"))

