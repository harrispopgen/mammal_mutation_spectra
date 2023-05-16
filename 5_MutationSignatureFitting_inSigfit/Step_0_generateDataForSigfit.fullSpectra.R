require(tidyverse)
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(sigfit)
require(spgs)  # for reverse complements
require(gtools) # for mixedsort
require(RColorBrewer)
require(scales)
require(sigfit)
set.seed(42)

separateCpGs="yes" # was always doing it with CpGs separated out. if you want to do it with CpGs not separated then you need to change aging signature because it has CpGs separated out.

speciesList=c('humans', 'Mus_musculus','Mus_spretus', 'polar_bear','brown_bear' ,'fin_whale', 'vaquita', 'Gorilla_gorilla', 'Pan_troglodytes', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus','wolves','ucla_wolves') #


# note we are NOT downsampling the data and it's based on 5 individuals' spectra per species/population that have been summed (shared alleles only counted once)

outdir="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/sigfit/"

## read in per-population mutation spectra
#updated on 11/29/2022 to be based on the sum of 5 inds that I'm using for rest of paper (generated using just those 5 inds and mutyper spectrum)
#were concatenated across all pops into one file

all_spectra_7mer <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/allspecies_mutyper_results_unified_SubsettingIndividuals/7_mer/concatenated_pop_level_spectra_allSpecies_basedOnSubsetOfIndividuals/masked_maskALL/including_ucla_wolves_comparison/allSpectra.PopulationLevel.BasedOn5individualsPerPop.notprojected.usethis.txt",header=T,sep="\t") # this is in my PCA directory because it's what I used for PCA 

head(all_spectra_7mer)

#This is what I am starting with. This is the 7-mer mutation spectra for  each population based on 5 individuals that were randomly sub-sampled from each species or population within a species.
#It has the per-individual 7-mer counts for 5 individuals per population: 
  
  
  #Note that the only processing that has occurred here is:
  
  #a)	If a mutation occurs in more than one individual it is only counted once toward the population spectrum

#b)	Counts were summed up across all chromosomes of the genome (were generated per-chromosome initially for parallel processing of the genomes)


#The sum across all individuals within a population represents the population-specific mutation spectrum, with no double-counting of mutations.

## collapse 7-mers into their central 3-mers 
#The data start as 24,000+ 7-mers but for SigFit we want 3-mers, so I am collapsing the 7-mers down in to 3-mers by contatenating basepairs 3-5 and 11-13 of #the 7-mer


all_spectra_7mer$mutation_type_3mer <- paste0(substr(all_spectra_7mer$mutation_type,3,5),".",substr(all_spectra_7mer$mutation_type,11,13))

# I then group by central 3mer and sum up the 7mer counts that are within that 3mer category (so that I have 3mer counts instead of 7mer counts)

all_spectra_3mer <- all_spectra_7mer %>%
  group_by(mutation_type_3mer,species,population,label) %>%
  summarise(total_mutations_3mer=sum(totalSites)) %>%
  ungroup()

# add ancestral 3mer info
all_spectra_3mer$ancestral3mer <- substr(all_spectra_3mer$mutation_type_3mer,1,3)
all_spectra_3mer$ancestral1mer <- substr(all_spectra_3mer$mutation_type_3mer,2,2)

# not adding 1mer info yet because want to modify it in below script

# to deal with whether or not cpgs are separated out need this special function:
make1merVersionOfSpectrum <- function(melted_3mer_spectrum_df,separateCpGs=separateCpGs){
  
  
  
  melted_3mer_spectrum_df$mutation_type_1mer <- paste0(substr(melted_3mer_spectrum_df$mutation_type_3mer,2,2),".",substr(melted_3mer_spectrum_df$mutation_type_3mer,6,6)) # get the central 1-mer.
  
  # add ancestral 1mer information (going to modify this if you separate cpgs)
  melted_3mer_spectrum_df$ancestral1mer <- substr(melted_3mer_spectrum_df$mutation_type_3mer,2,2) # if you don't want to separate out CpGs then this is just the regular center bp
  
  if(separateCpGs=="yes"){
    # add in cpg label: 
    print('separating CpG.TpG mutations')
    melted_3mer_spectrum_df[grepl("CG",melted_3mer_spectrum_df$ancestral3mer) & melted_3mer_spectrum_df$mutation_type_1mer=="C.T",]$mutation_type_1mer <- "CpG.TpG"
    
    # changing how I label ancestral 1mers for different mutation types. since we aren't separating out non C>T CpG mutations I want them to have the full "C" target size (CnonCpG+C_cpg) whereas noncpg C>T mtuations should just have nonCpG_Cs as their targets and cpg C>Ts should have CpGs as target size (NOTE this is different from how I've been doing it! -- i fixed the CpG-separated distance script to incorporate this as well)
    
    melted_3mer_spectrum_df[melted_3mer_spectrum_df$mutation_type_1mer=="C.T",]$ancestral1mer <- "C_nonCpG" # non-cpg C>Ts
    melted_3mer_spectrum_df[melted_3mer_spectrum_df$mutation_type_1mer=="CpG.TpG",]$ancestral1mer <- "C_CpG"
    # note that other C>N mutations will have "C" as ancestral 1mer because all kinds of Cs are valid (non cpg + cpg)
    
    
  } 
  
  melted_1mer_spectrum_df <- melted_3mer_spectrum_df %>%
    group_by(species,population,label,mutation_type_1mer,ancestral1mer) %>%
    summarise(total_mutations_1mer=sum(total_mutations_3mer)) %>%
    ungroup()
  
  return(melted_1mer_spectrum_df)
}

all_spectra_1mer =make1merVersionOfSpectrum(all_spectra_3mer,separateCpGs = "yes")

#### human genomic target correction ####
all3merTargets = data.frame()
all1merTargets = data.frame()

for(species in speciesList){
  # some 'samples' may contain multipops like mice and bears
  targetdir=paste0("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/allspecies_mutyper_results_unified_SubsettingIndividuals/7_mer/allspecies_summed_up_over_intervals_forTransfer/",species,"/mutyper_results_masked_maskALL/mutyper_target_files/")
  targetsfilename = paste0(targetdir,species,".summedup.mutyper.targets.SeeLogForFilters.maskALL.7mer.txt")
  
  targets=read.table(targetsfilename,header=T) # these are 7mer targets 
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  targets$species <- species
  # get central 3mer and 1mer:
  targets <- targets %>%
    mutate(target_3mer=substr(target_7mer,3,5),target_1mer= substr(target_7mer,4,4)) # get substrings 
  
  ######### relabel target_1mer with CpGs if you are separating out CpGs from other C>Ts ###########
  # 20220112 -- note that grepping central CpG is fine because there are no CGX ancestral 3mers because they would be revcomped to YCG so you're good here!  
  if(separateCpGs=="yes"){
    print("separating CpGs in 1mer targets as well")
    targets$CpGLabel <- ""
    targets[grepl("CG",targets$target_3mer),]$CpGLabel <- "_CpG"
    targets[targets$target_1mer=="C" & !grepl("CG",targets$target_3mer),]$CpGLabel <- "_nonCpG" # label the nonCpG C's as well.
    # so now nothing is labeled as plain "C"
    targets$target_1mer_CpGsNotLabeled <- targets$target_1mer # update with this 
    targets$target_1mer <- paste0(targets$target_1mer_CpGsNotLabeled,targets$CpGLabel)
    
    
    
  } # don't need an else statement, just keep going with targets target$target_1mer as is.
  
  
  # TARGETS per 3mer type #
  targets_3mer <- targets %>% group_by(species,target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup()
  # TARGETS per 1mer type (which may have separated CpGs from non CpGs)
  targets_1mer <- targets %>% group_by(species,target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup() 
  
  
  # combine  all targets:
  all3merTargets = bind_rows(all3merTargets,targets_3mer)
  all1merTargets = bind_rows(all1merTargets,targets_1mer)
  
  
}
# changing here. now the non CpG "C" are called C_nonCpG instead of just "C" since I don't want to use those for the C>A and C>G mutations (want to use CpG+nonCpG C sites)
head(all3merTargets)
head(all1merTargets)

sum(all3merTargets[all3merTargets$species=="humans",]$total_target_count)
sum(all1merTargets[all1merTargets$species=="humans",]$total_target_count)
# cool looks good. 
# these match. good. (~500Mb)

# a bit of ugly code to deal with the fact that I want Cnoncpg+c_cpg to be target size for C>G and C>A but not for C>T

#### get target proportions ####
# getting proportions is now tricky with the new "C" column! Hmmmmm... Gotta be super careful here with that. 
# so maybe shouldn't add C column until later? maybe not until things are getting merged? 
all1merTargets <- all1merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()
# want to build in some way to make sure proportions didn't take anything else into account 

# add an extra row that sums up the two types of C targets (cpg and non cpg) to be an extra row (MUST DO THIS *AFTER* getting target proportions above since don't want the Cs to be double counted toward proportions)

# only need to do this if CpGs are separated
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
##### note that at this point the total counts of C+A will be equal to C_CpG + C_nonCpG + A but you shouldn't SUM over all of these categories otherwise you'll double count Cs

all3merTargets <- all3merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()


# note that NOW C>A and C>G mutations have the full "C" set as their target size, whereas before I think they only had nonCpG "C" whcih isn't what we want. 

#### 1mer: add in target info
all_spectra_1mer_plusTargets <- merge(all_spectra_1mer, all1merTargets_plusCTotals,by.x=c("species","ancestral1mer"),by.y=c("species","target_1mer"))

head(all_spectra_1mer_plusTargets,n=25)


all_spectra_3mer_plusTargets <- merge(all_spectra_3mer,all3merTargets,by.x=c("species","ancestral3mer"),by.y=c("species","target_3mer"))

head(all_spectra_3mer_plusTargets,n=25)

#####  get human target info ####
humanTargetProportions_1mer <- all1merTargets_plusCTotals[all1merTargets_plusCTotals$species=="humans",]
humanTargetProportions_3mer <- all3merTargets[all3merTargets$species=="humans",] # must use the one that has C totals!

### 1mer
all_spectra_1mer_plusTargets_MERGEDWITHUMANTARGETPROPORTIONS <- merge(all_spectra_1mer_plusTargets,humanTargetProportions_1mer[,c("target_1mer","target_proportion")],by.x=c("ancestral1mer"),by.y=c("target_1mer"),suffixes=c("_thisSpecies","_HUMAN"))

### 3mer 
all_spectra_3mer_plusTargets_MERGEDWITHUMANTARGETPROPORTIONS <- merge(all_spectra_3mer_plusTargets,humanTargetProportions_3mer[,c("target_3mer","target_proportion")],by.x=c("ancestral3mer"),by.y=c("target_3mer"),suffixes=c("_thisSpecies","_HUMAN"))


#### rescale by human targets ####

all_spectra_1mer_plusTargets_RESCALEDBYHUMANTARGETS <- all_spectra_1mer_plusTargets_MERGEDWITHUMANTARGETPROPORTIONS %>%
    group_by(species,population,label) %>%
    mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = round(total_mutations_1mer* (target_proportion_HUMAN /target_proportion_thisSpecies))) %>% # adding "round" so you end up with integers
    ungroup() # adding ungroup on 20220727, previously didn't use # doing this based on non-rescaled by human version 
  # note that this leads to differing numbers of mutations due to rescaling of CpG targets etc.


all_spectra_3mer_plusTargets_RESCALEDBYHUMANTARGETS <- all_spectra_3mer_plusTargets_MERGEDWITHUMANTARGETPROPORTIONS %>%
  group_by(species,population,label) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = round(total_mutations_3mer* (target_proportion_HUMAN /target_proportion_thisSpecies))) %>% # adding "round" so you end up with integers
  ungroup() # adding ungroup on 20220727, previously didn't use # doing this based on non-rescaled by human version 
# note that this leads to differing numbers of mutations due to rescaling of CpG targets etc.

#### for 3mers, make labels match cosmic format (don't need to do for 1mers) ####
speciesOrder = c("Mus_spretus_Ms","Mus_musculus_Mmd","Mus_musculus_Mmm","Mus_musculus_Mmc","Pongo_abelii","Pongo_pygmaeus","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","humans_AFR","humans_EUR","humans_SAS","humans_EAS","humans_AMR","wolves","ucla_wolves","brown_bear_ABC","brown_bear_EUR","polar_bear_PB","fin_whale_GOC","fin_whale_ENP","vaquita")

revComp3mersForSigfit <- function(spectrum_3mer,speciesOrder){
  spectrum_3mer$mutation_label_for_sigfit <- spectrum_3mer$mutation_type_3mer
  # already have ancestral 3mer in there, now get derived 3mer:
  spectrum_3mer$derived3mer <- substr(spectrum_3mer$mutation_type_3mer,5,7) # XXX.YYY gets the YYY part using indices 5-7
  # need to revcomp any that have A as the ancestral 1mer
  spectrum_3mer[spectrum_3mer$ancestral1mer=="A",]$mutation_label_for_sigfit <- paste0(toupper(reverseComplement(spectrum_3mer[spectrum_3mer$ancestral1mer=="A",]$ancestral3mer)),".",toupper(reverseComplement(spectrum_3mer[spectrum_3mer$ancestral1mer=="A",]$derived3mer)))
  
  
  # and want to pivot wider 
  spectrum_3mer_wider <-  pivot_wider(spectrum_3mer,id_cols=c("label"),names_from=c("mutation_label_for_sigfit"),values_from =total_mutations_RESCALED_BY_HUMAN_TARGETS)
  
  # and want to swap . for >" (though when it gets read into sigfit you'll have to do this again)
  colnames(spectrum_3mer_wider) <- gsub("\\.",">",colnames(spectrum_3mer_wider))
  
  # add species labels as rowname (important! )
  spectrum_3mer_wider <- spectrum_3mer_wider %>%
    remove_rownames %>% # have to remove any preexisting rownames 
    column_to_rownames(var="label") # adds row names
  
  data("cosmic_signatures_v3.2") # this pulls from sigfit
  
  # reorder to match the order of the cosmic signatures
  spectrum_3mer_wider <- spectrum_3mer_wider[,colnames(cosmic_signatures_v3.2)]
  
  # reorder columns to match desired species order (make sure all are in there)
  if(sum(!rownames(spectrum_3mer_wider) %in% speciesOrder)!=0){
    stop("species missing from your species order vector")
  }
  spectrum_3mer_wider<- spectrum_3mer_wider[speciesOrder,] # order species in phylogenetic order 
  
  return(spectrum_3mer_wider)
  
  
}

all_spectra_3mer_plusTargets_RESCALEDBYHUMANTARGETS_READYFORSIGFIT_matrix <- revComp3mersForSigfit(all_spectra_3mer_plusTargets_RESCALEDBYHUMANTARGETS,speciesOrder)

write.table(all_spectra_3mer_plusTargets_RESCALEDBYHUMANTARGETS_READYFORSIGFIT_matrix,paste0(outdir,"all_spectra_3mer_Matrix_notdownsampled.notprojected.RESCALEDBYHUMANTARGETS_READYFORSIGFIT.txt"),row.names = T,col.names = T,quote=F,sep="\t")

######### format 1mers for sigfit ############
# note don't need to revcomp because not using cosmic signatures. 
format1mersForSigfit_dontneedtorevcomp <- function(spectrum_1mer,speciesOrder){
  # 1mers don't need revcomping or ">" subbed for . since they don't have to match cosmic. they do have to be in the same order as the aging signatures, however, but am doing that in the next script. 
  #  want to pivot wider 
  spectrum_1mer_wider <-  pivot_wider(spectrum_1mer,id_cols=c("label"),names_from=c("mutation_type_1mer"),values_from =total_mutations_RESCALED_BY_HUMAN_TARGETS)
  
  
  # add species labels as rowname (important! )
  spectrum_1mer_wider <- spectrum_1mer_wider %>%
    remove_rownames %>% # have to remove any preexisting rownames 
    column_to_rownames(var="label") # adds row names
  
  spectrum_1mer_wider <- spectrum_1mer_wider[,c("A.C","A.G","A.T","C.A","C.G","C.T","CpG.TpG")]
  
  # reorder columns to match desired species order (make sure all are in there)
  if(sum(!rownames(spectrum_1mer_wider) %in% speciesOrder)!=0){
    stop("species missing from your species order vector")
  }
  spectrum_1mer_wider<- spectrum_1mer_wider[speciesOrder,] # order species in phylogenetic order 
  
  
  
  return(spectrum_1mer_wider)
  
  
}

all_spectra_1mer_plusTargets_RESCALEDBYHUMANTARGETS_READYFORSIGFIT_matrix <- format1mersForSigfit_dontneedtorevcomp(all_spectra_1mer_plusTargets_RESCALEDBYHUMANTARGETS,speciesOrder)

write.table(all_spectra_1mer_plusTargets_RESCALEDBYHUMANTARGETS_READYFORSIGFIT_matrix,paste0(outdir,"all_spectra_1mer_Matrix_notdownsampled.notprojected.RESCALEDBYHUMANTARGETS_READYFORSIGFIT.txt"),row.names = T,col.names = T,quote=F,sep="\t")


# I checked that this matches the countsplit results and it does.
# okay cool so coded this all up correctly 