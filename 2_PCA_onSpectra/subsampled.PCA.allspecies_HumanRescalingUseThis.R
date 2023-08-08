############################# PCA on per individual ########################
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(gridExtra)
require(grid)
require(ggpubr)
require(stats)
require(ggfortify)
require(patchwork) # for putting plots together 
require(ggrepel)
require(compositions)

source("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/colors_and_labels.R") # get niceLabels df: 
#### IF YOU CHANGE SPECIES DESIGNATIONS NEED TO UPDATE THIS FILE ^^^^^
niceLabels
set.seed(42) # so same inds are sampled each time 
epsilon=1 # 
separateCpGs="no" # yes or no
todaysdate=format(Sys.Date(),"%Y%m%d")
maskLabel="maskALL" # if species get different mask labels this will be a pain.
wd="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/allspecies_mutyper_results_unified_SubsettingIndividuals/7_mer/" # using same individuals as were used for subset distances
flagOfAnalysis=paste0(todaysdate,".subsampledspectraPCA.epsilon.",epsilon,".sepCpGs.",separateCpGs,".masked",maskLabel,"_addingWolves_newPolForMsBBPB_HumanRescaling.SameIndsAsNonProjectedDists")  # a label that will be in title of plots and outdir 
#RColorBrewer::brewer.pal("Dark2",n=8)
#colorValuesPerSpecies=data.frame(humans=)

outdir=paste0("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/PCA_subsampled_individual_spectra/",flagOfAnalysis,"/")
dir.create(outdir,showWarnings = F,recursive = T)
kmersize="7"

# note: don't want to write anything to allspecies_summed_up_over_intervals_forTransfer because if I redownload things I don't want to get overwritten
filesuffix=".summedup.mutyper.spectrum.SeeLogForFilters.maskALL.7mer.SUBSETPRIORTOMUTYPERSPECTRUM.PERINDIVIDUAL.RANDOMIZEDWITHINPOPULATION.txt" # I manually gzipped them

#speciesList=c('humans', 'mice', 'bears', 'fin_whale', 'vaquita', 'Gorilla_gorilla', 'Pan_troglodytes', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus','wolves')
speciesList=c('humans', 'Mus_musculus','Mus_spretus' ,'polar_bear','brown_bear', 'fin_whale', 'vaquita', 'Gorilla_gorilla', 'Pan_troglodytes', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus','wolves')

popList=list(brown_bear=c("ABC","EUR"),polar_bear=c("PB"),fin_whale=c("ENP","GOC"),humans=c("AFR","EUR","SAS","EAS","AMR"),Mus_musculus=c("Mmc","Mmd","Mmm"),Mus_spretus=c("Ms"))

#THESE ARE NOW SET IN FILE THAT IS SOURCED niceLabels = data.frame(label=c("bears_ABC","bears_EUR","bears_PB","fin_whale_ENP","fin_whale_GOC","Gorilla_gorilla","humans_AFR","humans_AMR","humans_EAS","humans_EUR","humans_SAS","mice_Mmc","mice_Mmd","mice_Mmm","mice_Ms","Pan_paniscus","Pan_troglodytes", "Pongo_abelii","Pongo_pygmaeus","vaquita"),niceLabel=c("brown bear (ABC)","brown bear (EUR)","polar bear","fin whale (ENP)","fin whale (GOC)","gorilla","human (AFR)","human (AMR)","human (EAS)","human (EUR)","human (SAS)","mouse (Mmc)","mouse (Mmd)","mouse (Mmm)", "mouse (Ms)","bonobo","chimpanzee","orangutan (Sumatran)","orangutan (Bornean)","vaquita"),shortSpeciesLabel=c("brown bear","brown bear","polar bear","fin whale","fin whale","gorilla","human","human","human","human","human","house mouse","house mouse","house mouse","Algerian mouse","bonobo","chimpanzee","S. orangutan","B. orangutan","vaquita"))
#niceLabels
subsetcount=5 # 5 diploids per pop/species

all_subsampled_spectra <- data.frame()
# subsample as you go! 
for(species in speciesList) {
  print(species)
  indir=paste0(wd,"allspecies_summed_up_over_intervals_forTransfer/",species,"/mutyper_results_masked_",maskLabel,"/mutyper_spectrum_files/")
  pops=popList[[species]]
  if(is.null(pops)) {
    print("nopops")
    spectrum <- read.table(paste0(indir,species,filesuffix),header=T)
    spectrum$label <- species
    # subset the spectrum:
    all_samples_in_spectrum <- unique(spectrum$sample)
    
    randomsubsetofsamples <-sample(all_samples_in_spectrum,subsetcount)
    
    subset_spectrum <- spectrum[spectrum$sample %in% randomsubsetofsamples,]
    # add together:
    all_subsampled_spectra <- bind_rows(all_subsampled_spectra,subset_spectrum)
    
  } else {
    for(pop in pops){
      print(pop)
      spectrum <- read.table(paste0(indir,species,"_",pop,filesuffix),header=T)
      spectrum$label <- paste0(spectrum$species,"_",spectrum$pop)
      # subset the spectrum:
      all_samples_in_spectrum <- unique(spectrum$sample)
      
      randomsubsetofsamples <-sample(all_samples_in_spectrum,subsetcount)
      
      subset_spectrum <- spectrum[spectrum$sample %in% randomsubsetofsamples,]
      # add together:
      all_subsampled_spectra <- bind_rows(all_subsampled_spectra,subset_spectrum)
    }
  }

  
}

write.table(all_subsampled_spectra,paste0(outdir,"subsampled.individual.spectra.",subsetcount,".diploids.masked.",maskLabel,".txt"),row.names = F,quote=F,sep="\t")
# list of individuals per pop:
write.table(unique(all_subsampled_spectra[,c("species","population","sample")]),paste0(outdir,"subsampled.individualIDs.subsample.txt"),row.names = F,quote=F,sep="\t")


############## get 5,3,1 mer ################
##### 1mers are dealt with in the if statement below this (dealing with CpGs) ########
# get ancestral kmers:
all_subsampled_spectra$ancestral7mer <- substr(all_subsampled_spectra$mutation_type,1,7)
all_subsampled_spectra$ancestral5mer <- substr(all_subsampled_spectra$mutation_type,2,6)
all_subsampled_spectra$ancestral3mer <- substr(all_subsampled_spectra$mutation_type,3,5)

# get central mutation types:
all_subsampled_spectra$mutation_7mer <- all_subsampled_spectra$mutation_type
all_subsampled_spectra$mutation_5mer <- paste0(substr(all_subsampled_spectra$mutation_type,2,6),".",substr(all_subsampled_spectra$mutation_type,10,14))
all_subsampled_spectra$mutation_3mer <- paste0(substr(all_subsampled_spectra$mutation_type,3,5),".",substr(all_subsampled_spectra$mutation_type,11,13))

# want to get ancestral 1mers for all data frames so I can plot things with central mutattion type separated later:

#### SEPARATE OUT  -- OPTIONAL ############

if(separateCpGs=="yes"){
  #### label CpGs and specify them as ancestral target as well
  print("separating out CpG sites!")
  all_subsampled_spectra$mutation_1mer_CpGNotLabeled <- paste0(substr(all_subsampled_spectra$mutation_type,4,4),".",substr(all_subsampled_spectra$mutation_type,12,12))
  all_subsampled_spectra$ancestral1mer_CpGNotLabeled <- paste0(substr(all_subsampled_spectra$mutation_type,4,4)) # adding this on 20220112 
  #### label CpGs and specify them as ancestral target as well
  all_subsampled_spectra$CpGLabel <- ""
  # find all CpG sites:
  all_subsampled_spectra[grepl("CG",all_subsampled_spectra$ancestral3mer) & all_subsampled_spectra$ancestral1mer_CpGNotLabeled=="C",]$CpGLabel <- "_CpG" 
  all_subsampled_spectra$mutation_1mer <- paste0(all_subsampled_spectra$mutation_1mer_CpGNotLabeled,all_subsampled_spectra$CpGLabel)
  # label ancestral state with CpG as well
  all_subsampled_spectra$ancestral1mer <- all_subsampled_spectra$ancestral1mer_CpGNotLabeled
  all_subsampled_spectra[all_subsampled_spectra$CpGLabel=="_CpG",]$ancestral1mer <- "C_CpG"
  
  all_subsampled_spectra$mutation_1mer <- paste0(all_subsampled_spectra$mutation_1mer_CpGNotLabeled,all_subsampled_spectra$CpGLabel)
  
  all_subsampled_spectra$ancestral1mer_CpGNotLabeled <- substr(all_subsampled_spectra$mutation_type,4,4) # not yet CpG identified 
  
  
} else if(separateCpGs=="no"){
  print("note that CpGs are not separated out!")
  all_subsampled_spectra$mutation_1mer <- paste0(substr(all_subsampled_spectra$mutation_type,4,4),".",substr(all_subsampled_spectra$mutation_type,12,12)) # if you don't want to separate out CpGs then this is just the regular center bp
  # still get ancestral1mer (but don't specify if CpG or not)
  all_subsampled_spectra$ancestral1mer <- substr(all_subsampled_spectra$mutation_type,4,4) # if you don't want to separate out CpGs then this is just the regular center bp
  
  
  
} else {
  print("invalid choice for separateCpGs")
  break
}



# get spectrum summed per kmer type (7mers won't change)
all7merSpectraOnly <- all_subsampled_spectra %>%
  # adding groupbing by central 1mer but that doesn't actually change the sums since 7mers already group by those 
  # I just want to keep it in the dataframe 
  group_by(species,population,sample,label,ancestral7mer,mutation_7mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>% 
  ungroup() # this will be identical to original spectrum (7mer), just doing for consistency of formatting
# adding ungroup 

# sum up 5mers over 7mers:
all5merSpectraOnly <- all_subsampled_spectra %>%
  group_by(species,population,sample,label,ancestral5mer,mutation_5mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>%
  ungroup()
# check:
#sum(all_subsampled_spectra[all_subsampled_spectra$sample=="008_UARC_ABC_BGI_brownbear_43157" & all_subsampled_spectra$mutation_5mer=="AAAAA.AACAA",]$totalSites)
#head(all5merSpectraOnly)
all3merSpectraOnly <- all_subsampled_spectra %>%
  group_by(species,population,sample,label,ancestral3mer,mutation_3mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>%
  ungroup()

all1merSpectraOnly <- all_subsampled_spectra %>%
  group_by(species,population,sample,label,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations = sum(totalSites)) %>%
  ungroup()
# if you don't want to partition over CpGs you can use mutation_1mer_CpGNotLabeled instead
# okay NOTE HERE: "C" counts are nonCpG and then C_CpG are CpG ancestral targets

# okay do need to process the targets #
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
  
  ######### relabel target_1mer with CpGs if you are separating out CpGs ###########
  # 20220112 -- note that grepping central CpG is fine because there are no CGX ancestral 3mers because they would be revcomped to YCG so you're good here!  
  if(separateCpGs=="yes"){
    targets$CpGLabel <- ""
    targets[grepl("CG",targets$target_3mer),]$CpGLabel <- "_CpG"
    targets$target_1mer_CpGsNotLabeled <- targets$target_1mer # update with this 
    targets$target_1mer <- paste0(targets$target_1mer_CpGsNotLabeled,targets$CpGLabel)
    
  } # don't need an else statement, just keep going with targets target$target_1mer as is.
  
  
  
  # TARGETS per 3mer type #
  targets_7mer <- targets %>% group_by(species,target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) %>% ungroup() # original targets file, just doing to be consistent but it's identical to original targets; adding ungroup on 20220912
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

############## fill in missing mutation types (7mers) ########
# complete so that missing mutation types are filled in 
# okay this is kind of sinister, depending on how it was grouped above this may or may not work (but doesn't throw an error. Frustrating!)
# so put in a check
# keep an eye on this in other code. Don't think it's made things go wrong before but be alert.
all7merSpectraOnly_filledin <- all7merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_7mer,ancestral7mer,ancestral1mer,mutation_1mer),nesting(species,population,sample,label),fill=list(total_mutations=0))%>%
  mutate(mutation_label=mutation_7mer)

if(dim(all7merSpectraOnly_filledin)[1]==dim(all7merSpectraOnly)[1]){
  print("filling in didn't work!")
}
# s houldn't be needed for 5mer 3mer 1mer
all5merSpectraOnly_filledin <- all5merSpectraOnly %>%
  ungroup() %>% # need to pre-ungroup ; got grouped up above
  complete(nesting(mutation_5mer,ancestral5mer,ancestral1mer,mutation_1mer),nesting(species,population,label,sample),fill=list(total_mutations=0))%>%
  mutate(mutation_label=mutation_5mer) # adding mutation label for use in functions

all3merSpectraOnly_filledin <- all3merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_3mer,ancestral3mer,ancestral1mer,mutation_1mer),nesting(species,population,label,sample),fill=list(total_mutations=0))%>%
  mutate(mutation_label=mutation_3mer)

all1merSpectraOnly_filledin <- all1merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_1mer,ancestral1mer),nesting(species,population,label),fill=list(total_mutations=0,sample)) %>%
  mutate(mutation_label=mutation_1mer)
# fill in should only be needed for 7mer and sometimes 5mer if sparse

# so when I was merging before all species were combined there was a bug where missing mtuation types wouldn't get targets filled in which led to them being NA downstream 
# okay now fill in missing mutation types : as 0s PRIOR TO MERGING WITH TARGEST
#  note that targets are the same within species (mouse1 and mouse2 have same targets; bear 1 and bear2 have targets)


############## merge with targets ###############

all7merSpectra_plusTargets <- merge(all7merSpectraOnly_filledin, all7merTargets,by.x=c("species","ancestral7mer"),by.y=c("species","target_7mer"))

if(dim(all7merSpectra_plusTargets)[1] != dim(all7merSpectraOnly_filledin)[1]){
  print("something went wrong with merge!")
  
}


all5merSpectra_plusTargets <- merge(all5merSpectraOnly_filledin, all5merTargets,by.x=c("species","ancestral5mer"),by.y=c("species","target_5mer"))
all3merSpectra_plusTargets <- merge(all3merSpectraOnly_filledin, all3merTargets,by.x=c("species","ancestral3mer"),by.y=c("species","target_3mer"))
all1merSpectra_plusTargets <- merge(all1merSpectraOnly_filledin, all1merTargets,by.x=c("species","ancestral1mer"),by.y=c("species","target_1mer"))


############### >>>> RESCALE TO HUMAN TARGETS <<<<< #################
# rescaling all counts to be scaled by human targets total_mutations_of_typeX_in_SpeciesA* (target_proportion_HUMAN /target_proportion_speciesA)

humanTargetProportions_1mer <- all1merTargets[all1merTargets$species=="humans",]
humanTargetProportions_3mer <- all3merTargets[all3merTargets$species=="humans",]
humanTargetProportions_5mer <- all5merTargets[all5merTargets$species=="humans",]
humanTargetProportions_7mer <- all7merTargets[all7merTargets$species=="humans",]

# going to keep these names 
all1merSpectra_plusHumanInfoNotYetRescaled <- merge(all1merSpectra_plusTargets,humanTargetProportions_1mer[,c("target_1mer","total_target_count","target_proportion")],by.x=c("ancestral1mer"),by.y=c("target_1mer"),suffixes=c("_thisSpecies","_HUMAN"))

all3merSpectra_plusHumanInfoNotYetRescaled <- merge(all3merSpectra_plusTargets,humanTargetProportions_3mer[,c("target_3mer","total_target_count","target_proportion")],by.x=c("ancestral3mer"),by.y=c("target_3mer"),suffixes=c("_thisSpecies","_HUMAN"))

all5merSpectra_plusHumanInfoNotYetRescaled <- merge(all5merSpectra_plusTargets,humanTargetProportions_5mer[,c("target_5mer","total_target_count","target_proportion")],by.x=c("ancestral5mer"),by.y=c("target_5mer"),suffixes=c("_thisSpecies","_HUMAN"))

all7merSpectra_plusHumanInfoNotYetRescaled <- merge(all7merSpectra_plusTargets,humanTargetProportions_7mer[,c("target_7mer","total_target_count","target_proportion")],by.x=c("ancestral7mer"),by.y=c("target_7mer"),suffixes=c("_thisSpecies","_HUMAN"))

# okay you want to get spectrum and the mutation rate and want to add 1 to all counts to deal with 0 entries (assuming every mutation type would happen if you sampled enough -- though note this will have a bigger i)
print("rescaling to human targets")

# NOTE! need to be dividing by HUMAN target sizes for mutation rates!!! 
#### scale counts to human targets: 
########### HERE need to group by *SAMPLE* 
# note you don't call things 'total_mutations' since not projecting here, just using selected sample size of 5 inds
all1merSpectra_HumanRescaled <- all1merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label,sample) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

all3merSpectra_HumanRescaled <- all3merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label,sample) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

all5merSpectra_HumanRescaled <- all5merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label,sample) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

all7merSpectra_HumanRescaled <- all7merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label,sample) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>% 
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)

########### >>>> DOWNSAMPLING <<<<< #######################
print("downsampling!")

# do separately for every spectrum
downsampleFunc_RescaledByHumanTargets <- function(spectradf,plotdir){
  
  spectrumName <- deparse(substitute(spectradf)) # gets the name of the object that was put in
  ###### trying experiment 
  totalSegSites <- spectradf %>%
    group_by(species,label,sample) %>%
    summarise(totalSegSites=sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    ungroup()
  
  subsampleValue=min(totalSegSites$totalSegSites)
  
  subsampleSpecies=totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species
  subsampleSample=totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$sample
  
  print(paste0("downsampling to ",subsampleValue," sites (",subsampleSample," ",subsampleSpecies,") [ counts already rescaled by human targets ]"))
  
  
  multinomlogfile <- file(paste0(plotdir,"multinomialLogFile.",spectrumName,".Downsampling.RESCALED_BY_HUMAN_TARGETS.txt"))
  
  writeLines(paste0("downsampling to ",subsampleValue,"sites ",subsampleSample," ",subsampleSpecies," [after rescaling by human targets ])"),multinomlogfile)
  
  close(multinomlogfile)
  
  # first need to get frac segregating sites:
  spectradf <- spectradf %>% 
    group_by(species,population,label,sample) %>% # group by sample!
    mutate(fractionOfSegregatingSites_RESCALED_BY_HUMAN_TARGETS = total_mutations_RESCALED_BY_HUMAN_TARGETS / sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    mutate(total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS=as.numeric(rmultinom(n=1,size=subsampleValue,prob=fractionOfSegregatingSites_RESCALED_BY_HUMAN_TARGETS))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
    ungroup() %>% # ungroup just in case here 
    
    return(spectradf)
}

##### ok!! These have now been downsampled !!!!!!!!! ###########
all1merSpectra <- downsampleFunc_RescaledByHumanTargets(all1merSpectra_HumanRescaled,outdir)
all3merSpectra <- downsampleFunc_RescaledByHumanTargets(all3merSpectra_HumanRescaled,outdir)
all5merSpectra <- downsampleFunc_RescaledByHumanTargets(all5merSpectra_HumanRescaled,outdir)
all7merSpectra <- downsampleFunc_RescaledByHumanTargets(all7merSpectra_HumanRescaled,outdir)


######### merge with nice labels for plotting ##########
all1merSpectra <- merge(all1merSpectra,niceLabels,by="label")
all3merSpectra <- merge(all3merSpectra,niceLabels,by="label")
all5merSpectra <- merge(all5merSpectra,niceLabels,by="label")
all7merSpectra <- merge(all7merSpectra,niceLabels,by="label")

################### processing function ###############
# add multinomial variables to this function.
# note: nothing says 'projected' because these are per individual and not sample size projected


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
    group_by(species,population,label,sample) %>% # group by sample! 
    mutate(fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=(total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS/sum(total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS)),fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS=(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS/sum(total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS))) %>% 
    ungroup()
  # now that things are rescaled to human targets, I am thinking that frac seg sites may not have the same issues it used to have of baking in genome composition. 
  return(spectradf)
  
}
############# test: checking distances to see if wolf/mmd close at individual level ########
test=processSpectra(all1merSpectra,epsilon)
#test_wide <- pivot_wider(test[,c("mutation_label","mutation_rate_multinom_downsampled_normalized_plusEpsilon","label","sample")],id_cols=mutation_label,values_from=mutation_rate_multinom_downsampled_normalized_plusEpsilon,names_from=c(label,sample))

#dist <- tidy(dist(t(select(test_wide,-c(mutation_label)))))
# aha this is better now -- vaquita is lowest once again (24K sites )-- note that means VERY low density for 7mers.
#head(test)
#ggplot(test[test$species %in% c("wolves","bears","mice","vaquita"),],aes(x=mutation_1mer,y=mutation_rate_multinom_downsampled_normalized_plusEpsilon,fill=label,grouo=sample))+
#  geom_col(position="dodge")

######################## PCA ##################
# adding in desired pcs
identifyingMetadataColumnNames=unique(c("label","species","population","sample",names(niceLabels))) # these will be kept when pivoting and then removed prior to pca; can add more metadata columns here ; keeping 1mer mutation label
pca_function <- function(processeddf,kmersizelabel,variable,ilr_or_clr_or_none,plotloadings=T,plotdir,desiredPCs,identifyingMetadataColumnNames,colorList){
  outdir=paste0(plotdir,"pcaPlots/transform_",transform,"/",kmersizelabel,"/",variable,"/")
  dir.create(outdir,showWarnings = F,recursive = T)

  # specify if you want loadings (don't want for 5mer or 7mer -- just turns into a cloud) 
  # this will do pca and plot it from a processed --NOT PIVOTED-- df
  # for pca need a different pivot -- need df to be 'wide' with mutation types along the top (columns)
  # and species as the rows. 
  # need to make it so that mutaiton types are columns and species are rows: 
  df_WIDE <-  spread(data.frame(processeddf[,c(identifyingMetadataColumnNames,variable,"mutation_label")]),key = mutation_label,variable) # spread the variable you want.  # have to add "sample"
  # DO NOT use mutation 1 mer at this pivot point -- separates things incorrectly
  
  ### apply transform if it is specified: how that I've pivoted teh df to have mutation types along columns
  # I need to transform ROWS not colums (different from above)
  if(ilr_or_clr_or_none=="clr"){
    transformed_df <- data.frame(clr(select(df_WIDE,-all_of(identifyingMetadataColumnNames)))) # clr works ROWWISE which here is appropriate since species are ROWs. above the species were columns so instead I had to use sapply.
    # put info back in: 
    transformed_df <- bind_cols(transformed_df,df_WIDE[,identifyingMetadataColumnNames])
    #transformed_df$species <- df_WIDE$species
    #transformed_df$population <- df_WIDE$population
    #transformed_df$label <- df_WIDE$label
    #transformed_df$sample <- df_WIDE$sample
    
    
  } else if(ilr_or_clr_or_none=="ilr") {
    
    transformed_df <- data.frame(ilr(select(df_WIDE,select=-all_of(identifyingMetadataColumnNames)))) # clr works ROWWISE which here is appropriate since species are ROWs. above the species were columns so instead I had to use sapply.
    transformed_df <- bind_cols(transformed_df,df_WIDE[,identifyingMetadataColumnNames])
    
    #transformed_df$species <- df_WIDE$species
    #transformed_df$population <- df_WIDE$population
    #transformed_df$label <- df_WIDE$label
    #transformed_df$sample <- df_WIDE$sample
    
    # note with ilr you won't get mutation types as columns any more 
  } else if(ilr_or_clr_or_none=="none"){
    transformed_df = df_WIDE # don't transform it. 
  } else {
    print("invalid transform choice")
    break
  }
  
  # sub in ">" for "." in mutation labels:
  colnames(transformed_df) <- gsub("\\.",">",colnames(transformed_df))
  
  # only use colors in colorList that are the final df:
  colorListToUse = colorList[unique(transformed_df$niceLabel)]
  ############### carry out pca and plot ###############
  
  ####### sites that are missing in all individuals after downsampling (low freq) all end up with same freuency 1/total due to downsampling ; that means they have no variance and disrupt the pca. only really an issue for non-transformed pca with freq after downsampling. so could just skip those or could remove columns without variance
  # getting variances
  variances = apply(select(transformed_df, -all_of(identifyingMetadataColumnNames)),2,var) 
  columnsWith0Variance=names(variances[variances==0]) 
  if(length(columnsWith0Variance)>0){
  print(paste0("excluding mutation types with no variance (due to extremely low observations -- should only be relevant for non-transformed 7mers with frequency after multinomial downsampling: ",length(columnsWith0Variance)," mutation types:"))
    print(columnsWith0Variance)
    transformed_df <- select(transformed_df,-all_of(columnsWith0Variance))
  }
  
  pca <- prcomp(select(transformed_df, -all_of(identifyingMetadataColumnNames)),scale=T,center=T) # Michael says to scale/center; similar with and without
  
  # save pca as rds object:
  #saveRDS(pca, file = paste0(outdir,"pca.results.IfYouNeedToReplot.rds")) -- intead, I'm saving everything at the end. 
  
  # plot scree plot:
  #pdf(paste0(outdir,"pca.screePlot1.pdf"))
  #plot(pca)
  #dev.off()
  
  # PC combos to plot:
  PCCombos = combn(desiredPCs,2)
  allPCAPlotsAcrossPCCombos <- list()
  for(i in seq(1,dim(PCCombos)[2])) {
    xpc = PCCombos[1,i]
    ypc= PCCombos[2,i]
    # plot without loadings as well: 
    
    ######### autoplots #######
    pcaPlot_noLoadings <- autoplot(pca,x=xpc,y=ypc,data=transformed_df,colour="niceLabel",size=4,alpha=1,label.colour="black")+theme_classic()+ggtitle(paste0(kmersizelabel,"\n",flagOfAnalysis,"\n",variable,"\ntransform = ",ilr_or_clr_or_none))+ geom_text_repel(aes(color=niceLabel,label=niceLabel),size=3)+scale_color_manual(values=colorListToUse) # have color codes as part of it now
    #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".noLoadings.png"),pcaPlot_noLoadings,height=7,width=10)
    ggsave(paste0(outdir,"pca.Plot.",xpc,".",ypc,".noLoadings.WithLabels.png"),pcaPlot_noLoadings,height=7,width=10)
    
    # plot nicely for manuscript
    # taking out title ggtitle(paste0(kmersizelabel,"\n",flagOfAnalysis,"\n",variable,"\ntransform = ",ilr_or_clr_or_none))
    pcaPlot_forMS <- autoplot(pca,x=xpc,y=ypc,data=transformed_df,colour="niceLabel",size=4,label.colour="black")+theme_classic()+theme(legend.title = element_blank(),text = element_text(size=16),axis.text=element_text(size=12))+ geom_text_repel(aes(color=niceLabel,label=niceLabel),size=3)+scale_color_manual(values=colorListToUse)+guides(color=guide_legend(ncol=3)) # have color codes as part of it now
    #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".ForMS.png"),pcaPlot_forMS,height=7,width=10)
    ggsave(paste0(outdir,"pca.Plot.",xpc,".",ypc,".ForMS.pdf"),pcaPlot_forMS,height=7,width=16)
    
    # also plot without labels 
    # taking out title
    pcaPlot_forMS2 <- autoplot(pca,x=xpc,y=ypc,data=transformed_df,colour="niceLabel",size=4,label.colour="black")+theme_classic()+theme(legend.title = element_blank(),text = element_text(size=16),axis.text=element_text(size=12))+scale_color_manual(values=colorListToUse)+guides(color=guide_legend(ncol=3))+
      theme(legend.position="top")
    # have color codes as part of it now
    #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".ForMS.png"),pcaPlot_forMS,height=7,width=10)
    ggsave(paste0(outdir,"pca.Plot.",xpc,".",ypc,".ForMS.NOLABELS.COLORSONLY.pdf"),pcaPlot_forMS2,height=8,width=8)

    
    if(plotloadings==T) {
      pcaPlot_withLoadings <- autoplot(pca,x=xpc,y=ypc,data=transformed_df,colour="niceLabel",size=4,label.colour="black",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+geom_text(aes(color=niceLabel,label=niceLabel),size=3)+ggtitle(paste0(kmersizelabel,"\n",flagOfAnalysis,"\n",variable,"\ntransform = ",ilr_or_clr_or_none))+theme_classic()+theme(legend.title = element_blank())+scale_color_manual(values=colorListToUse) # have color codes as part of it now
      #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".withLoadings.png"),pcaPlot_withLoadings,height=7,width=10)
      ggsave(paste0(outdir,"pca.Plot.",xpc,".",ypc,".withLoadings.png"),pcaPlot_withLoadings,height=7,width=10)
      
      # just plot the loadings by themselves:
      ## this is using autoplot  ; make sure autoplot matches my manual plot
      pcaPlot_LoadingsOnly_autoplot <- autoplot(pca,x=xpc,y=ypc,data=transformed_df,colour="transparent",alpha=1,size=4,label.colour="black",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black",loadings.label.size=5)+theme_classic()+scale_color_manual(values=colorListToUse) # have color codes as part of it now # adding in legend to keep dimensions the same (?) #+ggtitle(paste0(kmersizelabel,"\n",flagOfAnalysis,"\n",variable,"\ntransform = ",ilr_or_clr_or_none))+theme_classic()
      
      #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".JustLoadings.png"),pcaPlot_LoadingsOnly,height=7,width=10)
      # don't need to return anything
      ## plot loadings colored by central bp 
      
      loadings=data.frame(pca$rotation[,c(paste0("PC",xpc),paste0("PC",ypc))]) # pull the columns by index of x and ypc
      loadings$mutation_label <- rownames(loadings) ## NEED TO MAKE THIS CENTRAL BP!!!
      ### get central bp: 
      if(nchar(loadings$mutation_label[1])==3){
        # then it's a 1mer
        loadings$centralBP <- loadings$mutation_label
      } else if(nchar(loadings$mutation_label[1])==7){
        # then it's a 3mer
        loadings$centralBP <- paste0(substr(loadings$mutation_label,2,2),">",substr(loadings$mutation_label,6,6))
      } else {
        "not 1mer or 3mer for loadings -- are you trying to do loadings for 5 or 7mers?"
        break
      }
      
      # https://clauswilke.com/blog/2020/09/07/pca-tidyverse-style/ <-- using this to plot
      # define the arrow style:
      arrow_style <- arrow(
        angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
      )
      # note: to keep the scales teh same between pca plot and loadings plot I am making the autoplot loadings transparent and overlaying new colored loadings! works great! 
      loadplot_manual <- autoplot(pca,x=xpc,y=ypc,data=transformed_df,colour="transparent",alpha=1,size=4,label.colour="black",loadings=T,loadings.label=F,loadings.colour="transparent",loadings.label.size=5)+theme_classic()+ # making everything transparents here so that I can keep scales the same as pca plot 
        geom_segment(data=loadings, aes_string(x=paste0("PC",xpc),y=paste0("PC",ypc),color="centralBP"),xend = 0, yend = 0, arrow = arrow_style,alpha=0.8)  +
        geom_text_repel(data=loadings,aes(label=mutation_label,color=centralBP),show.legend = FALSE,segment.color=NA,alpha=0.8,size=5) +
        theme_classic()+ theme(legend.title = element_blank(),legend.position="none",text=element_text(size=16),axis.text=element_text(size=12))+
        scale_color_brewer(palette = "Set2")
      
      
      loadplot_manual_NOLABELSFORLOADINGS <- autoplot(pca,x=xpc,y=ypc,data=transformed_df,colour="transparent",alpha=1,size=4,loadings=T,loadings.label=F,loadings.colour="transparent",loadings.label.size=5)+theme_classic()+ # making everything transparents here so that I can keep scales the same as pca plot 
        geom_segment(data=loadings, aes_string(x=paste0("PC",xpc),y=paste0("PC",ypc),color="centralBP"),xend = 0, yend = 0, arrow = arrow_style)  +
        # don't have labels: geom_text_repel(data=loadings,aes(label=mutation_label,color=centralBP),show.legend = FALSE,segment.color=NA,alpha=0.8,size=5) +
        theme_classic()+ theme(legend.title = element_blank(),legend.position=c(0.1,0.3),text=element_text(size=16),axis.text=element_text(size=12))+
        scale_color_brewer(palette = "Set2")+ # note these legend position coordinates are not rlative to x/y axis, they are relative to total plot so c(0.5,0.5 is middle) regardless of axis scales
        theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black")) # put box around legend 
      
      #loadplot_manual
      ### need to get grid at same size as other plot 
      
      grid1 <- ggarrange(pcaPlot_forMS,loadplot_manual,ncol=1,common.legend = T)
      #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".LoadingsBelow.png"),grid1,height=14,width=10) # simplifying names because changed dir structure:
      ggsave(paste0(outdir,"pca.Plot.",xpc,".",ypc,".LoadingsBelow.ColoredLoadings.pdf"),grid1,height=12,width=8)
      
      # save a version without labels (label manually in AI)
      grid2 <- ggarrange(pcaPlot_forMS2,loadplot_manual,ncol=1,common.legend = T)
      #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".LoadingsBelow.png"),grid1,height=14,width=10) # simplifying names because changed dir structure:
      ggsave(paste0(outdir,"pca.Plot.",xpc,".",ypc,".LoadingsBelow.ColoredLoadings.NOLABELS.pdf"),grid2,height=12,width=8)
      
      # save a version without loading labels
      grid3 <- ggarrange(pcaPlot_forMS2,loadplot_manual_NOLABELSFORLOADINGS,ncol=1,common.legend = F,heights=c(0.55,0.45))
      #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".LoadingsBelow.png"),grid1,height=14,width=10) # simplifying names because changed dir structure:
      ggsave(paste0(outdir,"pca.Plot.",xpc,".",ypc,".LoadingsBelow.ColoredLoadings.NOLABELS.NOLOADINGLABELS.pdf"),grid3,height=12,width=8)
      
      
      
    }
    
    
  }

}

# cannot rescale a constant/zero column to unit variance this error indicates a non-numeric column made it into the pca calc 
# test :
#pca_function(test,"3mer_spectrum","mutation_rate_multinom_downsampled_normalized_plusEpsilon",ilr_or_clr_or_none = "clr",plotloadings = T,plotdir=outdir,desiredPCs = c(1,2,3))
# taking out epsilon (not using any more)
combo_function_pca <- function(spectrumdf,kmsersizelabel,variable,ilr_or_clr_or_none,plotloadings=T,plotdir,speciesToExclude,desiredPCs,identifyingMetadataColumnNames,colorList,epsilon) {
  processSpectra(spectrumdf,epsilon) %>% # taking out epsilon (not using any more)
    filter(!label %in% speciesToExclude) %>%
    pca_function(.,kmsersizelabel,variable,ilr_or_clr_or_none,plotloadings=plotloadings,plotdir,desiredPCs,identifyingMetadataColumnNames,colorList)
  
}
############ lists for running combo function ###########
listOfDFs = list(spectra_1mer=all1merSpectra,spectra_3mer=all3merSpectra,spectra_5mer=all5merSpectra,spectra_7mer = all7merSpectra)

# lapply applies combo_function_phylo to each of the above spectra ^^ 
# if you set the names of a list, they get preserved in the output of lapply! cool!

#listOfVars = c("mutation_rate_normalized_plusEpsilon","fractionOfSegregatingSites_plusEpsilon","mutation_rate_multinom_downsampled_normalized_plusEpsilon","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon")

### NEW longer list of variables to test:
listOfVars=c("mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","mutation_rate_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","fractionOfSegregatingSites_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS","total_mutations_plusEpsilon_RESCALED_BY_HUMAN_TARGETS") # want to compare counts, fracs, and mutation rates

listOfTransforms=c("none", "clr") # skipping ilr for now 

# if you need to quickly rerun the main text figures:
#listOfVars=c("mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS")
#listOfTransforms=c("clr")
### note that excluding in this list will exclude from pca but won't exclude from downsampling (!!)
speciesToExclude=c("fin_whale_ENP","Mus_musculus_Mmc","Mus_musculus_Mmm","humans_EAS","humans_SAS","humans_AMR","humans_EUR","brown_bear_EUR") # excluding these to reduce sample sizes to make more comparable
#speciesToExclude=""
# want to also downsample humans mice to two each 
for(variable in listOfVars){
  #listOfPCAPlots = list()
  print(variable)
  for(transform in listOfTransforms){
    print(transform)
    # get for all dfs: 
    ########## pca ###############
    # hm okay need to deal with labeling each one 
    # maybe just run separately because 5mer/7mer loadings are crazy
    # not going to run on list of DFs -- too annoying. 
    # 1mer:
    combo_function_pca(all1merSpectra,"1mer_spectrum",variable,transform,plotloadings=T,outdir,speciesToExclude,desiredPCs = c(1,2,3),identifyingMetadataColumnNames,colorList,epsilon)
    # 3mer:
    combo_function_pca(all3merSpectra,"3mer_spectrum",variable,transform,plotloadings=T,outdir,speciesToExclude,desiredPCs = c(1,2,3),identifyingMetadataColumnNames,colorList,epsilon)
    # 5mer: (turn off loadings) -- too sparse? need to deal with sparsity error
    combo_function_pca(all5merSpectra,"5mer_spectrum",variable,transform,plotloadings=F,outdir,speciesToExclude,desiredPCs = c(1,2,3),identifyingMetadataColumnNames,colorList,epsilon)
    # 7mer: (turn off loadings) -- too sparse?
    combo_function_pca(all7merSpectra,"7mer_spectrum",variable,transform,plotloadings=F,outdir,speciesToExclude,desiredPCs = c(1,2,3),identifyingMetadataColumnNames,colorList,epsilon)
    #listOfPCAPlots=append(listOfPCAPlots,list(pcaplot1mer,pcaplot3mer,pcaplot5mer,pcaplot7mer))
    
  }
 
  
}


####### REVISIONS: NEW STUFF adding PCAs colored by study/sequencer/read length #############
# don't need color list or loadings
pca_function_ColorByMetadata <- function(processeddf,kmersizelabel,variable,ilr_or_clr_or_none,plotdir,desiredPCs,identifyingMetadataColumnNames){
  outdir=paste0(plotdir,"pcaPlots/transform_",transform,"/",kmersizelabel,"/",variable,"_COLORBYMETADATA/")
  dir.create(outdir,showWarnings = F,recursive = T)
  
  # specify if you want loadings (don't want for 5mer or 7mer -- just turns into a cloud) 
  # this will do pca and plot it from a processed --NOT PIVOTED-- df
  # for pca need a different pivot -- need df to be 'wide' with mutation types along the top (columns)
  # and species as the rows. 
  # need to make it so that mutaiton types are columns and species are rows: 
  df_WIDE <-  spread(data.frame(processeddf[,c(identifyingMetadataColumnNames,variable,"mutation_label")]),key = mutation_label,variable) # spread the variable you want.  # have to add "sample"
  # DO NOT use mutation 1 mer at this pivot point -- separates things incorrectly
  
  ### apply transform if it is specified: how that I've pivoted teh df to have mutation types along columns
  # I need to transform ROWS not colums (different from above)
  if(ilr_or_clr_or_none=="clr"){
    transformed_df <- data.frame(clr(select(df_WIDE,-all_of(identifyingMetadataColumnNames)))) # clr works ROWWISE which here is appropriate since species are ROWs. above the species were columns so instead I had to use sapply.
    # put info back in: 
    transformed_df <- bind_cols(transformed_df,df_WIDE[,identifyingMetadataColumnNames])
    #transformed_df$species <- df_WIDE$species
    #transformed_df$population <- df_WIDE$population
    #transformed_df$label <- df_WIDE$label
    #transformed_df$sample <- df_WIDE$sample
    
    
  } else if(ilr_or_clr_or_none=="ilr") {
    
    transformed_df <- data.frame(ilr(select(df_WIDE,select=-all_of(identifyingMetadataColumnNames)))) # clr works ROWWISE which here is appropriate since species are ROWs. above the species were columns so instead I had to use sapply.
    transformed_df <- bind_cols(transformed_df,df_WIDE[,identifyingMetadataColumnNames])
    
    #transformed_df$species <- df_WIDE$species
    #transformed_df$population <- df_WIDE$population
    #transformed_df$label <- df_WIDE$label
    #transformed_df$sample <- df_WIDE$sample
    
    # note with ilr you won't get mutation types as columns any more 
  } else if(ilr_or_clr_or_none=="none"){
    transformed_df = df_WIDE # don't transform it. 
  } else {
    print("invalid transform choice")
    break
  }
  
  # sub in ">" for "." in mutation labels:
  colnames(transformed_df) <- gsub("\\.",">",colnames(transformed_df))
  
  ############### carry out pca and plot ###############
  
  ####### sites that are missing in all individuals after downsampling (low freq) all end up with same freuency 1/total due to downsampling ; that means they have no variance and disrupt the pca. only really an issue for non-transformed pca with freq after downsampling. so could just skip those or could remove columns without variance
  # getting variances
  variances = apply(select(transformed_df, -all_of(identifyingMetadataColumnNames)),2,var)
  columnsWith0Variance=names(variances[variances==0])
  if(length(columnsWith0Variance)>0){
    print(paste0("excluding mutation types with no variance (due to extremely low observations -- should only be relevant for non-transformed 7mers with frequency after multinomial downsampling: ",length(columnsWith0Variance)," mutation types:"))
    print(columnsWith0Variance)
    transformed_df <- select(transformed_df,-all_of(columnsWith0Variance))
  }
  
  pca <- prcomp(select(transformed_df, -all_of(identifyingMetadataColumnNames)),scale=T,center=T) # Michael says to scale/center; similar with and without
  
  # save pca as rds object:
  #saveRDS(pca, file = paste0(outdir,"pca.results.IfYouNeedToReplot.rds")) -- intead, I'm saving everything at the end. 
  
  # plot scree plot:
  #pdf(paste0(outdir,"pca.screePlot1.pdf"))
  #plot(pca)
  #dev.off()
  
  # PC combos to plot:
  PCCombos = combn(desiredPCs,2)
  allPCAPlotsAcrossPCCombos <- list()
  for(i in seq(1,dim(PCCombos)[2])) {
    xpc = PCCombos[1,i]
    ypc= PCCombos[2,i]
    # plot without loadings as well: 
    
    ######### autoplots #######
    pcaPlot_noLoadings_ColorByMetadata <- autoplot(pca,x=xpc,y=ypc,data=transformed_df,colour="Sequencing_Platform",shape="Read_Length",size=4,alpha=1,label.colour="black")+theme_classic()+ggtitle(paste0(kmersizelabel,"\n",flagOfAnalysis,"\n",variable,"\ntransform = ",ilr_or_clr_or_none))+
    #geom_text_repel(aes(label=niceLabel),size=3)+
      scale_color_brewer(palette = "Accent")+
      scale_shape_manual(values=c(17,16))
    # +scale_color_manual(values=colorListToUse) # have color codes as part of it now
    #ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot.",xpc,".",ypc,".noLoadings.png"),pcaPlot_noLoadings,height=7,width=10)
    ggsave(paste0(outdir,"pca.Plot.",xpc,".",ypc,".noLoadings.WithLabels.COLORBYMETADATA.pdf"),pcaPlot_noLoadings_ColorByMetadata,height=7,width=10)
    

  }
  #return(pcaPlot_noLoadings_ColorByMetadata)
  
}

#### combo function to process (downsampling already occurred) and run pca 
combo_function_pca_COLORBYMETADATA <- function(spectrumdf,kmsersizelabel,variable,ilr_or_clr_or_none,plotdir,speciesToExclude,desiredPCs,identifyingMetadataColumnNames,epsilon) {
  processSpectra(spectrumdf,epsilon) %>% # taking out epsilon (not using any more)
    filter(!label %in% speciesToExclude) %>%
    pca_function_ColorByMetadata(.,kmsersizelabel,variable,ilr_or_clr_or_none,plotdir,desiredPCs,identifyingMetadataColumnNames)
  
}

############ add in sequencer/read len information for revisions ############
extraMetadata <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/information/revisions.confounders.SequencingPlatform.ReadLength.20230710.PERINDIVIDUALFORPCA.txt",sep="\t",header=T) # just has metadata for the species/pops that get run thru pca (not extra pops)

head(extraMetadata)
head(all1merSpectra)
# merge with extra metadata and exclude spp that aren't present
all1merSpectra_PlusExtraMetadata <- merge(all1merSpectra[!all1merSpectra$label %in% speciesToExclude,],extraMetadata,by=c("species","sample","label"),all=T)

all3merSpectra_PlusExtraMetadata <- merge(all3merSpectra[!all3merSpectra$label %in% speciesToExclude,],extraMetadata,by=c("species","sample","label"),all=T)

all5merSpectra_PlusExtraMetadata <- merge(all5merSpectra[!all5merSpectra$label %in% speciesToExclude,],extraMetadata,by=c("species","sample","label"),all=T)

all7merSpectra_PlusExtraMetadata <- merge(all7merSpectra[!all7merSpectra$label %in% speciesToExclude,],extraMetadata,by=c("species","sample","label"),all=T)


# add the new metadata in to this vector:
identifyingMetadataColumnNames_PlusExtraMetadata= c(identifyingMetadataColumnNames,c("Sequencing_Platform","Read_Length"))

dir.create(paste0(outdir,"/colorByMetadata/"),showWarnings = F)

# go through each one (just doing PCs 1,2 and one variable)
combo_function_pca_COLORBYMETADATA(all1merSpectra_PlusExtraMetadata, kmsersizelabel="1mer_spectrum",variable ="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",plotdir=paste0(outdir,"/colorByMetadata/"),desiredPCs = c(1,2),identifyingMetadataColumnNames=identifyingMetadataColumnNames_PlusExtraMetadata,speciesToExclude = speciesToExclude,epsilon = epsilon)

combo_function_pca_COLORBYMETADATA(all3merSpectra_PlusExtraMetadata, kmsersizelabel="3mer_spectrum",variable ="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",plotdir=paste0(outdir,"/colorByMetadata/"),desiredPCs = c(1,2),identifyingMetadataColumnNames=identifyingMetadataColumnNames_PlusExtraMetadata,speciesToExclude = speciesToExclude,epsilon = epsilon)

combo_function_pca_COLORBYMETADATA(all5merSpectra_PlusExtraMetadata, kmsersizelabel="5mer_spectrum",variable ="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",plotdir=paste0(outdir,"/colorByMetadata/"),desiredPCs = c(1,2),identifyingMetadataColumnNames=identifyingMetadataColumnNames_PlusExtraMetadata,speciesToExclude = speciesToExclude,epsilon = epsilon)

combo_function_pca_COLORBYMETADATA(all7merSpectra_PlusExtraMetadata, kmsersizelabel="7mer_spectrum",variable ="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",plotdir=paste0(outdir,"/colorByMetadata/"),desiredPCs = c(1,2),identifyingMetadataColumnNames=identifyingMetadataColumnNames_PlusExtraMetadata,speciesToExclude = speciesToExclude,epsilon = epsilon)



save.image(file = paste0(outdir,"RWORKSPACE.",todaysdate,".AtEndOfScript.RData"))


