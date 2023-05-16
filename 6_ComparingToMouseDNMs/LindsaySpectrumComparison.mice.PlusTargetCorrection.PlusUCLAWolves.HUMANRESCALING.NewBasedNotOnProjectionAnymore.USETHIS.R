####### Lindsay spectrum 
# there are lots of mutation types in this document, which ones do we want?
# going to start with all of them
require(spgs)
require(reshape2)
require(ggplot2)
require(dplyr)
require(tidyr)
require(broom)
require(compositions)
plotdir="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/lindsay_mice/"
lindsaySpectrum <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/information/lindsay_mice_naturecommunications_2019_doi.10.1038_s41467-019-12023-w/41467_2019_12023_MOESM5_ESM.txt",header=T,sep="\t")
head(lindsaySpectrum)
# rev comp to kh format 

### NOTE: on 20220616 I noticed a weird behavior of reverseComplement()
# but if you tried to do this with list of 1mers it would treat the whole thing as one string and revcomp backwards -- not good! scrambles labels. Once again, DID NOT HAPPEN IN MY 3MER REV COMP FOR SIGFIT (I checked) but be alert if using in future. Doesn't throw an error, just scrambles mutation assignment of 1mers when rev-comping. (if you make it a list it is better -- see spgs reversecomplement docs for more )
# if ref is already "C" or "A" then you can use non rev comped combo
lindsaySpectrum$mutation_type_raw <- paste0(lindsaySpectrum$ref,".",lindsaySpectrum$alt)
lindsaySpectrum$mutation_type_KHFormat <- lindsaySpectrum$mutation_type_raw # start with Raw then modify specific types:
lindsaySpectrum[lindsaySpectrum$mutation_type_raw=="G.T",]$mutation_type_KHFormat <- "C.A"
lindsaySpectrum[lindsaySpectrum$mutation_type_raw=="G.C",]$mutation_type_KHFormat <- "C.G"
lindsaySpectrum[lindsaySpectrum$mutation_type_raw=="G.A",]$mutation_type_KHFormat <- "C.T"
lindsaySpectrum[lindsaySpectrum$mutation_type_raw=="T.A",]$mutation_type_KHFormat <- "A.T"
lindsaySpectrum[lindsaySpectrum$mutation_type_raw=="T.C",]$mutation_type_KHFormat <- "A.G"
lindsaySpectrum[lindsaySpectrum$mutation_type_raw=="T.G",]$mutation_type_KHFormat <- "A.C"
# need to separate CpGs

lindsaySpectrum_Summarized <- lindsaySpectrum %>%
  group_by(mutation_type_KHFormat) %>%
  summarise(total_mutations=n()) %>%
  ungroup() #
# note these are similar fractions to what is in lindsay plot so I think using all the mutations was right move. 

# going to try target correction :
mouse_targets_7mer = read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/allspecies_mutyper_results_unified_SubsettingIndividuals/7_mer/allspecies_summed_up_over_intervals_forTransfer/Mus_musculus/mutyper_results_masked_maskALL/mutyper_target_files/Mus_musculus.summedup.mutyper.targets.SeeLogForFilters.maskALL.7mer.txt",header=T,sep="\t") # 20221220 this is latest set of target results. NOTE that this won't be the specific target for the lindsay DNMs but will hopefully allow general comparison

# note that these are mmd ancestral target for my polymorphism data. not sure how applicable these targets are to the lindsay data, but least are probably ballpark giving sense of composition of mouse genome compared ot human genome 
# these are 7mers, reduce to 1mers plx

mouse_targets_7mer$ancestralBP <- substr(mouse_targets_7mer$target,4,4)

mouse_targets_1mer <- mouse_targets_7mer %>%
  group_by(ancestralBP) %>%
  summarise(total_target_count=sum(countOverAllIntervals))  %>%
  ungroup() %>%
  mutate(target_proportion=total_target_count/sum(total_target_count))


mouse_targets_1mer
###### adjust lindsay spectrum 
head(lindsaySpectrum_Summarized)
lindsaySpectrum_Summarized$ancestralBP <- substr(lindsaySpectrum_Summarized$mutation_type_KHFormat,1,1)
# merge with target count:
lindsaySpectrum_Summarized_mergedWithABTargetCount <- merge(lindsaySpectrum_Summarized,mouse_targets_1mer,by="ancestralBP")

# want this to match myspectra so add columns:
lindsaySpectrum_Summarized_mergedWithABTargetCount$species <- "lindsay_mice"
lindsaySpectrum_Summarized_mergedWithABTargetCount$population <- "lindsay_mice"
lindsaySpectrum_Summarized_mergedWithABTargetCount$label <- "lindsay_mice"
lindsaySpectrum_Summarized_mergedWithABTargetCount$ancestral1mer <- lindsaySpectrum_Summarized_mergedWithABTargetCount$ancestralBP
lindsaySpectrum_Summarized_mergedWithABTargetCount$mutation_label <- lindsaySpectrum_Summarized_mergedWithABTargetCount$mutation_type_KHFormat
lindsaySpectrum_Summarized_mergedWithABTargetCount$mutation_1mer <- lindsaySpectrum_Summarized_mergedWithABTargetCount$mutation_type_KHFormat



# 20221129: now want to use 1mer spectra not based on projections but based on 5 individual spectra summed , so from: 
# 20221220 : this is based on 5 subset individuals and has ucla wolves (generated using makeDataFrameForLindsayAnalysis.includesuclawolves.R script)
myspectra <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/lindsay_mice/all1merSpectra_plusTargets.notyethumanrescaled.ForLindsayAnalysis.containsUCLAwolves.UsedSubsetIndividuals.notprojection.txt",header=T,sep="\t") # no longer projected!



# I want to downsample them. 
# and not add epsilon
# and add in Lindsay
myspectra_plusLindsay <- bind_rows(lindsaySpectrum_Summarized_mergedWithABTargetCount[names(myspectra)],myspectra) # bind rows together


############### >>>> RESCALE TO HUMAN TARGETS <<<<< #################
# rescaling all counts to be scaled by human targets total_mutations_of_typeX_in_SpeciesA* (target_proportion_HUMAN /target_proportion_speciesA)
# just need two targets (A and C)


humanTargetProportions_1mer <- unique(myspectra_plusLindsay[myspectra_plusLindsay$label=="humans_AFR",c("ancestral1mer","total_target_count","target_proportion")]) # gets the A and C target sizes for humans.

all1merSpectra_plusHumanInfoNotYetRescaled <- merge(myspectra_plusLindsay,humanTargetProportions_1mer[,c("ancestral1mer","total_target_count","target_proportion")],by="ancestral1mer",suffixes=c("_thisSpecies","_HUMAN"))

print("rescaling to human targets")

# NOTE! need to be dividing by HUMAN target sizes for mutation rates!!! 
#### scale counts to human targets: 

all1merSpectra_HumanRescaled <- all1merSpectra_plusHumanInfoNotYetRescaled %>%
  group_by(species,population,label) %>%
  mutate(total_mutations_RESCALED_BY_HUMAN_TARGETS = total_mutations* (target_proportion_HUMAN /target_proportion_thisSpecies)) %>%   
  ungroup() %>%
  rename(total_mutations_ORIGINAL=total_mutations)


# nice! this rescaled everything (including Lindsay) to human target sizes

##### don't want to downsample Lindsay but want to downsample everything else (?) ############
########### >>>> DOWNSAMPLING <<<<< #######################
print("downsampling!")

# do separately for every spectrum
# need to add plotdir to write files to 
downsampleFunc_RescaledByHumanTargets_ToVaquitaLevel_ExcludeLindsay <- function(spectradf_includingLindsay){
  
  spectrumName <- deparse(substitute(spectradf_includingLindsay)) # gets the name of the object that was put in
  ###### trying experiment 
  totalSegSites <- spectradf_includingLindsay %>%
    group_by(species,label) %>%
    summarise(totalSegSites=sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    ungroup()
  
  subsampleSpecies="vaquita"
  
  subsampleValue=totalSegSites[totalSegSites$label==subsampleSpecies,]$totalSegSites # setting to be vaquita
  
#  subsampleSpecies=totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species
  
  print(paste0("downsampling to ",subsampleValue," sites (",subsampleSpecies,") [ counts already rescaled by human targets ]"))
  
  
  #multinomlogfile <- file(paste0(plotdir,"multinomialLogFile.",spectrumName,".Downsampling.RESCALED_BY_HUMAN_TARGETS.txt"))
  
  #writeLines(paste0("downsampling to ",min(totalSegSites$totalSegSites)," sites (",totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species," [after rescaling by human targets ])"),multinomlogfile)
  
  #close(multinomlogfile)
  
  # want to EXCLUDE the lindsay DNMs from getting downsampled 
  print("NOT downsampling mouse DNMS")
  
  spectradf_excludingLindsay <- spectradf_includingLindsay[spectradf_includingLindsay$label!="lindsay_mice",]
  
  # first need to get frac segregating sites:
  spectradf_excludingLindsay <- spectradf_excludingLindsay %>% 
    group_by(species,population,label) %>%
    mutate(fractionOfSegregatingSites_RESCALED_BY_HUMAN_TARGETS = total_mutations_RESCALED_BY_HUMAN_TARGETS / sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
    mutate(total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled=as.numeric(rmultinom(n=1,size=subsampleValue,prob=fractionOfSegregatingSites_RESCALED_BY_HUMAN_TARGETS))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
    ungroup() # ungroup just in case here 
    
    # then want to add back in the lindsay mice:
    lindsayMice <- spectradf_includingLindsay[spectradf_includingLindsay$label=="lindsay_mice",]
    lindsayMice$total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled <- lindsayMice$total_mutations_RESCALED_BY_HUMAN_TARGETS
    # want to get the frac seg sites to match above:
    lindsayMice <- lindsayMice %>%
      group_by(species,population,label) %>% # only one pop but doing for safety 
      mutate(fractionOfSegregatingSites_RESCALED_BY_HUMAN_TARGETS = total_mutations_RESCALED_BY_HUMAN_TARGETS / sum(total_mutations_RESCALED_BY_HUMAN_TARGETS)) %>%
      ungroup()
    
    
    
    spectradf <- bind_rows(spectradf_excludingLindsay,lindsayMice) # 
    
    return(spectradf)
}

##### ok!! These have now been downsampled !!!!!!!!! ###########
all1merSpectra <- downsampleFunc_RescaledByHumanTargets_ToVaquitaLevel_ExcludeLindsay(all1merSpectra_HumanRescaled)


### 20221220 : changing to now just use proportions
# get proportions after havng rescaled ot human genome comp and 
# need to do it this way bc not downsampling all the way down to lindsay dnm counts, so they have much lower rate if you divide by target size than polymorphism data

all1merSpectra <- all1merSpectra %>%
  group_by(species,population,label) %>%
  mutate(fracSegSites_BASEDON_total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled=total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled/sum(total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled)) %>%
  ungroup()
# check:
# sum(all1merSpectra[all1merSpectra$label=="brown_bear_ABC",]$fracSegSites_BASEDON_total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled) # equals 1

pivotSpectra_perVariable <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("label"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}

# Frac Seg Sites: 
all1merSpectra_fracSegSites_BASEDON_total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_WIDE <- pivotSpectra_perVariable(all1merSpectra,"fracSegSites_BASEDON_total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled") # this is the variable I want to compare , don't need the plus epsilon version bc 1mers

# check all columsn are 1:
#colSums(select(all1merSpectra_fracSegSites_BASEDON_total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_WIDE,-c('mutation_label','variable')))

# clr transformation: 
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


clrTransform_processAllSpectraTogether_fracSegSites_BASEDON_total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_WIDE = clr_or_ilr_orNoTransform_calculations(all1merSpectra_fracSegSites_BASEDON_total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_WIDE,"clr")


# get distances but keep upper part of triangle (so that you have all:all distances to plot -- wouldn't be good for a scatterplot but works for this way I'm doing plotting )
euclideanDistance_KeepUpperPartOfTriangleForLindsayComparison <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- tidy(dist(t(pivotedspectradf_noIDVars),upper=T,diag=F)) # transpose for dist otherwise it goes along rows and is wrong
  colnames(distances) <- c("item1","item2","spectrum_distance")
  return(distances)
}


distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled <- euclideanDistance_KeepUpperPartOfTriangleForLindsayComparison(clrTransform_processAllSpectraTogether_fracSegSites_BASEDON_total_mutations_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_WIDE)


####### colors and nice labels ########
source("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/colors_and_labels.plusUCLAwolves.R") # get cols
# gives you niceLabels
head(niceLabels)
# add lindsay to nice labels
niceLabels = bind_rows(niceLabels, data.frame(label="lindsay_mice",niceLabel="Lindsay mouse DNMs",shortSpeciesLabel="Lindsay mice"))

# item1 nice labels
distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_NICELABELS <- merge(distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled,niceLabels,by.x="item1",by.y="label")
# item2 nice labels
distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_NICELABELS2 <- merge(distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_NICELABELS,niceLabels,by.x="item2",by.y="label",suffixes = c(".item1",".item2"))

head(distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_NICELABELS2)

# add to color list
colorList$`Lindsay mouse DNMs` <- "tomato"
colorList$`wolf (Broad)` <- "darkgrey" # updating to be dark gray 

colorList
# species to include in plot:
# 
# only having mice/wolves on x-axis, but want to show all the comparisons on y axis: 
niceLabelitem2sToPlot=c( "Lindsay mouse DNMs","mouse (Mmd)","mouse (Mmc)","mouse (Mmm)" ,"mouse (Ms)" , "orangutan (Sumatran)" ,"orangutan (Bornean)"  , "gorilla" ,"human (AFR)" , "bonobo" ,"chimpanzee","fin whale (GOC)" ,  "vaquita" ,"brown bear (ABC)" , "polar bear","wolf (Broad)"  ,"wolf (UCLA)") # human_AFR# switching to GOC fin whale for now due to ksfs weirdness
#niceLabelitem1sToPlot=c("Lindsay mouse DNMs","mouse (Mmd)","mouse (Mmc)","mouse (Mmm)" ,"mouse (Ms)","wolf (Broad)"  ,"wolf (UCLA)")

niceLabelitem1sToPlot=c("Lindsay mouse DNMs","mouse (Mmd)","mouse (Ms)","wolf (Broad)"  ,"wolf (UCLA)")

plot1 <-ggplot(distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_NICELABELS2[distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_NICELABELS2$niceLabel.item1 %in% niceLabelitem1sToPlot & distances_clrTransform_processAllSpectraTogether_wide_mutation_rate_multinom_downsampled_RESCALED_BY_HUMAN_TARGETS_noteLindsayMiceNotDownsampledButAreRescaled_NICELABELS2$niceLabel.item2 %in% niceLabelitem2sToPlot,],aes(x=niceLabel.item1,y=spectrum_distance,color=niceLabel.item2))+
  geom_text(aes(label=niceLabel.item2),size=3)+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("")+
  ylab("1-mer spectrum Aitchison Distance") +
  scale_color_manual(values=colorList)
plot1
ggsave(paste0(plotdir,"ComparingDistances.OurMice.LindsayMice.SpeciesDownsampledTOVaquitaLevel.LindsayNotDownsamp.RESCALEDTOHUMANTARGETS.UsingNonProjectedData.FracSegSites.usethis.Mmd.Ms.Only.pdf"),plot1,height=8,width=12)
# okay so based on my clr proof doesn't matter if we use frac seg sites or mut rate since they cancel out (even if not downsampling)
# cool
