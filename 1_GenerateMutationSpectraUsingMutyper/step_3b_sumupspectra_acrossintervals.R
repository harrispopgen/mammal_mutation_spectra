# want to run for each species separately, but do pops within the script
require(reshape2)
require(dplyr)

args = commandArgs(trailingOnly=TRUE) # supply args
# first arg will be species
# second arg will be ? 

indir=paste0(as.character(args[1]),"/")
outdir=paste0(as.character(args[2]) ,"/")
species=as.character(args[3])
intervalCount=as.character(args[4])
prepend0=as.character(args[5])
inputfilesuffix=as.character(args[6])

print(paste0("starting ", species))
print(paste0("interval count:", intervalCount))

######## going to get the pops from the dirs that are in the indir (if no dirs then it's no pops) #######


######## get the number of intervals ########
getIntervals <- function(intervalCount,prepend0) {
  
  if(intervalCount=="NA"){
    intervals="allautos"
  } else {
    
    if(prepend0=="FALSE") {
      intervals=as.character(c(seq(1,intervalCount))) # careful to make these characters when looping 
    } else if(prepend0=="TRUE") {
      intervals=as.character(sprintf("%02d", seq(1,intervalCount))) # if you need leading 0s 
    } else {
      stop("not a valid prepend0 option")
    }
  }
  return(intervals)
}

########## get the populations #######
getPopsFromIndir <- function(indir) {
  poplist=list.dirs(indir,full.names=F,recursive = F) 
  if(length(poplist)==0){
    print("no populations, using all inds")
    poplist=NA
  } else {
    print("using pops:")
    print(poplist)
  }
  return(poplist)
}


######### sum up spectra (two kinds of spectra need to be added up: per ind and per pop) ########
# note for species that don't have intervals, the interval count is set as "allautos" (vaquita) so this still works fine. 
sumupspectraacrossintervalsANDpopulations <- function(species,poplist,intervals,indir,inputfilesuffix,outdir) {
  if(!anyNA(poplist)) { # checks if poplist is NA (contains NA)
    for(pop in poplist){
      perIndSpectra_summedup = data.frame()
      perPopSpectra_summedup  = data.frame()
      print(paste0("starting ",pop))
      for(interval in intervals){
        print(interval)
        popdir=paste0(indir,"/",pop,"/") # this will differ for other species annoying
        perindspectrum=read.table(paste0(popdir,species,".",pop,".int_or_chr_",interval,inputfilesuffix,".SUBSETPRIORTOMUTYPERSPECTRUM.PERINDIVIDUAL.RANDOMIZEDWITHINPOPULATION.txt"),header=T) #  --randomize
        perpopspectrum=read.table(paste0(popdir,species,".",pop,".int_or_chr_",interval,inputfilesuffix,".SUBSETPRIORTOMUTYPERSPECTRUM.PERPOPULATION.txt"),header=T) # --population
        perindspectrum_melt <- melt(perindspectrum,id.vars = c("sample"))
        perpopspectrum_melt <- melt(perpopspectrum)
        colnames(perindspectrum_melt) <- c("sample","mutation_type","totalSites") # CHECK THIS!
        colnames(perpopspectrum_melt) <- c("mutation_type","totalSites")
        
        perpopspectrum_melt$population <- pop
        perpopspectrum_melt$species <- species
        perindspectrum_melt$population <- pop
        perindspectrum_melt$species <- species
        
        # combine with previous sum and then sum up again: 
        perIndSpectra_summedup=bind_rows(perIndSpectra_summedup,perindspectrum_melt)
        perPopSpectra_summedup=bind_rows(perPopSpectra_summedup,perpopspectrum_melt)
        
        perIndSpectra_summedup <- perIndSpectra_summedup %>% 
          group_by(mutation_type,population,species,sample) %>%  # need to group by populations
          summarise(totalSites=sum(totalSites)) %>%
          ungroup()
        
        perPopSpectra_summedup <- perPopSpectra_summedup %>% 
          group_by(mutation_type,population,species) %>%  # need to group by population
          summarise(totalSites=sum(totalSites)) %>%
          ungroup()
        
      }
      write.table(perIndSpectra_summedup,paste0(outdir,species,"_",pop,".summedup",inputfilesuffix,".SUBSETPRIORTOMUTYPERSPECTRUM.PERINDIVIDUAL.RANDOMIZEDWITHINPOPULATION.txt"),row.names = F,quote=F,sep="\t")
      write.table(perPopSpectra_summedup,paste0(outdir,species,"_",pop,".summedup",inputfilesuffix,".SUBSETPRIORTOMUTYPERSPECTRUM.PERPOPULATION.txt"),row.names = F,quote=F,sep="\t")
      
    }
  } else {
    perIndSpectra_summedup = data.frame()
    perPopSpectra_summedup  = data.frame()
    for(interval in intervals){
      print(interval)
      perindspectrum=read.table(paste0(indir,species,".int_or_chr_",interval,inputfilesuffix,".SUBSETPRIORTOMUTYPERSPECTRUM.PERINDIVIDUAL.RANDOMIZEDWITHINPOPULATION.txt"),header=T) #  --randomize
      perpopspectrum=read.table(paste0(indir,species,".int_or_chr_",interval,inputfilesuffix,".SUBSETPRIORTOMUTYPERSPECTRUM.PERPOPULATION.txt"),header=T) # --population
      perindspectrum_melt <- melt(perindspectrum,id.vars = c("sample"))
      #head(perindspectrum_melt)
      perpopspectrum_melt <- melt(perpopspectrum)
      #head(perpopspectrum_melt)
      colnames(perindspectrum_melt) <- c("sample","mutation_type","totalSites") # CHECK THIS!
      #print(head(perindspectrum_melt))
      colnames(perpopspectrum_melt) <- c("mutation_type","totalSites")
      #print(head(perpopspectrum_melt))
      
      perpopspectrum_melt$population <- species
      perpopspectrum_melt$species <- species
      perindspectrum_melt$population <- species
      perindspectrum_melt$species <- species
      
      # combine with previous sum and then sum up again: 
      perIndSpectra_summedup=bind_rows(perIndSpectra_summedup,perindspectrum_melt)
      perPopSpectra_summedup=bind_rows(perPopSpectra_summedup,perpopspectrum_melt)
      
      perIndSpectra_summedup <- perIndSpectra_summedup %>% 
        group_by(mutation_type,population,species,sample) %>%  # need to group by populations
        summarise(totalSites=sum(totalSites)) %>%
        ungroup()
      
      perPopSpectra_summedup <- perPopSpectra_summedup %>% 
        group_by(mutation_type,population,species) %>%  # need to group by population
        summarise(totalSites=sum(totalSites)) %>%
        ungroup()
      
    }
    write.table(perIndSpectra_summedup,paste0(outdir,species,".summedup",inputfilesuffix,".SUBSETPRIORTOMUTYPERSPECTRUM.PERINDIVIDUAL.RANDOMIZEDWITHINPOPULATION.txt"),row.names = F,quote=F,sep="\t")
    write.table(perPopSpectra_summedup,paste0(outdir,species,".summedup",inputfilesuffix,".SUBSETPRIORTOMUTYPERSPECTRUM.PERPOPULATION.txt"),row.names = F,quote=F,sep="\t")
    
  }
}
##### unified function ####
unifiedFunction_sumupspectra <- function(species,intervalCount,prepend0,indir,inputfilesuffix,outdir){
  intervals=getIntervals(intervalCount,prepend0)
  poplist=getPopsFromIndir(indir)
  # sum up and write out: 
  sumupspectraacrossintervalsANDpopulations(species,poplist,intervals,indir,inputfilesuffix,outdir)
}

##### run it (will write out) #######
unifiedFunction_sumupspectra(species,intervalCount,prepend0,indir,inputfilesuffix,outdir)


