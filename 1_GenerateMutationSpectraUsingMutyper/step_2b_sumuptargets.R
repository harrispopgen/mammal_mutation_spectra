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


######### sum up targets ########
# note for species that don't have intervals, the interval count is set as "allautos" (vaquita) so this still works fine. 

sumuptargets <- function(species,intervals,indir) {
  allTargets=data.frame() # must clear between species 
  for(interval in intervals){
    input=read.table(paste0(indir,species,".int_or_chr_",interval,inputfilesuffix))
    colnames(input) <- c("target","count")
    input$interval <- as.character(interval)
    input_melt <- melt(input,id.vars = c("target","interval"))
    input_melt$species <- species # safety net
    allTargets <- bind_rows(allTargets,data.frame(input_melt))
  }
  
  # sum up over chrs:
  allTargetsSummedOverIntervals <- allTargets %>%
    group_by(target,species) %>% # added insurance; should just be one species though
    summarise(countOverAllIntervals=sum(value)) %>%
    select(target,countOverAllIntervals) 
  
  return(allTargetsSummedOverIntervals)
}


##### unified function ####
unifiedFunction_sumuptargets <- function(species,intervalCount,prepend0,indir){
  intervals=getIntervals(intervalCount,prepend0)
  allTargetsSummedOverIntervals = sumuptargets(species,intervals,indir)
  return(allTargetsSummedOverIntervals)
}

##### run it #######
allTargetsSummedOverIntervals=unifiedFunction_sumuptargets(species,intervalCount,prepend0,indir)

##### write it out #######
write.table(allTargetsSummedOverIntervals,paste0(outdir,species,".summedup",inputfilesuffix),row.names = F,quote=F,sep="\t")
