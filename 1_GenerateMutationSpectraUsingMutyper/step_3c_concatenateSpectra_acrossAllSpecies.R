################## Concatenate all spectra ##############
# want to end up with format:
# species population label mutation_type totalSites_allFreqsSummed

require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)



kmersize="7"
maskLabel="maskALL" # if species get different mask labels this will be a pain.
wd=paste0("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/allspecies_mutyper_results_unified_SubsettingIndividuals/",kmersize,"_mer/")

outdir=paste0(wd,"/concatenated_pop_level_spectra_allSpecies_basedOnSubsetOfIndividuals/masked_",maskLabel,"/including_ucla_wolves_comparison/") ## note different outdir from normal version of script
dir.create(outdir,showWarnings = F,recursive = T)
# note: don't want to write anything to allspecies_summed_up_over_intervals_forTransfer because if I redownload things I don't want to get overwritten
filesuffix=".summedup.mutyper.spectrum.SeeLogForFilters.maskALL.7mer.SUBSETPRIORTOMUTYPERSPECTRUM.PERPOPULATION.txt" # I manually gzipped them
### will need to update these lists! # CHANGE THIS 
#speciesList=c('humans', 'mice', 'bears', 'fin_whale', 'vaquita', 'Gorilla_gorilla', 'Pan_troglodytes', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus','wolves')
speciesList=c('humans', 'Mus_musculus','Mus_spretus' ,'polar_bear','brown_bear', 'fin_whale', 'vaquita', 'Gorilla_gorilla', 'Pan_troglodytes', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus','wolves','ucla_wolves') # as of 20220718 the ucla wolves had new polarization

popList=list(brown_bear=c("ABC","EUR"),polar_bear=c("PB"),fin_whale=c("ENP","GOC"),humans=c("AFR","EUR","SAS","EAS","AMR"),Mus_musculus=c("Mmc","Mmd","Mmm"),Mus_spretus=c("Ms"))

all_spectra <- data.frame()
for(species in speciesList) {
  indir=paste0(wd,"allspecies_summed_up_over_intervals_forTransfer/",species,"/mutyper_results_masked_",maskLabel,"/mutyper_spectrum_files/")
  pops=popList[[species]]
  print(pops)
  if(is.null(pops)) {
    spectrum <- read.table(paste0(indir,species,filesuffix),header=T)
    spectrum$label <- species
    all_spectra <- bind_rows(all_spectra,spectrum)
  } else {
    for(pop in pops){
      print(pop)
      spectrum <- read.table(paste0(indir,species,"_",pop,filesuffix),header=T)
      spectrum$label <- paste0(spectrum$species,"_",spectrum$pop)
      all_spectra <- bind_rows(all_spectra,spectrum)
    }
  }
}

write.table(all_spectra,paste0(outdir,"allSpectra.PopulationLevel.BasedOn5individualsPerPop.notprojected.usethis.txt"),quote=F,sep="\t",row.names = F)
