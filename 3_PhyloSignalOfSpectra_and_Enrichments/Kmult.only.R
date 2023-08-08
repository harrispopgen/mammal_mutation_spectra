###### experimenting with kmult ######
#require(devtools)
#devtools::install_github("geomorphR/geomorph", ref = "Stable") # try without build vignettes (doesn't work if trying to build vignettes)

# had to install rgl using:
#remotes::install_github("dmurdoch/rgl") rather than cran rgl which caused aborts

require(geomorph)
require(ape)
require(dplyr)
require(ggplot2)
require(reshape2)
require(tidyr)
require(tidyverse)
require(ggrepel)

outdir="/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/kmult/"
niter=999 # for kmult - run into vector mem issues if go higher than this.
###### get tree ########
# process raxml tree here instead of in previous 
raxmlTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")

speciesToInclude_in_raxmlTree=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Canis_lupus") 

# want to get rid of the family/order info at end of tip. just keep species : 
raxmlTree_renamedTips <- raxmlTree
raxmlTree_renamedTips$tip.label <- paste0(unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",1)),"_",unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",2)))

# restrict to just the species to include
raxmlTree_renamedTips_subset <- keep.tip(raxmlTree_renamedTips,speciesToInclude_in_raxmlTree)


########### read in all spectra with all types of transformations ############
# human-scaled downsampled spectra:
all_spectra_all_conditions <- read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/results/compare_spectrumDist_toPhyloHammingDist/20221121.plots.SubsampledIndividualsNotProjection.epsilon.1.sepCpGs.no.maskedmaskALL.RESCALEDTOHUMANTARGETS.SimplifiedCode.9999999.MantelPermutations/AllProcessedSpectra.AllConditions.AllVariables.AllSpecies.RESCALED_BY_HUMAN_TARGETS.DOWNSAMPLED.PROCESSED.txt",header=T,sep="\t")

# need a matrix with rows = species, cols = traits (mutatio types)

# and need a phylo tree

# see variables:
unique(all_spectra_all_conditions$variable)

kmult_results <- list()
listOfIDs=c("1-mer spectrum","3-mer spectrum","5-mer spectrum","7-mer spectrum")

speciesToInclude_phylo=c("Mus_musculus_Mmd","Mus_spretus_Ms","brown_bear_ABC","polar_bear_PB","vaquita","fin_whale_GOC","humans_AFR","Gorilla_gorilla","Pan_paniscus","Pan_troglodytes","Pongo_abelii","Pongo_pygmaeus","wolves")


speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/multispecies_spectra_comparisons/scripts/multispecies_spectra/analyses/compare_spectrumDist_toPhyloHammingDist/speciesCodes.txt",header = T)
head(speciesCodes) # for getting phylo distances 

# restricting just to clr and mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS
for(id in listOfIDs){

print(id)
spectrum_mutrate_clr = all_spectra_all_conditions[all_spectra_all_conditions$id==id & all_spectra_all_conditions$variable=="mutation_rate_multinom_downsampled_plusEpsilon_RESCALED_BY_HUMAN_TARGETS" & all_spectra_all_conditions$transform_id=="clr",]

head(spectrum_mutrate_clr)
# print the dim out
dim(spectrum_mutrate_clr)

# reshape into a matrix where rows are species and cols are variables 
# row names must match tree tip names

spectrum_mutrate_clr_melt <- melt(spectrum_mutrate_clr,variable.name = "label")

spectrum_mutrate_clr_melt <- pivot_wider(spectrum_mutrate_clr_melt,values_from = value,id_cols =label,names_from = mutation_label )

# need to restrict to just species to include and then get their species names to match tree
spectrum_mutrate_clr_melt_subset <- spectrum_mutrate_clr_melt %>%
  filter(label %in% speciesToInclude_phylo)

# add in full spp names:
spectrum_mutrate_clr_melt_subset_codes <- merge(spectrum_mutrate_clr_melt_subset,speciesCodes[,c("code","species")],by.x="label",by.y="code")

# get rid of label and convert species to rownames 
spectrum_mutrate_clr_melt_subset_codes_MATRIXFORPHYSIGNAL <- spectrum_mutrate_clr_melt_subset_codes %>%
  column_to_rownames(var="species") %>%
  select(-"label") %>%
  as.matrix()

kmultResult <- physignal(spectrum_mutrate_clr_melt_subset_codes_MATRIXFORPHYSIGNAL, raxmlTree_renamedTips_subset, iter = niter, seed = NULL, print.progress = TRUE) # 

# note kmult performs lots of checks to make sure data matrix and tree are correctly formatted / contain right set of species
# https://github.com/geomorphR/geomorph/blob/Stable/R/physignal.r

kmult_results[[id]] <- kmultResult



}

kmult_results_summary=data.frame(id=listOfIDs,pvalue=NA,kmult=NA)


for(id in listOfIDs){
  
kmult_results_summary[kmult_results_summary$id==id,]$pvalue <- kmult_results[[id]]$pvalue
kmult_results_summary[kmult_results_summary$id==id,]$kmult <- kmult_results[[id]]$phy.signal

}

kmult_results_summary$significant <- "NS"
kmult_results_summary[kmult_results_summary$pvalue<0.01,]$significant <- "yes"
kmult_results_summary$pvalue_label <- paste0("p = ",kmult_results_summary$pvalue)
kmultplot <- ggplot(kmult_results_summary,aes(x=id,y=kmult))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  geom_text_repel(aes(label=pvalue_label),size=2)+
  ylab("Kmult")+
  xlab("")

ggsave(paste0(outdir,"kmultplot.niter.",niter,".pdf"),kmultplot,height=4,width=7)
ggsave(paste0(outdir,"kmultplot.niter.",niter,".png"),kmultplot,height=4,width=7)


write.table(kmult_results_summary,paste0(outdir,"kmultResults.",niter,".permutations.txt"),row.names = F,sep="\t",quote=F)

save.image(file = paste0(outdir,"RWORKSPACE.AtEndOfScript.image.RData"))



