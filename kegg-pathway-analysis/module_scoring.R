library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(biomaRt)
library(plyr)

options(future.globals.maxSize = 4000 * 1024^3)

updated_schwannv23_subset <- readRDS("~/glia-RDS-submission/schwann_cell_single-nuclei-atlas-12Jan2021.RDS")

#Calculate module score for all KEGG pathways
list_of_files <- read.csv('~/Schwann cell v23/expression_modules/index.txt',header=FALSE)
length(list_of_files[[1]])
list_of_files[[1]]

for (i in 1:length(list_of_files[[1]])){
  print(i)
  full_file <- paste('~/Schwann cell v23/expression_modules/',list_of_files[[1]][i],sep='')
  print(full_file)
  pathwaydata <- read.csv(full_file,header=TRUE)
  pathwaydata$Gene
  map_name <- unlist(pathwaydata$Pathway[[1]])
  map_name
  updated_schwannv23_subset <- AddModuleScore(object=updated_schwannv23_subset,features=list(pathwaydata$Gene),ctrl=5,nbin=50,name=map_name)
}
updated_schwannv23_subset

#Calculate average module score per cluster*
rm(report_results)
report_results <- list()
typeof(updated_schwannv23_subset)
list_of_maps <- read.csv('~/Schwann cell v23/expression_modules/index_name.txt',header=FALSE)
for (i in 1:length(list_of_maps[[1]])){
#for (i in 1:10){
  #print(i)
  map_name <- list_of_maps[[1]][i]
  print(map_name)
  mean_sig_clusters <-c() # define empty obj
  for (ident in levels(Idents(object = updated_schwannv23_subset))) {
    mean_sig_clusters <- c(mean_sig_clusters, mean(updated_schwannv23_subset@meta.data[[map_name]][which(colnames(x = updated_schwannv23_subset) %in% WhichCells(updated_schwannv23_subset,ident=ident))]))
  }
  #print(mean_sig_clusters)
  normalized_zscore <- scale(mean_sig_clusters)[,1]
  #print(normalized_zscore)
  print(c(map_name,mean_sig_clusters,normalized_zscore))
  report_results[i] <- list(c(map_name,mean_sig_clusters,normalized_zscore))
}

report_results[[1]]

#Individual test of pathways
mean_sig_clusters <- c() # define empty obj
for (ident in levels(Idents(object = updated_schwannv23_subset))) {
  print(ident)
  mean_sig_clusters <- c(mean_sig_clusters, mean(updated_schwannv23_subset@meta.data$glycolysis.gluconeogenesis1[which(colnames(x = updated_schwannv23_subset) %in%  WhichCells(updated_schwannv23_subset,ident=ident))]))
}
print(mean_sig_clusters)
normalized_zscore <- scale(mean_sig_clusters)[,1]
list(summary(normalized_zscore))[[1]]
print(normalized_zscore)

rm(report_results)
report_results <- list()
report_results[1] <- list(c('axon.guidance1',mean_sig_clusters,normalized_zscore))
report_results[2] <- list(c('axon.guidance2',mean_sig_clusters,normalized_zscore))
structure(report_results)

list(c('axon.guidance1',mean_sig_clusters,normalized_zscore))[[1]]

report_results

write.csv(report_results, file = "test.csv", quote=FALSE)

tweets.df = ldply(report_results, function(t) t$toDataFrame())

lapply(list(c('axon.guidance1',mean_sig_clusters,normalized_zscore)), function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))

report_df <- as.data.frame()

FeaturePlot(object = updated_schwannv23_subset, features = 'nicotinate.and.nicotinamide.metabolism1', reduction = "tsne", cols = c("steelblue","yellow","red"), min.cutoff = 0, pt.size =0.8)
FeaturePlot(object = updated_schwannv23_subset, features = 'fatty.acid.biosynthesis1', reduction = "tsne", cols = c("grey","red"),min.cutoff = 0, pt.size=1)

saveRDS(updated_schwannv23_subset,"~/utilities/seurat2-env/sciatic_v23/schwann_cell_v23_combined_aligned_finalized_PRESHINY_12Jan2021.RDS")
