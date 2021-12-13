library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(biomaRt)

options(future.globals.maxSize = 4000 * 1024^3)

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

sch1.data <- Read10X(data.dir = "sciatic-v3-gfp/mpz-s1r")
#sch1.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/VS/raw_matrices/sch1_raw_feature_bc_matrix/raw_feature_bc_matrix")
#sch1.data
sch1 <- CreateSeuratObject(counts=sch1.data, project = "sch1", assay="RNA", min.cells = 10, min.features = 100)
rm(sch1.data)
#sch1
sch1@meta.data$tech <- "sch1"
sch1.mito.genes <- grep(pattern = "^mt-", x = rownames(x = sch1@assays$RNA@data), value = TRUE)
sch1.percent.mito <- Matrix::colSums(sch1@assays$RNA@counts[sch1.mito.genes, ])/Matrix::colSums(sch1@assays$RNA@counts)
sch1 <- AddMetaData(object = sch1, metadata = sch1.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=sch1, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=sch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sch1 <- subset(sch1,subset=nFeature_RNA>1000 & nFeature_RNA<7500 & percent.mito<0.05 & percent.mito>0)
sch1 <- RenameCells(sch1,add.cell.id="s1r")
sch1 <- SCTransform(sch1,verbose=FALSE)
sch1

sch2.data <- Read10X(data.dir = "sciatic-v3-gfp/mpz-s2")
#sch2.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/VS/raw_matrices/sch2_raw_feature_bc_matrix/raw_feature_bc_matrix")
#sch2.data
sch2 <- CreateSeuratObject(counts=sch2.data, project = "sch2", assay="RNA", min.cells = 10, min.features = 100)
rm(sch2.data)
#sch2
sch2@meta.data$tech <- "sch2"
sch2.mito.genes <- grep(pattern = "^mt-", x = rownames(x = sch2@assays$RNA@data), value = TRUE)
sch2.percent.mito <- Matrix::colSums(sch2@assays$RNA@counts[sch2.mito.genes, ])/Matrix::colSums(sch2@assays$RNA@counts)
sch2 <- AddMetaData(object = sch2, metadata = sch2.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=sch2, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=sch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sch2 <- subset(sch2,subset=nFeature_RNA>1000 & nFeature_RNA<7500 & percent.mito<0.05 & percent.mito>0)
sch2 <- RenameCells(sch2,add.cell.id="s2")
sch2 <- SCTransform(sch2,verbose=FALSE)
sch2

addlist <- list("SCH1","SCH2")
ob.list <- list(sch1,sch2)
assay.list <- c("RNA","RNA")

for (i in 1:length(ob.list)) {
  ob.list[[i]] <- NormalizeData(ob.list[[i]], verbose = FALSE)
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]], selection.method = "vst",
                                             nfeatures = 2000, verbose = FALSE)
  #ob.list[[i]] <- SCTransform(ob.list[[i]])
}

rm(sch.anchors)
sch.anchors <- SelectIntegrationFeatures(object.list = ob.list, verbose=TRUE)
sch.anchors
'Prx' %in% sch.anchors

all_genes <- lapply(ob.list, row.names) %>% Reduce(intersect, .)

schwannoma.list <- PrepSCTIntegration(object.list = ob.list, anchor.features = sch.anchors, verbose = TRUE)
schwannoma.anchors <- FindIntegrationAnchors(object.list = schwannoma.list, normalization.method = "SCT", anchor.features = sch.anchors, verbose = TRUE)
schwannoma.integrated <- IntegrateData(anchorset = schwannoma.anchors, normalization.method = "SCT", features.to.integrate = all_genes, verbose = TRUE)
schwannoma.integrated <- RunPCA(schwannoma.integrated, verbose = FALSE)
schwannoma.integrated <- FindNeighbors(schwannoma.integrated, dims = 1:20)
schwannoma.integrated <- FindClusters(schwannoma.integrated, resolution = 0.5)
schwannoma.integrated <- RunUMAP(schwannoma.integrated, dims = 1:20)
schwannoma.integrated <- RunTSNE(schwannoma.integrated, dims.use = 1:30)

p1 <- DimPlot(schwannoma.integrated, reduction = "umap", group.by = "tech")
#p2 <- DimPlot(schwannoma.integrated,reduction = "umap", group.by = "Phase")
p2 <- DimPlot(schwannoma.integrated, reduction = "umap", label = TRUE,
              repel = TRUE) + NoLegend()
plot_grid(p1, p2)

DimPlot(schwannoma.integrated, reduction = "umap", group.by = "tech", pt.size = 0.3)
DimPlot(schwannoma.integrated, reduction = "umap", label = TRUE, pt.size = 0.3,
        repel = TRUE) + NoLegend()

DimPlot(schwannoma.integrated,reduction = "umap", group.by = "tech",pt.size=0.01)
immunegenes = c("Prx","Scn7a","Pdgfra","Cd34","Mki67","Ngfr","Fabp4")
immunegenes = c("CCR2","LYVE1","TREM2","MKI67")
FeaturePlot(schwannoma.integrated, features = immunegenes, cols = c("steelblue","yellow","red"), min.cutoff = 0, max.cutoff = 4,pt.size=0.01,reduction="umap")
p <- FeaturePlot(schwannoma.integrated, features = immunegenes, cols = c("steelblue","yellow","red"), min.cutoff = 0, max.cutoff = 4,pt.size=0.01,reduction="umap",combine=FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)

sciaticv3_neuralcrest <- subset(schwannoma.integrated, ident=c(0,1,2,3,4,8,7,6))
sciaticv3_neuralcrest <- RunPCA(sciaticv3_neuralcrest, verbose = FALSE)
sciaticv3_neuralcrest <- FindNeighbors(sciaticv3_neuralcrest, dims = 1:15)
sciaticv3_neuralcrest <- FindClusters(sciaticv3_neuralcrest, resolution = 0.2)
sciaticv3_neuralcrest <- RunUMAP(sciaticv3_neuralcrest, dims = 1:15)
#sciaticv3_neuralcrest <- RunTSNE(sciaticv3_neuralcrest, dims.use = 1:30)

cell.cycle.tirosh <- read.csv("http://genomedata.org/rnaseq-tutorial/scrna/CellCycleTiroshSymbol2ID.csv", header=TRUE); # read in the list of genes
cell.cycle.tirosh
s.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G1/S")]; # create a vector of S-phase genes
g2m.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G2/M")]; # create a vector of G2/M-phase genes
m.s.genes <- convertHumanGeneList(s.genes)
m.g2m.genes <- convertHumanGeneList(g2m.genes)
sciaticv3_neuralcrest <- CellCycleScoring(object=sciaticv3_neuralcrest, s.features=m.s.genes, g2m.features=m.g2m.genes, set.ident=FALSE)

glycolysis <- read.csv('sciatic-v3-gfp/seurat3_merged/pathwaySymbol.csv',header=TRUE)
glycolysis$Gene
sciaticv3_neuralcrest <- AddModuleScore(object=sciaticv3_neuralcrest,features=list(glycolysis$Gene),ctrl=5,nbin=50,name="CAMs")
names(x = sciaticv3_neuralcrest[[]])
FeaturePlot(object = sciaticv3_neuralcrest, features = 'CAMs1', cols = c("grey","red"), min.cutoff = 0, max.cutoff = 2 )

FindMarkers(sciaticv3_neuralcrest,ident.1="5")
sciaticv3_neuralcrest

p1 <- DimPlot(sciaticv3_neuralcrest, reduction = "umap", group.by = "tech")
p2 <- DimPlot(sciaticv3_neuralcrest,reduction = "umap", group.by = "Phase")
p3 <- DimPlot(sciaticv3_neuralcrest, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
p4 <- FeaturePlot(object = sciaticv3_neuralcrest, features = 'glycolysis')
plot_grid(p1,p2,p3)

immunegenes = c("Prx","Scn7a","Pdgfra","Cd34","Mki67","Ngfr","Egr2","Lama2","Plp1")
schwanngenes = c("Prx","Scn7a","Pmp2","Adamtsl1","Cldn14","Egr2","Lama2","Lmna","Aldh1b1")
FeaturePlot(sciaticv3_neuralcrest, features = immunegenes, cols = c("steelblue","yellow","red"), min.cutoff = 0, max.cutoff = 4,pt.size=0.01,reduction="umap")

sciaticv3_neuralcrest_v2 <- subset(sciaticv3_neuralcrest, ident=c(0,1,2,3,4,6))
p1 <- DimPlot(sciaticv3_neuralcrest_v2, reduction = "tsne", group.by = "rep")
p2 <- DimPlot(sciaticv3_neuralcrest_v2,reduction = "tsne", group.by = "Phase")
p3 <- DimPlot(sciaticv3_neuralcrest_v2, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
p4 <- FeaturePlot(object = sciaticv3_neuralcrest_v2, features = 'glycolysis')
DimPlot(sciaticv3_neuralcrest_v2, reduction = "tsne")
plot_grid(p1,p2,p3)

sciaticv3_neuralcrest_v2@misc$all_genes <- rownames(sciaticv3_neuralcrest_v2)
sciaticv3_neuralcrest_v2@misc$DE$all <- FindAllMarkers(object = sciaticv3_neuralcrest_v2,
                                        only.pos = FALSE,
                                        min.pct = 0.1)
sciaticv3_neuralcrest_v2@misc$DE$all %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% as.data.frame() -> sciaticv3_neuralcrest_v2@misc$DE$top10
sciaticv3_neuralcrest_v2@misc$DE$all %>% group_by(cluster) %>% top_n(30, avg_logFC) %>% as.data.frame() -> sciaticv3_neuralcrest_v2@misc$DE$top30
sciaticv3_neuralcrest_v2@misc$meta.info$Title <- "sciatic v3 seurat3 analysis"
DimPlot(sciaticv3_neuralcrest_v2, reduction = "tsne", group.by = "tech", pt.size=2)

saveRDS(sciaticv3_neuralcrest_v2,file="sciatic-v3-gfp/seurat3_merged/sciaticv3_neuralcrestOnly_visualize.RDS")
sciaticv3_neuralcrest_v2 <- readRDS("sciatic-v3-gfp/aligned_cc_sciatic_nerves_v3-endofib-schwann-finalized_GFP-16Dec2019.RData")
sciaticv3_neuralcrest_v2 <- UpdateSeuratObject(sciaticv3_neuralcrest_v2)

sciaticv3_neuralcrest_v2


FeatureScatter(object=sciaticv3_neuralcrest_v2, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=sciaticv3_neuralcrest_v2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object=sciaticv3_neuralcrest_v2, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by = "tech") + scale_y_continuous(limits = c(0,30000)) + scale_x_continuous(limits = c(0,7000))

VlnPlot(sciaticv3_neuralcrest_v2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "rep", pt.size=0)


VlnPlot(sciaticv3_neuralcrest_v2, features = "percent.mito", ncol = 3, group.by = "tech") + scale_y_continuous(limits = c(0,0.05))
VlnPlot(sciaticv3_neuralcrest_v2, features = "nFeature_RNA", ncol = 3, group.by = "tech") + scale_y_continuous(limits = c(0,7000))
VlnPlot(sciaticv3_neuralcrest_v2, features = "nCount_RNA", ncol = 3, group.by = "tech") + scale_y_continuous(limits = c(0,30000))


sciaticv3_neuralcrest_v2@reductions$umap
Embeddings(sciaticv3_neuralcrest_v2[["tsne"]])

sciaticv3_neuralcrest_v2 <- RunPCA(sciaticv3_neuralcrest_v2, verbose = FALSE)
sciaticv3_neuralcrest_v2 <- FindNeighbors(sciaticv3_neuralcrest_v2, dims = 1:15)
sciaticv3_neuralcrest_v2 <- FindClusters(sciaticv3_neuralcrest_v2, resolution = 0.2)
sciaticv3_neuralcrest_v2 <- RunUMAP(sciaticv3_neuralcrest_v2, dims = 1:15)
