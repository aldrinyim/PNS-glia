library(Rmagic)
library(ggplot2)
library(readr)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(plyr)

options(future.globals.maxSize = 4000 * 1024^3)

removing_genes <- function(nuc.data) {
  pcdh.index <- grep(pattern = "^Pcdh", x = rownames(nuc.data), value = FALSE) # Select row indices
  zfp.index <- grep(pattern = "^Zfp", x = rownames(nuc.data), value = FALSE) # Select row indices
  ugt.index <- grep(pattern = "^Ugt", x = rownames(nuc.data), value = FALSE) # Select row indices
  rps.index <- grep(pattern = "^Rps", x = rownames(nuc.data), value = FALSE) # Select row indices
  rpl.index <- grep(pattern = "^Rpl", x = rownames(nuc.data), value = FALSE) # Select row indices
  zc3h11a.index <- grep(pattern = "^Zc3h11a", x = rownames(nuc.data), value = FALSE)

  total.index <- unlist(list(pcdh.index,zfp.index,ugt.index,rps.index,rpl.index,zc3h11a.index))
  new.data <- nuc.data[-total.index, ]
  return(new.data)
}

#On server
sn1.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/SN1")
sn2.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/SN2")
sn3.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/SN3")
snc1.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/SN_C1")
snc2.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/SN_C2")
s1r.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/mpz-s1r")
s2.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/mpz-s2")
vagus1.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/vagus1/direct_expression")
vagus2.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/vagus2/direct_expression")
sural.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/sural1April/direct_expression")
peroneal.data <- Read10X(data.dir = "/scratch/jmlab/aldrinyim/schwanncell/Peroneal/direct_expression")

#sn1.data <- Read10X(data.dir = "sciatic-v2/SN1/")
sn1.data <- removing_genes(sn1.data)
sn1 <- CreateSeuratObject(counts=sn1.data, project = "SN1", assay="RNA", min.cells = 3, min.features = 50)
rm(sn1.data)
sn1
sn1@meta.data$rep <- "sn1"
sn1@meta.data$bio.rep <- "Sciatic"
sn1@meta.data$tech <- "v2"
sn1.mito.genes <- grep(pattern = "^mt-", x = rownames(x = sn1@assays$RNA@data), value = TRUE)
sn1.percent.mito <- Matrix::colSums(sn1@assays$RNA@counts[sn1.mito.genes, ])/Matrix::colSums(sn1@assays$RNA@counts)
sn1 <- AddMetaData(object = sn1, metadata = sn1.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=sn1, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=sn1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sn1 <- subset(sn1,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(sn1@meta.data$nFeature_RNA, 0.98)
count99
sn1 <- subset(sn1,nFeature_RNA < count99)
sn1 <- RenameCells(sn1,add.cell.id="sn1_")
sn1 <- SCTransform(sn1,verbose=FALSE, vars.to.regress = c("percent.mito"))
sn1

#sn2.data <- Read10X(data.dir = "sciatic-v2/SN2/")
sn2.data <- removing_genes(sn2.data)
sn2 <- CreateSeuratObject(counts=sn2.data, project = "sn2", assay="RNA", min.cells = 3, min.features = 50)
rm(sn2.data)
sn2
sn2@meta.data$rep <- "sn2"
sn2@meta.data$bio.rep <- "Sciatic"
sn2@meta.data$tech <- "v2"
sn2.mito.genes <- grep(pattern = "^mt-", x = rownames(x = sn2@assays$RNA@data), value = TRUE)
sn2.percent.mito <- Matrix::colSums(sn2@assays$RNA@counts[sn2.mito.genes, ])/Matrix::colSums(sn2@assays$RNA@counts)
sn2 <- AddMetaData(object = sn2, metadata = sn2.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=sn2, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=sn2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sn2 <- subset(sn2,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(sn2@meta.data$nFeature_RNA, 0.98)
count99
sn2 <- subset(sn2,nFeature_RNA < count99)
sn2 <- RenameCells(sn2,add.cell.id="sn2_")
sn2 <- SCTransform(sn2,verbose=FALSE, vars.to.regress = c("percent.mito"))
sn2

#sn3.data <- Read10X(data.dir = "sciatic-v2/SN3/")
sn3.data <- removing_genes(sn3.data)
sn3 <- CreateSeuratObject(counts=sn3.data, project = "sn3", assay="RNA", min.cells = 3, min.features = 50)
rm(sn3.data)
sn3
sn3@meta.data$rep <- "sn3"
sn3@meta.data$bio.rep <- "Sciatic"
sn3@meta.data$tech <- "v2"
sn3.mito.genes <- grep(pattern = "^mt-", x = rownames(x = sn3@assays$RNA@data), value = TRUE)
sn3.percent.mito <- Matrix::colSums(sn3@assays$RNA@counts[sn3.mito.genes, ])/Matrix::colSums(sn3@assays$RNA@counts)
sn3 <- AddMetaData(object = sn3, metadata = sn3.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=sn3, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=sn3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sn3 <- subset(sn3,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(sn3@meta.data$nFeature_RNA, 0.98)
count99
sn3 <- subset(sn3,nFeature_RNA < count99)
sn3 <- RenameCells(sn3,add.cell.id="sn3_")
sn3 <- SCTransform(sn3,verbose=FALSE, vars.to.regress = c("percent.mito"))
sn3

#snc1.data <- Read10X(data.dir = "sciatic-v2/SN_C1/")
snc1.data <- removing_genes(snc1.data)
snc1 <- CreateSeuratObject(counts=snc1.data, project = "snc1", assay="RNA", min.cells = 3, min.features = 50)
rm(snc1.data)
snc1
snc1@meta.data$rep <- "snc1"
snc1@meta.data$bio.rep <- "Sciatic"
snc1@meta.data$tech <- "v2"
snc1.mito.genes <- grep(pattern = "^mt-", x = rownames(x = snc1@assays$RNA@data), value = TRUE)
snc1.percent.mito <- Matrix::colSums(snc1@assays$RNA@counts[snc1.mito.genes, ])/Matrix::colSums(snc1@assays$RNA@counts)
snc1 <- AddMetaData(object = snc1, metadata = snc1.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=snc1, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=snc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
snc1 <- subset(snc1,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(snc1@meta.data$nFeature_RNA, 0.98)
count99
snc1 <- subset(snc1,nFeature_RNA < count99)
snc1 <- RenameCells(snc1,add.cell.id="snc1_")
snc1 <- SCTransform(snc1,verbose=FALSE, vars.to.regress = c("percent.mito"))
snc1

#snc2.data <- Read10X(data.dir = "sciatic-v2/SN_C2/")
snc2.data <- removing_genes(snc2.data)
snc2 <- CreateSeuratObject(counts=snc2.data, project = "snc2", assay="RNA", min.cells = 3, min.features = 50)
rm(snc2.data)
snc2
snc2@meta.data$rep <- "snc2"
snc2@meta.data$bio.rep <- "Sciatic"
snc2@meta.data$tech <- "v2"
snc2.mito.genes <- grep(pattern = "^mt-", x = rownames(x = snc2@assays$RNA@data), value = TRUE)
snc2.percent.mito <- Matrix::colSums(snc2@assays$RNA@counts[snc2.mito.genes, ])/Matrix::colSums(snc2@assays$RNA@counts)
snc2 <- AddMetaData(object = snc2, metadata = snc2.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=snc2, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=snc2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
snc2 <- subset(snc2,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(snc2@meta.data$nFeature_RNA, 0.98)
count99
snc2 <- subset(snc2,nFeature_RNA < count99)
snc2 <- RenameCells(snc2,add.cell.id="snc2_")
snc2 <- SCTransform(snc2,verbose=FALSE, vars.to.regress = c("percent.mito"))
snc2

#s1r.data <- Read10X(data.dir = "sciatic-v3-gfp/mpz-s1r/")
s1r.data <- removing_genes(s1r.data)
s1r <- CreateSeuratObject(counts=s1r.data, project = "s1r", assay="RNA", min.cells = 10, min.features = 50)
rm(s1r.data)
s1r
s1r@meta.data$rep <- "s1r"
s1r@meta.data$bio.rep <- "Sciatic"
s1r@meta.data$tech <- "v3"
s1r.mito.genes <- grep(pattern = "^mt-", x = rownames(x = s1r@assays$RNA@data), value = TRUE)
s1r.percent.mito <- Matrix::colSums(s1r@assays$RNA@counts[s1r.mito.genes, ])/Matrix::colSums(s1r@assays$RNA@counts)
s1r <- AddMetaData(object = s1r, metadata = s1r.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=s1r, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=s1r, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
s1r <- subset(s1r,subset=nFeature_RNA>500 & nFeature_RNA<8000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(s1r@meta.data$nFeature_RNA, 0.98)
count99
s1r <- subset(s1r,nFeature_RNA < count99)
s1r <- RenameCells(s1r,add.cell.id="s1r_")
s1r <- SCTransform(s1r,verbose=FALSE, vars.to.regress = c("percent.mito"))
s1r

#s2.data <- Read10X(data.dir = "sciatic-v3-gfp/mpz-s2/")
s2.data <- removing_genes(s2.data)
s2 <- CreateSeuratObject(counts=s2.data, project = "s2", assay="RNA", min.cells = 10, min.features = 50)
rm(s2.data)
s2
s2@meta.data$rep <- "s2"
s2@meta.data$bio.rep <- "Sciatic"
s2@meta.data$tech <- "v3"
s2.mito.genes <- grep(pattern = "^mt-", x = rownames(x = s2@assays$RNA@data), value = TRUE)
s2.percent.mito <- Matrix::colSums(s2@assays$RNA@counts[s2.mito.genes, ])/Matrix::colSums(s2@assays$RNA@counts)
s2 <- AddMetaData(object = s2, metadata = s2.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=s2, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=s2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
s2 <- subset(s2,subset=nFeature_RNA>500 & nFeature_RNA<8000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(s2@meta.data$nFeature_RNA, 0.98)
count99
s2 <- subset(s2,nFeature_RNA < count99)
s2 <- RenameCells(s2,add.cell.id="s2_")
s2 <- SCTransform(s2,verbose=FALSE, vars.to.regress = c("percent.mito"))
s2

vagus1.data <- Read10X(data.dir = "Vagus/vagus1/direct_expression")
vagus1.data <- removing_genes(vagus1.data)
vagus1 <- CreateSeuratObject(counts=vagus1.data, project = "vagus1", assay="RNA", min.cells = 3, min.features = 50)
rm(vagus1.data)
vagus1
vagus1@meta.data$rep <- "vagus1"
vagus1@meta.data$bio.rep <- "Sciatic"
vagus1@meta.data$tech <- "v2"
vagus1.mito.genes <- grep(pattern = "^mt-", x = rownames(x = vagus1@assays$RNA@data), value = TRUE)
vagus1.percent.mito <- Matrix::colSums(vagus1@assays$RNA@counts[vagus1.mito.genes, ])/Matrix::colSums(vagus1@assays$RNA@counts)
vagus1 <- AddMetaData(object = vagus1, metadata = vagus1.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=vagus1, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=vagus1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
vagus1 <- subset(vagus1,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(vagus1@meta.data$nFeature_RNA, 0.98)
count99
vagus1 <- subset(vagus1,nFeature_RNA < count99)
vagus1 <- RenameCells(vagus1,add.cell.id="vagus1_")
vagus1 <- SCTransform(vagus1,verbose=FALSE, vars.to.regress = c("percent.mito"))
vagus1

vagus2.data <- Read10X(data.dir = "Vagus/vagus2/direct_expression")
vagus2.data <- removing_genes(vagus2.data)
vagus2 <- CreateSeuratObject(counts=vagus2.data, project = "vagus2", assay="RNA", min.cells = 3, min.features = 50)
rm(vagus2.data)
vagus2
vagus2@meta.data$rep <- "vagus2"
vagus2@meta.data$bio.rep <- "Sciatic"
vagus2@meta.data$tech <- "v2"
vagus2.mito.genes <- grep(pattern = "^mt-", x = rownames(x = vagus2@assays$RNA@data), value = TRUE)
vagus2.percent.mito <- Matrix::colSums(vagus2@assays$RNA@counts[vagus2.mito.genes, ])/Matrix::colSums(vagus2@assays$RNA@counts)
vagus2 <- AddMetaData(object = vagus2, metadata = vagus2.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=vagus2, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=vagus2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
vagus2 <- subset(vagus2,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(vagus2@meta.data$nFeature_RNA, 0.98)
count99
vagus2 <- subset(vagus2,nFeature_RNA < count99)
vagus2 <- RenameCells(vagus2,add.cell.id="vagus2_")
vagus2 <- SCTransform(vagus2,verbose=FALSE, vars.to.regress = c("percent.mito"))
vagus2

#sural.data <- Read10X(data.dir = "sural1April/direct_expression")
sural.data <- removing_genes(sural.data)
sural <- CreateSeuratObject(counts=sural.data, project = "sural", assay="RNA", min.cells = 3, min.features = 50)
rm(sural.data)
sural
sural@meta.data$rep <- "sural"
sural@meta.data$bio.rep <- "Sciatic"
sural@meta.data$tech <- "v2"
sural.mito.genes <- grep(pattern = "^mt-", x = rownames(x = sural@assays$RNA@data), value = TRUE)
sural.percent.mito <- Matrix::colSums(sural@assays$RNA@counts[sural.mito.genes, ])/Matrix::colSums(sural@assays$RNA@counts)
sural <- AddMetaData(object = sural, metadata = sural.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=sural, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=sural, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sural <- subset(sural,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(sural@meta.data$nFeature_RNA, 0.98)
count99
sural <- subset(sural,nFeature_RNA < count99)
sural <- RenameCells(sural,add.cell.id="sural_")
sural <- SCTransform(sural,verbose=FALSE, vars.to.regress = c("percent.mito"))
sural

#peroneal.data <- Read10X(data.dir = "Peroneal/direct_expression")
peroneal.data <- removing_genes(peroneal.data)
peroneal <- CreateSeuratObject(counts=peroneal.data, project = "peroneal", assay="RNA", min.cells = 3, min.features = 50)
rm(peroneal.data)
peroneal
peroneal@meta.data$rep <- "peroneal"
peroneal@meta.data$bio.rep <- "Sciatic"
peroneal@meta.data$tech <- "v2"
peroneal.mito.genes <- grep(pattern = "^mt-", x = rownames(x = peroneal@assays$RNA@data), value = TRUE)
peroneal.percent.mito <- Matrix::colSums(peroneal@assays$RNA@counts[peroneal.mito.genes, ])/Matrix::colSums(peroneal@assays$RNA@counts)
peroneal <- AddMetaData(object = peroneal, metadata = peroneal.percent.mito, col.name = "percent.mito")
par(mfrow = c(1, 2))
FeatureScatter(object=peroneal, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object=peroneal, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
peroneal <- subset(peroneal,subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mito<0.05 & percent.mito>0)
count99 <- quantile(peroneal@meta.data$nFeature_RNA, 0.98)
count99
peroneal <- subset(peroneal,nFeature_RNA < count99)
peroneal <- RenameCells(peroneal,add.cell.id="peroneal_")
peroneal <- SCTransform(peroneal,verbose=FALSE, vars.to.regress = c("percent.mito"))
peroneal

allnerves.list <- list(sn1,sn2,sn3,snc1,snc2,s1r,s2,vagus1,vagus2,sural,peroneal)

all_genes <- lapply(allnerves.list, row.names) %>% Reduce(intersect, .)
all_genes

allnerves.features <- SelectIntegrationFeatures(object.list = allnerves.list, nfeatures = 3000)
allnerves.list <- PrepSCTIntegration(object.list = allnerves.list, anchor.features = allnerves.features,
                                   verbose = FALSE)
allnerves.anchors <- FindIntegrationAnchors(object.list = allnerves.list, normalization.method = "SCT", anchor.features = allnerves.features, verbose = TRUE)
allnerves.integrated <- IntegrateData(anchorset = allnerves.anchors, normalization.method = "SCT", features.to.integrate = all_genes, verbose = TRUE)

allnerves.integrated <- RunPCA(allnerves.integrated, verbose = FALSE)
allnerves.integrated <- FindNeighbors(allnerves.integrated, dims = 1:10)
allnerves.integrated <- FindClusters(allnerves.integrated, resolution = 0.5)
allnerves.integrated <- RunUMAP(allnerves.integrated, dims = 1:10)
allnerves.integrated <- RunTSNE(allnerves.integrated, dims.use = 1:8)

saveRDS(allnerves.integrated,"allnerves_v23_mpz-act-regressMito-merged_22Jan2021.RDS")

magic_sciatic.integrated <- magic(sciatic.integrated,genes=sciatic.features)

saveRDS(magic_sciatic.integrated,"magic_sciatic_v23_mpz-act-regressMito-merged_22Jan2021.RDS")
