### Take the data aligned with cell ranger and cleaned with cell bender,
### and remove doublets before continuing.

## Load packages, re/install as needed -----------------------------------------

library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(devtools)
library(hdf5r)
library(sctransform)
library(ggalluvial)
library(ggrepel)
BiocManager::install("plger/scDblFinder")
library(scDblFinder)

#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

## This subsection was written out by Mike DeBerardine -------------------------
# A function to do scDblFinder doublet detection on a Seurat object (x)
#
# Adds columns to Seurat metadata:
# - scDblFinder.class
# - scDblFinder.cxds_score
# - scDblFinder.score
# - scDblFinder.weighted
findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}
# ------------------------------------------------------------------------------

## For each one of the 20 samples:
# set the working directory
# read in the .h5 file from cell bender
# create a seurat object
# calculate the % mitochondrial
# create and print violin plots before doublet removal
# remove doublets and save the object

#Object 1/20: Infanticidal Mus 1 
setwd("Z:/Forrest/snRNAseq/rawreads2/Mus_Male_Infanticidal_1/")
mi1.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
mi1 <- CreateSeuratObject(counts = mi1.cb.data.h5, project = "MInfa1", min.cells = 10, min.features = 200)
mi1[["percent.mt"]] <- PercentageFeatureSet(mi1, pattern = "^mt-")
vlnplot.mi1 <- VlnPlot(mi1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_mi1.png", width = 1080, height = 1080, units = "px")
vlnplot.mi1
dev.off()
#mi1 <- subset(mi1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

mi1.1 <- findDoublets(mi1)

# to see how many doublets were detected
table(mi1.1$scDblFinder.class)

# to subset your seurat object to only singlets
mi1.1 <- subset(mi1.1, scDblFinder.class == "singlet")

mi1.1[["percent.mt"]] <- PercentageFeatureSet(mi1.1, pattern = "^mt-")
vlnplot.mi1.1 <- VlnPlot(mi1.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.mi1.1

#mi1.2 <- subset(mi1.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.mi1.2 <- VlnPlot(mi1.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.mi1.2

saveRDS(mi1.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/mi1.RDS")
rm(list = ls())

#Object 2/20: Infanticidal Mus 2
setwd("Z:/Forrest/snRNAseq/rawreads2/Mus_Male_Infanticidal_2/")
mi2.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
mi2 <- CreateSeuratObject(counts = mi2.cb.data.h5, project = "MInfa2", min.cells = 10, min.features = 200)
mi2[["percent.mt"]] <- PercentageFeatureSet(mi2, pattern = "^mt-")
vlnplot.mi2 <- VlnPlot(mi2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.mi2
#mi2 <- subset(mi2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
mi2.1 <- findDoublets(mi2)
# to see how many doublets were detected
table(mi2.1$scDblFinder.class)
# to subset your seurat object to only singlets
mi2.1 <- subset(mi2.1, scDblFinder.class == "singlet")
mi2.1[["percent.mt"]] <- PercentageFeatureSet(mi2.1, pattern = "^mt-")
vlnplot.mi2.1 <- VlnPlot(mi2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.mi2.1
#mi2.2 <- subset(mi2.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.mi2.2 <- VlnPlot(mi2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.mi2.2
saveRDS(mi2.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/mi2.RDS")
rm(list = ls())

#Object 3/20: Infanticidal Mus 3
setwd("Z:/Forrest/snRNAseq/rawreads2/Mus_Male_Infanticidal_3/")
mi3.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
mi3 <- CreateSeuratObject(counts = mi3.cb.data.h5, project = "MInfa3", min.cells = 10, min.features = 200)
mi3[["percent.mt"]] <- PercentageFeatureSet(mi3, pattern = "^mt-")
vlnplot.mi3 <- VlnPlot(mi3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.mi3
#mi3 <- subset(mi3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
mi3.1 <- findDoublets(mi3)
# to see how many doublets were detected
table(mi3.1$scDblFinder.class)
# to subset your seurat object to only singlets
mi3.1 <- subset(mi3.1, scDblFinder.class == "singlet")
mi3.1[["percent.mt"]] <- PercentageFeatureSet(mi3.1, pattern = "^mt-")
vlnplot.mi3.1 <- VlnPlot(mi3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.mi3.1
#mi3.2 <- subset(mi3.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.mi3.2 <- VlnPlot(mi3.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.mi3.2
saveRDS(mi3.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/mi3.RDS")
rm(list = ls())

#Object 4/20: Infanticidal Rhab 1
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Infanticidal_1/")
ri1.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
ri1 <- CreateSeuratObject(counts = ri1.cb.data.h5, project = "RInfa1", min.cells = 10, min.features = 200)
ri1[["percent.mt"]] <- PercentageFeatureSet(ri1, pattern = "^mt-")
vlnplot.ri1 <- VlnPlot(ri1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ri1
#ri1 <- subset(ri1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
ri1.1 <- findDoublets(ri1)
# to see how many doublets were detected
table(ri1.1$scDblFinder.class)
# to subset your seurat object to only singlets
ri1.1 <- subset(ri1.1, scDblFinder.class == "singlet")
ri1.1[["percent.mt"]] <- PercentageFeatureSet(ri1.1, pattern = "^mt-")
vlnplot.ri1.1 <- VlnPlot(ri1.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ri1.1
#ri1.2 <- subset(ri1.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.ri1.2 <- VlnPlot(ri1.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.ri1.2
saveRDS(ri1.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/ri1.RDS")
rm(list = ls())

#Object 5/20: Infanticidal Rhab 2
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Infanticidal_2/")
ri2.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
ri2 <- CreateSeuratObject(counts = ri2.cb.data.h5, project = "RInfa2", min.cells = 10, min.features = 200)
ri2[["percent.mt"]] <- PercentageFeatureSet(ri2, pattern = "^mt-")
vlnplot.ri2 <- VlnPlot(ri2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ri2
#ri2 <- subset(ri2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
ri2.1 <- findDoublets(ri2)
# to see how many doublets were detected
table(ri2.1$scDblFinder.class)
# to subset your seurat object to only singlets
ri2.1 <- subset(ri2.1, scDblFinder.class == "singlet")
ri2.1[["percent.mt"]] <- PercentageFeatureSet(ri2.1, pattern = "^mt-")
vlnplot.ri2.1 <- VlnPlot(ri2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ri2.1
#ri2.2 <- subset(ri2.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.ri2.2 <- VlnPlot(ri2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.ri2.2
saveRDS(ri2.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/ri2.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 6/20: Infanticidal Rhab 3
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Infanticidal_3/")
ri3.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
ri3 <- CreateSeuratObject(counts = ri3.cb.data.h5, project = "RInfa3", min.cells = 10, min.features = 200)
ri3[["percent.mt"]] <- PercentageFeatureSet(ri3, pattern = "^mt-")
vlnplot.ri3 <- VlnPlot(ri3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ri3
#ri3 <- subset(ri3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
ri3.1 <- findDoublets(ri3)
# to see how many doublets were detected
table(ri3.1$scDblFinder.class)
# to subset your seurat object to only singlets
ri3.1 <- subset(ri3.1, scDblFinder.class == "singlet")
ri3.1[["percent.mt"]] <- PercentageFeatureSet(ri3.1, pattern = "^mt-")
vlnplot.ri3.1 <- VlnPlot(ri3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ri3.1
#ri3.2 <- subset(ri3.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.ri3.2 <- VlnPlot(ri3.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.ri3.2
saveRDS(ri3.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/ri3.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 7/20: Control Rhab 1
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Control_1/")
rc1.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rc1 <- CreateSeuratObject(counts = rc1.cb.data.h5, project = "RCont1", min.cells = 10, min.features = 200)
rc1[["percent.mt"]] <- PercentageFeatureSet(rc1, pattern = "^mt-")
vlnplot.rc1 <- VlnPlot(rc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rc1
#rc1 <- subset(rc1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rc1.1 <- findDoublets(rc1)
# to see how many doublets were detected
table(rc1.1$scDblFinder.class)
# to subset your seurat object to only singlets
rc1.1 <- subset(rc1.1, scDblFinder.class == "singlet")
rc1.1[["percent.mt"]] <- PercentageFeatureSet(rc1.1, pattern = "^mt-")
vlnplot.rc1.1 <- VlnPlot(rc1.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rc1.1
#rc1.2 <- subset(rc1.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rc1.2 <- VlnPlot(rc1.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rc1.2
saveRDS(rc1.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rc1.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 8/20: Control Rhab 2
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Control_2/")
rc2.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rc2 <- CreateSeuratObject(counts = rc2.cb.data.h5, project = "RCont2", min.cells = 10, min.features = 200)
rc2[["percent.mt"]] <- PercentageFeatureSet(rc2, pattern = "^mt-")
vlnplot.rc2 <- VlnPlot(rc2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rc2
#rc2 <- subset(rc2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rc2.1 <- findDoublets(rc2)
# to see how many doublets were detected
table(rc2.1$scDblFinder.class)
# to subset your seurat object to only singlets
rc2.1 <- subset(rc2.1, scDblFinder.class == "singlet")
rc2.1[["percent.mt"]] <- PercentageFeatureSet(rc2.1, pattern = "^mt-")
vlnplot.rc2.1 <- VlnPlot(rc2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rc2.1
#rc2.2 <- subset(rc2.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rc2.2 <- VlnPlot(rc2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rc2.2
saveRDS(rc2.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rc2.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 9/20: Control Rhab 3
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Control_3/")
rc3.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rc3 <- CreateSeuratObject(counts = rc3.cb.data.h5, project = "RCont3", min.cells = 10, min.features = 200)
rc3[["percent.mt"]] <- PercentageFeatureSet(rc3, pattern = "^mt-")
vlnplot.rc3 <- VlnPlot(rc3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rc3
#rc3 <- subset(rc3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rc3.1 <- findDoublets(rc3)
# to see how many doublets were detected
table(rc3.1$scDblFinder.class)
# to subset your seurat object to only singlets
rc3.1 <- subset(rc3.1, scDblFinder.class == "singlet")
rc3.1[["percent.mt"]] <- PercentageFeatureSet(rc3.1, pattern = "^mt-")
vlnplot.rc3.1 <- VlnPlot(rc3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rc3.1
#rc3.2 <- subset(rc3.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rc3.2 <- VlnPlot(rc3.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rc3.2
saveRDS(rc3.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rc3.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 10/20: Allopaternal Rhab 1
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Allopat_1/")
ra1.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
ra1 <- CreateSeuratObject(counts = ra1.cb.data.h5, project = "RAllo1", min.cells = 10, min.features = 200)
ra1[["percent.mt"]] <- PercentageFeatureSet(ra1, pattern = "^mt-")
vlnplot.ra1 <- VlnPlot(ra1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ra1
#ra1 <- subset(ra1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
ra1.1 <- findDoublets(ra1)
# to see how many doublets were detected
table(ra1.1$scDblFinder.class)
# to subset your seurat object to only singlets
ra1.1 <- subset(ra1.1, scDblFinder.class == "singlet")
ra1.1[["percent.mt"]] <- PercentageFeatureSet(ra1.1, pattern = "^mt-")
vlnplot.ra1.1 <- VlnPlot(ra1.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ra1.1
#ra1.2 <- subset(ra1.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.ra1.2 <- VlnPlot(ra1.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.ra1.2
saveRDS(ra1.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/ra1.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 11/20: Allopaternal Rhab 1
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Allopat_2/")
ra2.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
ra2 <- CreateSeuratObject(counts = ra2.cb.data.h5, project = "RAllo2", min.cells = 10, min.features = 200)
ra2[["percent.mt"]] <- PercentageFeatureSet(ra2, pattern = "^mt-")
vlnplot.ra2 <- VlnPlot(ra2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ra2
#ra2 <- subset(ra2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
ra2.1 <- findDoublets(ra2)
# to see how many doublets were detected
table(ra2.1$scDblFinder.class)
# to subset your seurat object to only singlets
ra2.1 <- subset(ra2.1, scDblFinder.class == "singlet")
ra2.1[["percent.mt"]] <- PercentageFeatureSet(ra2.1, pattern = "^mt-")
vlnplot.ra2.1 <- VlnPlot(ra2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ra2.1
#ra2.2 <- subset(ra2.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.ra2.2 <- VlnPlot(ra2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.ra2.2
saveRDS(ra2.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/ra2.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 12/20: Allopaternal Rhab 1
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Allopat_3/")
ra3.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
ra3 <- CreateSeuratObject(counts = ra3.cb.data.h5, project = "RAllo3", min.cells = 10, min.features = 200)
ra3[["percent.mt"]] <- PercentageFeatureSet(ra3, pattern = "^mt-")
vlnplot.ra3 <- VlnPlot(ra3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ra3
#ra3 <- subset(ra3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
ra3.1 <- findDoublets(ra3)
# to see how many doublets were detected
table(ra3.1$scDblFinder.class)
# to subset your seurat object to only singlets
ra3.1 <- subset(ra3.1, scDblFinder.class == "singlet")
ra3.1[["percent.mt"]] <- PercentageFeatureSet(ra3.1, pattern = "^mt-")
vlnplot.ra3.1 <- VlnPlot(ra3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.ra3.1
#ra3.2 <- subset(ra3.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.ra3.2 <- VlnPlot(ra3.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.ra3.2
saveRDS(ra3.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/ra3.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 13/20: Sires Rhab 1
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Sire_1/")
rs1.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rs1 <- CreateSeuratObject(counts = rs1.cb.data.h5, project = "RSire1", min.cells = 10, min.features = 200)
rs1[["percent.mt"]] <- PercentageFeatureSet(rs1, pattern = "^mt-")
vlnplot.rs1 <- VlnPlot(rs1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rs1
#rs1 <- subset(rs1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rs1.1 <- findDoublets(rs1)
# to see how many doublets were detected
table(rs1.1$scDblFinder.class)
# to subset your seurat object to only singlets
rs1.1 <- subset(rs1.1, scDblFinder.class == "singlet")
rs1.1[["percent.mt"]] <- PercentageFeatureSet(rs1.1, pattern = "^mt-")
vlnplot.rs1.1 <- VlnPlot(rs1.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rs1.1
#rs1.2 <- subset(rs1.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rs1.2 <- VlnPlot(rs1.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rs1.2
saveRDS(rs1.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rs1.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 14/20: Sires Rhab 2
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Sire_2/")
rs2.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rs2 <- CreateSeuratObject(counts = rs2.cb.data.h5, project = "RSire2", min.cells = 10, min.features = 200)
rs2[["percent.mt"]] <- PercentageFeatureSet(rs2, pattern = "^mt-")
vlnplot.rs2 <- VlnPlot(rs2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rs2
#rs2 <- subset(rs2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rs2.1 <- findDoublets(rs2)
# to see how many doublets were detected
table(rs2.1$scDblFinder.class)
# to subset your seurat object to only singlets
rs2.1 <- subset(rs2.1, scDblFinder.class == "singlet")
rs2.1[["percent.mt"]] <- PercentageFeatureSet(rs2.1, pattern = "^mt-")
vlnplot.rs2.1 <- VlnPlot(rs2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rs2.1
#rs2.2 <- subset(rs2.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rs2.2 <- VlnPlot(rs2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rs2.2
saveRDS(rs2.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rs2.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 15/20: Sires Rhab 3
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Sire_3/")
rs3.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rs3 <- CreateSeuratObject(counts = rs3.cb.data.h5, project = "RSire3", min.cells = 10, min.features = 200)
rs3[["percent.mt"]] <- PercentageFeatureSet(rs3, pattern = "^mt-")
vlnplot.rs3 <- VlnPlot(rs3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rs3
#rs3 <- subset(rs3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rs3.1 <- findDoublets(rs3)
# to see how many doublets were detected
table(rs3.1$scDblFinder.class)
# to subset your seurat object to only singlets
rs3.1 <- subset(rs3.1, scDblFinder.class == "singlet")
rs3.1[["percent.mt"]] <- PercentageFeatureSet(rs3.1, pattern = "^mt-")
vlnplot.rs3.1 <- VlnPlot(rs3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rs3.1
#rs3.2 <- subset(rs3.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rs3.2 <- VlnPlot(rs3.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rs3.2
saveRDS(rs3.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rs3.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 16/20: Sires Rhab 4
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Male_Sire_4/")
rs4.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rs4 <- CreateSeuratObject(counts = rs4.cb.data.h5, project = "RSire4", min.cells = 10, min.features = 200)
rs4[["percent.mt"]] <- PercentageFeatureSet(rs4, pattern = "^mt-")
vlnplot.rs4 <- VlnPlot(rs4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rs4
#rs4 <- subset(rs4, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rs4.1 <- findDoublets(rs4)
# to see how many doublets were detected
table(rs4.1$scDblFinder.class)
# to subset your seurat object to only singlets
rs4.1 <- subset(rs4.1, scDblFinder.class == "singlet")
rs4.1[["percent.mt"]] <- PercentageFeatureSet(rs4.1, pattern = "^mt-")
vlnplot.rs4.1 <- VlnPlot(rs4.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rs4.1
#rs4.2 <- subset(rs4.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rs4.2 <- VlnPlot(rs4.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rs4.2
saveRDS(rs4.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rs4.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 17/20: Dam Rhab 1
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Female_Dam_1/")
rd1.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rd1 <- CreateSeuratObject(counts = rd1.cb.data.h5, project = "RDam1", min.cells = 10, min.features = 200)
rd1[["percent.mt"]] <- PercentageFeatureSet(rd1, pattern = "^mt-")
vlnplot.rd1 <- VlnPlot(rd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rd1
#rd1 <- subset(rd1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rd1.1 <- findDoublets(rd1)
# to see how many doublets were detected
table(rd1.1$scDblFinder.class)
# to subset your seurat object to only singlets
rd1.1 <- subset(rd1.1, scDblFinder.class == "singlet")
rd1.1[["percent.mt"]] <- PercentageFeatureSet(rd1.1, pattern = "^mt-")
vlnplot.rd1.1 <- VlnPlot(rd1.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rd1.1
#rd1.2 <- subset(rd1.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rd1.2 <- VlnPlot(rd1.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rd1.2
saveRDS(rd1.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rd1.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 18/20: Dam Rhab 2
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Female_Dam_2/")
rd2.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rd2 <- CreateSeuratObject(counts = rd2.cb.data.h5, project = "RDam2", min.cells = 10, min.features = 200)
rd2[["percent.mt"]] <- PercentageFeatureSet(rd2, pattern = "^mt-")
vlnplot.rd2 <- VlnPlot(rd2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rd2
#rd2 <- subset(rd2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rd2.1 <- findDoublets(rd2)
# to see how many doublets were detected
table(rd2.1$scDblFinder.class)
# to subset your seurat object to only singlets
rd2.1 <- subset(rd2.1, scDblFinder.class == "singlet")
rd2.1[["percent.mt"]] <- PercentageFeatureSet(rd2.1, pattern = "^mt-")
vlnplot.rd2.1 <- VlnPlot(rd2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rd2.1
#rd2.2 <- subset(rd2.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rd2.2 <- VlnPlot(rd2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rd2.2
saveRDS(rd2.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rd2.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 19/20: Dam Rhab 3
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Female_Dam_3/")
rd3.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rd3 <- CreateSeuratObject(counts = rd3.cb.data.h5, project = "RDam3", min.cells = 10, min.features = 200)
rd3[["percent.mt"]] <- PercentageFeatureSet(rd3, pattern = "^mt-")
vlnplot.rd3 <- VlnPlot(rd3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rd3
#rd3 <- subset(rd3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rd3.1 <- findDoublets(rd3)
# to see how many doublets were detected
table(rd3.1$scDblFinder.class)
# to subset your seurat object to only singlets
rd3.1 <- subset(rd3.1, scDblFinder.class == "singlet")
rd3.1[["percent.mt"]] <- PercentageFeatureSet(rd3.1, pattern = "^mt-")
vlnplot.rd3.1 <- VlnPlot(rd3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rd3.1
#rd3.2 <- subset(rd3.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rd3.2 <- VlnPlot(rd3.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rd3.2
saveRDS(rd3.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rd3.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}

#Object 20/20: Dam Rhab 4
setwd("Z:/Forrest/snRNAseq/rawreads2/Rhab_Female_Dam_4/")
rd4.cb.data.h5 <- Read10X_h5("CellBendFilteredOut.h5")
rd4 <- CreateSeuratObject(counts = rd4.cb.data.h5, project = "RDam4", min.cells = 10, min.features = 200)
rd4[["percent.mt"]] <- PercentageFeatureSet(rd4, pattern = "^mt-")
vlnplot.rd4 <- VlnPlot(rd4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rd4
#rd4 <- subset(rd4, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
rd4.1 <- findDoublets(rd4)
# to see how many doublets were detected
table(rd4.1$scDblFinder.class)
# to subset your seurat object to only singlets
rd4.1 <- subset(rd4.1, scDblFinder.class == "singlet")
rd4.1[["percent.mt"]] <- PercentageFeatureSet(rd4.1, pattern = "^mt-")
vlnplot.rd4.1 <- VlnPlot(rd4.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnplot.rd4.1
#rd4.2 <- subset(rd4.1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#vlnplot.rd4.2 <- VlnPlot(rd4.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vlnplot.rd4.2
saveRDS(rd4.1, file = "Z:/Forrest/snRNAseq/DoubletRemovedRDS/rd4.RDS")
rm(list = ls())

findDoublets <- function(x) {
  suppressMessages(
    sce <- scDblFinder(x[["RNA"]]$counts, verbose = FALSE)
  )
  AddMetaData(x, as.data.frame(colData(sce)))
}