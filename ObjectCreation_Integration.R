### Take the seurat objects with doublets removed and create 
### an integrated data set. 

## Load packages --------------------------------------------------------------
library(dplyr)
library(Seurat)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("harmony", version = "3.19")
library(harmony)
## Read in the individual objects ---------------------------------------------
setwd(dir = "/Volumes/Crucial X6/sequencing/DoubletRemovedRDS/")

#mi1 <- readRDS(file = "mi1.RDS")
#mi2 <- readRDS(file = "mi2.RDS")
#mi3 <- readRDS(file = "mi3.RDS")

ri1 <- readRDS(file = "ri1.RDS")
ri2 <- readRDS(file = "ri2.RDS")
ri3 <- readRDS(file = "ri3.RDS")
rc1 <- readRDS(file = "rc1.RDS")
rc2 <- readRDS(file = "rc2.RDS")
rc3 <- readRDS(file = "rc3.RDS")
ra1 <- readRDS(file = "ra1.RDS")
ra2 <- readRDS(file = "ra2.RDS")
ra3 <- readRDS(file = "ra3.RDS")
rs1 <- readRDS(file = "rs1.RDS")
rs2 <- readRDS(file = "rs2.RDS")
rs3 <- readRDS(file = "rs3.RDS")
rs4 <- readRDS(file = "rs4.RDS")
rd1 <- readRDS(file = "rd1.RDS")
rd2 <- readRDS(file = "rd2.RDS")
rd3 <- readRDS(file = "rd3.RDS")
rd4 <- readRDS(file = "rd4.RDS")

## Make violin plots of feature, count, and pct mito distributions ------------
vlnplot.mi1 <- VlnPlot(mi1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_mi1.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.mi1
dev.off()

vlnplot.mi2 <- VlnPlot(mi2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_mi2.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.mi2
dev.off()

vlnplot.mi3 <- VlnPlot(mi3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_mi3.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.mi3
dev.off()

vlnplot.ri1 <- VlnPlot(ri1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_ri1.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.ri1
dev.off()

vlnplot.ri2 <- VlnPlot(ri2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_ri2.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.ri2
dev.off()

vlnplot.ri3 <- VlnPlot(ri3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_ri3.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.ri3
dev.off()

vlnplot.rc1 <- VlnPlot(rc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rc1.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rc1
dev.off()

vlnplot.rc2 <- VlnPlot(rc2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rc2.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rc2
dev.off()

vlnplot.rc3 <- VlnPlot(rc3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rc3.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rc3
dev.off()

vlnplot.ra1 <- VlnPlot(ra1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_ra1.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.ra1
dev.off()

vlnplot.ra2 <- VlnPlot(ra2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_ra2.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.ra2
dev.off()

vlnplot.ra3 <- VlnPlot(ra3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_ra3.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.ra3
dev.off()

vlnplot.rs1 <- VlnPlot(rs1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rs1.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rs1
dev.off()

vlnplot.rs2 <- VlnPlot(rs2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rs2.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rs2
dev.off()

vlnplot.rs3 <- VlnPlot(rs3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rs3.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rs3
dev.off()

vlnplot.rs4 <- VlnPlot(rs4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rs4.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rs4
dev.off()

vlnplot.rd1 <- VlnPlot(rd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rd1.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rd1
dev.off()

vlnplot.rd2 <- VlnPlot(rd2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rd2.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rd2
dev.off()

vlnplot.rd3 <- VlnPlot(rd3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rd3.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rd3
dev.off()

vlnplot.rd4 <- VlnPlot(rd4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(filename = 
      "Z:/Forrest/snRNAseq/PostSequencingAnalysis/Plots/PostDoubletRemoval_ViolinPlots/violinplot_feature_count_pctmt_rd4.png", 
    width = 1080, height = 1080, units = "px")
vlnplot.rd4
dev.off()

## There are overlapping barcodes, so we need to rename the cells with ids ----
#mi1 <- RenameCells(object = mi1, add.cell.id = "mi1")
#mi2 <- RenameCells(object = mi2, add.cell.id = "mi2")
#mi3 <- RenameCells(object = mi3, add.cell.id = "mi3")
ri1 <- RenameCells(object = ri1, add.cell.id = "ri1")
ri2 <- RenameCells(object = ri2, add.cell.id = "ri2")
ri3 <- RenameCells(object = ri3, add.cell.id = "ri3")
rc1 <- RenameCells(object = rc1, add.cell.id = "rc1")
rc2 <- RenameCells(object = rc2, add.cell.id = "rc2")
rc3 <- RenameCells(object = rc3, add.cell.id = "rc3")
ra1 <- RenameCells(object = ra1, add.cell.id = "ra1")
ra2 <- RenameCells(object = ra2, add.cell.id = "ra2")
ra3 <- RenameCells(object = ra3, add.cell.id = "ra3")
rs1 <- RenameCells(object = rs1, add.cell.id = "rs1")
rs2 <- RenameCells(object = rs2, add.cell.id = "rs2")
rs3 <- RenameCells(object = rs3, add.cell.id = "rs3")
rs4 <- RenameCells(object = rs4, add.cell.id = "rs4")
rd1 <- RenameCells(object = rd1, add.cell.id = "rd1")
rd2 <- RenameCells(object = rd2, add.cell.id = "rd2")
rd3 <- RenameCells(object = rd3, add.cell.id = "rd3")
rd4 <- RenameCells(object = rd4, add.cell.id = "rd4")

##make a simple R list of the objects and normalize/find HVG for each ---------

parenting_list <- list()
#parenting_list[["mi1"]] <- mi1 # exclude the mus for the integrated sample
#parenting_list[["mi2"]] <- mi2 # exclude the mus for the integrated sample
#parenting_list[["mi3"]] <- mi3 # exclude the mus for the integrated sample
parenting_list[["ri1"]] <- ri1
parenting_list[["ri2"]] <- ri2
parenting_list[["ri3"]] <- ri3
parenting_list[["rc1"]] <- rc1
parenting_list[["rc2"]] <- rc2
parenting_list[["rc3"]] <- rc3
parenting_list[["ra1"]] <- ra1
parenting_list[["ra2"]] <- ra2
parenting_list[["ra3"]] <- ra3
parenting_list[["rs1"]] <- rs1
parenting_list[["rs2"]] <- rs2
parenting_list[["rs3"]] <- rs3
parenting_list[["rs4"]] <- rs4
parenting_list[["rd1"]] <- rd1
parenting_list[["rd2"]] <- rd2
parenting_list[["rd3"]] <- rd3
parenting_list[["rd4"]] <- rd4

parenting <- merge(parenting_list[[1]], parenting_list[-1])

parenting <- NormalizeData(parenting, verbose = FALSE)
parenting <- FindVariableFeatures(parenting, verbose = FALSE)
parenting <- ScaleData(parenting, verbose = FALSE)
parenting <- RunPCA(parenting, verbose = FALSE)

ElbowPlot(parenting, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.25,1.5), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = 35, colour = "navy", linewidth = .5, linetype = 3)

parenting <- RunHarmony(parenting, "orig.ident")
harmony.embeddings <- Embeddings(parenting, reduction = "harmony")
parenting <- JoinLayers(parenting)

parenting <- RunUMAP(parenting, dims = 1:35, reduction = "harmony")
parenting <- parenting %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.6) 

saveRDS(parenting, file="/Volumes/Crucial X6/sequencing//parenting_harmonyintegrated_v2.RDS")

png(filename = "/Volumes/Crucial X6/sequencing/nFeature_nCount_byType.png", 
    units = "mm", width = 360, height = 120, pointsize = 12, bg = "white", res = 600)
VlnPlot(parenting, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, log = T)
dev.off()
