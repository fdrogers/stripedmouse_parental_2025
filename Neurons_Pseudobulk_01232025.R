### Take the integrated data set that was clustered at the "neuronal"
###cell types and subset into main gaba-ergic and glutamatergic populations.
### Run pseudobulk anaysis with DESeq2. 

## Load Packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(harmony)
library(DESeq2)
library(tidyverse)
library(cowplot)
library(Matrix)
library(reshape2)
library(S4Vectors)


#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("SingleCellExperiment")

library(SingleCellExperiment)
library(png)
library(RColorBrewer)
library(data.table)
library(patchwork)
library(glmGamPoi)
library(scuttle)
## Load in Data
neurons <- readRDS(file="/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")
DimPlot(neurons)

neurons$contact <- NA
neurons$contact[neurons$orig.ident %in% c("RDam1","RDam2","RDam3", "RDam4")] <- "Dams"
neurons$contact[neurons$orig.ident %in% c("RAllo1", "RAllo2", "RAllo3", "RSire1", "RSire2")] <- "HighMale"
neurons$contact[neurons$orig.ident %in% c("RInfa1","RInfa2","RInfa3", "RSire3", "RSire4")] <- "LowMale"
neurons$contact[neurons$orig.ident %in% c("RCont1","RCont2","RCont3")] <- "NoContact"

DotPlot(neurons, features = c("Npas4", "Fos", "Egr1", "Arc", "Egr3", "Fosl1", "Fosb", "Jun", "Junb"), 
        group.by = "group", split.by = "orig.ident", cols = c("#960018", "#960018", "#960018",
                                                              "#DCB50B", "#DCB50B", "#DCB50B",
                                                              "#00B9B9", "#00B9B9", "#00B9B9",
                                                              "#008DC5", "#008DC5", "#008DC5", "#008DC5",
                                                              "#BB2A7F", "#BB2A7F", "#BB2A7F", "#BB2A7F")) 

## Analyze all neurons
new.cluster.ids <- c("1", "1", "1", "1", "1", "1", "1", "1", "1", "1",
                     "1", "1", "1", "1", "1", "1", "1")
names(new.cluster.ids) <- levels(neurons)
neurons <- RenameIdents(neurons, new.cluster.ids)
counts <- neurons@assays$RNA
metadata <- neurons@meta.data
metadata$cluster_id <- factor(neurons@active.ident)
metadata$sample_id <- factor(metadata$orig.ident)

neur_bulk_counts <- neurons %>% 
  SplitObject(split.by = "orig.ident") %>%  
  lapply(function(x) rowSums(x@assays$RNA$counts)) %>% 
  do.call(cbind, .)

neur_bulk_counts


write.csv(neur_bulk_counts, "/Volumes/Crucial X6/sequencing/bulk_counts.csv")

####

pseudo <- AggregateExpression(neurons, assays = "RNA", 
                              return.seurat = T, 
                              group.by = c("group", "orig.ident"))

tail(Cells(pseudo))


pseudo$id.group <- paste(pseudo$group, sep = "_")

Idents(pseudo) <- "id.group"

bulk.de <- FindMarkers(object = pseudo, 
                            ident.1 = "RInfa", 
                            ident.2 = "RAllo",
                            test.use = "DESeq2")

View(bulk.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(bulk.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(bulk.de, file = "/Volumes/Crucial X6/sequencing/neurons_bulk_deseq2.csv")
bulk.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_bulk_deseq2.csv", header = TRUE)
View(bulk.de)
bulk.de[1,1] <- "Agouti" #rename a to Agouti for clarity 

# add a column of NAs
bulk.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
bulk.de$Expression[bulk.de$avg_log2FC > 0.01 & bulk.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
bulk.de$Expression[bulk.de$avg_log2FC < -0.01 & bulk.de$p_val_adj < 0.1] <- "Up in Allo"

bulk.de$delabel <- NA
bulk.de$delabel[bulk.de$Expression != "Not Diff"] <- bulk.de$X[bulk.de$Expression != "Not Diff"]

options(ggrepel.max.overlaps = Inf)

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons.png",
    width = 250, height = 300, units = "mm", res = 600)
ggplot(data=bulk.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(bulk.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#00B9B9", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + xlim(-2,2) + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 30),
        axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40),
        text = element_text(size = 40, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### Now for GABAergic
rm(list = ls())

## Load in Data
neurons <- readRDS(file = "/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")
print(levels(Idents(neurons)))

## Analyze all neurons
gaba <- subset(x = neurons, idents = c("gab1", "gab2", "gab3", "gab4", "gab5", "gab6", 
                                       "gab7", "gab8", "gab9", "gab10", "gab11", "gab12"))

DimPlot(gaba) #confirm labels

####

pseudo_gaba <- AggregateExpression(gaba, assays = "RNA", 
                              return.seurat = T, 
                              group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba))


pseudo_gaba$id.group <- paste(pseudo_gaba$group, sep = "_")

Idents(pseudo_gaba) <- "id.group"

gaba.de <- FindMarkers(object = pseudo_gaba, 
                       ident.1 = "RInfa", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(gaba.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba_deseq2.csv")
gaba.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba_deseq2.csv", header = TRUE)
gaba.de[1,1] <- "Agouti" #rename a to Agouti for clarity 

# add a column of NAs
gaba.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba.de$Expression[gaba.de$avg_log2FC > 0.01 & gaba.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba.de$Expression[gaba.de$avg_log2FC < -0.01 & gaba.de$p_val_adj < 0.1] <- "Up in Allo"

gaba.de$delabel <- NA
gaba.de$delabel[gaba.de$Expression != "Not Diff"] <- gaba.de$X[gaba.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_GABANeurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#00B9B9", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + xlim(-2,2) + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### Now for GLUTAmatergic
rm(list = ls())

## Load in Data
neurons <- readRDS(file = "/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")
print(levels(Idents(neurons)))

## Analyze all neurons
gluta <- subset(x = neurons, idents = c("glu1", "glu2", "glu3", "glu4"))

DimPlot(gluta) #confirm labels

####

pseudo_gluta <- AggregateExpression(gluta, assays = "RNA", 
                                   return.seurat = T, 
                                   group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gluta))


pseudo_gluta$id.group <- paste(pseudo_gluta$group, sep = "_")

Idents(pseudo_gluta) <- "id.group"

gluta.de <- FindMarkers(object = pseudo_gluta, 
                       ident.1 = "RInfa", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(gluta.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gluta.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gluta.de, file = "/Volumes/Crucial X6/sequencing/neurons_gluta_deseq2.csv")
gluta.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gluta_deseq2.csv", header = TRUE)
gluta.de[1,1] <- "Agouti" #rename a to Agouti for clarity 

# add a column of NAs
gluta.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gluta.de$Expression[gluta.de$avg_log2FC > 0.01 & gluta.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gluta.de$Expression[gluta.de$avg_log2FC < -0.01 & gluta.de$p_val_adj < 0.1] <- "Up in Allo"

gluta.de$delabel <- NA
gluta.de$delabel[gluta.de$Expression != "Not Diff"] <- gluta.de$X[gluta.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_GLUTANeurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gluta.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gluta.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#00B9B9", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + xlim(-2,2) + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

#### Specific clusters

## gab1
gaba1 <- subset(x = neurons, idents = c("gab1"))

DimPlot(gaba1) #confirm labels

pseudo_gaba1 <- AggregateExpression(gaba1, assays = "RNA", 
                                   return.seurat = T, 
                                   group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba1))


pseudo_gaba1$id.group <- paste(pseudo_gaba1$group, sep = "_")

Idents(pseudo_gaba1) <- "id.group"

gaba1.de <- FindMarkers(object = pseudo_gaba1, 
                       ident.1 = "RInfa", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(gaba1.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba1.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba1.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba1_deseq2.csv")
gaba1.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba1_deseq2.csv", header = TRUE)
gaba1.de[1,1] <- "Agouti"

# add a column of NAs
gaba1.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba1.de$Expression[gaba1.de$avg_log2FC > 0.01 & gaba1.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba1.de$Expression[gaba1.de$avg_log2FC < -0.01 & gaba1.de$p_val_adj < 0.1] <- "Up in Allo"

gaba1.de$delabel <- NA
gaba1.de$delabel[gaba1.de$Expression != "Not Diff"] <- gaba1.de$X[gaba1.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_GABA1Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba1.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba1.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()



### GAB2

gaba2 <- subset(x = neurons, idents = c("gab2"))

DimPlot(gaba2) #confirm labels

pseudo_gaba2 <- AggregateExpression(gaba2, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba2))


pseudo_gaba2$id.group <- paste(pseudo_gaba2$group, sep = "_")

Idents(pseudo_gaba2) <- "id.group"

gaba2.de <- FindMarkers(object = pseudo_gaba2, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba2.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba2.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba2.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba2_deseq2.csv")
gaba2.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba2_deseq2.csv", header = TRUE)
gaba2.de[1,1] <- "Agouti"

# add a column of NAs
gaba2.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba2.de$Expression[gaba2.de$avg_log2FC > 0.01 & gaba2.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba2.de$Expression[gaba2.de$avg_log2FC < -0.01 & gaba2.de$p_val_adj < 0.1] <- "Up in Allo"

gaba2.de$delabel <- NA
gaba2.de$delabel[gaba2.de$Expression != "Not Diff"] <- gaba2.de$X[gaba2.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba2Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba2.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba2.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GAB3

gaba3 <- subset(x = neurons, idents = c("gab3"))

DimPlot(gaba3) #confirm labels

pseudo_gaba3 <- AggregateExpression(gaba3, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba3))


pseudo_gaba3$id.group <- paste(pseudo_gaba3$group, sep = "_")

Idents(pseudo_gaba3) <- "id.group"

gaba3.de <- FindMarkers(object = pseudo_gaba3, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba3.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba3.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba3.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba3_deseq2.csv")
gaba3.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba3_deseq2.csv", header = TRUE)
gaba3.de[1,1] <- "Agouti"

# add a column of NAs
gaba3.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba3.de$Expression[gaba3.de$avg_log2FC > 0.01 & gaba3.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba3.de$Expression[gaba3.de$avg_log2FC < -0.01 & gaba3.de$p_val_adj < 0.1] <- "Up in Allo"

gaba3.de$delabel <- NA
gaba3.de$delabel[gaba3.de$Expression != "Not Diff"] <- gaba3.de$X[gaba3.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba3Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba3.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba3.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#00B9B9", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GAB4

gaba4 <- subset(x = neurons, idents = c("gab4"))

DimPlot(gaba4) #confirm labels

pseudo_gaba4 <- AggregateExpression(gaba4, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba4))


pseudo_gaba4$id.group <- paste(pseudo_gaba4$group, sep = "_")

Idents(pseudo_gaba4) <- "id.group"

gaba4.de <- FindMarkers(object = pseudo_gaba4, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba4.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba4.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba4.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba4_deseq2.csv")
gaba4.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba4_deseq2.csv", header = TRUE)
gaba4.de[1,1] <- "Agouti"

# add a column of NAs
gaba4.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba4.de$Expression[gaba4.de$avg_log2FC > 0.01 & gaba4.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba4.de$Expression[gaba4.de$avg_log2FC < -0.01 & gaba4.de$p_val_adj < 0.1] <- "Up in Allo"

gaba4.de$delabel <- NA
gaba4.de$delabel[gaba4.de$Expression != "Not Diff"] <- gaba4.de$X[gaba4.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba4Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba4.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba4.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#00B9B9", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GAB5

gaba5 <- subset(x = neurons, idents = c("gab5"))

DimPlot(gaba5) #confirm labels

pseudo_gaba5 <- AggregateExpression(gaba5, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba5))


pseudo_gaba5$id.group <- paste(pseudo_gaba5$group, sep = "_")

Idents(pseudo_gaba5) <- "id.group"

gaba5.de <- FindMarkers(object = pseudo_gaba5, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba5.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba5.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba5.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba5_deseq2.csv")
gaba5.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba5_deseq2.csv", header = TRUE)
gaba5.de[1,1] <- "Agouti"

# add a column of NAs
gaba5.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba5.de$Expression[gaba5.de$avg_log2FC > 0.01 & gaba5.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba5.de$Expression[gaba5.de$avg_log2FC < -0.01 & gaba5.de$p_val_adj < 0.1] <- "Up in Allo"

gaba5.de$delabel <- NA
gaba5.de$delabel[gaba5.de$Expression != "Not Diff"] <- gaba5.de$X[gaba5.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba5Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba5.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba5.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GAB6

gaba6 <- subset(x = neurons, idents = c("gab6"))

DimPlot(gaba6) #confirm labels

pseudo_gaba6 <- AggregateExpression(gaba6, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba6))


pseudo_gaba6$id.group <- paste(pseudo_gaba6$group, sep = "_")

Idents(pseudo_gaba6) <- "id.group"

gaba6.de <- FindMarkers(object = pseudo_gaba6, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba6.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba6.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba6.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba6_deseq2.csv")
gaba6.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba6_deseq2.csv", header = TRUE)
gaba6.de[1,1] <- "Agouti"

# add a column of NAs
gaba6.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba6.de$Expression[gaba6.de$avg_log2FC > 0.01 & gaba6.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba6.de$Expression[gaba6.de$avg_log2FC < -0.01 & gaba6.de$p_val_adj < 0.1] <- "Up in Allo"

gaba6.de$delabel <- NA
gaba6.de$delabel[gaba6.de$Expression != "Not Diff"] <- gaba6.de$X[gaba6.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba6Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba6.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba6.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30))  
dev.off()

### GAB7

gaba7 <- subset(x = neurons, idents = c("gab7"))

DimPlot(gaba7) #confirm labels

pseudo_gaba7 <- AggregateExpression(gaba7, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba7))


pseudo_gaba7$id.group <- paste(pseudo_gaba7$group, sep = "_")

Idents(pseudo_gaba7) <- "id.group"

gaba7.de <- FindMarkers(object = pseudo_gaba7, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba7.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba7.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba7.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba7_deseq2.csv")
gaba7.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba7_deseq2.csv", header = TRUE)
gaba7.de[1,1] <- "Agouti"

# add a column of NAs
gaba7.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba7.de$Expression[gaba7.de$avg_log2FC > 0.01 & gaba7.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba7.de$Expression[gaba7.de$avg_log2FC < -0.01 & gaba7.de$p_val_adj < 0.1] <- "Up in Allo"

gaba7.de$delabel <- NA
gaba7.de$delabel[gaba7.de$Expression != "Not Diff"] <- gaba7.de$X[gaba7.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba7Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba7.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba7.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GAB8

gaba8 <- subset(x = neurons, idents = c("gab8"))

DimPlot(gaba8) #confirm labels

pseudo_gaba8 <- AggregateExpression(gaba8, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba8))


pseudo_gaba8$id.group <- paste(pseudo_gaba8$group, sep = "_")

Idents(pseudo_gaba8) <- "id.group"

gaba8.de <- FindMarkers(object = pseudo_gaba8, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba8.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba8.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba8.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba8_deseq2.csv")
gaba8.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba8_deseq2.csv", header = TRUE)


# add a column of NAs
gaba8.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba8.de$Expression[gaba8.de$avg_log2FC > 0.01 & gaba8.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba8.de$Expression[gaba8.de$avg_log2FC < -0.01 & gaba8.de$p_val_adj < 0.1] <- "Up in Allo"

gaba8.de$delabel <- NA
gaba8.de$delabel[gaba8.de$Expression != "Not Diff"] <- gaba8.de$X[gaba8.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba8Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba8.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba8.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#00B9B9", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()


### GAB9

gaba9 <- subset(x = neurons, idents = c("gab9"))

DimPlot(gaba9) #confirm labels

pseudo_gaba9 <- AggregateExpression(gaba9, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba9))


pseudo_gaba9$id.group <- paste(pseudo_gaba9$group, sep = "_")

Idents(pseudo_gaba9) <- "id.group"

gaba9.de <- FindMarkers(object = pseudo_gaba9, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba9.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba9.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba9.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba9_deseq2.csv")
gaba9.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba9_deseq2.csv", header = TRUE)
gaba9.de[1,1] <- "Agouti"

# add a column of NAs
gaba9.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba9.de$Expression[gaba9.de$avg_log2FC > 0.01 & gaba9.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba9.de$Expression[gaba9.de$avg_log2FC < -0.01 & gaba9.de$p_val_adj < 0.1] <- "Up in Allo"

gaba9.de$delabel <- NA
gaba9.de$delabel[gaba9.de$Expression != "Not Diff"] <- gaba9.de$X[gaba9.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba9Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba9.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba9.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GAB10

gaba10 <- subset(x = neurons, idents = c("gab10"))

DimPlot(gaba10) #confirm labels

pseudo_gaba10 <- AggregateExpression(gaba10, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba10))


pseudo_gaba10$id.group <- paste(pseudo_gaba10$group, sep = "_")

Idents(pseudo_gaba10) <- "id.group"

gaba10.de <- FindMarkers(object = pseudo_gaba10, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba10.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba10.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba10.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba10_deseq2.csv")
gaba10.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba10_deseq2.csv", header = TRUE)
gaba10.de[1,1] <- "Agouti"

# add a column of NAs
gaba10.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba10.de$Expression[gaba10.de$avg_log2FC > 0.01 & gaba10.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba10.de$Expression[gaba10.de$avg_log2FC < -0.01 & gaba10.de$p_val_adj < 0.1] <- "Up in Allo"

gaba10.de$delabel <- NA
gaba10.de$delabel[gaba10.de$Expression != "Not Diff"] <- gaba10.de$X[gaba10.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba10Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba10.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba10.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GAB11

gaba11 <- subset(x = neurons, idents = c("gab11"))

DimPlot(gaba11) #confirm labels

pseudo_gaba11 <- AggregateExpression(gaba11, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba11))


pseudo_gaba11$id.group <- paste(pseudo_gaba11$group, sep = "_")

Idents(pseudo_gaba11) <- "id.group"

gaba11.de <- FindMarkers(object = pseudo_gaba11, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba11.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba11.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba11.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba11_deseq2.csv")
gaba11.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba11_deseq2.csv", header = TRUE)
gaba11.de[1,1] <- "Agouti"

# add a column of NAs
gaba11.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba11.de$Expression[gaba11.de$avg_log2FC > 0.01 & gaba11.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba11.de$Expression[gaba11.de$avg_log2FC < -0.01 & gaba11.de$p_val_adj < 0.1] <- "Up in Allo"

gaba11.de$delabel <- NA
gaba11.de$delabel[gaba11.de$Expression != "Not Diff"] <- gaba11.de$X[gaba11.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba11Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba11.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba11.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GAB12

gaba12 <- subset(x = neurons, idents = c("gab12"))

DimPlot(gaba12) #confirm labels

pseudo_gaba12 <- AggregateExpression(gaba12, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_gaba12))


pseudo_gaba12$id.group <- paste(pseudo_gaba12$group, sep = "_")

Idents(pseudo_gaba12) <- "id.group"

gaba12.de <- FindMarkers(object = pseudo_gaba12, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(gaba12.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(gaba12.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(gaba12.de, file = "/Volumes/Crucial X6/sequencing/neurons_gaba12_deseq2.csv")
gaba12.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_gaba12_deseq2.csv", header = TRUE)
gaba12.de[1,1] <- "Agouti"

# add a column of NAs
gaba12.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
gaba12.de$Expression[gaba12.de$avg_log2FC > 0.01 & gaba12.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
gaba12.de$Expression[gaba12.de$avg_log2FC < -0.01 & gaba12.de$p_val_adj < 0.1] <- "Up in Allo"

gaba12.de$delabel <- NA
gaba12.de$delabel[gaba12.de$Expression != "Not Diff"] <- gaba12.de$X[gaba12.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_gaba12Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=gaba12.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(gaba12.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GLUTA1

glu1 <- subset(x = neurons, idents = c("glu1"))

DimPlot(glu1) #confirm labels

pseudo_glu1 <- AggregateExpression(glu1, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

tail(Cells(pseudo_glu1))


pseudo_glu1$id.group <- paste(pseudo_glu1$group, sep = "_")

Idents(pseudo_glu1) <- "id.group"

glu1.de <- FindMarkers(object = pseudo_glu1, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(glu1.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(glu1.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(glu1.de, file = "/Volumes/Crucial X6/sequencing/neurons_glu1_deseq2.csv")
glu1.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_glu1_deseq2.csv", header = TRUE)
glu1.de[1,1] <- "Agouti"

# add a column of NAs
glu1.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
glu1.de$Expression[glu1.de$avg_log2FC > 0.01 & glu1.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
glu1.de$Expression[glu1.de$avg_log2FC < -0.01 & glu1.de$p_val_adj < 0.1] <- "Up in Allo"

glu1.de$delabel <- NA
glu1.de$delabel[glu1.de$Expression != "Not Diff"] <- glu1.de$X[glu1.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_glu1Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=glu1.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(glu1.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GLUTA2

glu2 <- subset(x = neurons, idents = c("glu2"))

DimPlot(glu2) #confirm labels

pseudo_glu2 <- AggregateExpression(glu2, assays = "RNA", 
                                   return.seurat = T, 
                                   group.by = c("group", "orig.ident"))

tail(Cells(pseudo_glu2))


pseudo_glu2$id.group <- paste(pseudo_glu2$group, sep = "_")

Idents(pseudo_glu2) <- "id.group"

glu2.de <- FindMarkers(object = pseudo_glu2, 
                       ident.1 = "RInfa", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(glu2.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(glu2.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(glu2.de, file = "/Volumes/Crucial X6/sequencing/neurons_glu2_deseq2.csv")
glu2.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_glu2_deseq2.csv", header = TRUE)
glu2.de[1,1] <- "Agouti"

# add a column of NAs
glu2.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
glu2.de$Expression[glu2.de$avg_log2FC > 0.01 & glu2.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
glu2.de$Expression[glu2.de$avg_log2FC < -0.01 & glu2.de$p_val_adj < 0.1] <- "Up in Allo"

glu2.de$delabel <- NA
glu2.de$delabel[glu2.de$Expression != "Not Diff"] <- glu2.de$X[glu2.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_glu2Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=glu2.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(glu2.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GLUTA3

glu3 <- subset(x = neurons, idents = c("glu3"))

DimPlot(glu3) #confirm labels

pseudo_glu3 <- AggregateExpression(glu3, assays = "RNA", 
                                   return.seurat = T, 
                                   group.by = c("group", "orig.ident"))

tail(Cells(pseudo_glu3))


pseudo_glu3$id.group <- paste(pseudo_glu3$group, sep = "_")

Idents(pseudo_glu3) <- "id.group"

glu3.de <- FindMarkers(object = pseudo_glu3, 
                       ident.1 = "RInfa", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(glu3.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(glu3.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(glu3.de, file = "/Volumes/Crucial X6/sequencing/neurons_glu3_deseq2.csv")
glu3.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_glu3_deseq2.csv", header = TRUE)
glu3.de[1,1] <- "Agouti"

# add a column of NAs
glu3.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
glu3.de$Expression[glu3.de$avg_log2FC > 0.01 & glu3.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
glu3.de$Expression[glu3.de$avg_log2FC < -0.01 & glu3.de$p_val_adj < 0.1] <- "Up in Allo"

glu3.de$delabel <- NA
glu3.de$delabel[glu3.de$Expression != "Not Diff"] <- glu3.de$X[glu3.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_glu3Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=glu3.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(glu3.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### GLUTA4

glu4 <- subset(x = neurons, idents = c("glu4"))

DimPlot(glu4) #confirm labels

pseudo_glu4 <- AggregateExpression(glu4, assays = "RNA", 
                                   return.seurat = T, 
                                   group.by = c("group", "orig.ident"))

tail(Cells(pseudo_glu4))


pseudo_glu4$id.group <- paste(pseudo_glu4$group, sep = "_")

Idents(pseudo_glu4) <- "id.group"

glu4.de <- FindMarkers(object = pseudo_glu4, 
                       ident.1 = "RInfa", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(glu4.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(glu4.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(glu4.de, file = "/Volumes/Crucial X6/sequencing/neurons_glu4_deseq2.csv")
glu4.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_glu4_deseq2.csv", header = TRUE)
glu4.de[1,1] <- "Agouti"

# add a column of NAs
glu4.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
glu4.de$Expression[glu4.de$avg_log2FC > 0.01 & glu4.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
glu4.de$Expression[glu4.de$avg_log2FC < -0.01 & glu4.de$p_val_adj < 0.1] <- "Up in Allo"

glu4.de$delabel <- NA
glu4.de$delabel[glu4.de$Expression != "Not Diff"] <- glu4.de$X[glu4.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_glu4Neurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=glu4.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(glu4.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

### NPEP

npep <- subset(x = neurons, idents = c("npep"))

DimPlot(npep) #confirm labels

pseudo_npep <- AggregateExpression(npep, assays = "RNA", 
                                   return.seurat = T, 
                                   group.by = c("group", "orig.ident"))

tail(Cells(pseudo_npep))


pseudo_npep$id.group <- paste(pseudo_npep$group, sep = "_")

Idents(pseudo_npep) <- "id.group"

npep.de <- FindMarkers(object = pseudo_npep, 
                       ident.1 = "RInfa", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(npep.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(npep.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res
summary(sig_res)

write.csv(npep.de, file = "/Volumes/Crucial X6/sequencing/neurons_npep_deseq2.csv")
npep.de <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_npep_deseq2.csv", header = TRUE)


# add a column of NAs
npep.de$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
npep.de$Expression[npep.de$avg_log2FC > 0.01 & npep.de$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
npep.de$Expression[npep.de$avg_log2FC < -0.01 & npep.de$p_val_adj < 0.1] <- "Up in Allo"

npep.de$delabel <- NA
npep.de$delabel[npep.de$Expression != "Not Diff"] <- npep.de$X[npep.de$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_npepNeurons.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=npep.de, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = sqrt(-log10(npep.de$p_val))) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 50),
        axis.text.y = element_text(size = 50), axis.title.y = element_text(size = 50), axis.title.x = element_text(size = 50),
        text = element_text(size = 50, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()
#### Immediate Early Genes below

#####
### Fos Positive Subset
rm(list = ls())

## Load in Data
neurons <- readRDS(file = "/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")

fos.active <- subset(neurons, subset = Fos > 0)

FOS <- subset(x = fos.active)

DimPlot(FOS) #confirm labels

####

pseudo_fos <- AggregateExpression(fos.active, assays = "RNA", 
                                  return.seurat = T, 
                                  group.by = c("group", "orig.ident"))

tail(Cells(pseudo_fos))


pseudo_fos$id.group <- paste(pseudo_fos$group, sep = "_")

Idents(pseudo_fos) <- "id.group"

fos.de <- FindMarkers(object = pseudo_fos, 
                      ident.1 = "RInfa", 
                      ident.2 = "RAllo",
                      test.use = "DESeq2")

View(fos.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(fos.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(fos.de, file = "/Volumes/Crucial X6/sequencing/neurons_FOS_infaVallo_deseq2.csv")
FOSplot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_FOS_infaVallo_deseq2.csv", header = TRUE)
FOSplot[1,1] <- "Agouti"

# add a column of NAs
FOSplot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
FOSplot$Expression[FOSplot$avg_log2FC > 0.01 & FOSplot$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
FOSplot$Expression[FOSplot$avg_log2FC < -0.01 & FOSplot$p_val_adj < 0.1] <- "Up in Allo"

FOSplot$delabel <- NA
FOSplot$delabel[FOSplot$Expression != "Not Diff"] <- FOSplot$X[FOSplot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_FOS_infaVallo.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=FOSplot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = 1) + 
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 30),
        axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40),
        text = element_text(size = 40, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

# Change to Control vs. Allo

fos.de <- FindMarkers(object = pseudo_fos, 
                      ident.1 = "RCont", 
                      ident.2 = "RAllo",
                      test.use = "DESeq2")

#Cannot Run; Not enough Fos+ cells in RCont.


### Fosb Positive Subset

fosb.active <- subset(neurons, subset = Fosb > 0)

FOSB <- subset(x = fosb.active)

DimPlot(FOSB, group.by = "group") #confirm labels

####

pseudo_fosb <- AggregateExpression(fosb.active, assays = "RNA", 
                                  return.seurat = T, 
                                  group.by = c("group", "orig.ident"))

tail(Cells(pseudo_fosb))


pseudo_fosb$id.group <- paste(pseudo_fosb$group, sep = "_")

Idents(pseudo_fosb) <- "id.group"

fosb.de <- FindMarkers(object = pseudo_fosb, 
                      ident.1 = "RInfa", 
                      ident.2 = "RAllo",
                      test.use = "DESeq2")

View(fosb.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(fosb.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(fosb.de, file = "/Volumes/Crucial X6/sequencing/neurons_FOSB_infaVallo_deseq2.csv")
FOSBplot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_FOSB_infaVallo_deseq2.csv", header = TRUE)

# add a column of NAs
FOSBplot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
FOSBplot$Expression[FOSBplot$avg_log2FC > 0.01 & FOSBplot$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
FOSBplot$Expression[FOSBplot$avg_log2FC < -0.01 & FOSBplot$p_val_adj < 0.1] <- "Up in Allo"

FOSBplot$delabel <- NA
FOSBplot$delabel[FOSBplot$Expression != "Not Diff"] <- FOSBplot$X[FOSBplot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_FOSB_infaVallo.png",
    width = 150, height = 150, units = "mm", res = 300)
ggplot(data=FOSBplot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = -log2(FOSBplot$p_val_adj)*.05) + 
  theme_linedraw() +
  scale_color_manual(values=c("grey90", "#00B9B9", "#960018")) + 
  geom_text_repel(fontface = "bold", label.size = 1) + 
  xlim(-2.5,2.5) + ylim(-1, 40) + 
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + plot_annotation("FOSB+ Subset, Infa v. Allo") +
  theme(legend.position = "bottom", legend.justification = "center", legend.text = element_text(size = 12, face = "bold", family = "arial"),
        legend.direction = "horizontal", legend.location = "plot", legend.byrow = TRUE, 
        axis.title = element_text(face="bold", size = 20), axis.text = element_text(face = "bold", size = 20))
dev.off()

# Change to Control vs. Allo

fosb.de <- FindMarkers(object = pseudo_fosb, 
                      ident.1 = "RCont", 
                      ident.2 = "RAllo",
                      test.use = "DESeq2")

View(fosb.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(fosb.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(fosb.de, file = "/Volumes/Crucial X6/sequencing/neurons_FOSB_contVallo_deseq2.csv")
FOSBplot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_FOSB_contVallo_deseq2.csv", header = TRUE)

# add a column of NAs
FOSBplot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
FOSBplot$Expression[FOSBplot$avg_log2FC > 0.01 & FOSBplot$p_val_adj < 0.1] <- "Up in Control"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
FOSBplot$Expression[FOSBplot$avg_log2FC < -0.01 & FOSBplot$p_val_adj < 0.1] <- "Up in Allo"

FOSBplot$delabel <- NA
FOSBplot$delabel[FOSBplot$Expression != "Not Diff"] <- FOSBplot$X[FOSBplot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_FOSB_contVallo.png",
    width = 150, height = 150, units = "mm", res = 300)
ggplot(data=FOSBplot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = -log2(FOSBplot$p_val_adj)*.05) + 
  theme_linedraw() +
  scale_color_manual(values=c("grey90", "#00B9B9", "#DCB50B")) + 
  geom_text_repel(fontface = "bold", label.size = 1) + 
  xlim(-2.5,2.5) + ylim(-1, 40) + 
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + plot_annotation("FOSB+ Subset, Control v. Allo") +
  theme(legend.position = "bottom", legend.justification = "center", legend.text = element_text(size = 12, face = "bold", family = "arial"),
        legend.direction = "horizontal", legend.location = "plot", legend.byrow = TRUE, 
        axis.title = element_text(face="bold", size = 20), axis.text = element_text(face = "bold", size = 20))
dev.off()

### Fosl1 Positive Subset

fosl1.active <- subset(neurons, subset = Fosl1 > 0)

FOSL1 <- subset(x = fosl1.active)

DimPlot(FOSL1, group.by = "group") #confirm labels

####

pseudo_fosl1 <- AggregateExpression(fosl1.active, assays = "RNA", 
                                   return.seurat = T, 
                                   group.by = c("group", "orig.ident"))

pseudo_fosl1$id.group <- paste(pseudo_fosl1$group, sep = "_")

Idents(pseudo_fosl1) <- "id.group"

fosl1.de <- FindMarkers(object = pseudo_fosl1, 
                       ident.1 = "RInfa", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(fosl1.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(fosl1.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(fosl1.de, file = "/Volumes/Crucial X6/sequencing/neurons_FOSL1_infaVallo_deseq2.csv")
FOSL1plot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_FOSL1_infaVallo_deseq2.csv", header = TRUE)

# add a column of NAs
FOSL1plot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
FOSL1plot$Expression[FOSL1plot$avg_log2FC > 0.01 & FOSL1plot$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
FOSL1plot$Expression[FOSL1plot$avg_log2FC < -0.01 & FOSL1plot$p_val_adj < 0.1] <- "Up in Allo"

FOSL1plot$delabel <- NA
FOSL1plot$delabel[FOSL1plot$Expression != "Not Diff"] <- FOSL1plot$X[FOSL1plot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_FOSL1_infaVallo.png",
    width = 150, height = 150, units = "mm", res = 300)
ggplot(data=FOSL1plot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = -log2(FOSL1plot$p_val_adj)*.05) + 
  theme_linedraw() +
  scale_color_manual(values=c("grey90", "#00B9B9", "#960018")) + 
  geom_text_repel(fontface = "bold", label.size = 1) + 
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + plot_annotation("FOSL1+ Subset, Infa v. Allo") +
  theme(legend.position = "bottom", legend.justification = "center", legend.text = element_text(size = 12, face = "bold", family = "arial"),
        legend.direction = "horizontal", legend.location = "plot", legend.byrow = TRUE, 
        axis.title = element_text(face="bold", size = 20), axis.text = element_text(face = "bold", size = 20))
dev.off()

# Change to Control vs. Allo

fosl1.de <- FindMarkers(object = pseudo_fosl1, 
                       ident.1 = "RCont", 
                       ident.2 = "RAllo",
                       test.use = "DESeq2")

View(fosl1.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(fosl1.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(fosl1.de, file = "/Volumes/Crucial X6/sequencing/neurons_FOSL1_contVallo_deseq2.csv")
FOSL1plot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_FOSL1_contVallo_deseq2.csv", header = TRUE)

# add a column of NAs
FOSL1plot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
FOSL1plot$Expression[FOSL1plot$avg_log2FC > 0.01 & FOSL1plot$p_val_adj < 0.1] <- "Up in Control"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
FOSL1plot$Expression[FOSL1plot$avg_log2FC < -0.01 & FOSL1plot$p_val_adj < 0.1] <- "Up in Allo"

FOSL1plot$delabel <- NA
FOSL1plot$delabel[FOSL1plot$Expression != "Not Diff"] <- FOSL1plot$X[FOSL1plot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_FOSL1_contVallo.png",
    width = 150, height = 150, units = "mm", res = 300)
ggplot(data=FOSL1plot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = -log2(FOSL1plot$p_val_adj)*.05) + 
  theme_linedraw() +
  scale_color_manual(values=c("grey90", "#00B9B9", "#DCB50B")) + 
  geom_text_repel(fontface = "bold", label.size = 1) + 
  xlim(-2.5,2.5) + ylim(-1, 40) + 
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + plot_annotation("FOSL1+ Subset, Control v. Allo") +
  theme(legend.position = "bottom", legend.justification = "center", legend.text = element_text(size = 12, face = "bold", family = "arial"),
        legend.direction = "horizontal", legend.location = "plot", legend.byrow = TRUE, 
        axis.title = element_text(face="bold", size = 20), axis.text = element_text(face = "bold", size = 20))
dev.off()

### Npas4 Positive Subset

npas4.active <- subset(neurons, subset = Npas4 > 0)

NPAS4 <- subset(x = npas4.active)

DimPlot(NPAS4, group.by = "group") #confirm labels

####

pseudo_npas4 <- AggregateExpression(npas4.active, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

pseudo_npas4$id.group <- paste(pseudo_npas4$group, sep = "_")

Idents(pseudo_npas4) <- "id.group"

npas4.de <- FindMarkers(object = pseudo_npas4, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(npas4.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(npas4.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(npas4.de, file = "/Volumes/Crucial X6/sequencing/neurons_NPAS4_infaVallo_deseq2.csv")
NPAS4plot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_NPAS4_infaVallo_deseq2.csv", header = TRUE)
NPAS4plot[1,1] <- "Agouti"

# add a column of NAs
NPAS4plot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
NPAS4plot$Expression[NPAS4plot$avg_log2FC > 0.01 & NPAS4plot$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
NPAS4plot$Expression[NPAS4plot$avg_log2FC < -0.01 & NPAS4plot$p_val_adj < 0.1] <- "Up in Allo"

NPAS4plot$delabel <- NA
NPAS4plot$delabel[NPAS4plot$Expression != "Not Diff"] <- NPAS4plot$X[NPAS4plot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_NPAS4_infaVallo.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=NPAS4plot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = 1) + 
  theme_linedraw() +
  scale_color_manual(values=c("grey10", "#00B9B9", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 30),
        axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40),
        text = element_text(size = 40, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

# Change to Control vs. Allo

npas4.de <- FindMarkers(object = pseudo_npas4, 
                        ident.1 = "RCont", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(npas4.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(npas4.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(npas4.de, file = "/Volumes/Crucial X6/sequencing/neurons_NPAS4_contVallo_deseq2.csv")
NPAS4plot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_NPAS4_contVallo_deseq2.csv", header = TRUE)

# add a column of NAs
NPAS4plot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
NPAS4plot$Expression[NPAS4plot$avg_log2FC > 0.01 & NPAS4plot$p_val_adj < 0.1] <- "Up in Control"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
NPAS4plot$Expression[NPAS4plot$avg_log2FC < -0.01 & NPAS4plot$p_val_adj < 0.1] <- "Up in Allo"

NPAS4plot$delabel <- NA
NPAS4plot$delabel[NPAS4plot$Expression != "Not Diff"] <- NPAS4plot$X[NPAS4plot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_NPAS4_contVallo.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=NPAS4plot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = 1) + 
  theme_linedraw() +
  scale_color_manual(values=c("grey10","#DCB50B")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 30),
        axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40),
        text = element_text(size = 40, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()



### Egr1 Positive Subset

EGR1.active <- subset(neurons, subset = Egr1 > 0)

EGR1 <- subset(x = EGR1.active)

DimPlot(EGR1, group.by = "group") #confirm labels

####

pseudo_egr1 <- AggregateExpression(EGR1.active, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

pseudo_egr1$id.group <- paste(pseudo_egr1$group, sep = "_")

Idents(pseudo_egr1) <- "id.group"

egr1.de <- FindMarkers(object = pseudo_egr1, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

# Cannot run contrasts for Egr1: Insufficient Egr1+ neurons in both infa and control.

### Egr3 Positive Subset

Egr3.active <- subset(neurons, subset = Egr3 > 0)

EGR3 <- subset(x = Egr3.active)

DimPlot(EGR3, group.by = "group") #confirm labels

####

pseudo_egr3 <- AggregateExpression(Egr3.active, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

pseudo_egr3$id.group <- paste(pseudo_egr3$group, sep = "_")

Idents(pseudo_egr3) <- "id.group"

egr3.de <- FindMarkers(object = pseudo_egr3, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

#cannot draw comparisons, insufficient number of Egr3+ neurons in alloparents

### Arc Positive Subset

arc.active <- subset(neurons, subset = Arc > 0)

ARC <- subset(x = arc.active)

DimPlot(ARC, group.by = "group") #confirm labels

####

pseudo_arc <- AggregateExpression(arc.active, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

pseudo_arc$id.group <- paste(pseudo_arc$group, sep = "_")

Idents(pseudo_arc) <- "id.group"

arc.de <- FindMarkers(object = pseudo_arc, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

#Cannot run comparisons; insufficient Arc+ neurons in Alloparents

### Jun Positive Subset

jun.active <- subset(neurons, subset = Jun > 0)

JUN <- subset(x = jun.active)

DimPlot(JUN, group.by = "group") #confirm labels

####

pseudo_jun <- AggregateExpression(jun.active, assays = "RNA", 
                                    return.seurat = T, 
                                    group.by = c("group", "orig.ident"))

pseudo_jun$id.group <- paste(pseudo_jun$group, sep = "_")

Idents(pseudo_jun) <- "id.group"

jun.de <- FindMarkers(object = pseudo_jun, 
                        ident.1 = "RInfa", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

# There aren't enough Jun+ Neurons for analysis.

# Junb positive subset
junb.active <- subset(neurons, subset = Junb > 0)

JUNB <- subset(x = junb.active)

DimPlot(JUNB, group.by = "group") #confirm labels

####

pseudo_junb <- AggregateExpression(junb.active, assays = "RNA", 
                                  return.seurat = T, 
                                  group.by = c("group", "orig.ident"))

pseudo_junb$id.group <- paste(pseudo_junb$group, sep = "_")

Idents(pseudo_junb) <- "id.group"

junb.de <- FindMarkers(object = pseudo_junb, 
                      ident.1 = "RInfa", 
                      ident.2 = "RAllo",
                      test.use = "DESeq2")

View(junb.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(junb.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(junb.de, file = "/Volumes/Crucial X6/sequencing/neurons_JUNB_infaVallo_deseq2.csv")
JUNBplot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_JUNB_infaVallo_deseq2.csv", header = TRUE)

# add a column of NAs
JUNBplot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
JUNBplot$Expression[JUNBplot$avg_log2FC > 0.01 & JUNBplot$p_val_adj < 0.1] <- "Up in Infa"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
JUNBplot$Expression[JUNBplot$avg_log2FC < -0.01 & JUNBplot$p_val_adj < 0.1] <- "Up in Allo"

JUNBplot$delabel <- NA
JUNBplot$delabel[JUNBplot$Expression != "Not Diff"] <- JUNBplot$X[JUNBplot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_JUNB_infaVallo.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=JUNBplot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = 1) + 
  theme_linedraw() +
  scale_color_manual(values=c("grey10", "#960018")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 30),
        axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40),
        text = element_text(size = 40, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()



# Change to Control vs. Allo

junb.de <- FindMarkers(object = pseudo_junb, 
                        ident.1 = "RCont", 
                        ident.2 = "RAllo",
                        test.use = "DESeq2")

View(junb.de)

padj_cutoff <- .1

# Subset the significant results
sig_res <- dplyr::filter(junb.de, p_val_adj < padj_cutoff) %>%
  dplyr::arrange(p_val_adj)

# Check significant genes output
sig_res

write.csv(junb.de, file = "/Volumes/Crucial X6/sequencing/neurons_JUNB_contVallo_deseq2.csv")
JUNBplot <- read.csv( "/Volumes/Crucial X6/sequencing/neurons_JUNB_contVallo_deseq2.csv", header = TRUE)

# add a column of NAs
JUNBplot$Expression <- "Not Diff"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
JUNBplot$Expression[JUNBplot$avg_log2FC > 0.01 & JUNBplot$p_val_adj < 0.1] <- "Up in Control"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
JUNBplot$Expression[JUNBplot$avg_log2FC < -0.01 & JUNBplot$p_val_adj < 0.1] <- "Up in Allo"

JUNBplot$delabel <- NA
JUNBplot$delabel[JUNBplot$Expression != "Not Diff"] <- JUNBplot$X[JUNBplot$Expression != "Not Diff"]

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Pseudo_AllNeurons_JUNB_contVallo.png",
    width = 400, height = 400, units = "mm", res = 600)
ggplot(data=JUNBplot, aes(x=-(avg_log2FC), y=-log10(p_val), col=Expression, label=delabel)) + 
  geom_point(size = 1) + 
  theme_linedraw() +
  scale_color_manual(values=c("grey10", "#DCB50B")) + 
  geom_text_repel(size = 18, fontface = "bold", family = "arial") + scale_y_continuous(trans = "sqrt") +
  xlab("Average Log2 Fold Change") + ylab("-log10 of p-value") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 30),
        axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40),
        text = element_text(size = 40, face = "bold",), legend.position = "none", legend.title = element_blank(), 
        legend.text = element_text(size = 30)) 
dev.off()

neurons <- readRDS(file = "/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")

a1 <- length(Cells(subset(x = neurons, orig.ident == c("RAllo1"))))
a1a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RAllo1"))))
a1r.npas4 <- round((a1a.npas4 / a1)*100, digits = 1)
a2 <- length(Cells(subset(x = neurons, orig.ident == c("RAllo2"))))
a2a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RAllo2"))))
a2r.npas4 <- round((a2a.npas4 / a2)*100, digits = 1)
a3 <- length(Cells(subset(x = neurons, orig.ident == c("RAllo3"))))
a3a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RAllo3"))))
a3r.npas4 <- round((a3a.npas4 / a3)*100, digits = 1)

a1a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RAllo1"))))
a1r.fos <- round((a1a.fos / a1)*100, digits = 1)
a2a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RAllo2"))))
a2r.fos <- round((a2a.fos / a2)*100, digits = 1)
a3a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RAllo3"))))
a3r.fos <- round((a3a.fos / a3)*100, digits = 1)

a1a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RAllo1"))))
a1r.fosb <- round((a1a.fosb / a1)*100, digits = 1)
a2a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RAllo2"))))
a2r.fosb <- round((a2a.fosb / a2)*100, digits = 1)
a3a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RAllo3"))))
a3r.fosb <- round((a3a.fosb / a3)*100, digits = 1)

a1a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RAllo1"))))
a1r.fosl1 <- round((a1a.fosl1 / a1)*100, digits = 1)
a2a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RAllo2"))))
a2r.fosl1 <- round((a2a.fosl1 / a2)*100, digits = 1)
a3a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RAllo3"))))
a3r.fosl1 <- round((a3a.fosl1 / a3)*100, digits = 1)

a1a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RAllo1"))))
a1r.arc <- round((0 / a1)*100, digits = 1) #No Cells Found
a2a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RAllo2"))))
a2r.arc <- round((a2a.arc / a2)*100, digits = 1)
a3a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RAllo3"))))
a3r.arc <- round((a3a.arc / a3)*100, digits = 1)

a1a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RAllo1"))))
a1r.egr1 <- round((0 / a1)*100, digits = 1) #No Cells Found
a2a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RAllo2"))))
a2r.egr1 <- round((a2a.egr1 / a2)*100, digits = 1)
a3a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RAllo3"))))
a3r.egr1 <- round((a3a.egr1 / a3)*100, digits = 1)

a1a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RAllo1"))))
a1r.egr3 <- round((0 / a1)*100, digits = 1) #No Cells Found
a2a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RAllo2"))))
a2r.egr3 <- round((a2a.egr3 / a2)*100, digits = 1)
a3a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RAllo3"))))
a3r.egr3 <- round((a3a.egr3 / a3)*100, digits = 1)

a1a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RAllo1"))))
a1r.jun <- round((0 / a1)*100, digits = 1) #No Cells Found
a2a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RAllo2"))))
a2r.jun <- round((0 / a2)*100, digits = 1) #No Cells Found
a3a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RAllo3"))))
a3r.jun <- round((a3a.jun / a3)*100, digits = 1)

a1a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RAllo1"))))
a1r.junb <- round((a1a.junb / a1)*100, digits = 1) 
a2a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RAllo2"))))
a2r.junb <- round((a1a.junb / a2)*100, digits = 1) 
a3a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RAllo3"))))
a3r.junb <- round((a3a.junb / a3)*100, digits = 1)

c1 <- length(Cells(subset(x = neurons, orig.ident == c("RCont1"))))
c1a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RCont1"))))
c1r.npas4 <- round((c1a.npas4 / c1)*100, digits = 1)
c2 <- length(Cells(subset(x = neurons, orig.ident == c("RCont2"))))
c2a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RCont2"))))
c2r.npas4 <- round((c2a.npas4 / c2)*100, digits = 1)
c3 <- length(Cells(subset(x = neurons, orig.ident == c("RCont3"))))
c3a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RCont3"))))
c3r.npas4 <- round((c3a.npas4 / c3)*100, digits = 1)

c1a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RCont1"))))
c1r.fos <- round((c1a.fos / c1)*100, digits = 1)
c2a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RCont2"))))
c2r.fos <- round((0 / c2)*100, digits = 1) #No Cells Found
c3a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RCont3"))))
c3r.fos <- round((c3a.fos / c3)*100, digits = 1)

c1a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RCont1"))))
c1r.fosb <- round((c1a.fosb / c1)*100, digits = 1)
c2a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RCont2"))))
c2r.fosb <- round((c2a.fosb / c2)*100, digits = 1)
c3a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RCont3"))))
c3r.fosb <- round((c3a.fosb / c3)*100, digits = 1)

c1a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RCont1"))))
c1r.fosl1 <- round((c1a.fosl1 / c1)*100, digits = 1)
c2a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RCont2"))))
c2r.fosl1 <- round((c2a.fosl1 / c2)*100, digits = 1)
c3a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RCont3"))))
c3r.fosl1 <- round((c3a.fosl1 / c3)*100, digits = 1)

c1a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RCont1"))))
c1r.arc <- round((c1a.arc / c1)*100, digits = 1)
c2a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RCont2"))))
c2r.arc <- round((c2a.arc / c2)*100, digits = 1)
c3a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RCont3"))))
c3r.arc <- round((0 / c3)*100, digits = 1)  #No Cells Found

c1a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RCont1"))))
c1r.egr1 <- round((0 / c1)*100, digits = 1) #No Cells Found
c2a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RCont2"))))
c2r.egr1 <- round((0 / c2)*100, digits = 1)#No Cells Found
c3a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RCont3"))))
c3r.egr1 <- round((0 / c3)*100, digits = 1)#No Cells Found

c1a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RCont1"))))
c1r.egr3 <- round((c1a.egr3 / c1)*100, digits = 1) 
c2a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RCont2"))))
c2r.egr3 <- round((c2a.egr3 / c2)*100, digits = 1)
c3a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RCont3"))))
c3r.egr3 <- round((0 / c3)*100, digits = 1) #No Cells Found

c1a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RCont1"))))
c1r.jun <- round((c1a.jun / c1)*100, digits = 1) 
c2a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RCont2"))))
c2r.jun <- round((0 / c2)*100, digits = 1) #No Cells Found
c3a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RCont3"))))
c3r.jun <- round((c3a.jun / c3)*100, digits = 1)

c1a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RCont1"))))
c1r.junb <- round((c1a.junb / c1)*100, digits = 1) 
c2a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RCont2"))))
c2r.junb <- round((c1a.junb / c2)*100, digits = 1) 
c3a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RCont3"))))
c3r.junb <- round((c3a.junb / c3)*100, digits = 1)

i1 <- length(Cells(subset(x = neurons, orig.ident == c("RInfa1"))))
i1a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RInfa1"))))
i1r.npas4 <- round((i1a.npas4 / i1)*100, digits = 1)
i2 <- length(Cells(subset(x = neurons, orig.ident == c("RInfa2"))))
i2a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RInfa2"))))
i2r.npas4 <- round((i2a.npas4 / i2)*100, digits = 1)
i3 <- length(Cells(subset(x = neurons, orig.ident == c("RInfa3"))))
i3a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RInfa3"))))
i3r.npas4 <- round((i3a.npas4 / i3)*100, digits = 1)

i1a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RInfa1"))))
i1r.fos <- round((i1a.fos / i1)*100, digits = 1)
i2a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RInfa2"))))
i2r.fos <- round((i2a.fos / i2)*100, digits = 1)
i3a.fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RInfa3"))))
i3r.fos <- round((i3a.fos / i3)*100, digits = 1)

i1a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RInfa1"))))
i1r.fosb <- round((i1a.fosb / i1)*100, digits = 1)
i2a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RInfa2"))))
i2r.fosb <- round((i2a.fosb / i2)*100, digits = 1)
i3a.fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RInfa3"))))
i3r.fosb <- round((i3a.fosb / i3)*100, digits = 1)

i1a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RInfa1"))))
i1r.fosl1 <- round((i1a.fosl1 / i1)*100, digits = 1)
i2a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RInfa2"))))
i2r.fosl1 <- round((i2a.fosl1 / i2)*100, digits = 1)
i3a.fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RInfa3"))))
i3r.fosl1 <- round((i3a.fosl1 / i3)*100, digits = 1)

i1a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RInfa1"))))
i1r.arc <- round((i1a.arc / i1)*100, digits = 1) 
i2a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RInfa2"))))
i2r.arc <- round((i2a.arc / i2)*100, digits = 1)
i3a.arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RInfa3"))))
i3r.arc <- round((i3a.arc / i3)*100, digits = 1)

i1a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RInfa1"))))
i1r.egr1 <- round((i1a.egr1 / i1)*100, digits = 1) 
i2a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RInfa2"))))
i2r.egr1 <- round((0 / i2)*100, digits = 1) #No Cells Found
i3a.egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RInfa3"))))
i3r.egr1 <- round((i3a.egr1 / i3)*100, digits = 1)

i1a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RInfa1"))))
i1r.egr3 <- round((i1a.egr3 / i1)*100, digits = 1) 
i2a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RInfa2"))))
i2r.egr3 <- round((i2a.egr3 / i2)*100, digits = 1)
i3a.egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RInfa3"))))
i3r.egr3 <- round((i3a.egr3 / i3)*100, digits = 1)

i1a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RInfa1"))))
i1r.jun <- round((i1a.jun / i1)*100, digits = 1) 
i2a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RInfa2"))))
i2r.jun <- round((0 / i2)*100, digits = 1) #No Cells Found
i3a.jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RInfa3"))))
i3r.jun <- round((0 / i3)*100, digits = 1) #No Cells Found

i1a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RInfa1"))))
i1r.junb <- round((i1a.junb / i1)*100, digits = 1) 
i2a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RInfa2"))))
i2r.junb <- round((i1a.junb / i2)*100, digits = 1) 
i3a.junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RInfa3"))))
i3r.junb <- round((i3a.junb / i3)*100, digits = 1)

s1 <- length(Cells(subset(x = neurons, orig.ident == c("RSire1"))))
s1a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RSire1"))))
s1r.npas4 <- round((s1a.npas4 / s1)*100, digits = 1)
s2 <- length(Cells(subset(x = neurons, orig.ident == c("RSire2"))))
s2a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RSire2"))))
s2r.npas4 <- round((s2a.npas4 / s2)*100, digits = 1)
s3 <- length(Cells(subset(x = neurons, orig.ident == c("RSire3"))))
s3a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RSire3"))))
s3r.npas4 <- round((s3a.npas4 / s3)*100, digits = 1)
s4 <- length(Cells(subset(x = neurons, orig.ident == c("RSire4"))))
s4a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RSire4"))))
s4r.npas4 <- round((s4a.npas4 / s4)*100, digits = 1)


s1a.Fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RSire1"))))
s1r.Fos <- round((s1a.Fos / s1)*100, digits = 1)
s2a.Fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RSire2"))))
s2r.Fos <- round((s2a.Fos / s2)*100, digits = 1)
s3a.Fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RSire3"))))
s3r.Fos <- round((s3a.Fos / s3)*100, digits = 1)
s4a.Fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RSire4"))))
s4r.Fos <- round((s4a.Fos / s4)*100, digits = 1)

s1a.Fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RSire1"))))
s1r.Fosb <- round((0 / s1)*100, digits = 1) #No Cells Found
s2a.Fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RSire2"))))
s2r.Fosb <- round((s2a.Fosb / s2)*100, digits = 1)
s3a.Fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RSire3"))))
s3r.Fosb <- round((s3a.Fosb / s3)*100, digits = 1)
s4a.Fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RSire4"))))
s4r.Fosb <- round((s4a.Fosb / s4)*100, digits = 1)

s1a.Fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RSire1"))))
s1r.Fosl1 <- round((s1a.Fosl1 / s1)*100, digits = 1)
s2a.Fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RSire2"))))
s2r.Fosl1 <- round((s2a.Fosl1 / s2)*100, digits = 1)
s3a.Fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RSire3"))))
s3r.Fosl1 <- round((s3a.Fosl1 / s3)*100, digits = 1)
s4a.Fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RSire4"))))
s4r.Fosl1 <- round((s4a.Fosl1 / s4)*100, digits = 1)

s1a.Arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RSire1"))))
s1r.Arc <- round((s1a.Arc / s1)*100, digits = 1)
s2a.Arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RSire2"))))
s2r.Arc <- round((0 / s2)*100, digits = 1)  #No Cells Found
s3a.Arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RSire3"))))
s3r.Arc <- round((s3a.Arc / s3)*100, digits = 1)
s4a.Arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RSire4"))))
s4r.Arc <- round((s4a.Arc / s4)*100, digits = 1)

s1a.Egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RSire1"))))
s1r.Egr1 <- round((0 / s1)*100, digits = 1)  #No Cells Found
s2a.Egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RSire2"))))
s2r.Egr1 <- round((0 / s2)*100, digits = 1)  #No Cells Found
s3a.Egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RSire3"))))
s3r.Egr1 <- round((s3a.Egr1 / s3)*100, digits = 1)
s4a.Egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RSire4"))))
s4r.Egr1 <- round((s4a.Egr1 / s4)*100, digits = 1)

s1a.Egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RSire1"))))
s1r.Egr3 <- round((s1a.Egr3 / s1)*100, digits = 1)
s2a.Egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RSire2"))))
s2r.Egr3 <- round((0 / s2)*100, digits = 1)  #No Cells Found
s3a.Egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RSire3"))))
s3r.Egr3 <- round((s3a.Egr3 / s3)*100, digits = 1)
s4a.Egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RSire4"))))
s4r.Egr3 <- round((s4a.Egr3 / s4)*100, digits = 1)

s1a.Jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RSire1"))))
s1r.Jun <- round((0 / s1)*100, digits = 1)  #No Cells Found
s2a.Jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RSire2"))))
s2r.Jun <- round((0)*100, digits = 1)  #No Cells Found
s3a.Jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RSire3"))))
s3r.Jun <- round((0 / s3)*100, digits = 1)  #No Cells Found
s4a.Jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RSire4"))))
s4r.Jun <- round((0 / s4)*100, digits = 1)  #No Cells Found

s1a.Junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RSire1"))))
s1r.Junb <- round((s1a.Junb / s1)*100, digits = 1)
s2a.Junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RSire2"))))
s2r.Junb <- round((s2a.Junb / s2)*100, digits = 1)
s3a.Junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RSire3"))))
s3r.Junb <- round((s3a.Junb / s3)*100, digits = 1)
s4a.Junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RSire4"))))
s4r.Junb <- round((s4a.Junb / s4)*100, digits = 1)

d1 <- length(Cells(subset(x = neurons, orig.ident == c("RDam1"))))
d1a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RDam1"))))
d1r.npas4 <- round((d1a.npas4 / d1)*100, digits = 1)
d2 <- length(Cells(subset(x = neurons, orig.ident == c("RDam2"))))
d2a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RDam2"))))
d2r.npas4 <- round((d2a.npas4 / d2)*100, digits = 1)
d3 <- length(Cells(subset(x = neurons, orig.ident == c("RDam3"))))
d3a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RDam3"))))
d3r.npas4 <- round((d3a.npas4 / d3)*100, digits = 1)
d4 <- length(Cells(subset(x = neurons, orig.ident == c("RDam4"))))
d4a.npas4 <- length(Cells(subset(x = neurons, subset = Npas4 > 0 & orig.ident == c("RDam4"))))
d4r.npas4 <- round((d4a.npas4 / d4)*100, digits = 1)


d1a.Fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RDam1"))))
d1r.Fos <- round((0 / d1)*100, digits = 1) #No Cells Found
d2a.Fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RDam2"))))
d2r.Fos <- round((0 / d2)*100, digits = 1)  #No Cells Found
d3a.Fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RDam3"))))
d3r.Fos <- round((d3a.Fos / d3)*100, digits = 1)
d4a.Fos <- length(Cells(subset(x = neurons, subset = Fos > 0 & orig.ident == c("RDam4"))))
d4r.Fos <- round((d4a.Fos / d4)*100, digits = 1)

d1a.Fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RDam1"))))
d1r.Fosb <- round((d1a.Fosb / d1)*100, digits = 1)
d2a.Fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RDam2"))))
d2r.Fosb <- round((0 / d2)*100, digits = 1)  #No Cells Found
d3a.Fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RDam3"))))
d3r.Fosb <- round((d3a.Fosb / d3)*100, digits = 1)
d4a.Fosb <- length(Cells(subset(x = neurons, subset = Fosb > 0 & orig.ident == c("RDam4"))))
d4r.Fosb <- round((d4a.Fosb / d4)*100, digits = 1)

d1a.Fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RDam1"))))
d1r.Fosl1 <- round((d1a.Fosl1 / d1)*100, digits = 1)
d2a.Fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RDam2"))))
d2r.Fosl1 <- round((d2a.Fosl1 / d2)*100, digits = 1)
d3a.Fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RDam3"))))
d3r.Fosl1 <- round((d3a.Fosl1 / d3)*100, digits = 1)
d4a.Fosl1 <- length(Cells(subset(x = neurons, subset = Fosl1 > 0 & orig.ident == c("RDam4"))))
d4r.Fosl1 <- round((d4a.Fosl1 / d4)*100, digits = 1)

d1a.Arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RDam1"))))
d1r.Arc <- round((0 / d1)*100, digits = 1)  #No Cells Found
d2a.Arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RDam2"))))
d2r.Arc <- round((d2a.Arc / d2)*100, digits = 1)
d3a.Arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RDam3"))))
d3r.Arc <- round((d3a.Arc / d3)*100, digits = 1)
d4a.Arc <- length(Cells(subset(x = neurons, subset = Arc > 0 & orig.ident == c("RDam4"))))
d4r.Arc <- round((d4a.Arc / d4)*100, digits = 1)

d1a.Egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RDam1"))))
d1r.Egr1 <- round((0 / d1)*100, digits = 1)  #No Cells Found
d2a.Egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RDam2"))))
d2r.Egr1 <- round((0 / d2)*100, digits = 1)  #No Cells Found
d3a.Egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RDam3"))))
d3r.Egr1 <- round((d3a.Egr1 / d3)*100, digits = 1)
d4a.Egr1 <- length(Cells(subset(x = neurons, subset = Egr1 > 0 & orig.ident == c("RDam4"))))
d4r.Egr1 <- round((d4a.Egr1 / d4)*100, digits = 1)

d1a.Egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RDam1"))))
d1r.Egr3 <- round((0 / d1)*100, digits = 1)  #No Cells Found
d2a.Egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RDam2"))))
d2r.Egr3 <- round((d2a.Egr3 / d2)*100, digits = 1)
d3a.Egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RDam3"))))
d3r.Egr3 <- round((d3a.Egr3 / d3)*100, digits = 1)
d4a.Egr3 <- length(Cells(subset(x = neurons, subset = Egr3 > 0 & orig.ident == c("RDam4"))))
d4r.Egr3 <- round((d4a.Egr3 / d4)*100, digits = 1)

d1a.Jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RDam1"))))
d1r.Jun <- round((0 / d1)*100, digits = 1)  #No Cells Found
d2a.Jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RDam2"))))
d2r.Jun <- round((0 / d2)*100, digits = 1)  #No Cells Found
d3a.Jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RDam3"))))
d3r.Jun <- round((0 / d3)*100, digits = 1)  #No Cells Found
d4a.Jun <- length(Cells(subset(x = neurons, subset = Jun > 0 & orig.ident == c("RDam4"))))
d4r.Jun <- round((d4a.Jun / d4)*100, digits = 1)

d1a.Junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RDam1"))))
d1r.Junb <- round((0 / d1)*100, digits = 1)  #No Cells Found
d2a.Junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RDam2"))))
d2r.Junb <- round((d2a.Junb / d2)*100, digits = 1)
d3a.Junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RDam3"))))
d3r.Junb <- round((d3a.Junb / d3)*100, digits = 1)
d4a.Junb <- length(Cells(subset(x = neurons, subset = Junb > 0 & orig.ident == c("RDam4"))))
d4r.Junb <- round((d4a.Junb / d4)*100, digits = 1)

ids <- c("infa1", "infa2", "infa3", "con1", "con2", "con3", "allo1", "allo2", "allo3", "sire1", "sire2", "sire3", "sire4", "dam1", "dam2", "dam3", "dam4")
phenotype <- c("INFA", "INFA", "INFA", "CONT", "CONT", "CONT", "ALLO", "ALLO", "ALLO", "SIRE", "SIRE", "SIRE", "SIRE", "DAM", "DAM", "DAM", "DAM")
allneuro <- c(i1, i2, i3, c1, c2, c3, a1, a2, a3, s1, s2, s3, s4, d1, d2, d3, d4)
npas4 <- c(i1a.npas4, i2a.npas4, i3a.npas4, c1a.npas4, c2a.npas4, c3a.npas4, a1a.npas4, a2a.npas4, a3a.npas4, s1a.npas4, s2a.npas4, s3a.npas4, s4a.npas4, d1a.npas4, d2a.npas4, d3a.npas4, d4a.npas4)
ratio.npas4 <- c(i1r.npas4, i2r.npas4, i3r.npas4, c1r.npas4, c2r.npas4, c3r.npas4, a1r.npas4, a2r.npas4, a3r.npas4, s1r.npas4, s2r.npas4, s3r.npas4, s4r.npas4, d1r.npas4, d2r.npas4, d3r.npas4, d4r.npas4)
ratio.fos <- c(i1r.fos, i2r.fos, i3r.fos, c1r.fos, c2r.fos, c3r.fos, a1r.fos, a2r.fos, a3r.fos, s1r.Fos, s2r.Fos, s3r.Fos, s4r.Fos, d1r.Fos, d2r.Fos, d3r.Fos, d4r.Fos)
ratio.fosb <- c(i1r.fosb, i2r.fosb, i3r.fosb, c1r.fosb, c2r.fosb, c3r.fosb, a1r.fosb, a2r.fosb, a3r.fosb, s1r.Fosb, s2r.Fosb, s3r.Fosb, s4r.Fosb, d1r.Fosb, d2r.Fosb, d3r.Fosb, d4r.Fosb)
ratio.fosl1 <- c(i1r.fosl1, i2r.fosl1, i3r.fosl1, c1r.fosl1, c2r.fosl1, c3r.fosl1, a1r.fosl1, a2r.fosl1, a3r.fosl1, s1r.Fosl1, s2r.Fosl1, s3r.Fosl1, s4r.Fosl1, d1r.Fosl1, d2r.Fosl1, d3r.Fosl1, d4r.Fosl1)
ratio.arc <- c(i1r.arc, i2r.arc, i3r.arc, c1r.arc, c2r.arc, c3r.arc, a1r.arc, a2r.arc, a3r.arc, s1r.Arc, s2r.Arc, s3r.Arc, s4r.Arc, d1r.Arc, d2r.Arc, d3r.Arc, d4r.Arc)
ratio.egr1 <- c(i1r.egr1, i2r.egr1, i3r.egr1, c1r.egr1, c2r.egr1, c3r.egr1, a1r.egr1, a2r.egr1, a3r.egr1, s1r.Egr1, s2r.Egr1, s3r.Egr1, s4r.Egr1, d1r.Egr1, d2r.Egr1, d3r.Egr1, d4r.Egr1)
ratio.egr3 <- c(i1r.egr3, i2r.egr3, i3r.egr3, c1r.egr3, c2r.egr3, c3r.egr3, a1r.egr3, a2r.egr3, a3r.egr3, s1r.Egr3, s2r.Egr3, s3r.Egr3, s4r.Egr3, d1r.Egr3, d2r.Egr3, d3r.Egr3, d4r.Egr3)
ratio.jun <- c(i1r.jun, i2r.jun, i3r.jun, c1r.jun, c2r.jun, c3r.jun, a1r.jun, a2r.jun, a3r.jun, s1r.Jun, s2r.Jun, s3r.Jun, s4r.Jun, d1r.Jun, d2r.Jun, d3r.Jun, d4r.Jun)
ratio.junb <- c(i1r.junb, i2r.junb, i3r.junb, c1r.junb, c2r.junb, c3r.junb, a1r.junb, a2r.junb, a3r.junb, s1r.Junb, s2r.Junb, s3r.Junb, s4r.Junb, d1r.Junb, d2r.Junb, d3r.Junb, d4r.Junb)


dummycode <- c("-1", "-1", "-1", "0", "0", "0", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2")

IEGcounts <- data.frame(ids, phenotype, allneuro, ratio.npas4, ratio.fos, ratio.fosb, ratio.fosl1, ratio.arc, ratio.egr1, ratio.egr3, ratio.jun, ratio.junb, dummycode)
write.csv(IEGcounts, file = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/IEGcounts.csv")

#reformat in excel
ieg <- read.csv(file = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/IEGcounts_long.csv", header = T)

ieg$Group <- factor(ieg$Group, levels = c("INFA", "CONT", "ALLO", "SIRE", "DAM"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Percent_Expressing_IEG.png",
    width = 500, height = 300, units = "mm", res = 600, bg = "white")
ggplot(ieg, aes(x= Group, y= PCT, color = Group, fill = Group)) + facet_wrap(vars(IEG)) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .33, linewidth = 1) + 
  geom_point(size = 3) + 
  scale_y_continuous(transform = "sqrt", limits = c(0,11), breaks = c(0,2 ,5, 10)) + 
  theme_clean() + xlab("") + ylab("% neurons expressing > 0 counts") +
  scale_fill_manual(values = c("#C1A4A8", "#E4D9AA", "#9DD0D0", "#94C3D6", "#BF8AA9")) +
  scale_color_manual(values = c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) +
  theme(legend.position = "bottom", axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        text = element_text(face = "bold", family = "arial", size = 20),
        axis.text.y = element_text(face = "bold", family = "arial", size = 20),
        axis.title.y = element_text(face = "bold", family = "arial", size = 20), 
        legend.text = element_text(face = "bold", family = "arial", size = 20), 
        legend.title = element_blank())
dev.off()

library(car)

Anova(lm(ieg$PCT[which(ieg$IEG == "Arc")] ~ ieg$Group[which(ieg$IEG == "Arc")]), type = "II")
Anova(lm(ieg$PCT[which(ieg$IEG == "Egr1")] ~ ieg$Group[which(ieg$IEG == "Egr1")]), type = "II")
Anova(lm(ieg$PCT[which(ieg$IEG == "Egr3")] ~ ieg$Group[which(ieg$IEG == "Egr3")]), type = "II")
Anova(lm(ieg$PCT[which(ieg$IEG == "Fos")] ~ ieg$Group[which(ieg$IEG == "Fos")]), type = "II")
Anova(lm(ieg$PCT[which(ieg$IEG == "Fosb")] ~ ieg$Group[which(ieg$IEG == "Fosb")]), type = "II") #
Anova(lm(ieg$PCT[which(ieg$IEG == "Fosl")] ~ ieg$Group[which(ieg$IEG == "Fosl")]), type = "II")
Anova(lm(ieg$PCT[which(ieg$IEG == "Jun")] ~ ieg$Group[which(ieg$IEG == "Jun")]), type = "II")
Anova(lm(ieg$PCT[which(ieg$IEG == "Junb")] ~ ieg$Group[which(ieg$IEG == "Junb")]), type = "II")
Anova(lm(ieg$PCT[which(ieg$IEG == "Npas4")] ~ ieg$Group[which(ieg$IEG == "Npas4")]), type = "II")

Anova(lm(ieg$PCT ~ ieg$Group*ieg$IEG))
eta_squared(lm(ieg$PCT ~ ieg$Group*ieg$IEG), partial = TRUE) #eta-sq = 0.17 for Group; = 0.30 for IEG
TukeyHSD(aov(ieg$PCT ~ ieg$Group))
TukeyHSD(aov(ieg$PCT ~ ieg$IEG))

