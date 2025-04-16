### Take the integrated data set that was clustered at the 
### "higher-order" cell types-- i.e., neurons vs. glial sub-populations
### Subset and recluster the neuronal subset.

## Load Packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(harmony)
library(DESeq2)
library(ggthemes)

## Load in Data
parenting <- readRDS(file = '/Volumes/Crucial X6/sequencing/IntermediateRDS/parenting_harmonyintegrated_v2_reclustered.RDS')

parenting$group <- NA
parenting$group[parenting$orig.ident %in% c("MInfa1","MInfa2","MInfa3")] <- "MInfa"
parenting$group[parenting$orig.ident %in% c("RInfa1","RInfa2","RInfa3")] <- "RInfa"
parenting$group[parenting$orig.ident %in% c("RCont1","RCont2","RCont3")] <- "RCont"
parenting$group[parenting$orig.ident %in% c("RAllo1","RAllo2","RAllo3")] <- "RAllo"
parenting$group[parenting$orig.ident %in% c("RSire1","RSire2","RSire3", "RSire4")] <- "RSire"
parenting$group[parenting$orig.ident %in% c("RDam1","RDam2","RDam3", "RDam4")] <- "RDam"
parenting$group <- factor(parenting$group, levels = c("RInfa", "RCont", "RAllo", "RSire", "RDam"))

levels(parenting)

parenting$Agouti <- NA
parenting$Agouti[subset(x = parenting, )]

Astro <- subset(x = parenting, idents = c("Astrocytes"))
Micro <- subset(x = parenting, idents = c("Microglia"))
Oligo <- subset(x = parenting, idents = c("Pre_Oligo", "Myel_Oligo", "New_Oligo"))

FeaturePlot(object = parenting, features = "Npy", order = T, raster = F, label = T, repel = T)

Oligo <- NormalizeData(object = Oligo)
Oligo <- FindVariableFeatures(object = Oligo)
Oligo <- ScaleData(object = Oligo)
Oligo <- RunPCA(object = Oligo)

Oligo <- RunHarmony(Oligo, "group")
harmony.embeddings <- Embeddings(Oligo, reduction = "harmony")

ElbowPlot(Oligo, ndims = 50, reduction = "harmony")

Oligo <- RunUMAP(Oligo, dims = 1:14, reduction = "harmony")
Oligo <- Oligo %>%
  FindNeighbors(reduction = "harmony", dims = 1:14) %>%
  FindClusters(resolution = 0.18, algorithm = 1) 

DimPlot(object = Oligo, label = T)

markers <- FindAllMarkers(object = Oligo, only.pos = T, min.pct = .1, logfc.threshold = .25)
FeaturePlot(Oligo, features = c("Nkain2", "Mobp", "Tnr", "Diaph3", "Fth1", "Ablim2", "Psap", "Top2a"), 
            order = F, label = T, repel = T, ncol = 4)

DimPlot(object = Oligo, group.by = "group", split.by = "orig.ident")

  Astro <- NormalizeData(object = Astro)
Astro <- FindVariableFeatures(object = Astro)
Astro <- ScaleData(object = Astro)
Astro <- RunPCA(object = Astro)

Astro <- RunHarmony(Astro, "group")
harmony.embeddings <- Embeddings(Astro, reduction = "harmony")
#neurons <- JoinLayers(neurons)
levels(Astro)

ElbowPlot(Astro, ndims = 50, reduction = "harmony")

Astro <- RunUMAP(Astro, dims = 1:20, reduction = "harmony")
Astro <- Astro %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.18, algorithm = 1) 

DimPlot(Astro)

markers <- FindAllMarkers(object = Astro, only.pos = T, logfc.threshold = .2, min.pct = 0.10)
View(markers)

Astro1 <- Astro


new.cluster.ids <- c("1", "1", "1", "1", "1", "1", "1")
names(new.cluster.ids) <- levels(Astro1)
Astro1 <- RenameIdents(Astro1, new.cluster.ids)
counts <- Astro1@assays$RNA
metadata <- Astro1@meta.data
metadata$cluster_id <- factor(Astro1@active.ident)
metadata$sample_id <- factor(metadata$orig.ident)

Astro1_bulk_counts <- Astro1 %>% 
  SplitObject(split.by = "orig.ident") %>%  
  lapply(function(x) rowSums(x@assays$RNA$counts)) %>% 
  do.call(cbind, .)

Astro1_bulk_counts


write.csv(Astro1_bulk_counts, "/Volumes/Crucial X6/sequencing/Astro1_bulk_counts.csv")

####

pseudo <- AggregateExpression(Astro1, assays = "RNA", 
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

write.csv(bulk.de, file = "/Volumes/Crucial X6/sequencing/Astro1_bulk_deseq2.csv")
bulk.de <- read.csv( "/Volumes/Crucial X6/sequencing/Astro1_bulk_deseq2.csv", header = TRUE)
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

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Pseudo_Astro1.png",
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

FeaturePlot(object = Astro, features = c("a", "Syt1", "Snap25", "Rbfox3"))

FeaturePlot(object = Astro, features = c("a"), split.by = "group")

## Make subset of cells expressing FOXP3
geneAgouti <- subset(parenting, subset = a > 0)

## Get cell names
cellNames <- rownames(geneAgouti@meta.data)

## Mutate a column in original SO metadata
parenting$barcode <- rownames(parenting)
parenting@meta.data <- parenting@meta.data %>% mutate(AGOUTI = ifelse((parenting$barcode %in% cellNames), "Pos",  "Neg"))


