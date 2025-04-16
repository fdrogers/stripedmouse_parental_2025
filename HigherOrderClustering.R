### Take the integrated data set and recluster into appropriate clusters for
### "higher-order" cell types-- i.e., neurons vs. glial sub-populations.

## Load packages --------------------------------------------------------------
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)

## Import data ----------------------------------------------------------------
parenting <- readRDS(file = "/Volumes/Crucial X6/sequencing/parenting_harmonyintegrated_v1.RDS")

## Look at the original clusters after harmony integration, using a umap
## Resolution was originally set to 0.5

DimPlot(parenting, reduction = "umap", raster = F, pt.size = .001, 
        label = T, label.box = T) + NoLegend() + theme_dark()

elbow1 <- ElbowPlot(parenting, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.25,1.5), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = c(35, 45), colour = "navy", linewidth = .5, linetype = 3)

elbow2 <- ElbowPlot(parenting, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.25,1.5), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = c(35, 45), colour = "navy", linewidth = .5, linetype = 3) + xlim(34,50) + ylim(0,2.5)

elbow1 / elbow2 

parenting <- RunUMAP(parenting, dims = 1:45, reduction = "harmony")
parenting <- parenting %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.9) 

png(filename = "/Volumes/Crucial X6/sequencing/higherorder_umap_harmonyintegrated_v2.png", 
    width = 300, height = 300, units = "mm", res = 300, pointsize = 12, bg = "white")
DimPlot(parenting, label = T, raster = F) + NoLegend()
dev.off()

saveRDS(object = parenting, file = "/Volumes/Crucial X6/sequencing/parenting_harmonyintegrated_v2.RDS")

#https://doi-org.ezproxy.princeton.edu/10.1126/science.aau5324 Moffitt et al. marker genes
#https://doi.org/10.1038/s41467-020-17740-1 Muhl et al. marker genes
#https://doi.org/10.1016/j.celrep.2017.03.004 Chen et al. marker genes
#https://pmc.ncbi.nlm.nih.gov/articles/PMC5008443/ Enpp6 for NFO
#https://www-science-org.ezproxy.princeton.edu/doi/10.1126/science.aau0964 Macrophages Slco2b1

Astrocyte.Features <- c("Aqp4", "S100b", "Gfap", "Sox9", "Gja1", "Slc1a2", "Atp1a2", "Apoe")
Ependymal.Features <- c("Pltp", "Rarres2", "Dbi", "Ccdc153", "Tmem212", "Dynlrb2")
MyelinatedOligodendrocytes.Features <- c("Mobp", "Mog", "Plp1", "Mag", "Ermn", "Cldn11")
Fibroblast.Features <- c("Col1a1", "Col1a2", "Col5a1", "Lcxl1", "Lum", "Fbln1", "Fbln2")
Mural.Features <- c("Cspg4", "Notch3", "Des", "Cdl46", "Rgs5", "Mcam", "Tagln")
NewOligodendrocytes.Features <- c("Enpp6")
OPC.Features <- c("Pdgfra", "Serpine2", "Olig1", "Cspg5", "Lhfpl3")
Endothelial.Features <- c("Cldn5", "Vtn")
Microglia.Features <- c("Cx3cr1", "Aif1", "Spi1", "C3", "C1qb", "C1qa", "Ctss", "Hexb", "Tyrobp")
Macrophage.Features <- c("Pf4", "Lyz2", "Mrc1","Fcer1g", "C1qc")
Neuronal.Features <- c("Rbfox3", "Syt1", "Snap25", "Gad1", "Gad2", "Slc32a1", "Slc17a6", "Slc17a8", "Chat")
Tanycyte.Features <- c("Lhx2", "Ptn", "Prdx6", "Cst3", "Ntrk2", "Col23a1", "Nnat")
Epithelial.Features <- c("Timp3", "Slc9a3r2", "Igbp7", "Ifitm3", "Itm2a", "Slco1a4", "Bsg", "Abcb1a")

FeaturePlot(object = parenting, features = c(Astrocyte.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Ependymal.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(NewOligodendrocytes.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(MyelinatedOligodendrocytes.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Fibroblast.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Mural.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(OPC.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Endothelial.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Microglia.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Macrophage.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Neuronal.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Tanycyte.Features), order = T, label = T)
FeaturePlot(object = parenting, features = c(Epithelial.Features), order = T, label = T)

png(filename = "/Volumes/Crucial X6/sequencing/higherorder_clusters.png", 
    width = 3000, height = 3000, units = "px", res = 300, pointsize = .1, bg = "white")
DotPlot(object = parenting, 
        features = c(Astrocyte.Features, Ependymal.Features, NewOligodendrocytes.Features, MyelinatedOligodendrocytes.Features,
                     Fibroblast.Features, Mural.Features, OPC.Features, Endothelial.Features, Microglia.Features, Macrophage.Features,
                     Neuronal.Features, Tanycyte.Features, Epithelial.Features), cluster.idents = T) + coord_flip() + theme_linedraw()
dev.off()

parenting_small <- BuildClusterTree(object = parenting)
png(filename = "/Volumes/Crucial X6/sequencing/higherorder_clustertree.png", 
    width = 170, height = 170, units = "mm", res = 300, pointsize = 12, bg = "white")
PlotClusterTree(parenting_small, show.tip.label = T, type = "fan", use.edge.length = T)
dev.off()

higherordermarkers <- FindAllMarkers(object = parenting, min.pct = .51, test.use = "wilcox", only.pos = T)
View(higherordermarkers)
write.csv(x = higherordermarkers, file = "/Volumes/Crucial X6/sequencing/higherordermarkers.csv")

markers14v18 <- FindMarkers(object = parenting, ident.1 = 14, ident.2 = 18)
write.csv(x = markers14v18, file = "/Volumes/Crucial X6/sequencing/higherordermarkers_c14vc18.csv")
markers09v24 <- FindMarkers(object = parenting, ident.1 = 9, ident.2 = 24)
write.csv(x = markers09v24, file = "/Volumes/Crucial X6/sequencing/higherordermarkers_c09vc24.csv")
markers06v02 <- FindMarkers(object = parenting, ident.1 = 6, ident.2 = 2)
write.csv(x = markers06v02, file = "/Volumes/Crucial X6/sequencing/higherordermarkers_c06vc02.csv")
markers07v03  <- FindMarkers(object = parenting, ident.1 = 7, ident.2 = 3)
write.csv(x = markers07v03, file = "/Volumes/Crucial X6/sequencing/higherordermarkers_c07vc03.csv")
markers07v21 <- FindMarkers(object = parenting, ident.1 = 7, ident.2 = 21)
write.csv(x = markers07v21, file = "/Volumes/Crucial X6/sequencing/higherordermarkers_c07vc21.csv")
markers03v21 <- FindMarkers(object = parenting, ident.1 = 3, ident.2 = 21)
write.csv(x = markers03v21, file = "/Volumes/Crucial X6/sequencing/higherordermarkers_c03vc21.csv")
markers21v22  <- FindMarkers(object = parenting, ident.1 = 21, ident.2 = 22)
write.csv(x = markers21v22, file = "/Volumes/Crucial X6/sequencing/higherordermarkers_c21vc22.csv")
markers17v22  <- FindMarkers(object = parenting, ident.1 = 17, ident.2 = 22)
write.csv(x = markers17v22, file = "/Volumes/Crucial X6/sequencing/higherordermarkers_c17vc22.csv")
#### Plot selected features from above

png(filename = "/Volumes/Crucial X6/sequencing/higherorder_clusters_featureplot.png", 
    width = 300, height = 150, units = "mm", res = 300, pointsize = 12, bg = "white")
FeaturePlot(object = parenting, features = c("Aqp4", "Tmem212","Cx3cr1", "Cldn5", "Pdgfra", "Enpp6", "Mog", "Rbfox3"), 
            ncol = 4, order = T, cols = c("grey90","black"))
dev.off()

# rename clusters and replot
new_names <- c("Neurons", "Neurons", "Myel_Oligo","Astrocytes", "Neurons", "Neurons", "Myel_Oligo", "Astrocytes", "Neurons",
               "Pre_Oligo", "Neurons", "Neurons", "Neurons", "Neurons", "Microglia", "Neurons", "Endothelial", "Ependymal",
               "Microglia", "New_Oligo", "Neurons", "Astrocytes", "Astrocytes", "Microglia", "Microglia")

names(new_names) <- levels(parenting)

parenting <- RenameIdents(object = parenting, new_names)

saveRDS(object = parenting, 
        file = '/Volumes/Crucial X6/sequencing/IntermediateRDS/parenting_harmonyintegrated_v2_reclustered.RDS')

png(filename = "/Volumes/Crucial X6/sequencing/higherorder_umap_harmonyintegrated_v2_renamed.png", 
    width = 300, height = 300, units = "mm", res = 300, pointsize = 12, bg = "white")
DimPlot(parenting, reduction = "umap", label = T, label.box = T, label.color = "white" ,repel = T, raster = F, order = F,
        cols = c("firebrick", "darkorange3", "gold3", "orange", "turquoise4", "purple", "navy", "darkorange")) + 
  theme_minimal() + NoLegend()
dev.off()

png(filename = "/Volumes/Crucial X6/sequencing/higherorder_umap_harmonyintegrated_v2_renamed.png", 
    width = 150, height = 150, units = "mm", res = 300, pointsize = 12, bg = "white")
DimPlot(parenting, label = F, repel = T, label.box = T, label.color = "white", raster = FALSE, pt.size = .01, 
        cols = c("#000000", "#F0E442", "#56B4E9", "#D55E00", "#009e73", "#0072B2", "#CC79a7", "#E69F00"), order = F) + 
  NoLegend() + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=10))
dev.off()

png(filename = "/Volumes/Crucial X6/sequencing/higherorder_features_harmonyintegrated_v2_renamed.png", 
    width = 150, height = 150, units = "mm", res = 300, pointsize = 12, bg = "white")
FeaturePlot(object = parenting, features = c("Aqp4", "Tmem212","Cx3cr1", "Cldn5", "Pdgfra", "Enpp6", "Mog", "Rbfox3"), 
            ncol = 2, order = F, cols = c("#d6dce9","#5D0000")) +
  NoLegend() + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=10))
dev.off()

png(filename = "/Volumes/Crucial X6/sequencing/higherorder_umap_harmonyintegrated_v2_renamed_bygroup.png", 
    width = 150, height = 50, units = "mm", res = 300, pointsize = 12, bg = "white")
DimPlot(parenting, label = F, repel = T, label.box = T, label.color = "white", raster = FALSE, pt.size = .01, 
        group.by = "orig.ident", split.by = "group",
        cols = c("#40ffff", "#15ffff", "#00eaea", "#f6d649", "#f4cc20", "#dfb70b",
                 "#dc63aa", "#d44097", "#bf2b82", "#9c236a", "#ff153b", "#ea0025", "#bf001f",
                 "#40c9ff", "#15bdff", "#00a7ea", "#0089bf"), order = F) + 
  NoLegend() + theme(plot.title = element_blank())
dev.off()

DimPlot(parenting, reduction = "umap", group.by = "orig.ident", raster = "FALSE", pt.size = .01, split.by = "group",
        cols = c("#00B9B9", "#00B9B9", "#00B9B9", "#DCB50B", "#DCB50B", "#DCB50B",
                 "#BB2A7F", "#BB2A7F", "#BB2A7F", "#BB2A7F", "#960018", "#960018", "#960018",
                 "#008DC5", "#008DC5", "#008DC5", "#008DC5")) + 
  theme_minimal() + theme(legend.position = "bottom")

i1 <- length(Cells(subset(x = parenting, orig.ident == c("RInfa1"))))
i2 <- length(Cells(subset(x = parenting, orig.ident == c("RInfa2"))))
i3 <- length(Cells(subset(x = parenting, orig.ident == c("RInfa3"))))

c1 <- length(Cells(subset(x = parenting, orig.ident == c("RCont1"))))
c2 <- length(Cells(subset(x = parenting, orig.ident == c("RCont2"))))
c3 <- length(Cells(subset(x = parenting, orig.ident == c("RCont3"))))

a1 <- length(Cells(subset(x = parenting, orig.ident == c("RAllo1"))))
a2 <- length(Cells(subset(x = parenting, orig.ident == c("RAllo2"))))
a3 <- length(Cells(subset(x = parenting, orig.ident == c("RAllo3"))))

s1 <- length(Cells(subset(x = parenting, orig.ident == c("RSire1"))))
s2 <- length(Cells(subset(x = parenting, orig.ident == c("RSire2"))))
s3 <- length(Cells(subset(x = parenting, orig.ident == c("RSire3"))))
s4 <- length(Cells(subset(x = parenting, orig.ident == c("RSire4"))))

d1 <- length(Cells(subset(x = parenting, orig.ident == c("RDam1"))))
d2 <- length(Cells(subset(x = parenting, orig.ident == c("RDam2"))))
d3 <- length(Cells(subset(x = parenting, orig.ident == c("RDam3"))))
d4 <- length(Cells(subset(x = parenting, orig.ident == c("RDam4"))))

ids <- c("infa1", "infa2", "infa3", "con1", "con2", "con3", "allo1", "allo2", "allo3", "sire1", "sire2", "sire3", "sire4", "dam1", "dam2", "dam3", "dam4")
phenotype <- c("INFA", "INFA", "INFA", "CONT", "CONT", "CONT", "ALLO", "ALLO", "ALLO", "SIRE", "SIRE", "SIRE", "SIRE", "DAM", "DAM", "DAM", "DAM")
allcells <- c(i1, i2, i3, c1, c2, c3, a1, a2, a3, s1, s2, s3, s4, d1, d2, d3, d4)
sum(i1, i2, i3, c1, c2, c3, a1, a2, a3, s1, s2, s3, s4, d1, d2, d3, d4)

cell.counts <- data.frame(ids, phenotype, allcells)
summary(cell.counts$allcells)

cell.counts$phenotype <- factor(cell.counts$phenotype, levels = c("INFA", "CONT", "ALLO", "SIRE", "DAM"))

library(car)
a <- lm(cell.counts$allcells ~ cell.counts$phenotype)
car::Anova(a, type = "II")
effectsize::eta_squared(a, partial = T)
mean(cell.counts$allcells)

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Cell_Counts.png", 
                          width = 100, height = 80, units = "mm", res = 450)
    ggplot(cell.counts, aes(x=phenotype, y=(allcells)/1000, color = phenotype)) + ylim(0,18) +
      geom_hline(yintercept = c(0, 5,10,15), colour = "grey80", linewidth = .5, linetype = 2) +
      geom_violin() + geom_point(size = 2)  + theme_linedraw() + xlab("") +
      scale_color_manual(values=c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) + 
      ylab("Filtered Cell Count (Thousands)") + ggtitle("Total Sample: 164,507 Cells") + ggthemes::theme_few() +
      theme(axis.text = element_text(size = 12), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            text = element_text(size = 12, face = "bold")) + NoLegend()
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Cell_Counts_alt.png", 
    width = 300, height = 300, units = "mm", res = 450)
ggplot(cell.counts, aes(x=1, y=(allcells)/1000, color = phenotype, group = 1)) + 
  scale_y_continuous(breaks = c(0, 4, 8, 12, 16)) +
  geom_hline(yintercept = c(4, 8, 12), colour = "grey80", linewidth = .5, linetype = 2) +
  geom_violin(alpha = 40, fill = "grey70", linewidth = 3) + geom_jitter(width = .2, size = 20)  + theme_linedraw() + xlab("") +
  scale_color_manual(values=c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) + 
  ylab("") + ggtitle("Total Sample: 164,507 Nuclei") + ggthemes::theme_few() +
  theme(axis.text = element_text(size = 24), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 40), axis.title.y = element_text(size = 40),
        text = element_text(size = 40, face = "bold",), legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size = 33), plot.title = element_text(hjust = 0.5)) 
dev.off()

DimPlot(parenting)
