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

## Get Levels Names

levels(parenting)

neurons <- subset(x = parenting, idents = c("Neurons"))

neurons <- NormalizeData(object = neurons)
neurons <- FindVariableFeatures(object = neurons)
neurons <- ScaleData(object = neurons)
neurons <- RunPCA(object = neurons)

neurons <- RunHarmony(neurons, "group")
harmony.embeddings <- Embeddings(neurons, reduction = "harmony")
#neurons <- JoinLayers(neurons)
levels(neurons)

ElbowPlot(neurons, ndims = 50, reduction = "harmony") + geom_hline(yintercept = c(1.25,1.4))

elbow1 <- 
  ElbowPlot(neurons, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.25,1.1), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = c(35, 40, 50), colour = "navy", linewidth = .5, linetype = 3)

elbow2 <- ElbowPlot(parenting, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.5,1.1), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = c(35, 45,47, 50), colour = "navy", linewidth = .5, linetype = 3) + xlim(30,50) + ylim(0,2.5)

elbow1 / elbow2 

neurons <- RunUMAP(neurons, dims = 1:47, reduction = "harmony")
neurons <- neurons %>%
  FindNeighbors(reduction = "harmony", dims = 1:47) %>%
  FindClusters(resolution = 0.10, algorithm = 1) 

png(filename = "/Volumes/Crucial X6/sequencing/neurons_umap_harmonyintegrated_v2_ognames.png", 
    width = 300, height = 300, units = "mm", res = 300, pointsize = 12, bg = "white")
DimPlot(neurons, reduction = "umap", raster = F, label = T, 
        cols = c("cadetblue2", "cadetblue2", "cadetblue2", "cadetblue2", "firebrick1", "skyblue", "cadetblue2", "cadetblue2", "skyblue", "#FF5533",
                 "#FF5533", "firebrick1", "firebrick1", "firebrick1", "cadetblue2", "skyblue", "skyblue", "skyblue", "skyblue",
                 "skyblue", "gold1", "skyblue", "gold2", "gold3")) + theme_minimal() + NoLegend()
dev.off()

neurons_small <- BuildClusterTree(object = neurons, reduction = "harmony")
png(filename = "/Volumes/Crucial X6/sequencing/neurons_clustertree_v2_newnames.png", 
    width = 300, height = 300, units = "mm", res = 300, pointsize = 12, bg = "white")
PlotClusterTree(neurons_small, show.tip.label = T,  use.edge.length = T, direction = "down", type = "fan")
dev.off()

neurons <- neurons %>%
  FindNeighbors(reduction = "harmony", dims = 1:47) %>%
  FindClusters(resolution = 0.27, algorithm = 3) 

DimPlot(neurons, reduction = "umap", raster = F, label = T) 

markers <- FindAllMarkers(object = neurons, test.use = "roc", only.pos = T, logfc.threshold = .2, min.pct = 0.10)
View(markers)

FeaturePlot(neurons, features = c("Slc17a6", "Slc32a1", "Slc5a7", "Caprin2", "Gnrh1"), 
            label = T, order = F, reduction = "umap", pt.size = .1, ncol = 5, cols = c("pink", "firebrick4"))

DotPlot(object = neurons, features = c("Slc17a6", "Slc32a1", "Slc5a7", "Caprin2", "Gnrh1"), 
        cluster.idents = T, cols = c("navyblue", "red"), scale = F) + coord_flip()

neurons$group <- NA
#neurons$group[neurons$orig.ident %in% c("MInfa1","MInfa2","MInfa3")] <- "MInfa"
neurons$group[neurons$orig.ident %in% c("RInfa1","RInfa2","RInfa3")] <- "RInfa"
neurons$group[neurons$orig.ident %in% c("RCont1","RCont2","RCont3")] <- "RCont"
neurons$group[neurons$orig.ident %in% c("RAllo1","RAllo2","RAllo3")] <- "RAllo"
neurons$group[neurons$orig.ident %in% c("RSire1","RSire2","RSire3", "RSire4")] <- "RSire"
neurons$group[neurons$orig.ident %in% c("RDam1","RDam2","RDam3", "RDam4")] <- "RDam"

neurons$sex <- NA
#neurons$sex[neurons$orig.ident %in% c("MInfa1","MInfa2","MInfa3")] <- "Male"
neurons$sex[neurons$orig.ident %in% c("RInfa1","RInfa2","RInfa3")] <- "Male"
neurons$sex[neurons$orig.ident %in% c("RCont1","RCont2","RCont3")] <- "Male"
neurons$sex[neurons$orig.ident %in% c("RAllo1","RAllo2","RAllo3")] <- "Male"
neurons$sex[neurons$orig.ident %in% c("RSire1","RSire2","RSire3", "RSire4")] <- "Male"
neurons$sex[neurons$orig.ident %in% c("RDam1","RDam2","RDam3", "RDam4")] <- "Female"

saveRDS(neurons, file="/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v3_01202025.RDS")
neurons <- readRDS(file="/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v3_01202025.RDS")

DimPlot(neurons, reduction = "umap", label = T) + NoLegend()

VlnPlot(object = neurons, features = "nCount_RNA", pt.size = 0)
  #Interative inspection of individual clusters shows that C21 is absent across groups in 35% of samples. 
  #Remove and recluster. To improve signal, also remove Chat+ cluster (island identifiable in "parenting")

#Cluster 20 Astrocyte contamination. Remove.

neurons <- subset(x = neurons, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                                          "13", "14", "15", "16", "17", "18", "19", "20", "23"))


saveRDS(neurons, file="/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v4_01222025.RDS")
rm(list = ls())
neurons <- readRDS(file="/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v4_01222025.RDS")

neurons <- NormalizeData(object = neurons)
neurons <- FindVariableFeatures(object = neurons)
neurons <- ScaleData(object = neurons)
neurons <- RunPCA(object = neurons)

neurons <- RunHarmony(neurons, "orig.ident")
harmony.embeddings <- Embeddings(neurons, reduction = "harmony")
#neurons <- JoinLayers(neurons)


ElbowPlot(neurons, ndims = 50, reduction = "harmony") + geom_hline(yintercept = c(1.25,1.4)) + geom_vline(xintercept = c(45))

neurons <- RunUMAP(neurons, dims = 1:45, reduction = "harmony")
neurons <- neurons %>%
  FindNeighbors(reduction = "harmony", dims = 1:45) %>%
  FindClusters(resolution = 0.35, algorithm = 1) 

DimPlot(neurons, reduction = "umap", raster = F, label = T) + theme_minimal() + NoLegend()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Featureplot_Neurons_v2.png",
    width = 350, height = 350, units = "mm", res = 600)
FeaturePlot(neurons, features = c("Slc32a1", "Gad1", "Gad2", "Slc17a6", "Cacna2d1", "Cacna2d2", "Cacna2d3",
                                  "Zeb2", "Meis2", "Kcnab1", "Sox6", "Ebf1", "Trh","Ntng1", "Calcr", "Tcf7l2",
                                  "Rorb", "Gal", "Cdh23", "Foxp2", "Pbx3", "Prdm16", "Alk", "Tafa4", "Tac2", "Caprin2",
                                  "Sim1", "Avp", "Kcnh8", "a"), 
            raster = F, pt.size = .1, order = F, label = F,
            cols = c("grey95","#02341b"), keep.scale = "feature", ncol = 6) +  NoLegend()
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Featureplot_Neurons_v2_ordered.png",
    width = 350, height = 350, units = "mm", res = 600)
FeaturePlot(neurons, features = c("Slc32a1", "Gad1", "Gad2", "Slc17a6", "Cacna2d1", "Cacna2d2", "Cacna2d3",
                                  "Zeb2", "Meis2", "Kcnab1", "Sox6", "Ebf1", "Trh","Ntng1", "Calcr", "Tcf7l2",
                                  "Rorb", "Gal", "Cdh23", "Foxp2", "Pbx3", "Prdm16", "Alk", "Tafa4", "Tac2", "Caprin2",
                                  "Sim1", "Avp", "Kcnh8", "a"), 
            raster = F, pt.size = .1, order = T, label = F,
            cols = c("white","black"), keep.scale = "feature", ncol = 6) +  NoLegend()
dev.off()

markers <- FindAllMarkers(object = neurons, logfc.threshold = .1, only.pos = T, test.use = "roc", min.pct = .1)
View(markers)

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/DotPlot_Neurons_v2.png",
    width = 200, height = 150, units = "mm", res = 600)
DotPlot(neurons, features = c("Slc32a1", "Gad1", "Gad2", "Slc17a6", "Cacna2d1", "Cacna2d2", "Cacna2d3",
                              "Zeb2", "Meis2", "Kcnab1", "Sox6", "Ebf1", "Trh","Ntng1", "Calcr", "Tcf7l2",
                              "Rorb", "Gal", "Cdh23", "Foxp2", "Pbx3", "Prdm16", "Alk", "Tafa4", "Tac2", "Caprin2",
                              "Sim1", "Avp", "Kcnh8", "a"), cols = c("white", "#02341b"),
        scale = F,  scale.by = "radius", dot.min = 0, scale.max = 100, dot.scale = 5, cluster.idents = F) + theme_linedraw() +
  theme(axis.text = element_text(size = 9, face = c("bold")), 
        text = element_blank(), 
        legend.position = "bottom", legend.location = "plot", legend.justification = "center", 
        legend.text =  element_text(size = 10), legend.title = element_text(size = 10, face = "bold")) + coord_flip()
dev.off()

neurons <- RenameIdents(object = neurons, `7` = "merge7_12")
neurons <- RenameIdents(object = neurons, `12` = "merge7_12")

neurons <- RenameIdents(object = neurons, `5` = "merge5_15")
neurons <- RenameIdents(object = neurons, `15` = "merge5_15")

neurons <- RenameIdents(object = neurons, `11` = "merge11_14")
neurons <- RenameIdents(object = neurons, `14` = "merge11_14")

neurons <- RenameIdents(object = neurons, `4` = "merge4_18")
neurons <- RenameIdents(object = neurons, `18` = "merge4_18")

neurons <- RenameIdents(object = neurons, `1` = "merge1_23")
neurons <- RenameIdents(object = neurons, `23` = "merge1_23")

neurons <- RenameIdents(object = neurons, `0` = "merge0_8_16")
neurons <- RenameIdents(object = neurons, `8` = "merge0_8_16")
neurons <- RenameIdents(object = neurons, `16` = "merge0_8_16")

neurons <- RenameIdents(object = neurons, `6` = "merge6_25")
neurons <- RenameIdents(object = neurons, `25` = "merge6_25")

neurons <- RenameIdents(object = neurons, `merge0_8_16` = "gab1") #Not Specific Beyond GABA
neurons <- RenameIdents(object = neurons, `merge6_25` = "gab2") # Zeb2
neurons <- RenameIdents(object = neurons, `merge1_23` = "gab3") # Kcnab1
neurons <- RenameIdents(object = neurons, `merge4_18` = "gab1") #Ptprt not well defined
neurons <- RenameIdents(object = neurons, `merge11_14` = "gab5") #Sox6
neurons <- RenameIdents(object = neurons, `merge5_15` = "glu1") #Cacna2d1
neurons <- RenameIdents(object = neurons, `merge7_12` = "glu2") #Ebf1 and Trh
neurons <- RenameIdents(object = neurons, `2` = "gab6") #Ntng1
neurons <- RenameIdents(object = neurons, `3` = "glu3") # Calcr 
neurons <- RenameIdents(object = neurons, `9` = "glu4") #Tcf7l2
neurons <- RenameIdents(object = neurons, `10` = "gab7") #Rorb and Gal
neurons <- RenameIdents(object = neurons, `13` = "gab8") # Cdh23
neurons <- RenameIdents(object = neurons, `17` = "gab9") #Meis2 and Foxp2 and Pbx3
neurons <- RenameIdents(object = neurons, `19` = "gab10") #Prdm16
neurons <- RenameIdents(object = neurons, `20` = "gab11") #Alk and Tafa4
neurons <- RenameIdents(object = neurons, `21` = "gab12") #Tac2
neurons <- RenameIdents(object = neurons, `22` = "npep") #Caprin2 and Sim1 and Avp
neurons <- RenameIdents(object = neurons, `24` = "gab13") #Kcnh8

#Correct numbering issue
neurons <- RenameIdents(object = neurons, `gab5` = "gab4")
neurons <- RenameIdents(object = neurons, `gab6` = "gab5")
neurons <- RenameIdents(object = neurons, `gab7` = "gab6")
neurons <- RenameIdents(object = neurons, `gab8` = "gab7")
neurons <- RenameIdents(object = neurons, `gab9` = "gab8")
neurons <- RenameIdents(object = neurons, `gab10` = "gab9")
neurons <- RenameIdents(object = neurons, `gab11` = "gab10")
neurons <- RenameIdents(object = neurons, `gab12` = "gab11")
neurons <- RenameIdents(object = neurons, `gab13` = "gab12")

#Reorder
neurons@meta.data$seurat_clusters <- factor(neurons@active.ident, levels = c("gab1", "gab2", "gab3", "gab4", "gab5", "gab6", "gab7", "gab8",
                                                                "gab9", "gab10", "gab11", "gab12", 
                                                                "glu1", "glu2", "glu3", "glu4", "npep"))



#Make updated Dimplot
png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Dimplot_Neurons_v2.png",
    width = 300, height = 300, units = "mm", res = 600)
DimPlot(neurons, label = F, repel = T, label.box = T, label.color = "white", raster = FALSE, pt.size = .01, 
        cols = c("#afafff", "#8f8fff", "#7070ff", "#5050ff", "#3030ff",
                                     "#1010ff", "#0000ef", "#0000cf",
                                     "#0000af", "#00008f", "#000070", "#000050",
                                     "#ff8094", "#ff1a3e", "#e60025", "#b3001d", "#CE910E"), order = F) + NoLegend()
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Dimplot_Neurons_v2.png",
    width = 150, height = 150, units = "mm", res = 600)
DimPlot(neurons, label = F, repel = T, label.box = T, 
        label.color = "white", raster = FALSE, pt.size = .01, order = F, 
        cols = c("magenta", "magenta2", "magenta3", "cyan", 
                 "yellow1", "yellow3", "cyan2", "magenta4", 
                 "cyan3", "cyan4", "yellow4", "darkcyan", 
                 "green1", "green2", "green3", "green4", "darkgreen")) + NoLegend() + xlim(-10,15) + ylim(-13, 13)
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/FeaturePlot_Neurons_v2.png",
    width = 150, height = 150, units = "mm", res = 600)
FeaturePlot(object = neurons, features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2","Caprin2", "Cacna2d1", "Cacna2d2", "Cacna2d3"), 
            ncol = 2, order = F, cols = c("#d6dce9","#5D0000")) +
  NoLegend() + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=10))
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/FeaturePlot_a_Neurons_v2.png",
    width = 500, height = 100, units = "mm", res = 450)
FeaturePlot(object = neurons, features = c("a"), split.by = "group", pt.size = 1,
            ncol = 1, order = F, cols = c("#d6dce9","#5D0000")) +
  NoLegend() + 
  guides(color = guide_legend(override.aes = list(size=2), ncol=10))
dev.off()

saveRDS(neurons, file="/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")
neurons <- readRDS(file = "/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")

markers <- FindAllMarkers(object = neurons, logfc.threshold = .1, only.pos = T, test.use = "roc", min.pct = .1)
View(markers)
write.csv(markers, file = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/neuronal_markers_01232025.csv")

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Dotplot_NeuronClusters_groupedhormones.png",
    width = 400, height = 150, units = "mm", res = 600)
DotPlot(neurons, features = c("Esr1", "Esr2", "Ar", "Prlr", "Oxt", "Oxtr", "Avp", "Avpr1a", "Avpr1b", "Calcr", "Gal", "Trh", "Ghr"), 
        cols = c("#b3e9ff","#5D0000"),
        scale = F,  scale.by = "radius", dot.min = 0, dot.scale = 10, group.by = "group", cluster.idents = F) + theme_linedraw() +
  theme(axis.text = element_text(size = 22),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        text = element_text(size = 20, face = "bold"), 
        legend.text = element_text(size = 10), legend.position = "bottom") 
dev.off()

FeaturePlot(neurons, features = "Avp")

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Dotplot_NeuronClusters_groupedIEG.png",
    width = 400, height = 150, units = "mm", res = 600)
DotPlot(neurons, features = c("Arc", "Egr1", "Egr3", "Fos", "Fosb", "Fosl1", "Jun", "Junb","Npas4"), 
        cols = c("#b3e9ff","#5D0000"),
        scale = F,  scale.by = "radius", dot.min = 0, dot.scale = 10, group.by = "group", cluster.idents = F) + theme_linedraw() +
  theme(axis.text = element_text(size = 22),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        text = element_text(size = 20, face = "bold"), 
        legend.text = element_text(size = 10), legend.position = "bottom") 
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Dotplot_NeuronClusters_groupedIEG_c1only.png",
    width = 175, height = 175, units = "mm", res = 600)
DotPlot(neurons, features = c("Fos", "Fosb","Npas4", "Egr1", "Jun", "Junb"), cols = c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F"), idents = 'npep',
        scale = F,  scale.by = "radius", dot.min = 0, dot.scale = 10, split.by = "group") + theme_linedraw() +
  theme(axis.text = element_text(size = 5), 
        text = element_blank(), 
        legend.position = "bottom", legend.location = "plot", legend.justification = "center") + coord_flip()
dev.off()

##Moffitt Kiss1, Fam19a1, and Gsc were not found
png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Dotplot_MoffittClusters.png",
    width = 500, height = 225, units = "mm", res = 600)
DotPlot(neurons, features = c("Gal", "Calcr","Slc32a1", "Slc17a6", "Amigo2", "Th", "Tac2", "Moxd1", "Pmaip1", "Ucn3", 
                              "Rxfp1", "Esr1", "Prlr", "Ar", "Pgr", "Ddc", "Slc18a2", "Lhx8", "Trh", "Cxcl14",
                              "Adcyap1", "Bdnf", "Otp", "Fezf1", "Eomes", "Tmem163", "C1ql1", "Sncg", "Tacr3",
                              "Angpt1", "Ghrh"), cols = c("white", "magenta4"),
        scale = F,  scale.by = "radius", dot.min = 0, dot.scale = 5) + theme_linedraw() +
  theme(axis.text = element_text(size = 5), 
        text = element_blank(), 
        legend.position = "bottom", legend.location = "plot", legend.justification = "center")
dev.off()

write.csv(table(neurons@meta.data$seurat_clusters, neurons@meta.data$orig.ident), file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_v3.csv")
p <- read.csv(file = "/Volumes/Crucial X6/sequencing/percent_by_cluster_v3.csv", header = T)
View(p)
p$Group <- factor(p$Group, levels = c("RInfa", "RCont", "RAllo", "RSire", "RDam"))
p$Cluster <- as.factor(p$Cluster)
perc.plot <- ggplot(p, aes(x=as.factor(Group), y=Percent, color = Group)) + facet_wrap(. ~ Cluster) +
  geom_boxplot() + ylim(0,40) + 
  xlab("") + ylab("Percent of Nuclei by Cluster") +
  scale_color_manual(values = c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) +
  theme_linedraw() + NoLegend()

perc.plot2 <- ggplot(p, aes(x=totalcontact, y=Percent)) + facet_wrap(. ~ Cluster) +
  geom_point(aes(color = Group), size = 4) + scale_color_manual(values = c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) + 
  geom_smooth(method = "lm", color = "black") + 
  xlab("Total Contact Time (s)") + ylab("Percent of Nuclei by Cluster") +
  theme_linedraw() + theme(legend.position = "bottom")

png(filename = "/Volumes/Crucial X6/sequencing/perc.plots.png", 
    width = 8.5, height = 11, units = "in", res = 400)
perc.plot / perc.plot2
dev.off() 

p$Cluster <- as.factor(p$Cluster)
anova(lm(p$Percent ~ p$high_low*p$Cluster))

perc.plot <- ggplot(p, aes(x=as.factor(high_low), y=Percent, color = high_low)) + facet_wrap(. ~ Cluster) +
  geom_boxplot() + ylim(0,40) + 
  xlab("") + ylab("Percent of Nuclei by Cluster") +
  scale_color_manual(values = c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) +
  theme_linedraw() + NoLegend()

neurons <- SetIdent(neurons, value = "seurat_clusters")
newMarkers <- FindAllMarkers(object = neurons, only.pos = T, min.pct = .51, logfc.threshold = 1)
View(newMarkers)

neurons <- SetIdent(neurons, value = "seurat_clusters")
spec.markers <- FindMarkers(object = neurons, ident.1 = "0", ident.2 = "3", logfc.threshold = 1)
View(spec.markers)

neurons <- SetIdent(neurons, value = "seurat_clusters")
plot0 <- FeaturePlot(object = neurons, features = c('Fos', 'Junb', 'Arc'), 
            order = T, split.by = "group", label = T, repel = T, max.cutoff = 10,
            cols = c("gold1", "magenta4"), alpha = .8, keep.scale = "feature")

neurons <- SetIdent(neurons, value = "group")
plot1 <- DotPlot(object = neurons, features = c('Fos', 'Fosb','Npas4', 'Arc', "Jun", "Junb", "Egr1"), 
                 group.by = "seurat_clusters", scale = F, split.by = "group", 
        idents = c("RInfa"), scale.min = 0, scale.max = 10,
        cols = c("#960018")) + theme_linedraw() + 
  coord_flip()

plot2 <- DotPlot(object = neurons, features = c('Fos', 'Fosb','Npas4', 'Arc', "Jun", "Junb", "Egr1"), 
                 group.by = "seurat_clusters", scale = F, split.by = "group", 
                 idents = c("RCont"), scale.min = 0, scale.max = 10,
                 cols = c("#DCB50B")) + theme_linedraw() + 
  coord_flip()

plot3 <- DotPlot(object = neurons, features = c('Fos', 'Fosb','Npas4', 'Arc', "Jun", "Junb", "Egr1"), 
                 group.by = "seurat_clusters", scale = F, split.by = "group", 
                 idents = c("RAllo"), scale.min = 0, scale.max = 10,
                 cols = c("#00B9B9")) + theme_linedraw() + 
  coord_flip()

plot4 <- DotPlot(object = neurons, features = c('Fos', 'Fosb','Npas4', 'Arc', "Jun", "Junb", "Egr1"), 
                 group.by = "seurat_clusters", scale = F, split.by = "group", 
                 idents = c("RSire"), scale.min = 0, scale.max = 10,
                 cols = c("#008BD5")) + theme_linedraw() + 
  coord_flip()

plot5 <- DotPlot(object = neurons, features = c('Fos', 'Fosb','Npas4', 'Arc', "Jun", "Junb", "Egr1"), 
                 group.by = "seurat_clusters", scale = F, split.by = "group", 
                 idents = c("RDam"), scale.min = 0, scale.max = 10,
                 cols = c("#BB2A7F")) + theme_linedraw() + 
  coord_flip()

png(filename = "/Volumes/Crucial X6/sequencing/dotplots_stacked.png", 
    width = 22, height = 16, units = "in", res = 600)
plot0 / (plot1 / plot2 / plot3 / plot4 / plot5)
dev.off()

## add "high_low"

neurons$contact <- NA
neurons$contact[neurons$orig.ident %in% c("RInfa1","RInfa2","RInfa3", "RCont1","RCont2","RCont3",
                                          "RSire3","RSire4")] <- "Low"
neurons$contact[neurons$orig.ident %in% c("RAllo1","RAllo2","RAllo3", "RSire1", "RSire2",
                                          "RDam1","RDam2","RDam3", "RDam4")] <- "High"

neurons <- SetIdent(neurons, value = "seurat_clusters")
plot0.1 <- FeaturePlot(object = neurons, features = c('Fos', 'Junb', 'Arc'), 
                     order = T, split.by = "contact", label = T, repel = T, max.cutoff = 10,
                     cols = c("gold1", "magenta4"), alpha = .8, keep.scale = "feature")

neurons <- SetIdent(neurons, value = "contact")
plot1.1 <- DotPlot(object = neurons, features = c('Fos', 'Fosb','Npas4', 'Arc', "Jun", "Junb", "Egr1"), 
                 group.by = "seurat_clusters", scale = F, split.by = "contact", 
                 idents = c("High"), scale.min = 0, scale.max = 10,
                 cols = c("#960018")) + theme_linedraw() + 
  coord_flip()

plot2.1 <- DotPlot(object = neurons, features = c('Fos', 'Fosb','Npas4', 'Arc', "Jun", "Junb", "Egr1"), 
                 group.by = "seurat_clusters", scale = F, split.by = "contact", 
                 idents = c("Low"), scale.min = 0, scale.max = 10,
                 cols = c("#DCB50B")) + theme_linedraw() + 
  coord_flip()


png(filename = "/Volumes/Crucial X6/sequencing/dotplots_stacked.png", 
    width = 22, height = 16, units = "in", res = 600)
plot0.1 / (plot1.1 / plot2.1)
dev.off()

high_low <- FindAllMarkers(object = neurons, min.diff.pct = .1, only.pos = T, return.thresh = .01)
View(high_low)
write.csv(high_low, file = "/Volumes/Crucial X6/sequencing/high_low.csv")

DimPlot(object = neurons)
threeversuszero <- FindMarkers(object = neurons, ident.1 = '0', ident.2 = '3', only.pos = F, min.dif.pct = .1)
View(threeversuszero)
write.csv(threeversuszero, file = "/Volumes/Crucial X6/sequencing/cluster3versuscluster0.csv")

DimPlot(object = neurons, group.by = "contact", 
        cols = c("magenta", "orchid4"), 
        alpha = .7, order = T, split.by = "contact")

FeaturePlot(neurons, features = c('a'), order = T, split.by = "contact", 
            cols = c("grey90","skyblue", "blue", "navy"))

#Cluster0 appears increased by releative percgroup.by = #Cluster0 appears increased by releative percentage in Dams, Allopaternal Males, and some sires. 
#Cluster 6 appears increased by relative percentage in Infanticidal and Control males. 

markers.cluster0v6 <- FindMarkers(object = neurons, ident.1 = '0', ident.2 = '6', min.cells.feature = 1,
                                  test.use = 'LR')
View(markers.cluster0v6)
summary(lm(p$Percent[which(p$Cluster == "0")] ~ p$Group[which(p$Cluster == "0")]))


feat_found2 <- FindAllMarkers(neurons, min.pct = .51, only.pos = T, return.thresh = .1, logfc.threshold = .1, 
                             test.use = "roc")
View(feat_found2)
write.csv(feat_found2, file = "/Volumes/Crucial X6/sequencing/features_found2.csv")
FeaturePlot(neurons, features = c("Slc17a6", "Slc32a1", "Chat","Syt1"), order = T, label = T, split.by = "group")

write.csv(table(neurons@meta.data$seurat_clusters, neurons@meta.data$orig.ident), file = "/Volumes/Crucial X6/sequencing/counts_by_cluster.csv")

DotPlot(object = neurons, features = c("Slc17a6", "Slc32a1", "Chat", "Avp"), scale = F) + theme_linedraw()

DotPlot(object = neurons, features = c("Avp"), scale = F, group.by = "orig.ident") + theme_linedraw()
#RAllo1 over expresses Avp; indentify cluster; either cluster 21, 24, or 25. 
VlnPlot(neurons, features = c("Avp"), split.by = "seurat_clusters", group.by = "group", split.plot = T)
#Perhaps both 21 and 25
FeaturePlot(neurons, features = c("Avp"), order = T, label = T, split.by = "sex")

RidgePlot(neurons, features = c("Avp", "Esr1"), 
          group.by = "sex", 
          same.y.lims = F, 
          log = T, layer = "data", sort = T) + NoLegend() 
