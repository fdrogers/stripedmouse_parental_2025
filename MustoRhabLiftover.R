library(Seurat)
library(Matrix)
library(irlba)  # For PCA
library(dplyr)
library(ggplot2)
library(ggthemes)

#####
RhabE <- readRDS(, file = "/Volumes/Crucial X6/sequencing/excitatoryneurons_harmonyintegrated_REV2025.RDS")
MusE <- readRDS("~/Desktop/Black6POA/GSE280964_e16_to_p65_C57Bl6j_snRNAseq_excitatory.RDS")

excit.p65 <- subset(x = MusE, subset = MusE@meta.data$age == "p65")
excit.p65

MusE <- excit.p65

Mus.reference <- NormalizeData(MusE)
Mus.reference <- FindVariableFeatures(Mus.reference)
Mus.reference <- ScaleData(Mus.reference)

Rhab.query <- NormalizeData(RhabE)
Rhab.query <- FindVariableFeatures(Rhab.query)
Rhab.query <- ScaleData(Rhab.query)

# find anchors
anchors <- FindTransferAnchors(reference = Mus.reference, query = Rhab.query)

# transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = Mus.reference$cell.type
)
Rhab.query <- AddMetaData(object = Rhab.query, metadata = predictions)
}

Rhab.query$prediction.match <- Rhab.query$predicted.id == Rhab.query$predicted.id
table(Rhab.query$prediction.match)

table(Rhab.query$predicted.id)
VlnPlot(Rhab.query, c("Calcr"), group.by = "predicted.id")

DimPlot(object = Rhab.query, group.by = "predicted.id", label = T) + NoLegend()

Excitatory <- RunHarmony(Rhab.query, c("sex", "group"))
harmony.embeddings <- Embeddings(Excitatory, reduction = "harmony")
#neurons <- JoinLayers(neurons)
levels(Excitatory)

ElbowPlot(Excitatory, ndims = 50, reduction = "harmony") + geom_hline(yintercept = c(1.25,1.4))

elbow1 <- 
  ElbowPlot(Excitatory, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.25,1.1), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = c(35, 40, 50), colour = "navy", linewidth = .5, linetype = 3)

elbow2 <- ElbowPlot(Excitatory, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.5,1.1), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = c(35, 45,47, 50), colour = "navy", linewidth = .5, linetype = 3) + xlim(30,50) + ylim(0,2.5)

elbow1 / elbow2 

#Excitatory; maximum clusters (0-24)
Excitatory <- RunUMAP(Excitatory, dims = 1:47, reduction = "harmony", min.dist = .001) #DecreaseMinDist
Excitatory <- Excitatory %>%
  FindNeighbors(reduction = "harmony", dims = 1:47) %>%
  FindClusters(resolution = 0.7, algorithm = 1) 

a <- DimPlot(Excitatory, reduction = "umap", raster = F, label = T, repel = T) + theme_minimal() + NoLegend()

b <- DimPlot(Excitatory, reduction = "umap", raster = F, label = T, group.by = "predicted.id", repel = T) + theme_minimal() + NoLegend()

a + b

#ExcitatoryMarkers <- FindAllMarkers(object = Excitatory, logfc.threshold = 0.2, only.pos = T, min.pct = .2)
#View(ExcitatoryMarkers)
#write.csv(x = ExcitatoryMarkers, file = "/Volumes/Crucial X6/sequencing/ExcitatoryMarkersMaximumClusters.csv")

excitatoryneurons_large <- BuildClusterTree(object = Excitatory, reduction = "harmony")

##PlotClusterTree(excitatoryneurons_large, show.tip.label = T,  use.edge.length = T, direction = "down", type = "phylogram")

#Excitatory: reduce n of clusters ()
Excitatory <- RunUMAP(Excitatory, dims = 1:47, reduction = "harmony", min.dist = .001) #DecreaseMinDist
Excitatory <- Excitatory %>%
  FindNeighbors(reduction = "harmony", dims = 1:47) %>%
  FindClusters(resolution = 0.475, algorithm = 1) 

a <- DimPlot(Excitatory, reduction = "umap", raster = F, label = T, repel = T, group.by = "seurat_clusters",
             cols = c("red1", "red4", "orange1", "orange4", "yellow", "gold", "gold3", "green1", "green2", "green3", 
                      "cyan4", "cyan2", "green4", "cyan4", "magenta", "purple", "purple3", "pink", "pink3", "orchid")) + theme_minimal() + NoLegend()

b <- DimPlot(Excitatory, reduction = "umap", raster = F, label = T, group.by = "predicted.id", repel = T) + theme_minimal() + NoLegend()

a + b

#ExcitatoryMarkersMed <- FindAllMarkers(object = Excitatory, logfc.threshold = 0.2, only.pos = T, min.pct = .2)
#View(ExcitatoryMarkersMed)
#write.csv(x = ExcitatoryMarkersMed, file = "/Volumes/Crucial X6/sequencing/ExcitatoryMarkersMedClusters.csv")

#excitatoryneurons_med <- BuildClusterTree(object = Excitatory, reduction = "harmony")

#PlotClusterTree(excitatoryneurons_med, show.tip.label = T,  use.edge.length = T, direction = "down", type = "phylogram")

Excitatory <- RenameIdents(object = Excitatory, `0` = "e-C1")
Excitatory <- RenameIdents(object = Excitatory, `1` = "e-N1")
Excitatory <- RenameIdents(object = Excitatory, `2` = "e-M1")
Excitatory <- RenameIdents(object = Excitatory, `3` = "e-L1")
Excitatory <- RenameIdents(object = Excitatory, `4` = "e-H1")
Excitatory <- RenameIdents(object = Excitatory, `5` = "e-M2")
Excitatory <- RenameIdents(object = Excitatory, `6` = "e-M3")
Excitatory <- RenameIdents(object = Excitatory, `7` = "e-L2")
Excitatory <- RenameIdents(object = Excitatory, `8` = "e-H2")
Excitatory <- RenameIdents(object = Excitatory, `9` = "e-A1")
Excitatory <- RenameIdents(object = Excitatory, `10` = "e-A2")
Excitatory <- RenameIdents(object = Excitatory, `11` = "e-H3")
Excitatory <- RenameIdents(object = Excitatory, `12` = "e-N2")
Excitatory <- RenameIdents(object = Excitatory, `13` = "e-A3")
Excitatory <- RenameIdents(object = Excitatory, `14` = "e-C2")
Excitatory <- RenameIdents(object = Excitatory, `15` = "e-T1")
Excitatory <- RenameIdents(object = Excitatory, `16` = "e-H4")
Excitatory <- RenameIdents(object = Excitatory, `17` = "e-L3")
Excitatory <- RenameIdents(object = Excitatory, `18` = "e-Misc")
Excitatory <- RenameIdents(object = Excitatory, `19` = "e-P1")

Excitatory <- RenameIdents(object = Excitatory, `e-C1` = "GLUT1")
Excitatory <- RenameIdents(object = Excitatory, `e-C2` = "GLUT2")
Excitatory <- RenameIdents(object = Excitatory, `e-N1` = "GLUT3")
Excitatory <- RenameIdents(object = Excitatory, `e-N2` = "GLUT4")
Excitatory <- RenameIdents(object = Excitatory, `e-M1` = "GLUT5")
Excitatory <- RenameIdents(object = Excitatory, `e-M2` = "GLUT6")
Excitatory <- RenameIdents(object = Excitatory, `e-M3` = "GLUT7")
Excitatory <- RenameIdents(object = Excitatory, `e-L1` = "GLUT8")
Excitatory <- RenameIdents(object = Excitatory, `e-L2` = "GLUT9")
Excitatory <- RenameIdents(object = Excitatory, `e-L3` = "GLUT10")
Excitatory <- RenameIdents(object = Excitatory, `e-H1` = "GLUT11")
Excitatory <- RenameIdents(object = Excitatory, `e-H2` = "GLUT12")
Excitatory <- RenameIdents(object = Excitatory, `e-H3` = "GLUT13")
Excitatory <- RenameIdents(object = Excitatory, `e-H4` = "GLUT14")
Excitatory <- RenameIdents(object = Excitatory, `e-A1` = "GLUT15")
Excitatory <- RenameIdents(object = Excitatory, `e-A2` = "GLUT16")
Excitatory <- RenameIdents(object = Excitatory, `e-A3` = "GLUT17")
Excitatory <- RenameIdents(object = Excitatory, `e-T1` = "GLUT18")
Excitatory <- RenameIdents(object = Excitatory, `e-Misc` = "GLUT19")
Excitatory <- RenameIdents(object = Excitatory, `e-P1` = "GLUT20")

Excitatory@active.ident <- factor(Excitatory@active.ident, 
                                            levels = c("GLUT1", "GLUT2", "GLUT3", "GLUT4", "GLUT5", "GLUT6", "GLUT7", "GLUT8", "GLUT9", "GLUT10",
                                                       "GLUT11", "GLUT12", "GLUT13", "GLUT14", "GLUT15", "GLUT16", "GLUT17", "GLUT18", "GLUT19", "GLUT20"))

FeaturePlot(object = Excitatory, features = c("a", "Mc3r", "Mc4r"), order = T, 
            cols = c("grey90","darkred"), split.by = "group", pt.size = .1)

VlnPlot(object = Excitatory, features = c("a"), log = F, pt.size = 0, group.by = "group", split.by = "seurat_clusters", 
        cols = c("red1", "red4", "orange1", "orange4", "yellow", "gold", "gold3", "green1", "green2", "green3", 
                 "cyan4", "cyan2", "green4", "cyan4", "magenta", "purple", "purple3", "pink", "pink3", "orchid"))


saveRDS(object = Excitatory, file = "/Volumes/Crucial X6/sequencing/ExcitatoryNeuronsReLabeled.RDS")
Excitatory <- readRDS(file = "/Volumes/Crucial X6/sequencing/ExcitatoryNeuronsReLabeled.RDS")

a <- DimPlot(Excitatory, reduction = "umap", raster = F, label = T, group.by = "seurat_clusters") +  NoLegend()

b <- DimPlot(Excitatory, reduction = "umap", raster = F, label = T, group.by = "predicted.id", repel = T) + NoLegend()

a + b

ExcitatoryMarkers <- FindAllMarkers(object = Excitatory, logfc.threshold = 1, only.pos = T, min.pct = .51)
View(ExcitatoryMarkers)
write.csv(x = ExcitatoryMarkers, file = "/Volumes/Crucial X6/sequencing/ExcitatoryNeuronMarkers.csv")

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/ExcitatoryNeuronMarkers.png",
    width = 900, height = 325, units = "mm", res = 600)
DotPlot(object = Excitatory, features = c("Ebf1", "Npas3",
                                          "Ebf2", "Unc13c",
                                          "Trh", "Nkain3", 
                                          "Adarb2", "Sim1",
                                          "Trps1", "Prlr",
                                          "Calcr", "Cpne4", 
                                          "Rorb", "Hs3st5",
                                          "Slit1", "Cntnap4",
                                          "Dchs2", "Ammecr1",
                                          "Prkd1", "Grik1",
                                          "Il1rapl2", "Nxph1",
                                          "Myo16", "Meis2",
                                          "Erbb4", "Tox",
                                          "Rarb", "Lypd6",
                                          "Nrg1", "Zfp804b", 
                                          "Reln", "Dpyd",
                                          "Adcy2", "Onecut2",
                                          "Epha7", "Arpp21",
                                          "Scd2", "Plp1",
                                          "Pde3a", "Ebf3"),
        scale = F, cols = c("white", "navy"), scale.max = 100, cluster.idents = F) +  theme_linedraw() +
theme(legend.position = "top", 
      text = element_text(family = "arial", size = 30, face = "bold"), 
      axis.text.x = element_text(angle = 90, vjust = 0.90), axis.title.x = element_blank(), axis.title.y = element_blank(),
      legend.text = element_text(size = 30), legend.title = element_text(size = 30))
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/ExcitatoryNeuronswithLabel.png",
    width = 250, height = 200, units = "mm", res = 600)
DimPlot(Excitatory, label = T, repel = T, label.box = T, 
        label.color = "white", raster = FALSE, pt.size = .01, order = F, 
        cols = c("red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4",
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4",
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4", 
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4")) + NoLegend()
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/ExcitatoryNeuronsNoLabel.png",
    width = 250, height = 200, units = "mm", res = 600)
DimPlot(Excitatory, label = F, repel = T, label.box = T, 
        label.color = "white", raster = FALSE, pt.size = .01, order = F, 
        cols = c("red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4",
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4",
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4", 
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4")) + NoLegend()
dev.off()

#####

RhabI <- readRDS(file = "/Volumes/Crucial X6/sequencing/InhibitoryNeuronsSubset.RDS")
MusI <- readRDS("~/Desktop/Black6POA/GSE280964_e16_to_p65_C57Bl6j_snRNAseq_inhibitory.RDS")

inhib.p65 <- subset(x = MusI, subset = MusI@meta.data$age == "p65")
inhib.p65

MusI <- inhib.p65

Mus.reference <- NormalizeData(MusI)
Mus.reference <- FindVariableFeatures(Mus.reference)
Mus.reference <- ScaleData(Mus.reference)

Rhab.query <- NormalizeData(RhabI)
Rhab.query <- FindVariableFeatures(Rhab.query)
Rhab.query <- ScaleData(Rhab.query)

# find anchors
anchors <- FindTransferAnchors(reference = Mus.reference, query = Rhab.query)

# transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = Mus.reference$cell.type)
Rhab.query <- AddMetaData(object = Rhab.query, metadata = predictions)

Rhab.query$prediction.match <- Rhab.query$predicted.id == Rhab.query$predicted.id
table(Rhab.query$prediction.match)

table(Rhab.query$predicted.id)
VlnPlot(Rhab.query, c("Calcr"), group.by = "predicted.id")

DimPlot(object = Rhab.query, group.by = "predicted.id", label = T) + NoLegend()

Inhibitory <- RunHarmony(Rhab.query, c("sex", "group"))
harmony.embeddings <- Embeddings(Inhibitory, reduction = "harmony")
#neurons <- JoinLayers(neurons)
levels(Inhibitory)

ElbowPlot(Inhibitory, ndims = 50, reduction = "harmony") + geom_hline(yintercept = c(1.25,1.4))

elbow1 <- 
  ElbowPlot(Inhibitory, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.25,1.1), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = c(35, 40, 50), colour = "navy", linewidth = .5, linetype = 3)

elbow2 <- ElbowPlot(Inhibitory, ndims = 50, reduction = "harmony") + 
  geom_hline(yintercept = c(1.5,1.1), colour = "grey30", linewidth = .5, linetype = 2) + 
  geom_vline(xintercept = c(35, 45,47, 50), colour = "navy", linewidth = .5, linetype = 3) + xlim(30,50) + ylim(0,2.5)

elbow1 / elbow2 

Inhibitory <- RunUMAP(Inhibitory, dims = 1:50, reduction = "harmony", min.dist = .0005) #DecreaseMinDist
Inhibitory <- Inhibitory %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.85, algorithm = 1) 

a <- DimPlot(Inhibitory, reduction = "umap", raster = F, label = T) + theme_minimal() + NoLegend()

b <- DimPlot(Inhibitory, reduction = "umap", raster = F, label = T, group.by = "predicted.id", repel = T) + theme_minimal() + NoLegend()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/ExtendedDataFig5_REV25_umaps_inhibitory.png",
    width = 500, height = 300, units = "mm", res = 600, bg = "white")
a + b 
dev.off()

neurons_small <- BuildClusterTree(object = Inhibitory, reduction = "harmony")

PlotClusterTree(neurons_small, show.tip.label = T,  use.edge.length = T, direction = "down", type = "fan")

DotPlot(object = Inhibitory, features = c("Chat", "Caprin2"), idents = c("27", "28", "29", "30", "31", "32"))
#Ach cluster is cluster 28; Npep cluster is cluster 30

Inhibitory <- RenameIdents(object = Inhibitory, `0` = "GABA7")
Inhibitory <- RenameIdents(object = Inhibitory, `1` = "GABA27")
Inhibitory <- RenameIdents(object = Inhibitory, `2` = "GABA31")
Inhibitory <- RenameIdents(object = Inhibitory, `3` = "GABA14")
Inhibitory <- RenameIdents(object = Inhibitory, `4` = "GABA15")
Inhibitory <- RenameIdents(object = Inhibitory, `5` = "GABA1")
Inhibitory <- RenameIdents(object = Inhibitory, `6` = "GABA28")
Inhibitory <- RenameIdents(object = Inhibitory, `7` = "GABA18")
Inhibitory <- RenameIdents(object = Inhibitory, `8` = "GABA5")
Inhibitory <- RenameIdents(object = Inhibitory, `9` = "GABA6")
Inhibitory <- RenameIdents(object = Inhibitory, `10` = "GABA19")
Inhibitory <- RenameIdents(object = Inhibitory, `11` = "GABA22")
Inhibitory <- RenameIdents(object = Inhibitory, `12` = "GABA23")
Inhibitory <- RenameIdents(object = Inhibitory, `13` = "GABA10")
Inhibitory <- RenameIdents(object = Inhibitory, `14` = "GABA11")
Inhibitory <- RenameIdents(object = Inhibitory, `15` = "GABA12")
Inhibitory <- RenameIdents(object = Inhibitory, `16` = "GABA16")
Inhibitory <- RenameIdents(object = Inhibitory, `17` = "GABA17")
Inhibitory <- RenameIdents(object = Inhibitory, `18` = "GABA26")
Inhibitory <- RenameIdents(object = Inhibitory, `19` = "GABA8")
Inhibitory <- RenameIdents(object = Inhibitory, `20` = "GABA20")
Inhibitory <- RenameIdents(object = Inhibitory, `21` = "GABA30")
Inhibitory <- RenameIdents(object = Inhibitory, `22` = "GABA33")
Inhibitory <- RenameIdents(object = Inhibitory, `23` = "GABA21")
Inhibitory <- RenameIdents(object = Inhibitory, `24` = "GABA2")
Inhibitory <- RenameIdents(object = Inhibitory, `25` = "GABA25")
Inhibitory <- RenameIdents(object = Inhibitory, `26` = "GABA34")
Inhibitory <- RenameIdents(object = Inhibitory, `27` = "GABA4")
Inhibitory <- RenameIdents(object = Inhibitory, `28` = "ACH")
Inhibitory <- RenameIdents(object = Inhibitory, `29` = "GABA24")
Inhibitory <- RenameIdents(object = Inhibitory, `30` = "NPEP")
Inhibitory <- RenameIdents(object = Inhibitory, `31` = "GABA9")
Inhibitory <- RenameIdents(object = Inhibitory, `32` = "GABA32")
Inhibitory <- RenameIdents(object = Inhibitory, `33` = "GABA3")
Inhibitory <- RenameIdents(object = Inhibitory, `34` = "GABA13")
Inhibitory <- RenameIdents(object = Inhibitory, `35` = "GABA29")

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/InhibitoryNeuronswithLabel.png",
    width = 250, height = 200, units = "mm", res = 600)
DimPlot(Inhibitory, label = T, repel = T, label.box = T, 
        label.color = "white", raster = FALSE, pt.size = .01, order = F, 
        cols = c("red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4",
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4",
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4", 
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4")) + NoLegend()
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/InhibitoryNeuronsNoLabel.png",
    width = 250, height = 200, units = "mm", res = 600)
DimPlot(Inhibitory, label = F, repel = T, label.box = T, 
        label.color = "white", raster = FALSE, pt.size = .01, order = F, 
        cols = c("red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4",
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4",
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4", 
                 "red", "blue", "gold", "red2", "blue2", "gold2", "darkred", "darkblue", "gold4")) + NoLegend()
dev.off()

saveRDS(object = Inhibitory, file = "/Volumes/Crucial X6/sequencing/InhibitoryNeuronsLabeled.RDS")
Inhibitory <- readRDS(file = "/Volumes/Crucial X6/sequencing/InhibitoryNeuronsLabeled.RDS")

InhibitoryMarkers <- FindAllMarkers(object = Inhibitory, logfc.threshold = .7, only.pos = T, min.pct = .4)
View(InhibitoryMarkers)
write.csv(x = InhibitoryMarkers, file = "/Volumes/Crucial X6/sequencing/InhibitoryNeuronMarkers.csv")

Inhibitory@active.ident <- factor(Inhibitory@active.ident, levels = c("GABA1", "GABA2", "GABA3", "GABA4", "GABA5", "GABA6", "GABA7", "GABA8", "GABA9", "GABA10",
                                                                      "GABA11", "GABA12", "GABA13", "GABA14", "GABA15", "GABA16", "GABA17", "GABA18", "GABA19",
                                                                      "GABA20", "GABA21", "GABA22", "GABA23", "GABA24", "GABA25", "GABA26", "GABA27", "GABA28", "GABA29", 
                                                                      "GABA30", "GABA31", "GABA32", "GABA33", "GABA34", "NPEP", "ACH"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/InhibitoryNeuronMarkers.png",
    width = 900, height = 525, units = "mm", res = 600)
DotPlot(object = Inhibitory, features = c("Ntng1", "Pdzrn4", "Satb2", "Tac2",
                                          "Pou6f2", "Grm8", "Pcdh15", "Ephb1",
                                          "Calcr", "Rorb", "Zfpm2", "Ecel1",
                                          "Tmem130", "Slc22a17", "Nrg1", "Kctd8",
                                          "Ghr", "Maml2", "Sh3rf3", "Vwc2l",
                                          "Nfia", "Nfib", "Col25a1", "Sdk1",
                                          "Arl2bp", "Impad1", "Nkain3", "Nxph1",
                                          "Ptprt", "Sorcs1", "Sox5",
                                          "Pard3b", "Rasgef1b","Mast4", "Tafa2",
                                          "Gfra1", "Prkca","Egfem1", "Plcxd3",
                                          "L3mbtl4", "Kcnh8","Gbx1", "Pde3a",
                                          "Sox6", "Stxbp6","Npy", "Sst",
                                          "Zeb2", "Foxp2","Tshz1", "Slc4a4",
                                          "Meis2", "Pbx3", "Rarb", "Fbn2",
                                          "Tmem132c", "Scn5a", "Gpc5", "Gabra1",
                                          "Prr16", "Thsd7b", "Ccbe1", "Zfp536",
                                          "Tafa4", "Adamtsl1", "Drd3", "Fstl1",
                                          "Chat", "Slc5a7", "Arhgap15", "Avp", "Caprin2"),
        scale = F, cols = c("white", "navy"), scale.max = 100, cluster.idents = F) + theme_linedraw() +
  theme(legend.position = "top", 
        text = element_text(family = "arial", size = 30, face = "bold"), 
        axis.text.x = element_text(angle = 90, vjust = 0.90), axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.text = element_text(size = 30), 
        legend.title = element_text(size = 30))
dev.off()


