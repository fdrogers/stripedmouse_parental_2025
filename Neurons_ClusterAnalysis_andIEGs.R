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

## Load in Data
neurons <- readRDS(file="/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")

group_colors <- c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")
orig.ident_colors <- c("#960018", "#960018", "#960018",
                       "#DCB50B", "#DCB50B", "#DCB50B",
                       "#00B9B9", "#00B9B9", "#00B9B9",
                       "#008DC5", "#008DC5", "#008DC5", "#008DC5",
                       "#BB2A7F", "#BB2A7F", "#BB2A7F", "#BB2A7F")

#Idents(neurons) <- neurons$seurat_clusters
by_id <- table(neurons@active.ident, neurons@meta.data$orig.ident)
write.csv(by_id, file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_01292025.csv")
#Reformat in Excel
d <- read.csv(file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_01292025_reformatted.csv", header = T)
View(d)

d$Group <- factor(d$Group, levels = c("INFA", "CONT", "ALLO", "SIRE", "DAM"))

d$Cluster <- factor(d$Cluster, levels = c("gab1", "gab2", "gab3", "gab4", "gab5", "gab6", "gab7", "gab8", "gab9", "gab10",
                               "gab11", "gab12", "glu1", "glu2", "glu3", "glu4", "npep"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/ExtendedDataFig5_raw.png",
    width = 500, height = 300, units = "mm", res = 600, bg = "white")
ggplot(data=d, aes(x=Group, y=PCT, colour = Group, fill = Group)) + 
facet_wrap(facets = vars(Cluster), ncol = 6) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .33, linewidth = 1) + 
  geom_point(size = 3) + 
  scale_y_continuous(transform = "sqrt", limits = c(0,50), breaks = c(0,5,25,50)) + 
  theme_clean() + xlab("") + ylab("% of sample") +
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
Anova(lm(d$PCT[which(d$Cluster == "gab1")] ~ d$Group[which(d$Cluster == "gab1")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab2")] ~ d$Group[which(d$Cluster == "gab2")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab3")] ~ d$Group[which(d$Cluster == "gab3")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab4")] ~ d$Group[which(d$Cluster == "gab4")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab5")] ~ d$Group[which(d$Cluster == "gab5")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab6")] ~ d$Group[which(d$Cluster == "gab6")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab7")] ~ d$Group[which(d$Cluster == "gab7")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab8")] ~ d$Group[which(d$Cluster == "gab8")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab9")] ~ d$Group[which(d$Cluster == "gab9")]), type = "II")
  clustergab9 <- lm(d$PCT[which(d$Cluster == "gab9")] ~ d$Group[which(d$Cluster == "gab9")])  
eta_squared(clustergab9, partial = TRUE) #eta-sq = 0.64 for cohort
TukeyHSD(aov(d$PCT[which(d$Cluster == "gab9")] ~ d$Group[which(d$Cluster == "gab9")]))
#Dams and Allo had marginally smaller differences in PCT compared to control group: 
  #Dams - 0.9344%, padj 0.0269
  #Allo - 0.9779%, padj 0.0306

Anova(lm(d$PCT[which(d$Cluster == "gab10")] ~ d$Group[which(d$Cluster == "gab10")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab11")] ~ d$Group[which(d$Cluster == "gab11")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab12")] ~ d$Group[which(d$Cluster == "gab12")]), type = "II")

Anova(lm(d$PCT[which(d$Cluster == "glu1")] ~ d$Group[which(d$Cluster == "glu1")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "glu2")] ~ d$Group[which(d$Cluster == "glu2")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "glu3")] ~ d$Group[which(d$Cluster == "glu3")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "glu4")] ~ d$Group[which(d$Cluster == "glu4")]), type = "II")
  clusterglu4 <- lm(d$PCT[which(d$Cluster == "glu4")] ~ d$Group[which(d$Cluster == "glu4")])  
eta_squared(clusterglu4, partial = TRUE) #eta-sq = 0.65 for group
TukeyHSD(aov(d$PCT[which(d$Cluster == "glu4")] ~ d$Group[which(d$Cluster == "glu4")]))
#Sires had marginally smaller differences in PCT compared to allo: 
#Sire - 2.7325%, padj 0.00469

Anova(lm(d$PCT[which(d$Cluster == "npep")] ~ d$Group[which(d$Cluster == "npep")]), type = "II")

#### Subset by active state
neurons.active <- subset(neurons, subset = Npas4 > 0 | Fos > 0 | Egr1 > 0 | Arc > 0 | Egr3 > 0 | Fosl1 > 0 | Fosb > 0 | Jun > 0 | Junb > 0)

#Idents(neurons) <- neurons$seurat_clusters
by_id.active <- table(neurons.active@active.ident, neurons.active@meta.data$orig.ident)
write.csv(by_id.active, file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_active_01292025.csv")
#Reformat in Excel
d <- read.csv(file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_active_01292025_reformatted.csv", header = T)
View(d)

d$Group <- factor(d$Group, levels = c("INFA", "CONT", "ALLO", "SIRE", "DAM"))

d$Cluster <- factor(d$Cluster, levels = c("gab1", "gab2", "gab3", "gab4", "gab5", "gab6", "gab7", "gab8", "gab9", "gab10",
                                          "gab11", "gab12", "glu1", "glu2", "glu3", "glu4", "npep"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/ExtendedDataFig6_raw.png",
    width = 500, height = 300, units = "mm", res = 600, bg = "white")
ggplot(data=d, aes(x=Group, y=PCT_ACTIVE, colour = Group, fill = Group)) + 
  facet_wrap(facets = vars(Cluster), ncol = 6) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .33, linewidth = 1) + 
  geom_point(size = 3) + 
  scale_y_continuous(transform = "sqrt", limits = c(0,100), breaks = c(0,5,25,50,100)) + 
  theme_clean() + xlab("") + ylab("% of cluster active") +
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
Anova(lm(d$PCT[which(d$Cluster == "gab1")] ~ d$Group[which(d$Cluster == "gab1")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab2")] ~ d$Group[which(d$Cluster == "gab2")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab3")] ~ d$Group[which(d$Cluster == "gab3")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab4")] ~ d$Group[which(d$Cluster == "gab4")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab5")] ~ d$Group[which(d$Cluster == "gab5")]), type = "II")#
  clustergab5 <- lm(d$PCT[which(d$Cluster == "gab5")] ~ d$Group[which(d$Cluster == "gab5")])  
eta_squared(clustergab5, partial = TRUE) #eta-sq = 0.65 for group
TukeyHSD(aov(d$PCT[which(d$Cluster == "gab5")] ~ d$Group[which(d$Cluster == "gab5")]))
# Alloparents had greater activity than infa up 8.652%, p = 0.0493
# Alloparents had greater activity than control up 11.2957% p = 0.00925
# Alloparents had greater activity than sires up 10.22% p = 0.0117

Anova(lm(d$PCT[which(d$Cluster == "gab6")] ~ d$Group[which(d$Cluster == "gab6")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab7")] ~ d$Group[which(d$Cluster == "gab7")]), type = "II")#
clustergab7 <- lm(d$PCT[which(d$Cluster == "gab7")] ~ d$Group[which(d$Cluster == "gab7")])  
eta_squared(clustergab7, partial = TRUE) #eta-sq = 0.54 for group
TukeyHSD(aov(d$PCT[which(d$Cluster == "gab7")] ~ d$Group[which(d$Cluster == "gab7")]))
# Alloparents had greater activity than control up 8.395%, p = 0.0334

Anova(lm(d$PCT[which(d$Cluster == "gab8")] ~ d$Group[which(d$Cluster == "gab8")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab9")] ~ d$Group[which(d$Cluster == "gab9")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab10")] ~ d$Group[which(d$Cluster == "gab10")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab11")] ~ d$Group[which(d$Cluster == "gab11")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "gab12")] ~ d$Group[which(d$Cluster == "gab12")]), type = "II")

Anova(lm(d$PCT[which(d$Cluster == "glu1")] ~ d$Group[which(d$Cluster == "glu1")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "glu2")] ~ d$Group[which(d$Cluster == "glu2")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "glu3")] ~ d$Group[which(d$Cluster == "glu3")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "glu4")] ~ d$Group[which(d$Cluster == "glu4")]), type = "II")

Anova(lm(d$PCT[which(d$Cluster == "npep")] ~ d$Group[which(d$Cluster == "npep")]), type = "II")

###

DotPlot(neurons, features = c("Avp", "Sox6"), scale = T, cols = c("cyan", "magenta"))

DimPlot(neurons, cols = c("magenta", "magenta", "magenta", "cyan", 
                          "yellow2", "yellow2", "cyan", "magenta", 
                          "cyan", "cyan", "yellow2", "cyan", 
                          "green", "green", "green", "green", "green"), label = T) 

IEGlist <- list(c("Arc", "Egr1", "Egr3", "Fos", "Fosb", "Fosl1", "Jun", "Junb","Npas4"))
neurons <- AddModuleScore(object = neurons, features = IEGlist, name = "IEG_composit")
FeaturePlot(object = neurons, features = "IEG_composit1", order = T)


## Supplemental Figure 1
png(filename = "/Volumes/Macintosh HD/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/SupplementalFigure1.png", 
    units = "mm", width = 350, height = 400, res = 600, bg = "white")
DotPlot(neurons, features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2", "Caprin2","Cacna2d1", "Cacna2d2", "Cacna2d3", 
                              "Zeb2", "Grm7", "Foxp2",
                              "Celf2", "Meis2", "Rbfox1",
                              "Sox6", "Bcl11a", "Adarb2",
                              "Ntng1", "Pdzrn4",
                              "Rorb", "Zpf804b", "Zfhx4",
                              "Unc5d", "Cdh18", "Zfhx3",
                              "Erbb4",
                              "Prkca", "Megf11", "Prdm16",
                              "Alk", "Pbx3", "Adcy2",
                              "Tac2", "Grik1", "Cadps2",
                              "Nrg3", "Gria4",
                              "Ebf1", "Gpc6", "Ralyl",
                              "Calcr", "Il1rapl1",
                              "Nxph1", "Galntl6",
                              "Rbms3", "Pde4b", 
                              "Esr1", "Esr2", "Ar", "Prlr", "Oxt", "Oxtr", "Avp", "Avpr1a", "Avpr1b", "Gal", "Trh", "Ghr", "Brs3"),
        cols = c("#b3e9ff","#5D0000"), scale = F, scale.min = 0, scale.max = 100, cluster.idents = F) + 
  coord_flip() + theme_linedraw() + theme(axis.text = element_text(size = 16, family = "arial", face = "italic"), 
                                          axis.title = element_blank(), legend.position = "right", 
                                          legend.text = element_text(size = 16, family = "arial"), 
                                          legend.title = element_text(size = 16, family = "arial", face = "bold"))
dev.off()

png(filename = "/Volumes/Macintosh HD/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/SupplementalFigure2a.png", 
    units = "mm", width = 375, height = 100, res = 600, bg = "white")
DotPlot(neurons, features = c("Arc", "Egr1", "Egr3", "Fos", "Fosb", "Fosl1", "Jun", "Junb","Npas4"),
        cols = c("#b3e9ff","#5D0000"), scale = F, cluster.idents = F) + 
  coord_flip() + theme_linedraw() + theme(axis.text = element_text(size = 16, family = "arial", face = "italic"), 
                                          axis.title = element_blank(), legend.position = "right", 
                                          legend.text = element_text(size = 16, family = "arial"), 
                                          legend.title = element_text(size = 16, family = "arial", face = "bold"))
dev.off()

png(filename = "/Volumes/Macintosh HD/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/SupplementalFigure2b.png", 
    units = "mm", width = 375, height = 100, res = 600, bg = "white")
DotPlot(neurons, features = c("Arc", "Egr1", "Egr3", "Fos", "Fosb", "Fosl1", "Jun", "Junb","Npas4"), group.by = "group",
        cols = c("#b3e9ff","#5D0000"), scale = F, cluster.idents = F) + 
  coord_flip() + theme_linedraw() + theme(axis.text = element_text(size = 16, family = "arial", face = "italic"), 
                                          axis.title = element_blank(), legend.position = "right", 
                                          legend.text = element_text(size = 16, family = "arial"), 
                                          legend.title = element_text(size = 16, family = "arial", face = "bold"))
dev.off()

#split into allo vs infa
allo <- subset(x = neurons, orig.ident == c("RAllo1", "RAllo2", "RAllo3"))
infa <- subset(x = neurons, orig.ident == c("RInfa1", "RInfa2", "RInfa3"))
cont <- subset(x = neurons, orig.ident == c("RCont1", "RCont2", "RCont3"))
sire <- subset(x = neurons, orig.ident == c("RSire1", "RSire2", "RSire3", "RSire4"))
dams <- subset(x = neurons, orig.ident == c("RDam1", "RDam2", "RDam3", "RDam4"))

png(filename = "/Volumes/Macintosh HD/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/SupplementalFigure2allo.png", 
    units = "mm", width = 500, height = 100, res = 600, bg = "white")
DotPlot(allo, features = c("IEG_composit1"), scale.min = 0, scale.max = 10, 
        cols = c("#b3e9ff","#5D0000"), scale = F, cluster.idents = F) + 
  coord_flip() + theme_linedraw() + theme(axis.text = element_text(size = 16, family = "arial", face = "italic"), 
                                          axis.title = element_blank(), legend.position = "right", 
                                          legend.text = element_text(size = 16, family = "arial"), 
                                          legend.title = element_text(size = 16, family = "arial", face = "bold"), 
                                          axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

png(filename = "/Volumes/Macintosh HD/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/SupplementalFigure2infa.png", 
    units = "mm", width = 500, height = 100, res = 600, bg = "white")
DotPlot(infa, features = c("IEG_composit1"), scale.min = 0, scale.max = 10, 
        cols = c("#b3e9ff","#5D0000"), scale = F, cluster.idents = F) + 
  coord_flip() + theme_linedraw() + theme(axis.text = element_text(size = 16, family = "arial", face = "italic"), 
                                          axis.title = element_blank(), legend.position = "right", 
                                          legend.text = element_text(size = 16, family = "arial"), 
                                          legend.title = element_text(size = 16, family = "arial", face = "bold"), 
                                          axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

png(filename = "/Volumes/Macintosh HD/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/SupplementalFigure2cont.png", 
    units = "mm", width = 500, height = 100, res = 600, bg = "white")
DotPlot(cont, features = c("IEG_composit1"), scale.min = 0, scale.max = 10, 
        cols = c("#b3e9ff","#5D0000"), scale = F, cluster.idents = F) + 
  coord_flip() + theme_linedraw() + theme(axis.text = element_text(size = 16, family = "arial", face = "italic"), 
                                          axis.title = element_blank(), legend.position = "right", 
                                          legend.text = element_text(size = 16, family = "arial"), 
                                          legend.title = element_text(size = 16, family = "arial", face = "bold"), 
                                          axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

png(filename = "/Volumes/Macintosh HD/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/SupplementalFigure2sire.png", 
    units = "mm", width = 500, height = 100, res = 600, bg = "white")
DotPlot(sire, features = c("IEG_composit1"), scale.min = 0, scale.max = 10, 
        cols = c("#b3e9ff","#5D0000"), scale = F, cluster.idents = F) + 
  coord_flip() + theme_linedraw() + theme(axis.text = element_text(size = 16, family = "arial", face = "italic"), 
                                          axis.title = element_blank(), legend.position = "right", 
                                          legend.text = element_text(size = 16, family = "arial"), 
                                          legend.title = element_text(size = 16, family = "arial", face = "bold"), 
                                          axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

png(filename = "/Volumes/Macintosh HD/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/SupplementalFigure2dams.png", 
    units = "mm", width = 500, height = 100, res = 600, bg = "white")
DotPlot(dams, features = c("IEG_composit1"), scale.min = 0, scale.max = 10, 
        cols = c("#b3e9ff","#5D0000"), scale = F, cluster.idents = F) + 
  coord_flip() + theme_linedraw() + theme(axis.text = element_text(size = 16, family = "arial", face = "italic"), 
                                          axis.title = element_blank(), legend.position = "right", 
                                          legend.text = element_text(size = 16, family = "arial"), 
                                          legend.title = element_text(size = 16, family = "arial", face = "bold"), 
                                          axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()