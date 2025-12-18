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
Inhibitory <- readRDS(file = "/Volumes/Crucial X6/sequencing/InhibitoryNeuronsLabeled.RDS")

group_colors <- c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")
orig.ident_colors <- c("#960018", "#960018", "#960018",
                       "#DCB50B", "#DCB50B", "#DCB50B",
                       "#00B9B9", "#00B9B9", "#00B9B9",
                       "#008DC5", "#008DC5", "#008DC5", "#008DC5",
                       "#BB2A7F", "#BB2A7F", "#BB2A7F", "#BB2A7F")

#Idents(neurons) <- neurons$seurat_clusters
by_id <- table(Inhibitory@active.ident, Inhibitory@meta.data$orig.ident)
write.csv(by_id, file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_inhibitory_REV2025.csv")
#Reformat in Excel
d <- read.csv(file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_inhibitory_REV2025_reformatted.csv", header = T)
View(d)

d$Group <- factor(d$Group, levels = c("INFA", "CONT", "ALLO", "SIRE", "DAM"))

d$Cluster <- factor(d$Cluster, levels = c("GABA1", "GABA2", "GABA3", "GABA4", "GABA5", "GABA6", "GABA7", "GABA8", "GABA9", "GABA10",
                               "GABA11", "GABA12", "GABA13", "GABA14", "GABA15", "GABA16", "GABA17", "GABA18", "GABA19",
                               "GABA20", "GABA21", "GABA22", "GABA23", "GABA24", "GABA25", "GABA26", "GABA27", "GABA28", "GABA29", 
                               "GABA30", "GABA31", "GABA32", "GABA33", "GABA34", "NPEP", "ACH"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/ExtendedDataFig5_REV25_raw.png",
    width = 500, height = 300, units = "mm", res = 600, bg = "white")
ggplot(data=d, aes(x=Group, y=PCT, colour = Group, fill = Group)) + 
facet_wrap(facets = vars(Cluster), ncol = 4) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .33, linewidth = 1) + 
  geom_point(size = 3) + 
  scale_y_continuous(limits = c(0,24), breaks = c(0,12,24)) + 
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
Anova(lm(d$PCT[which(d$Cluster == "GABA1")] ~ d$Group[which(d$Cluster == "GABA1")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA2")] ~ d$Group[which(d$Cluster == "GABA2")]), type = "II") #p = 0.07
Anova(lm(d$PCT[which(d$Cluster == "GABA3")] ~ d$Group[which(d$Cluster == "GABA3")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA4")] ~ d$Group[which(d$Cluster == "GABA4")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA5")] ~ d$Group[which(d$Cluster == "GABA5")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA6")] ~ d$Group[which(d$Cluster == "GABA6")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA7")] ~ d$Group[which(d$Cluster == "GABA7")]), type = "II") #p = 0.07
Anova(lm(d$PCT[which(d$Cluster == "GABA8")] ~ d$Group[which(d$Cluster == "GABA8")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA9")] ~ d$Group[which(d$Cluster == "GABA9")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA10")] ~ d$Group[which(d$Cluster == "GABA10")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA11")] ~ d$Group[which(d$Cluster == "GABA11")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA12")] ~ d$Group[which(d$Cluster == "GABA12")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA13")] ~ d$Group[which(d$Cluster == "GABA13")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA14")] ~ d$Group[which(d$Cluster == "GABA14")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA15")] ~ d$Group[which(d$Cluster == "GABA15")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA16")] ~ d$Group[which(d$Cluster == "GABA16")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA17")] ~ d$Group[which(d$Cluster == "GABA17")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA18")] ~ d$Group[which(d$Cluster == "GABA18")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA19")] ~ d$Group[which(d$Cluster == "GABA19")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA20")] ~ d$Group[which(d$Cluster == "GABA20")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA21")] ~ d$Group[which(d$Cluster == "GABA21")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA22")] ~ d$Group[which(d$Cluster == "GABA22")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA23")] ~ d$Group[which(d$Cluster == "GABA23")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA24")] ~ d$Group[which(d$Cluster == "GABA24")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA25")] ~ d$Group[which(d$Cluster == "GABA25")]), type = "II") #p = 0.099
Anova(lm(d$PCT[which(d$Cluster == "GABA26")] ~ d$Group[which(d$Cluster == "GABA26")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA27")] ~ d$Group[which(d$Cluster == "GABA27")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA28")] ~ d$Group[which(d$Cluster == "GABA28")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA29")] ~ d$Group[which(d$Cluster == "GABA29")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA30")] ~ d$Group[which(d$Cluster == "GABA30")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA31")] ~ d$Group[which(d$Cluster == "GABA31")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA32")] ~ d$Group[which(d$Cluster == "GABA32")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA33")] ~ d$Group[which(d$Cluster == "GABA33")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GABA34")] ~ d$Group[which(d$Cluster == "GABA34")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "NPEP")] ~ d$Group[which(d$Cluster == "NPEP")]), type = "II") #p = 0.068
Anova(lm(d$PCT[which(d$Cluster == "ACH")] ~ d$Group[which(d$Cluster == "ACH")]), type = "II")

### Excitatory

## Load in Data
Excitatory <- readRDS(file = "/Volumes/Crucial X6/sequencing/ExcitatoryNeuronsReLabeled.RDS")

group_colors <- c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")
orig.ident_colors <- c("#960018", "#960018", "#960018",
                       "#DCB50B", "#DCB50B", "#DCB50B",
                       "#00B9B9", "#00B9B9", "#00B9B9",
                       "#008DC5", "#008DC5", "#008DC5", "#008DC5",
                       "#BB2A7F", "#BB2A7F", "#BB2A7F", "#BB2A7F")

#Idents(neurons) <- neurons$seurat_clusters
by_id <- table(Excitatory@meta.data$seurat_clusters, Excitatory@meta.data$orig.ident)
write.csv(by_id, file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_excitatory_REV2025.csv")
#Reformat in Excel
d <- read.csv(file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_excitatory_REV2025_reformatted.csv", header = T)
View(d)

d$Group <- factor(d$Group, levels = c("INFA", "CONT", "ALLO", "SIRE", "DAM"))

d$Cluster <- factor(d$Cluster, levels = c("GLUT1", "GLUT2", "GLUT3", "GLUT4", "GLUT5", "GLUT6", "GLUT7", "GLUT8", "GLUT9", "GLUT10",
                                          "GLUT11", "GLUT12", "GLUT13", "GLUT14", "GLUT15", "GLUT16", "GLUT17", "GLUT18", "GLUT19",
                                          "GLUT20"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/excitatory_ExtendedDataFig5_REV25_raw.png",
    width = 500, height = 166.67, units = "mm", res = 600, bg = "white")
ggplot(data=d, aes(x=Group, y=PCT, colour = Group, fill = Group)) + 
  facet_wrap(facets = vars(Cluster), ncol = 4) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .33, linewidth = 1) + 
  geom_point(size = 3) + 
  scale_y_continuous(limits = c(0,26.5), breaks = c(0,12,24)) + 
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
Anova(lm(d$PCT[which(d$Cluster == "GLUT1")] ~ d$Group[which(d$Cluster == "GLUT1")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT2")] ~ d$Group[which(d$Cluster == "GLUT2")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT3")] ~ d$Group[which(d$Cluster == "GLUT3")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT4")] ~ d$Group[which(d$Cluster == "GLUT4")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT5")] ~ d$Group[which(d$Cluster == "GLUT5")]), type = "II") #p = 0.05677
Anova(lm(d$PCT[which(d$Cluster == "GLUT6")] ~ d$Group[which(d$Cluster == "GLUT6")]), type = "II") #p = 0.0254
clusterglu6 <- lm(d$PCT[which(d$Cluster == "GLUT6")] ~ d$Group[which(d$Cluster == "GLUT6")])  
eta_squared(clusterglu6, partial = TRUE) #eta-sq = 0.65 for group
TukeyHSD(aov(d$PCT[which(d$Cluster == "GLUT6")] ~ d$Group[which(d$Cluster == "GLUT6")]))

Anova(lm(d$PCT[which(d$Cluster == "GLUT7")] ~ d$Group[which(d$Cluster == "GLUT7")]), type = "II") #p = 0.009
clusterglu7 <- lm(d$PCT[which(d$Cluster == "GLUT7")] ~ d$Group[which(d$Cluster == "GLUT7")])  
eta_squared(clusterglu7, partial = TRUE) #eta-sq = 0.65 for group
TukeyHSD(aov(d$PCT[which(d$Cluster == "GLUT7")] ~ d$Group[which(d$Cluster == "GLUT7")]))

Anova(lm(d$PCT[which(d$Cluster == "GLUT8")] ~ d$Group[which(d$Cluster == "GLUT8")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT9")] ~ d$Group[which(d$Cluster == "GLUT9")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT10")] ~ d$Group[which(d$Cluster == "GLUT10")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT11")] ~ d$Group[which(d$Cluster == "GLUT11")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT12")] ~ d$Group[which(d$Cluster == "GLUT12")]), type = "II") #p = 0.099
Anova(lm(d$PCT[which(d$Cluster == "GLUT13")] ~ d$Group[which(d$Cluster == "GLUT13")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT14")] ~ d$Group[which(d$Cluster == "GLUT14")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT15")] ~ d$Group[which(d$Cluster == "GLUT15")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT16")] ~ d$Group[which(d$Cluster == "GLUT16")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT17")] ~ d$Group[which(d$Cluster == "GLUT17")]), type = "II") #p = 0.05669
Anova(lm(d$PCT[which(d$Cluster == "GLUT18")] ~ d$Group[which(d$Cluster == "GLUT18")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT19")] ~ d$Group[which(d$Cluster == "GLUT19")]), type = "II")
Anova(lm(d$PCT[which(d$Cluster == "GLUT20")] ~ d$Group[which(d$Cluster == "GLUT20")]), type = "II")



#### Subset by active state
inhibitory.active <- subset(Inhibitory, subset = Npas4 > .9 | Fos >.9 | Egr1 > .9 | Arc > .9 | Egr3 > .9 | Fosl1 > .9 | Fosb > .9 | Jun > .9 | Junb > .9)

a <- DimPlot(object = Inhibitory, label = F) + NoLegend()
b <- DimPlot(object = inhibitory.active, label = F) + NoLegend()
a + b

#Idents(neurons) <- neurons$seurat_clusters
by_id.active <- table(inhibitory.active@active.ident, inhibitory.active@meta.data$orig.ident)
write.csv(by_id.active, file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_inhibitoryneuronsactive_rev25.csv")
#Reformat in Excel
d <- read.csv(file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_inhibitory_REV2025_reformatted.csv", header = T)
View(d)

d$Group <- factor(d$Group, levels = c("INFA", "CONT", "ALLO", "SIRE", "DAM"))

d$Cluster <- factor(d$Cluster, levels = c("GABA1", "GABA2", "GABA3", "GABA4", "GABA5", "GABA6", "GABA7", "GABA8", "GABA9", "GABA10",
                                          "GABA11", "GABA12", "GABA13", "GABA14", "GABA15", "GABA16", "GABA17", "GABA18", "GABA19",
                                          "GABA20", "GABA21", "GABA22", "GABA23", "GABA24", "GABA25", "GABA26", "GABA27", "GABA28", "GABA29", 
                                          "GABA30", "GABA31", "GABA32", "GABA33", "GABA34", "NPEP", "ACH"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/inhibitoryExtendedDataFig7_REV25_raw.png",
    width = 500, height = 300, units = "mm", res = 600, bg = "white")
ggplot(data=d, aes(x=Group, y=Active.Pct, colour = Group, fill = Group)) + 
  facet_wrap(facets = vars(Cluster), ncol = 4) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .33, linewidth = 1) + 
  geom_point(size = 3) + 
  scale_y_continuous(limits = c(0,75), breaks = c(0,25,50,75)) + 
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
library(effectsize)
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA1")] ~ d$Group[which(d$Cluster == "GABA1")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA2")] ~ d$Group[which(d$Cluster == "GABA2")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA3")] ~ d$Group[which(d$Cluster == "GABA3")]), type = "II") #p = 0.086

clusterGABA3a <- lm(d$Active.Pct[which(d$Cluster == "GABA3")] ~ d$Group[which(d$Cluster == "GABA3")])  
eta_squared(clusterGABA3a, partial = TRUE) #eta-sq = 0.65 for group
TukeyHSD(aov(d$Active.Pct[which(d$Cluster == "GABA3")] ~ d$Group[which(d$Cluster == "GABA3")]))

Anova(lm(d$Active.Pct[which(d$Cluster == "GABA4")] ~ d$Group[which(d$Cluster == "GABA4")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA5")] ~ d$Group[which(d$Cluster == "GABA5")]), type = "II") #p = 0.0009

clusterGABA5a <- lm(d$Active.Pct[which(d$Cluster == "GABA5")] ~ d$Group[which(d$Cluster == "GABA5")])  
eta_squared(clusterGABA5a, partial = TRUE) #eta-sq = 0.65 for group
TukeyHSD(aov(d$Active.Pct[which(d$Cluster == "GABA5")] ~ d$Group[which(d$Cluster == "GABA5")]))
#Allo > Infa (p = 0.005), Allo > Cont (p = 0.0008), Allo > Sire (p = 0.002), Allo > Dam (p = 0.009)

Anova(lm(d$Active.Pct[which(d$Cluster == "GABA6")] ~ d$Group[which(d$Cluster == "GABA6")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA7")] ~ d$Group[which(d$Cluster == "GABA7")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA8")] ~ d$Group[which(d$Cluster == "GABA8")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA9")] ~ d$Group[which(d$Cluster == "GABA9")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA10")] ~ d$Group[which(d$Cluster == "GABA10")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA11")] ~ d$Group[which(d$Cluster == "GABA11")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA12")] ~ d$Group[which(d$Cluster == "GABA12")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA13")] ~ d$Group[which(d$Cluster == "GABA13")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA14")] ~ d$Group[which(d$Cluster == "GABA14")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA15")] ~ d$Group[which(d$Cluster == "GABA15")]), type = "II") #p = 0.03975

clusterGABA15a <- lm(d$Active.Pct[which(d$Cluster == "GABA15")] ~ d$Group[which(d$Cluster == "GABA15")])  
eta_squared(clusterGABA15a, partial = TRUE) #eta-sq = 0.65 for group
TukeyHSD(aov(d$Active.Pct[which(d$Cluster == "GABA15")] ~ d$Group[which(d$Cluster == "GABA15")]))

Anova(lm(d$Active.Pct[which(d$Cluster == "GABA16")] ~ d$Group[which(d$Cluster == "GABA16")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA17")] ~ d$Group[which(d$Cluster == "GABA17")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA18")] ~ d$Group[which(d$Cluster == "GABA18")]), type = "II") #p = 0.01918

clusterGABA18a <- lm(d$Active.Pct[which(d$Cluster == "GABA18")] ~ d$Group[which(d$Cluster == "GABA18")])  
eta_squared(clusterGABA18a, partial = TRUE) #eta-sq = 0.65 for group
TukeyHSD(aov(d$Active.Pct[which(d$Cluster == "GABA18")] ~ d$Group[which(d$Cluster == "GABA18")]))
#Allo > Cont (p = 0.0119)

Anova(lm(d$Active.Pct[which(d$Cluster == "GABA19")] ~ d$Group[which(d$Cluster == "GABA19")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA20")] ~ d$Group[which(d$Cluster == "GABA20")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA21")] ~ d$Group[which(d$Cluster == "GABA21")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA22")] ~ d$Group[which(d$Cluster == "GABA22")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA23")] ~ d$Group[which(d$Cluster == "GABA23")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA24")] ~ d$Group[which(d$Cluster == "GABA24")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA25")] ~ d$Group[which(d$Cluster == "GABA25")]), type = "II") 
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA26")] ~ d$Group[which(d$Cluster == "GABA26")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA27")] ~ d$Group[which(d$Cluster == "GABA27")]), type = "II") # p = 0.0564
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA28")] ~ d$Group[which(d$Cluster == "GABA28")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA29")] ~ d$Group[which(d$Cluster == "GABA29")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA30")] ~ d$Group[which(d$Cluster == "GABA30")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA31")] ~ d$Group[which(d$Cluster == "GABA31")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA32")] ~ d$Group[which(d$Cluster == "GABA32")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA33")] ~ d$Group[which(d$Cluster == "GABA33")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GABA34")] ~ d$Group[which(d$Cluster == "GABA34")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "NPEP")] ~ d$Group[which(d$Cluster == "NPEP")]), type = "II") 
Anova(lm(d$Active.Pct[which(d$Cluster == "ACH")] ~ d$Group[which(d$Cluster == "ACH")]), type = "II")


#Active Excitatory

excitatory.active <- subset(Excitatory, subset = Npas4 > .9 | Fos >.9 | Egr1 > .9 | Arc > .9 | Egr3 > .9 | Fosl1 > .9 | Fosb > .9 | Jun > .9 | Junb > .9)

a <- DimPlot(object = Excitatory, label = F) + NoLegend()
b <- DimPlot(object = excitatory.active, label = F) + NoLegend()
a + b

#Idents(neurons) <- neurons$seurat_clusters
by_id.active <- table(excitatory.active@active.ident, excitatory.active@meta.data$orig.ident)
write.csv(by_id.active, file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_excitatoryneuronsactive_rev25.csv")

#Reformat in Excel
d <- read.csv(file = "/Volumes/Crucial X6/sequencing/counts_by_cluster_excitatory_REV2025_reformatted.csv", header = T)
View(d)

d$Group <- factor(d$Group, levels = c("INFA", "CONT", "ALLO", "SIRE", "DAM"))

d$Cluster <- factor(d$Cluster, levels = c("GLUT1", "GLUT2", "GLUT3", "GLUT4", "GLUT5", "GLUT6", "GLUT7", "GLUT8", "GLUT9", "GLUT10",
                                          "GLUT11", "GLUT12", "GLUT13", "GLUT14", "GLUT15", "GLUT16", "GLUT17", "GLUT18", "GLUT19",
                                          "GLUT20"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/excitatory_ExtendedDataFig7_REV25_raw.png",
    width = 500, height = 166.67, units = "mm", res = 600, bg = "white")
ggplot(data=d, aes(x=Group, y=Active.Pct, colour = Group, fill = Group)) + 
  facet_wrap(facets = vars(Cluster), ncol = 4) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .33, linewidth = 1) + 
  geom_point(size = 3) + 
  scale_y_continuous(limits = c(0,75), breaks = c(0,25,50,75)) + 
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
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT1")] ~ d$Group[which(d$Cluster == "GLUT1")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT2")] ~ d$Group[which(d$Cluster == "GLUT2")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT3")] ~ d$Group[which(d$Cluster == "GLUT3")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT4")] ~ d$Group[which(d$Cluster == "GLUT4")]), type = "II") #p = 0.091
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT5")] ~ d$Group[which(d$Cluster == "GLUT5")]), type = "II") #p = 0.0068

clusterGLUT5a <- lm(d$Active.Pct[which(d$Cluster == "GLUT5")] ~ d$Group[which(d$Cluster == "GLUT5")])  
eta_squared(clusterGLUT5a, partial = TRUE) 
TukeyHSD(aov(d$Active.Pct[which(d$Cluster == "GLUT5")] ~ d$Group[which(d$Cluster == "GLUT5")]))
#Infa > Cont (p = 0.040), Infa > Dam (p = 0.041), Allo > Cont (p = 0.0457), Allo > Dam (p = 0.0479)

Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT6")] ~ d$Group[which(d$Cluster == "GLUT6")]), type = "II") 
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT7")] ~ d$Group[which(d$Cluster == "GLUT7")]), type = "II") #p = 0.034

clusterGLUT7a <- lm(d$Active.Pct[which(d$Cluster == "GLUT7")] ~ d$Group[which(d$Cluster == "GLUT7")])  
eta_squared(clusterGLUT7a, partial = TRUE) 
TukeyHSD(aov(d$Active.Pct[which(d$Cluster == "GLUT7")] ~ d$Group[which(d$Cluster == "GLUT7")]))
#Allo > Cont (p = 0.0119)

Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT8")] ~ d$Group[which(d$Cluster == "GLUT8")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT9")] ~ d$Group[which(d$Cluster == "GLUT9")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT10")] ~ d$Group[which(d$Cluster == "GLUT10")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT11")] ~ d$Group[which(d$Cluster == "GLUT11")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT12")] ~ d$Group[which(d$Cluster == "GLUT12")]), type = "II") 
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT13")] ~ d$Group[which(d$Cluster == "GLUT13")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT14")] ~ d$Group[which(d$Cluster == "GLUT14")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT15")] ~ d$Group[which(d$Cluster == "GLUT15")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT16")] ~ d$Group[which(d$Cluster == "GLUT16")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT17")] ~ d$Group[which(d$Cluster == "GLUT17")]), type = "II") 
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT18")] ~ d$Group[which(d$Cluster == "GLUT18")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT19")] ~ d$Group[which(d$Cluster == "GLUT19")]), type = "II")
Anova(lm(d$Active.Pct[which(d$Cluster == "GLUT20")] ~ d$Group[which(d$Cluster == "GLUT20")]), type = "II")

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

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/inhibitory_agouti.png",
    units = "mm", width = 100, height = 100, res = 950, bg = "white")
FeaturePlot(object = Inhibitory, features = "a", cols = c("grey90", "#5D0000"), order = T, pt.size = .1, alpha = .75) +
  theme(text = element_blank(), legend.text = element_blank())
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/excitatory_agouti.png",
    units = "mm", width = 100, height = 100, res = 950, bg = "white")
FeaturePlot(object = Excitatory, features = "a", cols = c("grey90", "#5D0000"), order = T, pt.size = .5, alpha = .75) + 
  theme(text = element_blank(), legend.text = element_blank())
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/inhibitory_agouti_dots.png",
    units = "mm", width = 25, height = 100, res = 900, bg = "white")
DotPlot(object = Inhibitory, features = c("a"), scale.min = 5, scale.max = 50, scale = F, group.by = "group", cols = c("navy", "red")) + 
  theme(text = element_blank(), legend.text = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank()) + NoLegend()
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/excitatory_agouti_dots.png",
    units = "mm", width = 25, height = 100, res = 900, bg = "white")
DotPlot(object = Excitatory, features = c("a"), scale.min = 5, scale.max = 50, scale = F, group.by = "group", cols = c("navy", "red")) + 
  theme(text = element_blank(), legend.text = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank()) + NoLegend()
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/excitatory_agouti_dots_forscale.png",
    units = "mm", width = 50, height = 100, res = 900, bg = "white")
DotPlot(object = Excitatory, features = c("a"), scale.min = 5, scale.max = 50, scale = F, group.by = "group", cols = c("navy", "red")) + 
  theme(text = element_blank(), legend.text = element_blank()) 
dev.off()


png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/inhibitory_Mc4r.png",
    units = "mm", width = 100, height = 100, res = 950, bg = "white")
FeaturePlot(object = Inhibitory, features = "Mc4r", cols = c("grey90", "darkgreen"), order = T, pt.size = .5, alpha = .75) + 
  theme(text = element_blank(), legend.text = element_blank())
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/excitatory_Mc4r.png",
    units = "mm", width = 100, height = 100, res = 950, bg = "white")
FeaturePlot(object = Excitatory, features = "Mc4r", cols = c("grey90", "darkgreen"), order = T, pt.size = .5, alpha = .75) + 
  theme(text = element_blank(), legend.text = element_blank())
dev.off()

Inhibitory@active.ident <- factor(Inhibitory@active.ident, levels = c("GABA1", "GABA2", "GABA3", "GABA4", "GABA5", "GABA6", "GABA7", "GABA8", "GABA9", "GABA10",
                                                                      "GABA11", "GABA12", "GABA13", "GABA14", "GABA15", "GABA16", "GABA17", "GABA18", "GABA19",
                                                                      "GABA20", "GABA21", "GABA22", "GABA23", "GABA24", "GABA25", "GABA26", "GABA27", "GABA28", "GABA29", 
                                                                      "GABA30", "GABA31", "GABA32", "GABA33", "GABA34", "NPEP", "ACH"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/inhibitory_Mcr_dots.png",
    units = "mm", width = 100, height = 250, res = 600, bg = "white")
DotPlot(object = Inhibitory, features = c("Mc3r","Mc4r"), scale = F, cols = c("navy", "red")) + theme_linedraw() + theme(text = element_text(family = "arial"))
dev.off()

Excitatory@active.ident <- factor(Excitatory@active.ident, levels = c("GLUT1", "GLUT2", "GLUT3", "GLUT4", "GLUT5", "GLUT6", "GLUT7", "GLUT8", "GLUT9", "GLUT10",
                                                                      "GLUT11", "GLUT12", "GLUT13", "GLUT14", "GLUT15", "GLUT16", "GLUT17", "GLUT18", "GLUT19",
                                                                      "GLUT20"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/ReviseResubmitJune2025/Bioinformatics exports/excitatory_Mcr_dots.png",
    units = "mm", width = 100, height = 250, res = 600, bg = "white")
DotPlot(object = Excitatory, features = c("Mc3r","Mc4r"), scale = F, cols = c("navy", "red")) + theme_linedraw() + theme(text = element_text(family = "arial"))
dev.off()