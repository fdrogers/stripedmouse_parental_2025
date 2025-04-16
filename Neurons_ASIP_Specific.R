### Take the integrated data set that was clustered at the "neuronal"
### cell types and look at ASIP specifically.

## Load Packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)

## Load in Data
neurons <- readRDS(file="/Volumes/Crucial X6/sequencing/neurons_harmonyintegrated_v5_01222025.RDS")

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/FeaturePlot_Asip.png",
    width = 250, height = 300, units = "mm", res = 600)
FeaturePlot(object = neurons, features = c("a"),  pt.size = .5,
            ncol = 1, order = T, cols = c("#d6dce9","#5D0000")) + ggthemes::theme_few() +
  theme(text = element_text(size = 30), plot.title = element_blank(), legend.position = "right",
        axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold", size = 40)) 
dev.off()


#split into allo vs infa
allo <- subset(x = neurons, orig.ident == c("RAllo1", "RAllo2", "RAllo3"))
infa <- subset(x = neurons, orig.ident == c("RInfa1", "RInfa2", "RInfa3"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/FeaturePlot_Asip_Allo.png",
    width = 250, height = 150, units = "mm", res = 600)
FeaturePlot(object = allo, features = c("a"),  pt.size = 2,
            ncol = 1, order = F, cols = c("#d6dce9","#5D0000")) + ggthemes::theme_few() +
  theme(text = element_text(size = 30), plot.title = element_blank(), legend.position = "right",
        axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold", size = 40)) 
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/FeaturePlot_Asip_Infa.png",
    width = 250, height = 150, units = "mm", res = 600)
FeaturePlot(object = infa, features = c("a"),  pt.size = 2, 
            ncol = 1, order = F, cols = c("#d6dce9","#5D0000")) + ggthemes::theme_few() +
  theme(text = element_text(size = 30), plot.title = element_blank(), legend.position = "right",
        axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold", size = 40)) 
dev.off()

a1 <- length(Cells(subset(x = neurons, orig.ident == c("RAllo1"))))
a1a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RAllo1"))))
a1r <- round((a1a / a1)*100, digits = 1)
a2 <- length(Cells(subset(x = neurons, orig.ident == c("RAllo2"))))
a2a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RAllo2"))))
a2r <- round((a2a / a2)*100, digits = 1)
a3 <- length(Cells(subset(x = neurons, orig.ident == c("RAllo3"))))
a3a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RAllo3"))))
a3r <- round((a3a / a3)*100, digits = 1)

c1 <- length(Cells(subset(x = neurons, orig.ident == c("RCont1"))))
c1a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RCont1"))))
c1r <- round((c1a / c1)*100, digits = 1)
c2 <- length(Cells(subset(x = neurons, orig.ident == c("RCont2"))))
c2a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RCont2"))))
c2r <- round((c2a / c2)*100, digits = 1)
c3 <- length(Cells(subset(x = neurons, orig.ident == c("RCont3"))))
c3a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RCont3"))))
c3r <- round((c3a / c3)*100, digits = 1)

i1 <- length(Cells(subset(x = neurons, orig.ident == c("RInfa1"))))
i1a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RInfa1"))))
i1r <- round((i1a / i1)*100, digits = 1)
i2 <- length(Cells(subset(x = neurons, orig.ident == c("RInfa2"))))
i2a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RInfa2"))))
i2r <- round((i2a / i2)*100, digits = 1)
i3 <- length(Cells(subset(x = neurons, orig.ident == c("RInfa3"))))
i3a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RInfa3"))))
i3r <- round((i3a / i3)*100, digits = 1)

s1 <- length(Cells(subset(x = neurons, orig.ident == c("RSire1"))))
s1a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RSire1"))))
s1r <- round((s1a / s1)*100, digits = 1)
s2 <- length(Cells(subset(x = neurons, orig.ident == c("RSire2"))))
s2a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RSire2"))))
s2r <- round((s2a / s2)*100, digits = 1)
s3 <- length(Cells(subset(x = neurons, orig.ident == c("RSire3"))))
s3a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RSire3"))))
s3r <- round((s3a / s3)*100, digits = 1)
s4 <- length(Cells(subset(x = neurons, orig.ident == c("RSire4"))))
s4a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RSire4"))))
s4r <- round((s4a / s4)*100, digits = 1)

d1 <- length(Cells(subset(x = neurons, orig.ident == c("RDam1"))))
d1a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RDam1"))))
d1r <- round((d1a / d1)*100, digits = 1)
d2 <- length(Cells(subset(x = neurons, orig.ident == c("RDam2"))))
d2a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RDam2"))))
d2r <- round((d2a / d2)*100, digits = 1)
d3 <- length(Cells(subset(x = neurons, orig.ident == c("RDam3"))))
d3a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RDam3"))))
d3r <- round((d3a / d3)*100, digits = 1)
d4 <- length(Cells(subset(x = neurons, orig.ident == c("RDam4"))))
d4a <- length(Cells(subset(x = neurons, subset = a > 0 & orig.ident == c("RDam4"))))
d4r <- round((d4a / d4)*100, digits = 1)

ids <- c("infa1", "infa2", "infa3", "con1", "con2", "con3", "allo1", "allo2", "allo3", "sire1", "sire2", "sire3", "sire4", "dam1", "dam2", "dam3", "dam4")
phenotype <- c("infa", "infa", "infa", "ambicon", "ambicon", "ambicon", "allo", "allo", "allo", "sire", "sire", "sire", "sire", "dam", "dam", "dam", "dam")
allneuro <- c(i1, i2, i3, c1, c2, c3, a1, a2, a3, s1, s2, s3, s4, d1, d2, d3, d4)
asipneuro <- c(i1a, i2a, i3a, c1a, c2a, c3a, a1a, a2a, a3a, s1a, s2a, s3a, s4a, d1a, d2a, d3a, d4a)
ratio <- c(i1r, i2r, i3r, c1r, c2r, c3r, a1r, a2r, a3r, s1r, s2r, s3r, s4r, d1r, d2r, d3r, d4r)
dummycode <- c("-1", "-1", "-1", "0", "0", "0", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2")

agouticounts <- data.frame(ids, phenotype, allneuro, asipneuro, ratio, dummycode)
agouticounts$dummycode <- as.numeric(agouticounts$dummycode)
summary(agouticounts$dummycode)

agouticounts$phenotype <- factor(agouticounts$phenotype, levels = c("infa", "ambicon", "allo", "sire", "dam"))

head(agouticounts)
library(car)
a <- lm(agouticounts$ratio ~ agouticounts$phenotype)
Anova(a, type = "II")
effectsize::eta_squared(a, partial = T)
TukeyHSD(aov(agouticounts$ratio ~ agouticounts$phenotype))

png(filemicrobenchmarkpng(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Percent_Expressing_ASIP.png",
   width = 200, height = 200, units = "mm", res = 600)
ggplot(agouticounts, aes(x=phenotype, y=ratio, color = phenotype)) + 
  geom_boxplot(size = 1) + geom_point(size = 6)  + theme_linedraw() + 
  scale_color_manual(values=c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) + 
  ylim(0,100) + ylab("% Neurons Expressing > 0 counts ASIP") + 
  theme(axis.text = element_text(size = 20), 
        text = element_text(size = 20, face = "bold")) + NoLegend() 
dev.off()

b <- read.csv(file = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/ASIP_Behavior_Connection.csv", header = T)
View(b)

b$Group <- factor(b$Group, levels = c("infa", "ambicon", "allo", "sire", "dam"))

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Behavior.png",
    width = 175, height = 175, units = "mm", res = 600)
ggplot(b, aes(x=Group, y=Care/3, color = Group)) + 
  geom_boxplot() + geom_point()  + theme_linedraw() + 
  scale_color_manual(values=c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) + 
  ylab("Percent of Time in Care Behavior") + ylim(0,100) + xlab("phenotype") +
  theme(axis.text = element_text(size = 10), 
        text = element_text(size = 12, face = "bold")) + NoLegend() 
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Normalized_ASIP.png",
    width = 175, height = 175, units = "mm", res = 600)
ggplot(b, aes(x=Group, y=AvgASIP, color = Group)) + 
  geom_boxplot() + geom_point()  + theme_linedraw() + 
  scale_color_manual(values=c("#960018", "#DCB50B", "#00B9B9", "#008DC5", "#BB2A7F")) + 
  ylab("Aggregate Counts ASIP ÷ # Neurons ")  + xlab("phenotype") +
  theme(axis.text = element_text(size = 10), 
        text = element_text(size = 12, face = "bold")) + NoLegend() 
dev.off()

library(car)
agg.lm <- lm(b$ASIP ~ b$Group)
Anova(agg.lm, type = "III")

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/asip_behavior_correlation.png",
    width = 150, height = 200, units = "mm", res = 600)
ggplot(b, aes(x= Care/3, y = ASIP, color = Group)) + 
  geom_point(size = 6) + geom_smooth(method = "lm", color = "black", linewidth = 2, se = T) + 
  xlab("") + ylab("") + 
  scale_color_manual(values=c("#00B9B9", "#DCB50B","#BB2A7F", "#960018", "#008DC5")) + 
  theme(axis.text = element_text(size = 20), legend.position = "none",
        text = element_text(size = 30, face = "bold"))
dev.off()
cor.test(x = (b$Care), y = b$ASIP, method = "pearson")

c <- read.csv(file = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/Exemplars_qPCR_Behavior_Correlation.csv", header = T)
head(c)

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure3/asip_behavior_correlation_confirmation.png",
    width = 200, height = 200, units = "mm", res = 600)
ggplot(c, aes(x= pct_contact, y = NormalizedASIP, color = BroadPheno)) + 
  geom_point(size = 6) + geom_smooth(method = "lm", color = "black", linewidth = 3, se = F) +
  xlab("Percent of Time in Care Behavior") + ylab("Pseudocounts ASIP") + 
  scale_color_manual(values=c("#DCB50B", "#960018", "#00B9B9")) + 
  theme(axis.text = element_text(size = 20), legend.position = "bottom",
        text = element_text(size = 20, face = "bold"))
dev.off()

ggplot(c, aes(x= pct_contact, y = LEPTIN, color = BroadPheno)) + 
  geom_point(size = 6) + geom_smooth(method = "lm", color = "black", linewidth = 3, se = F) +
  xlab("Percent of Time in Care Behavior") + ylab("Pseudocounts ASIP") + 
  scale_color_manual(values=c("#DCB50B", "#960018", "#00B9B9")) + 
  theme(axis.text = element_text(size = 20), legend.position = "bottom",
        text = element_text(size = 20, face = "bold"))

ggplot(c, aes(x= BroadPheno, y = LEPTIN, color = BroadPheno)) + 
  geom_point(size = 6) + geom_violin() +
  xlab("Percent of Time in Care Behavior") + ylab("LEPTIN") + 
  scale_color_manual(values=c("#DCB50B", "#960018", "#00B9B9")) + 
  theme(axis.text = element_text(size = 20), legend.position = "bottom",
        text = element_text(size = 20, face = "bold"))

cor.test(x = c$pct_contact, y = c$NormalizedASIP, method = "pearson")

cor.test(x = c$LEPTIN, y = c$pct_contact, method = "pearson")


  