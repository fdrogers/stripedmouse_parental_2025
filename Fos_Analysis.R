### This R script was created for the analysis of normalized FOS+ neuron densities

## Load Packages
library(car)
library(ggplot2)
library(lmtest)
library(lme4)
library(afex)
library(sjPlot)
library(effectsize)
library(corrplot)
library(patchwork)
library(agricolae)
library(ggthemes)

## Import Data

d <- read.csv(file = "Dropbox/Princeton/African striped mouse/cFOS_count_data/Consolidated_FOS.csv", header = T)
View(d)

e <- read.csv(file = "Dropbox/Princeton/African striped mouse/cFOS_count_data/Consolidated_FOS_notnormalized.csv", header = T)
View(e)

# Check for and remove outliers

e$MPOA[13] <- NA

Boxplot(e$AON, id.method="y") #No Outliers
Boxplot(e$CP, id.method="y")
e$CP[13] <- NA
e$CP[17] <- NA
e$CP[20] <- NA
e$CP[36] <- NA

Boxplot(e$NAcc, id.method="y")
e$NAcc[2] <- NA
e$NAcc[5] <- NA
e$NAcc[13] <- NA
e$NAcc[20] <- NA
e$NAcc[30] <- NA

Boxplot(e$NAcsh, id.method="y")
e$NAcsh[7] <- NA
e$NAcsh[13] <- NA
e$NAcsh[20] <- NA
e$NAcsh[23] <- NA
e$NAcsh[30] <- NA

Boxplot(e$vBST, id.method="y") #No Outliers
Boxplot(e$MPOA, id.method="y") #1 outlier
e$MPOA[13] <- NA

Boxplot(e$PVT, id.method="y")
e$PVT[36] <- NA

Boxplot(e$PFC, id.method="y") #No Outliers
Boxplot(e$VP, id.method="y")
e$VP[23] <- NA
e$VP[36] <- NA

Boxplot(e$PAG, id.method="y") #No Outliers
Boxplot(e$VTA ~ e$Cohort, id.method="y")
e$VTA[23] <- NA
e$VTA[21] <- NA
e$VTA[36] <- NA
e$VTA[30] <- NA

Boxplot(e$BLA, id.method="y")
e$BLA[22] <- NA
e$BLA[31] <- NA

Boxplot(e$LS, id.method="y")
e$LS[13] <- NA
e$LS[20] <- NA
e$LS[27] <- NA
e$LS[36] <- NA

Boxplot(e$LSd, id.method="y")
e$LSd[13] <- NA
e$LSd[20] <- NA
e$LSd[27] <- NA
e$LSd[36] <- NA

### By Housing

aon.plot.housing <- ggplot(data = e, aes(x = Housing, y = AON, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

cp.plot.housing <- ggplot(data = e, aes(x = Housing, y = CP, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

pfc.plot.housing <- ggplot(data = e, aes(x = Housing, y = PFC, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

nacc.plot.housing <- ggplot(data = e, aes(x = Housing, y = NAcc, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

nacsh.plot.housing <- ggplot(data = e, aes(x = Housing, y = NAcsh, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

vbst.plot.housing <- ggplot(data = e, aes(x = Housing, y = vBST, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

mpoa.plot.housing <- ggplot(data = e, aes(x = Housing, y = MPOA, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

pvt.plot.housing <- ggplot(data = e, aes(x = Housing, y = PVT, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

vp.plot.housing <- ggplot(data = e, aes(x = Housing, y = VP, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

pag.plot.housing <- ggplot(data = e, aes(x = Housing, y = PAG, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

vta.plot.housing <- ggplot(data = e, aes(x = Housing, y = VTA, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

bla.plot.housing <- ggplot(data = e, aes(x = Housing, y = BLA, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

ls.plot.housing <- ggplot(data = e, aes(x = Housing, y = LS, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

pvt.plot.housing <- ggplot(data = e, aes(x = Housing, y = PVT, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

lsd.plot.housing <- ggplot(data = e, aes(x = Housing, y = LSd, fill = Housing, color = Housing)) + 
  facet_grid(. ~ Cohort, ) +
  stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("white", "grey50")) + 
  scale_color_manual(values = c("black", "grey")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())


png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2a.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
aon.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2b.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
bla.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2c.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
cp.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2d.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
ls.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2e.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
lsd.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2f.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
mpoa.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2g.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
nacc.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2h.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
nacsh.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2i.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
pag.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2j.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
pfc.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2k.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
pvt.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2l.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
vbst.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2m.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
vp.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2n.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
vta.plot.housing
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure2.png", 
    width = 240, height = 500, units = "mm", res = 600, bg = "white")
(aon.plot.housing | bla.plot.housing) / (cp.plot.housing | ls.plot.housing) / 
  (lsd.plot.housing | mpoa.plot.housing) / (nacc.plot.housing | nacsh.plot.housing) / 
  (pag.plot.housing | pfc.plot.housing ) / (pvt.plot.housing | vbst.plot.housing ) /  (vp.plot.housing | vta.plot.housing)
dev.off()

aon.housing <- lm(e$AON ~ e$Housing*e$Cohort)
Anova(aon.housing, type = "III") #F(1, 34) = 6.5318, p = .01524 Pup Effect, No Housing
eta_squared(aon.housing, partial = TRUE) #eta-sq = 0.34 for Cohort
marginal.aon.housing = emmeans(aon.housing, ~ Cohort)
pairs(marginal.aon.housing, adjust="fdr")
#contrast               estimate   SE df t.ratio p.value
#control - dysparental    -282.2 73.1 35  -3.859  0.0014
#control - parental       -265.0 86.5 35  -3.062  0.0063
#dysparental - parental     17.2 88.6 35   0.195  0.8469

bla.housing <- lm(e$BLA ~ e$Housing*e$Cohort)
Anova(bla.housing, type = "III") # ns
eta_squared(bla.housing, partial = TRUE) #eta-sq = 0.34 for Cohort

cp.housing <- lm(e$CP ~ e$Housing*e$Cohort)
Anova(cp.housing, type = "III") #ns 
eta_squared(cp.housing, partial = TRUE) 

ls.housing <- lm(e$LS ~ e$Housing*e$Cohort)
Anova(ls.housing, type = "III") #ns
eta_squared(ls.housing, partial = TRUE) 

lsd.housing <- lm(e$LSd ~ e$Housing*e$Cohort)
Anova(lsd.housing, type = "III") #ns 
eta_squared(lsd.housing, partial = TRUE)

mpoa.housing <- lm(e$MPOA ~ e$Housing*e$Cohort)
Anova(mpoa.housing, type = "III") # For cohort F(1,31) = 4.3280, p = .045845
eta_squared(mpoa.housing, partial = TRUE) #eta-sq = 0.29 for cohort
marginal.mpoa.housing = emmeans(mpoa.housing, ~ Cohort)
pairs(marginal.mpoa.housing, adjust="fdr") # fdr corrected cohort 2 > cohort 3 p = .0011

nacc.housing <- lm(e$NAcc ~ e$Housing*e$Cohort)
Anova(nacc.housing, type = "III") #ns 
eta_squared(nacc.housing, partial = TRUE) 

nacsh.housing <- lm(e$NAcsh ~ e$Housing*e$Cohort)
Anova(nacsh.housing, type = "III") #ns 
eta_squared(nacsh.housing, partial = TRUE) 

pag.housing <- lm(e$PAG ~ e$Housing*e$Cohort)
Anova(pag.housing, type = "III") #ns 
eta_squared(pag.housing, partial = TRUE)

pfc.housing <- lm(e$PFC ~ e$Housing*e$Cohort)
Anova(pfc.housing, type = "III") #F(1,37) = 12.9154, p = .0009445 for cohort
eta_squared(pfc.housing, partial = TRUE) #eta-sq = 0.49 for cohort
marginal.pfc.housing = emmeans(pfc.housing, ~ Cohort)
pairs(marginal.pfc.housing, adjust="fdr") # fdr corrected cohort 2 > cohort 3 p < .0001

pvt.housing <- lm(e$PVT ~ e$Housing*e$Cohort)
Anova(pvt.housing, type = "III") #ns 
eta_squared(pvt.housing, partial = TRUE) 

vbst.housing <- lm(e$vBST ~ e$Housing*e$Cohort)
Anova(vbst.housing, type = "III") #ns 
eta_squared(vbst.housing, partial = TRUE) 

vp.housing <- lm(e$VP ~ e$Housing*e$Cohort)
Anova(vp.housing, type = "III") #ns 
eta_squared(vp.housing, partial = TRUE) 

vta.housing <- lm(e$VTA ~ e$Housing*e$Cohort)
Anova(vta.housing, type = "III") #ns
eta_squared(vta.housing, partial = TRUE)


### By MetaGroup

aon.plot <- ggplot(data = e, aes(x = MetaGroup, y = AON, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

cp.plot <- ggplot(data = e, aes(x = MetaGroup, y = CP, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

pfc.plot <- ggplot(data = e, aes(x = MetaGroup, y = PFC, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

nacc.plot <- ggplot(data = e, aes(x = MetaGroup, y = NAcc, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

nacsh.plot <- ggplot(data = e, aes(x = MetaGroup, y = NAcsh, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

vbst.plot <- ggplot(data = e, aes(x = MetaGroup, y = vBST, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

mpoa.plot <- ggplot(data = e, aes(x = MetaGroup, y = MPOA, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())
##For main figure 2
plot0 <- ggplot(data = e, aes(x = MetaGroup, y = MPOA, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 10) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous() + 
  scale_fill_manual(values = c("grey", "orange", "#9DD0D0")) + 
  scale_color_manual(values = c("grey20", "darkorange3", "#00B9B9")) + theme_few() + xlab("Phenotype") +
  ylab("cFOS+ Neurons per mm^2") + 
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank(), 
        title = element_text(face = "bold", family = "arial", size = 20)) + ylim(-2,400)

pvt.plot <- ggplot(data = e, aes(x = MetaGroup, y = PVT, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

vp.plot <- ggplot(data = e, aes(x = MetaGroup, y = VP, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

pag.plot <- ggplot(data = e, aes(x = MetaGroup, y = PAG, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

vta.plot <- ggplot(data = e, aes(x = MetaGroup, y = VTA, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

bla.plot <- ggplot(data = e, aes(x = MetaGroup, y = BLA, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

ls.plot <- ggplot(data = e, aes(x = MetaGroup, y = LS, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

pvt.plot <- ggplot(data = e, aes(x = MetaGroup, y = PVT, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

lsd.plot <- ggplot(data = e, aes(x = MetaGroup, y = LSd, fill = MetaGroup, color = MetaGroup)) + stat_summary(fun = mean, geom = "col") +
  geom_point(size = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75) + scale_y_continuous(limits = c(0,1000)) + 
  scale_fill_manual(values = c("cyan", "magenta", "green")) + 
  scale_color_manual(values = c("darkcyan", "darkmagenta", "darkgreen")) + theme_minimal() + xlab("") +
  ylab("cFOS+ Neurons per mm2") + ggtitle("") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_text(family = "arial", size = 16), strip.text = element_blank())

(aon.plot | bla.plot | cp.plot | ls.plot | lsd.plot | mpoa.plot | nacc.plot ) / 
  (nacsh.plot | pag.plot | pfc.plot | pvt.plot | vbst.plot | vp.plot | vta.plot)

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3raw.png", 
    width = 240, height = 500, units = "mm", res = 600, bg = "white")
(aon.plot | bla.plot | cp.plot) / (ls.plot | lsd.plot | mpoa.plot) / (nacc.plot | nacsh.plot | pag.plot | pfc.plot ) / (pvt.plot | vbst.plot | vp.plot | vta.plot )
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3a.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
aon.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3b.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
bla.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3c.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
cp.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3d.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
ls.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3e.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
lsd.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3f.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
mpoa.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3g.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
nacc.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3h.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
nacsh.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3i.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
pag.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3j.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
pfc.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3k.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
pvt.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3l.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
vbst.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3m.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
vp.plot
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure3n.png", 
    width = 90, height = 90, units = "mm", res = 600, bg = "white")
vta.plot
dev.off()

aon <- lm(e$AON ~ e$MetaGroup)
Anova(aon, type = "III") #F(2, 35) = 8.849, p = .0007759
eta_squared(aon, partial = TRUE) #eta-sq = 0.34
marginal.aon = emmeans(aon, ~ MetaGroup)
pairs(marginal.aon, adjust="fdr")
#contrast               estimate   SE df t.ratio p.value
#control - dysparental    -282.2 73.1 35  -3.859  0.0014
#control - parental       -265.0 86.5 35  -3.062  0.0063
#dysparental - parental     17.2 88.6 35   0.195  0.8469

bla <- lm(e$BLA ~ e$MetaGroup)
Anova(bla, type = "III") # F(2, 25) = 5.1391, p = .0135
eta_squared(bla, partial = TRUE) #eta-sq = 0.29
marginal.bla = emmeans(bla, ~ MetaGroup)
pairs(marginal.bla, adjust="fdr")
#contrast               estimate   SE df t.ratio p.value
#control - dysparental     -81.1 39.3 25  -2.061  0.0748
#control - parental       -137.2 43.1 25  -3.183  0.0116
#dysparental - parental    -56.1 39.3 25  -1.427  0.1661

cp <- lm(e$CP ~ e$MetaGroup)
Anova(cp, type = "III") #ns F(2,34) = 0.2294, p = .7962
eta_squared(cp, partial = TRUE) #eta-sq = 0.01

ls <- lm(e$LS ~ e$MetaGroup)
Anova(ls, type = "III") #ns F(2,35) = 0.1009, p = .9043
eta_squared(ls, partial = TRUE) #eta-sq = 5.73e-03

lsd <- lm(e$LSd ~ e$MetaGroup)
Anova(lsd, type = "III") #ns F(2,35) = 0.2676, p = .7668
eta_squared(lsd, partial = TRUE) #eta-sq = 0.02

mpoa <- lm(e$MPOA ~ e$MetaGroup)
Anova(mpoa, type = "III") # F(2,32) = 11.1892, p = .0002068
eta_squared(mpoa, partial = TRUE) #eta-sq = 0.41
marginal.mpoa = emmeans(mpoa, ~ MetaGroup)
pairs(marginal.mpoa, adjust="fdr")
#contrast               estimate   SE df t.ratio p.value
#control - dysparental     -75.4 32.0 32  -2.353  0.0249
#control - parental       -173.7 36.9 32  -4.712  0.0001
#dysparental - parental    -98.3 37.4 32  -2.630  0.0195

nacc <- lm(e$NAcc ~ e$MetaGroup)
Anova(nacc, type = "III") #ns F(2,29) = 0.9963, p = .38153
eta_squared(nacc, partial = TRUE) #eta-sq = 0.06

nacsh <- lm(e$NAcsh ~ e$MetaGroup)
Anova(nacsh, type = "III") #ns F(2,30) = 1.0186, p = .373241
eta_squared(nacsh, partial = TRUE) #eta-sq = 0.06

pag <- lm(e$PAG ~ e$MetaGroup)
Anova(pag, type = "III") #ns F(2,31) = 1.3128, p = .2836
eta_squared(pag, partial = TRUE) #eta-sq = 0.08

pfc <- lm(e$PFC ~ e$MetaGroup)
Anova(pfc, type = "III") #F(2,38) = 18.378, p = 2.609e-06
eta_squared(pfc, partial = TRUE) #eta-sq = 0.49
marginal.pfc = emmeans(pfc, ~ MetaGroup)
pairs(marginal.pfc, adjust="fdr")
#contrast               estimate   SE df t.ratio p.value
#control - dysparental    -384.5 67.4 38  -5.705  <.0001
#control - parental       -317.9 80.7 38  -3.942  0.0005
#dysparental - parental     66.6 84.8 38   0.785  0.4372

pvt <- lm(e$PVT ~ e$MetaGroup)
Anova(pvt, type = "III") #ns F(2,28) = 2.3109, p = .1178
eta_squared(pvt, partial = TRUE) #eta-sq = 0.14

vbst <- lm(e$vBST ~ e$MetaGroup)
Anova(vbst, type = "III") #ns F(2,39) = 0.234, p = 0.7925
eta_squared(vbst, partial = TRUE) #eta-sq = 0.01

vp <- lm(e$VP ~ e$MetaGroup)
Anova(vp, type = "III") #ns F(2,36) = 0.1507, p = .8606
eta_squared(vp, partial = TRUE) #eta-sq = 8.30e-03

vta <- lm(e$VTA ~ e$MetaGroup)
Anova(vta, type = "III") #ns F(2,30) = 5.5943, p = 0.008613
eta_squared(vta, partial = TRUE) #eta-sq = 0.11
marginal.vta = emmeans(vta, ~ MetaGroup)
pairs(marginal.vta, adjust="fdr")
#contrast               estimate   SE df t.ratio p.value
#control - dysparental     -68.2 26.0 30  -2.621  0.0204
#control - parental        -92.8 31.1 30  -2.984  0.0168
#dysparental - parental    -24.6 31.1 30  -0.791  0.4353

##### Correlation based analyses

#Check for normality
mpoa.lm <- lm(d$MPOA ~ d$Total_Contact, data=d, na.action = na.exclude) 
mpoa.res <- resid(mpoa.lm, na.rm = TRUE) 
shapiro.test(mpoa.res) #Data is not normally distributed. Consider using spearman instaed of pearson

#Run correlation tests
cor.test(d$Total_Contact, d$AON, method = "spearman", exact = F) #p = .7913, rho = -0.0599
cor.test(d$Total_Contact, d$CP, method = "spearman", exact = F) #p = .6971, rho = .0903
cor.test(d$Total_Contact, d$NAcc, method = "spearman", exact = F) #p = .864, rho = .0421
cor.test(d$Total_Contact, d$NAcsh, method = "spearman", exact = F) #p = .6178, rho = .1188
cor.test(d$Total_Contact, d$vBST, method = "spearman", exact = F) #p = .7568, rho = .0700
cor.test(d$Total_Contact, d$MPOA, method = "spearman", exact = F) #p = .0232, rho = .4818
cor.test(d$Total_Contact, d$PVT, method = "spearman", exact = F) #p = .4123, rho = .1941
cor.test(d$Total_Contact, d$PFC, method = "spearman", exact = F) #p = .4034, rho = .1875
cor.test(d$Total_Contact, d$VP, method = "spearman", exact = F) #p = .2135, rho = .2832
cor.test(d$Total_Contact, d$PAG, method = "spearman", exact = F) #p = .423, rho = - .2013
cor.test(d$Total_Contact, d$VTA, method = "spearman", exact = F) #p = .5612, rho = .1344
cor.test(d$Total_Contact, d$BLA, method = "spearman", exact = F) #p = .1391, rho = .3339
cor.test(d$Total_Contact, d$LS, method = "spearman", exact = F) #p = .631, rho = - .1084
cor.test(d$Total_Contact, d$LSd, method = "spearman", exact = F) #p = .8457, rho = - .0441

cor.test(d$Full_Contact, d$AON, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$CP, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$NAcc, method = "spearman", exact = F)
cor.test(d$Full_Contact, d$NAcsh, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$vBST, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$MPOA, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$PVT, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$PFC, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$VP, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$PAG, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$VTA, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$BLA, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$LS, method = "spearman", exact = F) 
cor.test(d$Full_Contact, d$LSd, method = "spearman", exact = F) 

cor.test(d$MPOA, d$AON, method = "spearman", exact= F) #p = .1133
cor.test(d$MPOA, d$PFC, method = "spearman", exact= F) #p = .007777 rho = .5517
cor.test(d$MPOA, d$BLA, method = "spearman", exact = F) #p = .3661
cor.test(d$MPOA, d$PAG, method = "spearman", exact = F) #p = .7416
cor.test(d$MPOA, d$PVT, method = "spearman", exact = F) #p = .03759 rho = 0.4676692
cor.test(d$MPOA, d$LS, method = "spearman", exact = F) #p = .06653 rho = 0.3980802
cor.test(d$MPOA, d$LSd, method = "spearman", exact = F) #p = .01975 rho = .4929418
cor.test(d$MPOA, d$CP, method = "spearman", exact = F) #p = .003297 rho = .6103896
cor.test(d$MPOA, d$VP, method = "spearman", exact = F) #p = 1.507e-06, rho = .8441558
cor.test(d$MPOA, d$NAcc, method = "spearman", exact = F) #p = .0004292 rho = 0.7263158
cor.test(d$MPOA, d$NAcsh, method = "spearman", exact = F) #p = .0002635 rho = .7293233
cor.test(d$MPOA, d$vBST, method = "spearman", exact = F) #p = .0007002 rho = .6668549
cor.test(d$MPOA, d$VTA, method = "spearman", exact = F) #p = .02341 rho = .4922078

cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$AON[which(d$Phenotype == "paternal")], method = "spearman", exact= F) 
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$PFC[which(d$Phenotype == "paternal")], method = "spearman", exact= F) 
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$BLA[which(d$Phenotype == "paternal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$PAG[which(d$Phenotype == "paternal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$PVT[which(d$Phenotype == "paternal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$LS[which(d$Phenotype == "paternal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$LSd[which(d$Phenotype == "paternal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$CP[which(d$Phenotype == "paternal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$VP[which(d$Phenotype == "paternal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$NAcc[which(d$Phenotype == "paternal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$NAcsh[which(d$Phenotype == "paternal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$vBST[which(d$Phenotype == "paternal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "paternal")], d$VTA[which(d$Phenotype == "paternal")], method = "spearman", exact = F) 

cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$AON[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact= F) 
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$PFC[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact= F) 
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$BLA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$PAG[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$PVT[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$LS[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$LSd[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$CP[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$VP[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$NAcc[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$NAcsh[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F)
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$vBST[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F) 
cor.test(d$MPOA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], d$VTA[which(d$Phenotype == "ambivalent" | d$Phenotype == "infanticidal")], method = "spearman", exact = F) 


#Make chord diagrams

library(circlize)
chord <- read.csv(file = "Dropbox/Princeton/African striped mouse/cFOS_count_data/Chord_all_cFOS.csv", header = T)

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure4a.png", 
    width = 300, height = 300, units = "mm", res = 600, bg = "white")
grid.col = c(AON = "red",
             PFC = "orange",
             BLA = "yellow",
             PAG = "green",
             PVT = "blue",
             LS = "purple",
             LSd = "pink",
             CP = "gold",
             VP = "lightyellow",
             NAcc = "lightgreen",
             NAcsh = "lightblue",
             vBST = "lavender",
             VTA = "firebrick", 
             MPOA = "black") 
circos.clear()
circos.par(start.degree = 0, clock.wise = TRUE)
par(cex = 2.5, mar = c(0, 0, 0, 0))
chordDiagram(chord, scale = FALSE, preAllocateTracks = list(track.height = 0.01), big.gap = 5, 
             directional = 2, small.gap = 5,
             grid.border = 3, grid.col = grid.col, 
             link.border = "black", symmetric = FALSE)
dev.off()

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure4b.png", 
    width = 300, height = 300, units = "mm", res = 600, bg = "white")
chord2 <- read.csv(file = "Dropbox/Princeton/African striped mouse/cFOS_count_data/Chord_paternal_cFOS.csv", header = T)

grid.col = c(AON = "red",
             PFC = "orange",
             BLA = "yellow",
             PAG = "green",
             PVT = "blue",
             LS = "purple",
             LSd = "pink",
             CP = "gold",
             VP = "lightyellow",
             NAcc = "lightgreen",
             NAcsh = "lightblue",
             vBST = "lavender",
             VTA = "firebrick", 
             MPOA = "black") 
circos.clear()
circos.par(start.degree = 0, clock.wise = TRUE)
par(cex = 2.5, mar = c(0, 0, 0, 0))
chordDiagram(chord2, scale = FALSE, preAllocateTracks = list(track.height = 0.01), big.gap = 5, directional = 2, small.gap = 5,
             grid.border = 3, grid.col = grid.col, link.border = "black", symmetric = FALSE)
dev.off()

chord3 <- read.csv(file = "Dropbox/Princeton/African striped mouse/cFOS_count_data/Chord_dyspaternal_cFOS.csv", header = T)

png(filename = "/Users/forrestrogers/Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/ExtendedDataFigure4c.png", 
    width = 300, height = 300, units = "mm", res = 600, bg = "white")
grid.col = c(AON = "red",
             PFC = "orange",
             BLA = "yellow",
             PAG = "green",
             PVT = "blue",
             LS = "purple",
             LSd = "pink",
             CP = "gold",
             VP = "lightyellow",
             NAcc = "lightgreen",
             NAcsh = "lightblue",
             vBST = "lavender",
             VTA = "firebrick", 
             MPOA = "black") 
circos.clear()
circos.par(start.degree = 0, clock.wise = TRUE)
par(cex = 2.5, mar = c(0, 0, 0, 0))
chordDiagram(chord3, scale = FALSE, preAllocateTracks = list(track.height = 0.01), big.gap = 5, directional = 2, small.gap = 5,
             grid.border = 3, grid.col = grid.col, link.border = "black", symmetric = FALSE)
dev.off()

#Make correlation matrix
mat <- read.csv(file = "Dropbox/Princeton/African striped mouse/cFOS_count_data/fos_matrix.csv", header = T)

#subset data for dysparental
e.dysp <- subset(x = e, subset = MetaGroup == "dysparental") #subset
e.dysp <- e.dysp[,8:21] #select relevant columns
corr_mat_dys=cor(e.dysp,method="spearman", use = "pairwise.complete.obs")#generate correlation matrix

corrplot(corr_mat_dys, diag = F, order = "hclust", type = "full", 
         sig.level = 0.05, insig = c("p-value"), 
         col = c("firebrick4","firebrick","pink","skyblue","blue3","navy"), 
         tl.col = "black", tl.srt = 45, tl.offset = 1,)

png(filename = "Dropbox/Princeton/African striped mouse/cFOS_count_data/fos_matrix_dysparental.png", 
    width = 225, height = 200, units = "mm", bg = "white", res = 600)
corrplot.mixed(corr_mat_dys, lower = "number", upper = "circle",
               tl.pos = "lt", tl.srt = 45, tl.offset = .5, tl.col = "black",
               order = "alphabet", bg = "black", addgrid.col = "white") 
dev.off()

#subset data for allopaternal
e.allo <- subset(x = e, subset = MetaGroup == "parental") #subset
e.allo <- e.allo[,8:21] #select relevant columns
corr_mat_allo=cor(e.allo,method="spearman", use = "pairwise.complete.obs") #generate correlation matrix
corrplot(corr_mat_allo, diag = F, order = "hclust", type = "full", 
         sig.level = 0.05, insig = c("p-value"), 
         col = c("firebrick4","firebrick","pink","skyblue","blue3","navy"), 
         tl.col = "black", tl.srt = 45, tl.offset = 1,)
png(filename = "Dropbox/Princeton/African striped mouse/cFOS_count_data/fos_matrix_allopaternal.png", 
    width = 225, height = 200, units = "mm", bg = "white", res = 600)
corrplot.mixed(corr_mat_allo, lower = "number", upper = "circle",
               tl.pos = "lt", tl.srt = 45, tl.offset = .5, tl.col = "black",
               order = "alphabet", bg = "black", addgrid.col = "white")
dev.off()

corr_mat=cor(mat,method="spearman", use = "pairwise.complete.obs")
corr_mat[1:14,1:14]
corrplot(corr_mat, diag = F, order = "hclust", type = "full", 
         sig.level = 0.05, insig = c("p-value"), 
         col = c("firebrick4","firebrick","pink","skyblue","blue3","navy"), 
         tl.col = "black", tl.srt = 45, tl.offset = 1,)


corrplot.mixed(corr_mat, lower = "number", upper = "circle",
               tl.pos = "lt", tl.srt = 45, tl.offset = .5, tl.col = "black",
               order = "hclust", bg = "black", addgrid.col = "white")

mat <- read.csv(file = "Dropbox/Princeton/African striped mouse/cFOS_count_data/fos_matrix.csv", header = T)

#Change order of phenotypes
d$Phenotype <- factor(d$Phenotype, levels = c("infanticidal", "ambivalent", "paternal"))

#Make correlations plots
plot1 <- ggplot(data = d, aes(x = (Total_Contact/1200)*100, y = MPOA, color = Phenotype)) + 
  geom_point(size = 30) + geom_smooth(method = "lm", se = F, color = "black", linewidth = 3) + 
  ylab("Normalized cFOS+ Density") + xlab("% Time in Any Contact") + 
  scale_color_manual(values = c("#960018", "#DCB50B", "#00B9B9")) + theme_few() +
  theme(legend.position = "none", legend.text = element_text(size = 14, face = "bold"), legend.title = element_blank(), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        title = element_blank())

plot2 <- ggplot(data = d, aes(x = (Full_Contact/1200)*100, y = MPOA, color = Phenotype)) + 
  geom_point(size = 10) + geom_smooth(method = "lm", se = F, color = "black") + 
  ylab("Normalized cFOS+ Density") + xlab("% Time Huddling ") + xlim(0,80) +
  scale_color_manual(values = c("#960018", "#DCB50B", "#00B9B9")) + theme_few() +
  theme(legend.position = "none", legend.text = element_text(size = 14, face = "bold"),legend.title = element_blank(),
        axis.title = element_text(face = "bold", size = 40, family = "arial"),
        axis.text = element_text(face = "bold", size = 20, family = "arial"), 
        title = element_text(face = "bold", family = "arial", size = 20))

plot1 + plot2

f <- read.csv(file = "Dropbox/Princeton/African striped mouse/cFOS_count_data/Consolidated_FOS_restandardized.csv", header = T)
View(f)

png(filename = "Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/Consolidated_FOS_restandardized2.png", 
    width = 400, height = 150, units = "mm", res = 600, bg = "white")
ggplot(data = f, aes(x = Phenotype, y = Density, fill = Phenotype, color = Phenotype, group = ROI)) + 
  geom_hline(yintercept = c(-1, -.5, .5, 1), colour = "black", linewidth = .5, linetype = 20, alpha = .5) +
  geom_hline(yintercept = c(0), colour = "black", linewidth = .75, linetype = 20, alpha = 1) +
  stat_summary(fun = mean, geom = "point", size = 5) + stat_summary(fun.data = mean_se, geom = "errorbar", width = .75, size = .75) +
  facet_wrap(facets = vars(ROI), ncol = 14, strip.position = "bottom") +
  scale_y_continuous(breaks = c(-1, -.5, 0, .5, 1)) + 
  scale_fill_manual(values = c("orange", "#9DD0D0")) + 
  scale_color_manual(values = c("darkorange3", "#00B9B9")) + theme_blank() + xlab("") +
  ylab("Standardized cFos Density") + ggtitle("") +
  theme(legend.position = "none", 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.title = element_blank(), strip.text = element_blank())
dev.off()

stat_summary(fun.data = mean_se, geom = "errorbar", width = .25, size = .75)

png(filename = "Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/MPOA_Fos_PlotB.png", 
    width = 450, height = 450, units = "mm", res = 600, bg = "white")
plot0
dev.off()

png(filename = "Dropbox/Princeton/African striped mouse/Manuscript_ASIP/Figure2/MPOA_Fos_PlotC.png", 
    width = 450, height = 450, units = "mm", res = 600, bg = "white")
plot1
dev.off()
