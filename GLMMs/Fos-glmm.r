#!/bin/env Rscript

set.seed(1234)
library("glmmTMB")
library("marginaleffects")
library("stringr")

inpath <- "GLMMs/Fos-GLMMs/"
flavor <- read.table(paste(inpath, "flavor.csv", sep = ""), sep = ",", col.names = "flavor")
phase <- read.table(paste(inpath, "phase.csv", sep = ""), sep = ",", col.names = "phase")
batch <- read.table(paste(inpath, "batch.csv", sep = ""), sep = ",", col.names = "batch")
sex <- read.table(paste(inpath, "sex.csv", sep = ""), sep = ",", col.names = "sex")
totalcounts <- read.table(paste(inpath, "offset.csv", sep = ""), sep = ",", col.names = "totalcounts")
pbncounts <- read.table(paste(inpath, "pbn.csv", sep = ""), sep = ",", col.names = "pbncounts")
weights <- read.table(paste(inpath, "weights.csv", sep = ""), sep = ",", col.names = "weights")

for (i in 1:201) {

    cat(paste("Modeling brain region #", i, sep = ""), "\n")
    counts <- read.table(paste(inpath, "counts_region_", str_pad(as.character(i), 4, pad = "0"), ".csv", sep = ""), sep = ",", col.names = "counts")
    DF <- data.frame(counts = counts, flavor = flavor, phase = phase, batch = batch, sex = sex, totalcounts = totalcounts, pbncounts = pbncounts, weights = weights)
    DF$weights <- DF$weights/sum(DF$weights)*length(DF$weights)
    DF_main<- subset(DF, phase == "Consumption" | phase == "Malaise" | phase == "Retrieval")
    DF_main$weights <- DF_main$weights/sum(DF_main$weights)*length(DF_main$weights)

    ##### Equation 2 #####
    glmm <- glmmTMB(counts ~ flavor * phase + sex + (1 | batch) + offset(log(totalcounts)), data = DF_main, family = "nbinom2", weights = weights)
    outpath1 <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "0"), "_Eq2_coefficients.csv", sep = "")
    file.remove(outpath1)
    capture.output(summary(glmm)$coefficients$cond, file = outpath1, append = TRUE)

    glmmx0 <- update(glmm, . ~ . - flavor * phase)
    anv0 <- anova(glmmx0, glmm)
    glmmx1 <- update(glmm, . ~ . - flavor - flavor:phase)
    anv1 <- anova(glmmx1, glmm)
    glmmx2 <- update(glmm, . ~ . - phase - flavor:phase)
    anv2 <- anova(glmmx2, glmm)
    glmmx3 <- update(glmm, . ~ . - flavor:phase)
    anv3 <- anova(glmmx3, glmm)
    glmmx4 <- update(glmm, . ~ . - sex)
    anv4 <- anova(glmmx4, glmm)
    outpath2 <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "0"), "_Eq2_reducedmodels.csv", sep = "")
    file.remove(outpath2)
    capture.output(anv0$BIC, file = outpath2, append = TRUE)
    capture.output(anv0$logLik, file = outpath2, append = TRUE)
    capture.output(anv0$"Pr(>Chisq)", file = outpath2, append = TRUE)
    capture.output(anv1$"Pr(>Chisq)", file = outpath2, append = TRUE)
    capture.output(anv2$"Pr(>Chisq)", file = outpath2, append = TRUE)
    capture.output(anv3$"Pr(>Chisq)", file = outpath2, append = TRUE)
    capture.output(anv4$"Pr(>Chisq)", file = outpath2, append = TRUE)
    capture.output(anv0$"Chisq", file = outpath2, append = TRUE)
    capture.output(anv1$"Chisq", file = outpath2, append = TRUE)
    capture.output(anv2$"Chisq", file = outpath2, append = TRUE)
    capture.output(anv3$"Chisq", file = outpath2, append = TRUE)
    capture.output(anv4$"Chisq", file = outpath2, append = TRUE)

    mfx1 <- comparisons(glmm, variables = "flavor", by = "phase", type = "response", newdata = DF_main, transform_pre = "difference", wts = DF_main$weights)
    mfx2 <- predictions(glmm, by = "phase", type = "response", newdata = DF_main, wts = DF_main$weights)
    mfx3 <- predictions(glmm, by = "phase", type = "response", newdata = DF_main, wts = DF_main$weights, hypothesis = "pairwise")
    outpath3 <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "0"), "_Eq2_marginaleffects.csv", sep = "")
    file.remove(outpath3)
    capture.output(mfx1$"estimate"[1:4], file = outpath3, append = TRUE)
    capture.output(mfx1$"std.error"[1:4], file = outpath3, append = TRUE)
    capture.output(mfx1$"p.value"[1:4], file = outpath3, append = TRUE)
    capture.output(mfx2$"estimate"[1:4], file = outpath3, append = TRUE)
    capture.output(mfx2$"std.error"[1:4], file = outpath3, append = TRUE)
    capture.output(mfx3$"p.value"[1:4], file = outpath3, append = TRUE)

    ##### Equation 3 #####
    DF_malaise <- DF
    DF_malaise$weights <- DF_malaise$weights/sum(DF_malaise$weights)*length(DF_malaise$weights)
    glmm <- glmmTMB(counts ~ phase + sex + (1 | batch) + offset(log(pbncounts)), data = DF_malaise, family = "nbinom2", weights = weights)
    mfx5 <- predictions(glmm, by = "phase", type = "response", newdata = DF_malaise, wts = DF_malaise$weights)
    mfx6 <- predictions(glmm, by = "phase", type = "response", newdata = DF_malaise, wts = DF_malaise$weights, hypothesis = "pairwise")
    outpath4 <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "0"), "_Eq3_coefficients.csv", sep = "")
    file.remove(outpath4)
    capture.output(summary(glmm)$coefficients$cond, file = outpath4, append = TRUE)

    ##### Equation 4 #####
    DF_malaise <- subset(DF, phase == "Malaise" | phase == "PBNCGRP")
    DF_malaise$weights <- DF_malaise$weights/sum(DF_malaise$weights)*length(DF_malaise$weights)
    glmm <- glmmTMB(counts ~ flavor * phase + sex + (1 | batch) + offset(log(pbncounts)), data = DF_malaise, family = "nbinom2", weights = weights)
    mfx4 <- comparisons(glmm, variables = "flavor", by = "phase", type = "response", newdata = DF_malaise, transform_pre = "difference", wts = DF_malaise$weights)
    outpath5 <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "0"), "_Eq4_marginaleffects.csv", sep = "")
    file.remove(outpath5)
    capture.output(mfx4$"estimate"[1:4], file = outpath5, append = TRUE)
    capture.output(mfx4$"std.error"[1:4], file = outpath5, append = TRUE)
    capture.output(mfx4$"p.value"[1:4], file = outpath5, append = TRUE)
    capture.output(mfx5$"estimate"[1:4], file = outpath5, append = TRUE)
    capture.output(mfx5$"std.error"[1:4], file = outpath5, append = TRUE)
    capture.output(mfx6$"p.value"[1:6], file = outpath5, append = TRUE)

}

##### Equation 5 #####
DF <- read.table(paste(inpath, "LS-hM3D.csv", sep = ""), sep = ",", header=TRUE)
glmm_LS <- glmmTMB(LS.Fos ~ Phase + Sex + (1 | Batch.code) + offset(log(Total.Fos)), data = DF, family = "nbinom2")
summary(glmm_LS)
glmm_CEA <- glmmTMB(CEA.Fos ~ Phase + Sex + (1 | Batch.code) + offset(log(Total.Fos)), data = DF, family = "nbinom2")
summary(glmm_CEA)