#!/bin/env Rscript

set.seed(1234)
library("glmmTMB")
library("marginaleffects")
library("stringr")

inpath <- "GLMMs/Fos-GLMMs/"
flavor <- read.table(paste(inpath, "flavor.csv", sep = ""), sep = ",", col.names = "flavor")
timepoint <- read.table(paste(inpath, "timepoint.csv", sep = ""), sep = ",", col.names = "timepoint")
batch <- read.table(paste(inpath, "batch.csv", sep = ""), sep = ",", col.names = "batch")
sex <- read.table(paste(inpath, "sex.csv", sep = ""), sep = ",", col.names = "sex")
totalcounts <- read.table(paste(inpath, "totalcounts.csv", sep = ""), sep = ",", col.names = "totalcounts")
pbcounts <- read.table(paste(inpath, "pbcounts.csv", sep = ""), sep = ",", col.names = "pbcounts")
weights <- read.table(paste(inpath, "weights.csv", sep = ""), sep = ",", col.names = "weights")

for (i in 1:221) {

    cat(paste("Modeling brain region #", i, sep = ""), "\n")
    counts <- read.table(paste(inpath, "counts_region_", str_pad(as.character(i), 4, pad = "2"), ".csv", sep = ""), sep = ",", col.names = "counts")
    DF <- data.frame(counts = counts, flavor = flavor, timepoint = timepoint, batch = batch, sex = sex, totalcounts = totalcounts, pbcounts = pbcounts, weights = weights)
    DF$weights <- DF$weights/sum(DF$weights)*length(DF$weights)

    ##### Equation 2 #####
    DF_main<- subset(DF, timepoint == "Consumption" | timepoint == "Malaise" | timepoint == "Retrieval")
    DF_main$weights <- DF_main$weights/sum(DF_main$weights)*length(DF_main$weights)
    glmm2 <- glmmTMB(counts ~ flavor * timepoint + sex + (1 | batch) + offset(log(totalcounts)), data = DF_main, family = "nbinom2", weights = weights)
    outpath2a <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "2"), "_Eq2_coefficients.csv", sep = "")
    file.remove(outpath2a)
    capture.output(summary(glmm2)$coefficients$cond, file = outpath2a, append = TRUE)

    glmm2x <- update(glmm2, . ~ . - flavor * timepoint)
    anv2 <- anova(glmm2x, glmm2)
    outpath2b <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "2"), "_Eq2_modelstats.csv", sep = "")
    file.remove(outpath2b)
    capture.output(anv2$"Chisq", file = outpath2b, append = TRUE)
    capture.output(anv2$"Pr(>Chisq)", file = outpath2b, append = TRUE)

    mfx2 <- comparisons(glmm2, variables = "flavor", by = "timepoint", type = "response", newdata = DF_main, transform_pre = "difference", wts = DF_main$weights)
    outpath2c <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "2"), "_Eq2_flavor.csv", sep = "")
    file.remove(outpath2c)
    capture.output(mfx2$"estimate"[1:4], file = outpath2c, append = TRUE)
    capture.output(mfx2$"std.error"[1:4], file = outpath2c, append = TRUE)
    capture.output(mfx2$"p.value"[1:4], file = outpath2c, append = TRUE)

    ##### Equation 3 #####
    glmm3 <- glmmTMB(counts ~ timepoint + sex + (1 | batch) + offset(log(pbcounts)), data = DF, family = "nbinom2", weights = weights)
    outpath3 <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "2"), "_Eq3_coefficients.csv", sep = "")
    file.remove(outpath3)
    capture.output(summary(glmm3)$coefficients$cond, file = outpath3, append = TRUE)

    ##### Equation 4 #####
    DF_malaise <- subset(DF, timepoint == "Malaise" | timepoint == "PBNCGRP")
    DF_malaise$weights <- DF_malaise$weights/sum(DF_malaise$weights)*length(DF_malaise$weights)
    glmm4 <- glmmTMB(counts ~ flavor * timepoint + sex + (1 | batch) + offset(log(pbcounts)), data = DF_malaise, family = "nbinom2", weights = weights)
    mfx4 <- comparisons(glmm4, variables = "flavor", by = "timepoint", type = "response", newdata = DF_malaise, transform_pre = "difference", wts = DF_malaise$weights)
    outpath4 <- paste(inpath, "GLMM_output_region_", str_pad(as.character(i), 4, pad = "2"), "_Eq4_flavor.csv", sep = "")
    file.remove(outpath4)
    capture.output(mfx4$"estimate"[1:4], file = outpath4, append = TRUE)
    capture.output(mfx4$"std.error"[1:4], file = outpath4, append = TRUE)
    capture.output(mfx4$"p.value"[1:4], file = outpath4, append = TRUE)

}

##### Equation 5 #####
DF <- read.table(paste(inpath, "LS-hM3D.csv", sep = ""), sep = ",", header=TRUE)
glmm5_LS <- glmmTMB(LS.Fos ~ Timepoint + Sex + (1 | Batch) + offset(log(Total.Fos)), data = DF, family = "nbinom2")
summary(glmm5_LS)
glmm5_CEA <- glmmTMB(CEA.Fos ~ Timepoint + Sex + (1 | Batch) + offset(log(Total.Fos)), data = DF, family = "nbinom2")
summary(glmm5_CEA)