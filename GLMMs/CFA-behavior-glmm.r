#!/bin/env Rscript

set.seed(1234)
library("glmmTMB")
library("marginaleffects")
library("stringr")

DF <- read.table(paste(inpath, "GLMMs/CFA-behavior-GLMM-data.csv", sep = ""), sep = ",", header=TRUE)

##### Equation 1 #####
glmm <- glmmTMB(Preference ~ Novel*LiCl*Day + Sex + (1|Subject), data = DF, family = "gaussian")

DF_Novel = subset(DF, Novel==1)
mfx_Novel <- comparisons(glmm, variables = "Injection", by = "Day", type = "response", newdata = DF_Novel, transform_pre = "difference")
summary(mfx_Novel)

DF_Familiar = subset(DF, Novel==0)
mfx_Familiar <- comparisons(glmm, variables = "Injection", by = "Day", type = "response", newdata = DF_Familiar, transform_pre = "difference")
summary(mfx_Familiar)