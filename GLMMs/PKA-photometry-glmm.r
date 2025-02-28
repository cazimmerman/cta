#!/bin/env Rscript

set.seed(1234)
library("glmmTMB")
library("marginaleffects")
library("stringr")

DF <- read.table(paste(inpath, "GLMMs/PKA-photometry-GLMM-data.csv", sep = ""), sep = ",", header=TRUE)

##### Equation 6 #####
glmm <- glmmTMB(Response ~ Port*Day + (1|Subject), data = DF, family = "gaussian")

mfx <- comparisons(glmm, variables = "Port", by = "Day", type = "response", newdata = DF, transform_pre = "difference")
summary(mfx)