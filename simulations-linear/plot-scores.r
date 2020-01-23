setwd("~/HVLF")
source("plotting-functions.r")
resultsDir="simulations-linear"
data.models = c("gauss")
plotScores(resultsDir, data.models)
#plotSims(resultsDir, "gauss", 1)
#plotSims(resultsDir, "gauss", 2)
rm(list=ls())