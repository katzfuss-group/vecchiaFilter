setwd("~/vecchiaFilter")
source("data-application/settings.r")
source("data-application/particle-filtering.r")
source("aux-functions.r")
source("data-application/data/process-data.r")
source("data-application/plot-results.r")
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidyverse))
registerDoParallel(cores = NCORES)


init_covparms = c(SIG_02, RANGE, SMOOTH)
lik.params = list(data.model = DATA_MODEL, alpha = ALPHA)


X_UPPER_LIMIT = -21.76952
LOC_QUANTILE = 1.0

TPW = readr::read_csv("TPW_10k.csv")
#TPW = readr::read_csv("data-application/data/TPW.csv") %>%
#    dplyr::filter( x < X_UPPER_LIMIT ) %>%
#    dplyr::filter( x < quantile(x, LOC_QUANTILE), y < quantile(y, LOC_QUANTILE) ) %>%
#    sample_n(size = 10000)
#readr::write_csv(TPW, sprintf("TPW_%fk.csv", nrow(TPW)/1000))


locs = TPW %>% select(x, y) %>% data.matrix()
nx = length(unique(locs[,1]))
ny = length(unique(locs[,2]))

listOfDataColumns = split.default(TPW, colnames(TPW))[as.character(1:TMAX)]

preProcessColumn = function(column){
    column = column %>% pull()
    column[column==0] = 1
    column = as.numeric(remove.mean(column))
}


Y = lapply(listOfDataColumns, preProcessColumn)

       

## filter ---------------------------------------
mra = GPvecchia::vecchia_specify(locs, COND_SET_SIZE, conditioning = 'mra')
predsMRA = filter(mra, Y, N_PARTS, lik.params, init_covparms, saveUQ = "L")


YY = list(x = Y, y = Y)
#true_params = list(c = 0.8)
plotFields(YY, predsMRA$preds, locs)
#plotParamPaths(predsMRA$particles, predsMRA$resampled.indices, DATA_MODEL, true_params)
#plot1dLogLikelihood(predsMRA$particles, predsMRA$logliks, true_params)
#plotMarginalDistribution(predsMRA$particles, predsMRA$resampled.indices, true_params)
par(mfcol = c(4, 5))
for (t in 1:TMAX) {
    logliks = predsMRA$logliks[[t]]
    logliks = logliks - max(logliks)
    plot(predsMRA$particles[[t]][, "c"], exp(logliks), pch = 16, main = t, ylab = "likelihood", type = "l", xaxt = "n", xlab = "")
    plot(predsMRA$particles[[t]][, "c"], logliks, pch = 16, xlab = "c", ylab = "log-likelihood", type = "l")
}
