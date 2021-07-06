setwd("/home/marcin/vecchiaFilter")
source("getMatCov.r")
source("data-application/settings.r")
#source("data-application/particle-filtering-recompute.r")
source("data-application/particle-filtering.r")
source("aux-functions.r")
source("data-application/data/process-data.r")
source("data-application/plot-results.r")
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidyverse))

registerDoParallel(cores = NCORES)


init_covparms = c(SIG_02, RANGE, SMOOTH)
lik.params = list(data.model = DATA_MODEL)



TPW = readr::read_csv("TPW_43.264000k.csv")
## TPW_full = readr::read_csv("data-application/data/TPW.csv") %>%
##     dplyr::filter( x < X_UPPER_LIMIT ) %>%
##     dplyr::filter( x < quantile(x, LOC_QUANTILE), y > quantile(y, 1-LOC_QUANTILE) ) #%>%
## TPW = TPW_full %>% sample_n(size = round(SUB_SAMPLE * nrow(TPW_full)))
## readr::write_csv(TPW, sprintf("TPW_%fk.csv", nrow(TPW)/1000))

locs = TPW %>% select(x, y) %>% data.matrix()
nx = length(unique(locs[,1]))
ny = length(unique(locs[,2]))

listOfDataColumns = split.default(TPW, colnames(TPW))[as.character(2:(TMAX+1))]

preProcessColumn = function(column){
    column = column %>% pull()
    column = column / 1000
    #column[column==0] = 1
    #column = as.numeric(remove.mean(column))
}
Y = lapply(listOfDataColumns, preProcessColumn)

       

## filter ---------------------------------------
mra = GPvecchia::vecchia_specify(locs, COND_SET_SIZE, conditioning = 'mra')
#predsMRA = filter_lean(mra, Y, N_PARTS, lik.params, prior_mean = rep(MEAN, nrow(locs)), init_covparms, saveUQ = "L")
predsMRA = filter(mra, Y, N_PARTS, lik.params, init_covparms, saveUQ = "L")




true_params = list(c = C)
YY = list(x = Y, y = Y)
plotFields(YY, predsMRA, locs)
plotParamPaths(predsMRA$particles, predsMRA$resampled.indices, DATA_MODEL, true_params)
plot1dLogLikelihood(predsMRA$particles, predsMRA$logliks, true_params)
plotMarginalDistribution(predsMRA$particles, predsMRA$resampled.indices, true_params)
par(mfcol = c(2, TMAX))
for (t in 1:(TMAX)) {
    logliks = predsMRA$logliks[[t]]
    finite.inds = which(is.finite(logliks))
    logliks = logliks[finite.inds]
    particles = predsMRA$particles[[t]][finite.inds, "c"]   
    logliks = logliks - max(logliks)
    o = order(predsMRA$particles[[t]][, "c"])
    plot(predsMRA$particles[[t]][o, "c"], exp(logliks[o]), pch = 16, main = t, ylab = "likelihood", type = "l", xaxt = "n", xlab = "")
    plot(predsMRA$particles[[t]][o, "c"], logliks[o], pch = 16, xlab = "c", ylab = "log-likelihood", type = "l")
}
