## imports ------------------------------------
setwd("~/vecchiaFilter")
source("data-application/test-filter-settings.r")
#source("data-application/particle-filtering.r")
#source("data-application/particle-filtering-recompute.r")
source("data-application/particle-filtering.r")
source("aux-functions.r")
source("data-application/plot-results.r")
library(doParallel)
registerDoParallel(cores = NCORES)


init_covparms = c(SIG_02, RANGE, SMOOTH)
lik.params = list(data.model = DATA_MODEL, alpha = ALPHA, sigma = sqrt(ME_VAR))


## simulate data --------------------------------
grid.oneside = seq(0, 0.9282838, length = round(sqrt(N)))
locs = as.matrix(expand.grid(grid.oneside, grid.oneside))
`%>%` = dplyr::`%>%`
locs = readr::read_csv("locs.csv", col_types = readr::cols()) %>% data.matrix()
N = nrow(locs)
mu0 = matrix(rep(0, N), ncol = 1)
mu0 = diffAdvVec2d(sqrt(N), ny=sqrt(N), height=1, rnge=10)
XY = simulate.xy(mu0, SIG_02, RANGE, SMOOTH, evolFun, NULL, FRAC_OBS, lik.params, TMAX, sig2 = SIG2, smooth = SMOOTH, range = RANGE, locs = locs)


## filter ---------------------------------------
mra = GPvecchia::vecchia_specify(locs, COND_SET_SIZE, conditioning = 'mra')
predsMRA = filter(mra, XY$y, N_PARTS, lik.params, init_covparms, saveUQ = "L")


## Plots ----------------------------------------
true_params = list(c = C, sig2 = SIG2, range = RANGE)
preds = list()
for(t in 1:TMAX){
    print(t)
    ind = predsMRA$resampled.indices[[t]][1]
    preds[[t]] = predsMRA$preds[[t]][[ind]]
}
plotFields(XY, preds, locs)
#plotParamPaths(predsMRA$particles, predsMRA$resampled.indices, DATA_MODEL, true_params)
#plotMarginalDistribution(predsMRA$particles, predsMRA$resampled.indices, true_params)
#plot1dLogLikelihood(predsMRA$particles, predsMRA$logliks, true_params)
#par(mfcol = c(2, TMAX))
#for (t in 1:TMAX) {
#    logliks = predsMRA$logliks[[t]]
#    logliks = logliks - max(logliks)
#    o = order(predsMRA$particles[[t]][, "c"])
#    plot(predsMRA$particles[[t]][o, "c"], exp(logliks[o]), pch = 16, main = t, ylab = "likelihood", type = "l", xaxt = "n", xlab = "")
#    plot(predsMRA$particles[[t]][o, "c"], logliks[o], pch = 16, xlab = "c", ylab = "log-likelihood", type = "l")
#}
