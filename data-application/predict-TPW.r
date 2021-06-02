cat(sprintf("%s Script started\n", date()))

setwd("~/vecchiaFilter")
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
lik.params = list(data.model = DATA_MODEL)#, alpha = ALPHA)


#TPW = readr::read_csv("TPW_10k.csv")
TPW = readr::read_csv("data-application/data/TPW.csv", col_types = cols()) %>%
    dplyr::filter( x < X_UPPER_LIMIT ) %>%
    dplyr::filter( x > quantile(x, 1-LOC_QUANTILE), y > quantile(y, 1-LOC_QUANTILE) ) #%>%
    #sample_n(size = 100000)
#readr::write_csv(TPW, sprintf("TPW_%fk.csv", nrow(TPW)/1000))


locs = TPW %>% select(x, y) %>% data.matrix()
nx = length(unique(locs[,1]))
ny = length(unique(locs[,2]))

listOfDataColumns = split.default(TPW, colnames(TPW))[as.character(1:TMAX)]

preProcessColumn = function(column, truncate=TRUE){
    column = column %>% pull()
    column = column / 1000
    if (truncate) {
        inds.obs = sample(1:length(column), round(length(column) * FRAC_OBS))
        column[-inds.obs] = NA
    }
    return(column)
    #column[column==0] = 1
    #column = as.numeric(remove.mean(column))
}
Y = lapply(listOfDataColumns, preProcessColumn)
Yfull = lapply(listOfDataColumns, function(c) preProcessColumn(c, truncate = FALSE))
       

## filter ---------------------------------------
mra = GPvecchia::vecchia_specify(locs, COND_SET_SIZE, conditioning = 'mra')
lrf = GPvecchia::vecchia_specify(locs, COND_SET_SIZE, conditioning = 'firstm')
#predsMRA = filter_lean(mra, Y, N_PARTS, lik.params, prior_mean = rep(MEAN, nrow(locs)), init_covparms, saveUQ = "L")
predsMRA = filter(mra, Y, N_PARTS, lik.params, init_covparms, saveUQ = "L")
predsLRF = filter(lrf, Y, N_PARTS, lik.params, init_covparms, saveUQ = "L")



T = TMAX
xMRF = zMRF = list()
xLRF = zLRF = list()
for (t in 1:T) {
    obs = Y[[t]]
    xMRF[[t]] = predsMRA$preds[[t]][[1]]$state
    xLRF[[t]] = predsLRF$preds[[t]][[1]]$state
    na.inds = which(is.na(obs))
    zMRF[[t]] = zLRF[[t]] = obs
    zMRF[[t]][na.inds] = rgamma(nrow(locs[na.inds,]), ALPHAS[t], ALPHAS[t]*exp(-xMRF[[t]][na.inds]))
    zLRF[[t]][na.inds] = rgamma(nrow(locs[na.inds,]), ALPHAS[t], ALPHAS[t]*exp(-xLRF[[t]][na.inds]))
}
zlim = c(1, 1)
xlim = c(MEAN)
for (t in 1:T) {
    zlim = range(c(zlim, Y[[t]], zMRF[[t]], zLRF[[t]], Yfull[[t]]), na.rm = TRUE)
    xlim = range(c(xlim, xMRF[[t]], xLRF[[t]]))
}
## conditional simulation
MSE_HV = rep(0, T)
MSE_LR = rep(0, T)
n = nrow(locs)
for(t in 1:T) {
    pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=18, height=9)
    oldpar = par(mfrow = c(2, 3))
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(Y[[t]]), zlim = zlim, nx = nx, ny = ny, main = sprintf("data %d", t))
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(zMRF[[t]]), zlim = zlim, nx = nx, ny = ny, main = "reconstruction HV")
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(zLRF[[t]]), zlim = zlim, nx = nx, ny = ny, main = "reconstruction LR")
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(Yfull[[t]]), zlim = zlim, nx = nx, ny = ny, main = sprintf("truth %d", t))
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(xMRF[[t]]), zlim = xlim, nx = nx, ny = ny, main = sprintf("HV %d", t))
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(xLRF[[t]]), zlim = xlim, nx = nx, ny = ny, main = sprintf("LR %d", t))
    par(oldpar)
    dev.off()
    MSE_HV[t] = sqrt(mean((Yfull[[t]] - zMRF[[t]])**2, na.rm = TRUE))
    MSE_LR[t] = sqrt(mean((Yfull[[t]] - zLRF[[t]])**2, na.rm = TRUE))
}
print(mean(MSE_HV))
print(mean(MSE_LR))
