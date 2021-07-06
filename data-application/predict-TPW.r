cat(sprintf("%s Script started\n", date()))

setwd("/home/marcin/vecchiaFilter")
source("getMatCov.r")
source("aux-functions.r")
source("data-application/settings.r")
#source("data-application/particle-filtering-recompute.r")
source("data-application/particle-filtering.r")
source("data-application/data/process-data.r")
source("data-application/plot-results.r")
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidyverse))

registerDoParallel(cores = NCORES)


init_covparms = c(SIG_02, RANGE, SMOOTH)
lik.params = list(data.model = DATA_MODEL)#, alpha = ALPHA)


#TPW = readr::read_csv("TPW_10.816000k.csv")
TPW = readr::read_csv("data-application/data/TPW.csv", col_types = cols()) %>%
    dplyr::filter( x < X_UPPER_LIMIT ) %>%
    dplyr::filter( x > quantile(x, 1-LOC_QUANTILE), y > quantile(y, 1-LOC_QUANTILE) ) #%>%
    #sample_n(size = 100000)
readr::write_csv(TPW, sprintf("TPW_%fk.csv", nrow(TPW)/1000))


locs = TPW %>% select(x, y) %>% data.matrix()
nx = length(unique(locs[,1]))
ny = length(unique(locs[,2]))

x_min = quantile(locs[, 1], 0.5 - sqrt(FRAC_OBS)/2)
x_max = quantile(locs[, 1], 0.5 + sqrt(FRAC_OBS)/2)
y_min = quantile(locs[, 2], 0.5 - sqrt(FRAC_OBS)/2)
y_max = quantile(locs[, 2], 0.5 + sqrt(FRAC_OBS)/2)
inds.nobs = which(locs[, 1] > x_min & locs[, 1] < x_max & locs[, 2] > y_min & locs[, 2] < y_max)

listOfDataColumns = split.default(TPW, colnames(TPW))[as.character(1:TMAX)]
preProcessColumn = function(column, truncate=TRUE){    
    column = column %>% pull()
    column = column / 1000
    if (truncate) {
        inds.nobs = sample(1:length(column), round(length(column) * (1 - FRAC_OBS)))
        column[inds.nobs] = NA
    }
    return(column)
    #column[column==0] = 1
    #column = as.numeric(remove.mean(column))
}
Y = lapply(listOfDataColumns, preProcessColumn)
Yfull = lapply(listOfDataColumns, function(c) preProcessColumn(c, truncate = FALSE))
       

## filter ---------------------------------------
mra = GPvecchia::vecchia_specify(locs, COND_SET_SIZE, conditioning = 'mra', mra.options = MRA_OPTIONS, verbose = TRUE)
#plot(NULL, xlim = range(locs[, 1]), ylim = range(locs[, 2]), xlab = "lat", ylab = "lon")
#text(mra$locsord[1:30, ], labels = 1:30)
lrf = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray)-1, conditioning = 'firstm')
predsMRA = filter(mra, Y, N_PARTS, lik.params, init_covparms, saveUQ = "L")
predsLRF = filter(lrf, Y, N_PARTS, lik.params, init_covparms, saveUQ = "L")
#predsLRF = predsMRA


T = TMAX
xMRF = zMRF = list()
xLRF = zLRF = list()
for (t in 1:T) {
    obs = Y[[t]]
    xMRF[[t]] = predsMRA$preds[[t]][[1]]$state
    vMRF = predsMRA$vars[[t]][[1]]
    xLRF[[t]] = predsLRF$preds[[t]][[1]]$state
    vLRF = predsLRF$vars[[t]][[1]]
    mu = MEAN_COEFS[t, 1] + MEAN_COEFS[t, 2] * mean(locs[, 2])
    zMRF[[t]] = exp(xMRF[[t]] + vMRF/2 - mu)
    zLRF[[t]] = exp(xLRF[[t]] + vLRF/2 - mu)
}

zlim = c(1, 1)
xlim = c(MEAN)
for (t in 1:T) {
    zlim = range(c(zlim, Y[[t]], zMRF[[t]], zLRF[[t]], Yfull[[t]]), na.rm = TRUE)
    xlim = range(c(xlim, xMRF[[t]], xLRF[[t]]))
}

## conditional simulation
MSE_HV_train = MSE_HV = rep(0, T)
MSE_LR_train = MSE_LR = rep(0, T)
n = nrow(locs)
for(t in 1:T) {
    #pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=5, height=20)
    pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=10, height=10)
    oldpar = par(mfrow = c(2, 2))
    #oldpar = par(mfrow = c(4, 1), mar = c(1, 1, 1, 1), cex = 1.2)
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(Y[[t]]), zlim = zlim, nx = nx, ny = ny, main = sprintf("observations, t=%d", t))
    par(cex = 1.2)
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(zMRF[[t]]), zlim = zlim, nx = nx, ny = ny, main = "predictions HV")
    par(cex = 1.2)
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(zLRF[[t]]), zlim = zlim, nx = nx, ny = ny, main = "predictions LR")
    par(cex = 1.2)
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(Yfull[[t]]), zlim = zlim, nx = nx, ny = ny, main = sprintf("full data, t=%d", t))
    #fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(xMRF[[t]]), zlim = xlim, nx = nx, ny = ny, main = sprintf("HV %d", t))
    #fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(xLRF[[t]]), zlim = xlim, nx = nx, ny = ny, main = sprintf("LR %d", t))
    par(oldpar)
    dev.off()
    ind.obs = which(!is.na(Y[[t]]))
    ind.pred = setdiff(which(!is.na(Yfull[[t]])), ind.obs)
    MSE_HV[t] = sqrt(mean((Yfull[[t]][ind.pred] - zMRF[[t]][ind.pred])**2, na.rm = TRUE))
    MSE_LR[t] = sqrt(mean((Yfull[[t]][ind.pred] - zLRF[[t]][ind.pred])**2, na.rm = TRUE))
    MSE_HV_train[t] = sqrt(mean((Yfull[[t]][ind.obs] - zMRF[[t]][ind.obs])**2, na.rm = TRUE))
    MSE_LR_train[t] = sqrt(mean((Yfull[[t]][ind.obs] - zLRF[[t]][ind.obs])**2, na.rm = TRUE))
}

ylim = range(MSE_HV, MSE_LR, 0)
defpar = par(mfrow = c(1, 2))
plot(1:TMAX, MSE_HV, type = "l", col = "red", ylim = ylim, ylab = "RMSE", xlab = "time", lwd = 2, main = "Test set error")
lines(1:TMAX, MSE_LR, type = "l", col = "blue", lty = 2, lwd = 2)
abline(h=0, lwd = 2)
plot(1:TMAX, MSE_HV_train, type = "l", col = "red", ylim = ylim, ylab = "RMSE", xlab = "time", lwd = 2, main = "Training set error")
lines(1:TMAX, MSE_LR_train, type = "l", col = "blue", lty = 2, lwd = 2)
abline(h=0, lwd = 2)
par(defpar)
