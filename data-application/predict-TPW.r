cat(sprintf("%s Script started\n", date()))

setwd("/home/marcin/vecchiaFilter")
source("getMatCov.r")
source("data-application/aux-functions.r")
source("data-application/settings.r")
source("data-application/filtering.r")
suppressPackageStartupMessages(library(tidyverse))


init_covparms = c(SIG_02, RANGE, SMOOTH)


## load data -----------------------------------
TPW = readr::read_csv(FILENAME, col_types = cols()) 

locs = TPW %>% select(lon, lat) %>% data.matrix()
lat = sort(unique(locs[, "lat"]))
lon = sort(unique(locs[, "lon"]))[1:length(lat)]
locs = as_tibble(expand.grid(lon, lat))
colnames(locs) = c("lon", "lat")

Y = list()
Yf = list()
for (p in unique(TPW$day)) {
    if (p <= TMAX) {
        Yf[[p]] = TPW %>% dplyr::filter(day==p) %>% select(-day)
        Yf[[p]] = Yf[[p]] %>% right_join(locs, by = c("lon", "lat"))
        Yf[[p]] = Yf[[p]] %>%  arrange(lon, lat) %>% select(value) %>% pull()
        inds.obs = sample(1:length(Yf[[p]]), size = FRAC_OBS * length(Yf[[p]]))
        Y[[p]] = Yf[[p]]
        Y[[p]][-inds.obs] = NA
    }
}
nx = length(unique(locs$lon))
ny = length(unique(locs$lat))


## filter ---------------------------------------
mra = GPvecchia::vecchia_specify(locs %>% data.matrix(), COND_SET_SIZE, conditioning = 'mra', mra.options = MRA_OPTIONS, verbose = TRUE)
predsMRA = filter(mra, Y, lik.params, init_covparms)
lrf = GPvecchia::vecchia_specify(locs %>% data.matrix(), ncol(mra$U.prep$revNNarray)-1, conditioning = 'firstm')
predsLRF = filter(lrf, Y, lik.params, init_covparms)


## process and plot results --------------------
vMRF = vLRF = xLRF = zLRF = xMRF = zMRF = list()
for (t in 1:TMAX) {
    obs = Y[[t]]
    xMRF[[t]] = as.numeric(predsMRA$preds[[t]][[1]]$state)
    vMRF[[t]] = as.numeric(predsMRA$vars[[t]][[1]])
    xLRF[[t]] = as.numeric(predsLRF$preds[[t]][[1]]$state)
    vLRF[[t]] = as.numeric(predsLRF$vars[[t]][[1]])

    if (lik.params[["data.model"]]=="gamma") {
        mu = as.numeric(MEAN_COEFS[t, 1] + MEAN_COEFS[t, 2] * mean(locs$lat))
        zMRF[[t]] = as.numeric(exp(xMRF[[t]] + vMRF[[t]]/2 - mu))
        zLRF[[t]] = as.numeric(exp(xLRF[[t]] + vLRF[[t]]/2 - mu))
    } else {
        zMRF = xMRF
        zLRF = xLRF
    }
}



## plot filtering means ------------------------
zlim = c(1, 1)
xlim = c(1)
vlim = c(0, 0)
for (t in 1:TMAX) {
    zlim = range(c(zlim, Y[[t]], zMRF[[t]], zLRF[[t]], Yf[[t]]), na.rm = TRUE)
    xlim = range(c(xlim, xMRF[[t]], xLRF[[t]]))
    vlim = range(c(vlim, sqrt(vMRF[[t]]), sqrt(vLRF[[t]])))
}

for (t in 1:TMAX) {
    #pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=5, height=20)
    pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=15, height=10)
    oldpar = par(mfcol = c(2, 3))
    fields::quilt.plot(locs$lon, locs$lat, Yf[[t]], zlim = zlim, nx = nx, ny = ny, main = sprintf("full data, t=%d", t))
    fields::quilt.plot(locs$lon, locs$lat, Y[[t]], zlim = zlim, nx = nx, ny = ny, main = sprintf("observations, t=%d", t))
    fields::quilt.plot(locs$lon, locs$lat, zMRF[[t]], zlim = zlim, nx = nx, ny = ny, main = "predictions HV")
    fields::quilt.plot(locs$lon, locs$lat, zLRF[[t]], zlim = zlim, nx = nx, ny = ny, main = "predictions LRF")
    fields::quilt.plot(locs$lon, locs$lat, sqrt(vLRF[[t]]), zlim = vlim, nx = nx, ny = ny, main = "sd LRF")
    fields::quilt.plot(locs$lon, locs$lat, sqrt(vMRF[[t]]), zlim = vlim, nx = nx, ny = ny, main = "sd HV")
    par(oldpar)
    dev.off()
}


MSE_LR_train = MSE_LR = MSE_HV_train = MSE_HV = rep(0, TMAX)
for(t in 1:TMAX) {
    ind.obs = which(!is.na(Y[[t]]))
    ind.pred = setdiff(which(!is.na(Yf[[t]])), ind.obs)
    MSE_HV[t] = sqrt(mean((Yf[[t]][ind.pred] - zMRF[[t]][ind.pred])**2, na.rm = TRUE))
    MSE_LR[t] = sqrt(mean((Yf[[t]][ind.pred] - zLRF[[t]][ind.pred])**2, na.rm = TRUE))
    MSE_HV_train[t] = sqrt(mean((Yf[[t]][ind.obs] - zMRF[[t]][ind.obs])**2, na.rm = TRUE))
    MSE_LR_train[t] = sqrt(mean((Yf[[t]][ind.obs] - zLRF[[t]][ind.obs])**2, na.rm = TRUE))
}

ylim = range(c(MSE_HV, MSE_LR))
plot(MSE_HV[1:9], col = "red", type = "l", ylim = ylim, ylab = "RMSE", xlab = "time", lwd = 2)
lines(MSE_LR[1:9], col = "blue", lwd = 2)
legend("topleft", legend = c("HV", "LR"), col = c("red", "blue"), lty = c(1, 1))
