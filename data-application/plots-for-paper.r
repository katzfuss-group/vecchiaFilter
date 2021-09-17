cat(sprintf("%s Script started\n", date()))
setwd("/home/marcin/vecchiaFilter")
suppressPackageStartupMessages(library(tidyverse))


## init_covparms = c(SIG_02, RANGE, SMOOTH)
load("~/vecchiaFilter/predict-TPW-MIRS.RData")
results = tibble()

TMAX = 11

for (t in 1:TMAX) {
    xMRF_tibble = tibble(locs, value = xMRF[[t]], val_type = "xMRF", day = t)
    vMRF_tibble = tibble(locs, value = vMRF[[t]], val_type = "vMRF", day = t)
    xLRF_tibble = tibble(locs, value = xLRF[[t]], val_type = "xLRF", day = t)
    vLRF_tibble = tibble(locs, value = vLRF[[t]], val_type = "vLRF", day = t)
    Y_tibble = tibble(locs, value = Y[[t]], val_type = "data", day = t)
    Yf_tibble = tibble(locs, value = Yf[[t]], val_type = "full data", day = t)

    results = bind_rows(results, xMRF_tibble, vMRF_tibble, xLRF_tibble, vLRF_tibble, Y_tibble, Yf_tibble)
}
rm(list = c("xMRF", "vMRF", "xLRF", "vLRF"))




all_lon = results %>% select(lon) %>% unique() %>% pull()
all_lat = results %>% select(lat) %>% unique() %>% pull()
border = c(0.03, 0.97)
border_lat = quantile(all_lat, border)
border_lon = quantile(all_lon, border)
nx = length(which(all_lon >= border_lon[1] & all_lon <= border_lon[2]))
ny = length(which(all_lat >= border_lat[1] & all_lat <= border_lat[2]))


results_no_border = results %>%
    dplyr::filter(lon >= quantile(unique(lon), border[1]), lon <= quantile(unique(lon), border[2])) %>%
    dplyr::filter(lat >= quantile(unique(lat), border[1]), lat <= quantile(unique(lat), border[2]))


## plot filtering means ------------------------
xlim = results_no_border %>% dplyr::filter(val_type %in% c("xMRF", "xLRF", "full data")) %>%
    select(value) %>% summarize(range(value, na.rm = TRUE)) %>% pull()
vlim = results_no_border %>% dplyr::filter(val_type %in% c("vMRF", "vLRF")) %>%
    select(value) %>% summarize(range(sqrt(value), na.rm = TRUE)) %>% pull()


MSE_LR_train = MSE_LR = MSE_HV_train = MSE_HV = rep(0, TMAX)
for (t in 1:TMAX) {
    #pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=5, height=20)
    pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width = 5, height = 13)
    oldpar = par(mfcol = c(3, 1))


    res_day = dplyr::filter(results_no_border, day == t)
    long = dplyr::filter(res_day, val_type == "full data") %>% dplyr::select(lon) %>% pull()
    lati = dplyr::filter(res_day, val_type == "full data") %>% dplyr::select(lat) %>% pull()
    Yfull = dplyr::filter(res_day, val_type == "full data") %>% select(value) %>% pull()
    Y = dplyr::filter(res_day, val_type == "data") %>% select(value) %>% pull()
    xMRF = dplyr::filter(res_day, val_type == "xMRF") %>% select(value) %>% pull()
    xLRF = dplyr::filter(res_day, val_type == "xLRF") %>% select(value) %>% pull()
    sdMRF = dplyr::filter(res_day, val_type == "vMRF") %>% select(value) %>% mutate(value = sqrt(value)) %>% pull()
    sdLRF = dplyr::filter(res_day, val_type == "vLRF") %>% select(value) %>% mutate(value = sqrt(value)) %>% pull()

    
    ind.obs = which(!is.na(Y))
    ind.pred = setdiff(which(!is.na(Yfull)), ind.obs)
    MSE_HV[t] = sqrt(mean((Yfull[ind.pred] - xMRF[ind.pred])**2, na.rm = TRUE))
    MSE_LR[t] = sqrt(mean((Yfull[ind.pred] - xLRF[ind.pred])**2, na.rm = TRUE))
    MSE_HV_train[t] = sqrt(mean((Yfull[ind.obs] - xMRF[ind.obs])**2, na.rm = TRUE))
    MSE_LR_train[t] = sqrt(mean((Yfull[ind.obs] - xLRF[ind.obs])**2, na.rm = TRUE))
    
    
    par(cex = 1.3)
    fields::quilt.plot(long, lati, Yfull, zlim = xlim, nx = nx, ny = ny)#, main = sprintf("full data, t=%d", t))
    par(cex = 1.3)
    fields::quilt.plot(long, lati, xMRF,  zlim = xlim, nx = nx, ny = ny)#, main = "predictions HV")
    par(cex = 1.3)
    fields::quilt.plot(long, lati, xLRF,  zlim = xlim, nx = nx, ny = ny)#, main = "predictions LRF")
    par(oldpar)
    dev.off()

    pdf(sprintf("data-application/tests/test-err-%d.pdf", t), width = 5, height = 10)
    oldpar = par(mfcol = c(2, 1))
    par(cex = 1.3)#, mar = rep(2, 4))
    fields::quilt.plot(long, lati, sdLRF, zlim = vlim, nx = nx, ny = ny)#, main = "sd LRF")
    par(cex = 1.3)#, mar = rep(2, 4))
    fields::quilt.plot(long, lati, sdMRF, zlim = vlim, nx = nx, ny = ny)#, main = "sd HV")
    par(oldpar)
    dev.off()
}



ylim = range(c(0, MSE_HV, MSE_LR))
xlim = c(0, TMAX+1)
pdf("data-application/RMSEscoresDataApp.pdf", width = 8, height = 4)
plot(1:TMAX, MSE_HV, col = "red", type = "l", ylim = ylim, ylab = "RMSPE", xlab = "time", lwd = 2)
lines(1:TMAX, MSE_LR, col = "blue", lwd = 2)
legend("topleft", legend = c("HV", "LR"), col = c("red", "blue"), lty = c(1, 1), lwd = c(2, 2))
dev.off()
