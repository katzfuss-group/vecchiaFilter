## imports ------------------------------------
setwd("~/vecchiaFilter")
source("data-application/particle-filtering.r")
source("aux-functions.r")
source("data-application/settings.r")
library(doParallel)
registerDoParallel(cores = NCORES)


init_covparms = c(SIG_02, RANGE, SMOOTH)
evolFun = function(X) C * X
lik.params = list(data.model = DATA_MODEL, alpha = ALPHA, sigma = sqrt(ME_VAR))


## simulate data --------------------------------
grid.oneside = seq(0, 1, length = round(sqrt(N)))
locs = as.matrix(expand.grid(grid.oneside, grid.oneside))
mu0 = matrix(rep(0, N), ncol = 1)
XY = simulate.xy(mu0, SIG_02, RANGE, SMOOTH, evolFun, NULL, FRAC_OBS, lik.params, TMAX, sig2 = SIG2, smooth = SMOOTH, range = RANGE, locs = locs)


## filter ---------------------------------------
mra = GPvecchia::vecchia_specify(locs, COND_SET_SIZE, conditioning = 'mra')
predsMRA = filter(mra, XY$y, N_PARTS, lik.params, init_covparms, saveUQ = "L")


## Plots ---------------------------------------- 
# latent field + likelihood
zlim_x = range(sapply(1:TMAX, function(t) range(XY$x[[t]])), XY$x0)
zlim_y = range(sapply(1:TMAX, function(t) range(XY$y[[t]])))
zlim_p = range(sapply(1:TMAX, function(t) range(predsMRA$preds[[t]]$state)))


pdf(sprintf("data-application/tests/init-field.pdf"))
fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(XY$x0), zlim = zlim_x, nx = sqrt(N), ny = sqrt(N), main = "init")
dev.off()

for (t in 1:TMAX) {
    pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=18, height=6)
    oldpar = par(mfrow = c(1, 3))
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(XY$x[[t]]), zlim = zlim_x, nx = sqrt(N), ny = sqrt(N), main = sprintf("truth %d", t))
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(XY$y[[t]]), zlim = zlim_y, nx = sqrt(N), ny = sqrt(N), main = sprintf("data %d", t))
    fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(predsMRA$preds[[t]]$state), zlim = zlim_p, nx = sqrt(N), ny = sqrt(N), main = "reconstruction")
    par(oldpar)
    dev.off()
}
par(oldpar)


## parameters
pdf(sprintf("data-application/tests/%s-particles-over-time.pdf", DATA_MODEL), width = 6, height = 4)
oldpar = par(mfrow = c(2, 3))

for (par.name in colnames(predsMRA$particles[[1]])) {    
    time = 1:TMAX
    
    uq = sapply(predsMRA$particles, function(t) as.numeric(quantile(t[, par.name ], 0.9)))
    means = sapply(predsMRA$particles, function(t) as.numeric(mean(t[, par.name ])))
    lq = sapply(predsMRA$particles, function(t) as.numeric(quantile(t[, par.name ], 0.1)))

    if (all(uq == lq)) {
       cat(sprintf("%s was not sampled\n", par.name))
       next()
    }
    
    ylim = range(c(uq, lq))
    ylim[2] = 1.1 * ylim[2]
    plot(time, uq, type = "l", lty = 2, ylim = ylim, ylab = par.name, main = par.name)
    lines(time, means, type = "l")
    lines(time, lq, type = "l", lty = 2)    
}
par(oldpar)
dev.off()


## likelihood
param = "sig2"
truth = list(c = C, sig2 = SIG2, range = RANGE)
params = predsMRA$particles[[TMAX]][, param]
ord = order(params)
plot(params[ord], predsMRA$logliks[[1]][ord], type = "l", xlab = param, ylab = "likelihood")
abline(v = truth[[param]], col = "red")
