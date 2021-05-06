## imports ------------------------------------
setwd("~/vecchiaFilter")
source("data-application/test-recompute-settings.r")
source("data-application/particle-filtering.r")
source("data-application/particle-filtering-recompute.r")
source("aux-functions.r")
source("data-application/plot-results.r")
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
preds = filter(mra, XY$y, N_PARTS, lik.params, init_covparms, saveUQ = "L", seed = SEED)

preds_r = filter_lean(mra, XY$y, N_PARTS, lik.params, init_covparms, saveUQ = "L", seed = SEED)

success = TRUE
for (t in 1:TMAX) {
    if (length(preds$preds[[t]]) != length(preds_r$preds[[t]])) {
        success = FALSE
        stop(sprintf("Problem at t = %d\n", t))
    }
    for (l in 1:length(preds$preds[[t]])) {
        d = preds$preds[[t]][[l]]$state - preds_r$preds[[t]][[l]]$state
        if (sum(abs(d))>0) {
            success = FALSE
            stop(sprintf("Problem at l = %d and t = %d\n", l, t))
        }
    }
}
cat(sprintf("All matches up: %s\n", success))
    
