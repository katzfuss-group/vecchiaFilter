setwd("~/vecchiaFilter")
source("aux-functions.r")
source("simulations-lorenz/EnKF-filtering.r")
source("simulations-lorenz/aux-functions-Lorenz.r")
source("scores.r")
resultsDir = "simulations-lorenz"
#library(VEnKF)
library(Matrix)
#library(GPvecchia)
#library(foreach)
#library(iterators)
#library(parallel)
library(doParallel)
registerDoParallel(cores=4)



######### set parameters #########
set.seed(1988)
n = 960
m = 50
frac.obs = 0.1
Tmax = 2

## evolution function ##
Force = 10
K = 32
dt = 0.005
M = 5
b = 0.1
evolFun = function(X) b*VEnKF::Lorenz04M2Sim(as.numeric(X)/b, Force, K, dt, M, iter = 1, burn = 0)
max.iter = 1



## covariance function
sig2 = 0.01; range = .15; smooth = 0.5;
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)



## likelihood settings
me.var = 0.1;
data.model = "gauss"
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)



## generate grid of pred.locs
locs = matrix( seq(0, 1, length = round(n)), ncol = 1)


## set initial state
cat("Loading the moments of the long-run Lorenz\n")
moments = getLRMuCovariance(n, Force, dt, K)
Sig0 = (b**2)*moments[["Sigma"]] + diag(1e-10, n)
mu = b*moments[["mu"]]
#x0 = b*getX0(n, Force, K, dt)
x0 = t(chol(Sig0)) %*% matrix(rnorm(n), ncol=1) + mu
#Sigt = sig2*Sig0
Sigt = covfun(locs)



## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', ordering = 'maxmin')
exact = GPvecchia::vecchia_specify(locs, nrow(locs) - 1, ordering = 'maxmin', conditioning = 'firstm')
approximations = list(mra = mra, exact = exact)


scores = foreach( iter=1:max.iter) %dopar% {
#for (iter in 1:max.iter) {

    XY = simulate.xy(x0, evolFun, Sigt, frac.obs, lik.params, Tmax)

    predsEnKF  = filterEnKF(m, XY, saveUQ="W")
    predsE = filterLorenz('exact', XY)
    predsMRA = filterLorenz('mra', XY)    

    RRMSPE = calculateRRMSPE(predsMRA, predsEnKF, predsE, XY$x)
    LogSc = calculateLSs(predsMRA, predsEnKF, predsE, XY$x)
    #write.csv(RRMSPE, file = paste(resultsDir, "/", data.model, "/RRMSPE.", iter, sep = ""))
    #write.csv(LogSc, file = paste(resultsDir, "/", data.model, "/LogSc.", iter, sep = ""))

    if (iter == 1) {
        plotResults(XY, predsE, predsMRA, predsEnKF, resultsDir)
    }
    list(RRMSPE, LogSc)
}
save(scores, file="EnKF.scores.lorenz")
