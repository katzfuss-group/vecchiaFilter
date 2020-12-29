setwd("~/vecchiaFilter")
source("filtering.r")
source("filtering-EnKF.r")
source('aux-functions.r')
source('getMatCov.r')
source('plotting.r')
Rcpp::sourceCpp('src/getMatCovFromFactor.cpp')
source('scores.r')
resultsDir = "simulations-linear"
library(doParallel)
library(Matrix)
registerDoParallel(cores = 10)


######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 34**2
m = 50
frac.obs = 0.1
Tmax = 20
diffusion = 0.00004
advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
max.iter = 80


## covariance parameters
sig2 = 0.5; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.05;
data.model = 'gauss'
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## set initial state
x0 = matrix(RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range, var=sig2),
                                          x = locs[,1], y = locs[,2], spConform = FALSE), ncol = 1)


## define Vecchia approximation
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm', verbose = TRUE)
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)
#low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
#approximations = list(mra = mra, low.rank = low.rank, exact = exact)
approximations = list(exact = exact, mra=mra)


RMSPE = list(); LogSc = list()
scores = foreach( iter=1:max.iter) %dopar% {
#for (iter in 1:max.iter) {  

    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)

    cat(paste("MRA", "\n", sep = ""))
    predsMRA = filter('mra', XY, saveUQ="W")    
    
    cat(paste("EnKF", "\n", sep = ""))
    predsEnKF  = filterEnKF(m, XY)
    
    cat(paste("KF", "\n", sep = ""))
    predsE = filter('exact', XY, saveUQ="W")
    
    RRMSPE = calculateRRMSPE(predsMRA, predsEnKF, predsE, XY$x)
    LogSc = calculateLSs(predsMRA, predsEnKF, predsE, XY$x)

    if(iter==1){
        plotResults2d(XY, predsE, predsEnKF, data.model, resultsDir, "Kalman", "EnKF")        
    }

    list(RRMSPE, LogSc)
}
save(scores, file="EnKF.scores")
