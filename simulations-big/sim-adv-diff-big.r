setwd("~/vecchiaFilter")
source("simulations-linear/filtering.r")
source('aux-functions.r')
source('getMatCov.r')
source('plotting.r')
Rcpp::sourceCpp('src/getMatCovFromFactor.cpp')
source('scores.r')
resultsDir = "simulations-linear"
library(doParallel)
library(Matrix)
registerDoParallel(cores = 6)


######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 300**2
m = 50
frac.obs = 0.1
Tmax = 2
diffusion = 0.0000001 #diffusion = 0.00004
advection = 0.001 #advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
max.iter = 2

## covariance parameters
sig2 = 0.2; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.05;
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  if (!(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))) {
    stop("One of the models has to be passed as argument")
  } else {
    data.model = args[1]
  }
} else {
  data.model = "gauss"  
}
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 



## set initial state
x0 = matrix(0.5*RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range),
                                          x = locs[,1], y = locs[,2], spConform = FALSE), ncol = 1)


## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
approximations = list(mra = mra, low.rank = low.rank)


RMSPE = list(); LogSc = list()
#foreach( iter=1:max.iter) %dopar% {
for (iter in 1:max.iter) {  

    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
    
    cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
    predsMRA = filter('mra', XY)
    cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
    predsLR  = filter('low.rank', XY)
    
    RMSPE = calculateRMSPE(predsMRA, predsLR, XY$x)

    write.csv(RMSPE, file = paste(resultsDir, "/", data.model, "/RMSPE.", iter, sep = ""))
    
    print(RMSPE)

    if(iter==1){
        plotResults2d(XY, predsMRA, predsLR, data.model, resultsDir, "MRA", "LR")        
    }

    
}
