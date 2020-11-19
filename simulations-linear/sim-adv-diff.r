setwd("~/vecchiaFilter")
source("filtering.r")
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
n = 20**2
m = 40
frac.obs = 0.1
Tmax = 2
diffusion = 0.00004
advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
max.iter = 2
Nparticles = 2

## covariance parameters
sig2 = 0.5; range = .1; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.05;
#args = commandArgs(trailingOnly = TRUE)
#if (length(args) == 1) {
#  if (!(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))) {
#    stop("One of the models has to be passed as argument")
#  } else {
#    data.model = args[1]
#  }
#} else {
  data.model = "gauss"  
#}
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)



## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 



## set initial state
x0 = matrix(0.5*RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range),
                                          x = locs[,1], y = locs[,2], spConform = FALSE), ncol = 1)


## define Vecchia approximation
exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm', verbose = TRUE)
#mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)
#low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
#approximations = list(mra = mra, low.rank = low.rank, exact = exact)
approximations = list(exact = exact)

#RMSPE = list(); LogSc = list()
#scores = foreach( iter=1:max.iter) %dopar% {
for (iter in 1:max.iter) {  

    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)

    cat(paste("iteration: ", iter, ", exact", "\n", sep = ""))
    predsE = filter('exact', XY, saveUQ="L")
    #cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
    #predsMRA = filter('mra', XY, saveUQ="W")
    #cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
    #predsLR  = filter('low.rank', XY, saveUQ="W")
    #cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
    #predsEnKF  = filterEnKF(XY, saveUQ="W")


    
    #RRMSPE = calculateRRMSPE(predsMRA, predsLR, predsE, XY$x)
    #LogSc = calculateLSs(predsMRA, predsLR, predsE, XY$x)

    #write.csv(RRMSPE, file = paste(resultsDir, "/", data.model, "/RRMSPE.", iter, sep=""))
    #write.csv(LogSc, file = paste(resultsDir, "/", data.model, "/LogSc.", iter, sep=""))

    #if(iter==1){
    #    plotResults2d(XY, predsMRA, predsLR, data.model, resultsDir, "MRA", "LR")        
    #}

    #list(RRMSPE, LogSc)
}


#avgRRMSPE = Reduce("+", lapply(scores, `[[`, 1))/length(scores)
#avgdLS = Reduce("+", lapply(scores, `[[`, 2))/length(scores)
#resultsAsString = list("===== avg. RRMSPE: ====\n")
#resultsAsString = c(resultsAsString, capture.output(print(avgRRMSPE)))
#resultsAsString = c(resultsAsString, "\n==== avg. dLS ====\n")
#resultsAsString = c(resultsAsString, capture.output(print(avgdLS)))
#resultsAsString = c(resultsAsString, "\n")
#resultsAsString = paste(resultsAsString, collapse="\n")
#cat(resultsAsString)
#output = paste(c(AllParamsAsString, resultsAsString), sep="\n")
#writeLines(output, paste(resultsDir, "/logs/", hash, sep=""))
