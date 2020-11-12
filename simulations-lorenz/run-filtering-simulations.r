setwd("~/vecchiaFilter")
rm(list = ls())
source("aux-functions.r")
source("simulations-lorenz/aux-functions-Lorenz.r")
source("scores.r")
resultsDir = "simulations-lorenz"
library(VEnKF)
library(Matrix)
library(GPvecchia)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
registerDoParallel(cores=5)



AllParamsAsString = list()



######### set parameters #########
set.seed(1988)
n = 96
m = 20
frac.obs = 0.1
Tmax = 20

## evolution function ##
Force = 10
K = 10
dt = 0.005
M = 5
b = 0.1
evolFun = function(X) b*Lorenz04M2Sim(as.numeric(X)/b, Force, K, dt, M, iter = 1, burn = 0)
max.iter = 1#getDoParWorkers()



## covariance function
sig2 = 0.01; range = .15; smooth = 0.5;
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)


## likelihood settings
me.var = 0.1;
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


AllParamsAsString = list("======= Lorenz simulation ========\n",
                         "started at: ", format(Sys.time(), "%H:%M:%S, %a, %b %d"), "\n",
                         "iterations: ", max.iter, "\n",
                         "T = ", Tmax, ", frac. obs = ", frac.obs, ", n = ", n, ", m = ", m, "\n",
                         "\n",
                         "+++ likelihood +++\n",
                         "data model: ", data.model, "\n")
if(data.model == "gauss") {
    AllParamsAsString = c(AllParamsAsString, list("sd = ", lik.params[["sigma"]], "\n"))
} else if(data.model == "gamma") {
    AllParamsAsString = c(AllParamsAsString, list("alpha = ", lik.params[["alpha"]], "\n"))
}
AllParamsAsString = c(AllParamsAsString,
                    list("\n+++ Lorenz settings +++\n",
                         "Force = ", Force, ", K = ", K, ", dt = ", dt, ", M = ", M, ", b = ", b, "\n",
                         "\n",
                         "+++ model error covariance +++\n",
                         "sig2 = ", sig2, ", range = ", range, ", smoothness = ", smooth, "\n"
                         ))

AllParamsAsString = paste(AllParamsAsString, collapse="")
hash = digest::digest(AllParamsAsString, serialize=FALSE)
cat(AllParamsAsString)

## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(n))
locs = matrix(grid.oneside, ncol = 1)


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
cat("Calculating the approximations\n")
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', ordering = 'maxmin')
exact = GPvecchia::vecchia_specify(locs, nrow(locs) - 1, ordering = 'maxmin', conditioning = 'firstm')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, ordering = 'maxmin', conditioning = 'firstm')
approximations = list(mra = mra, low.rank = low.rank, exact = exact)


scores = foreach( iter=1:max.iter) %dopar% {
#for (iter in 1:max.iter) {

    cat("Simulating data\n")
    XY = simulate.xy(x0, evolFun, Sigt, frac.obs, lik.params, Tmax)

    cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
    start = proc.time()
    predsMRA = filter('mra', XY)
    d = as.numeric(proc.time() - start)
    cat(paste("MRA filtering took ", d[3], "s\n", sep = ""))
    
    cat(paste("iteration: ", iter, ", exact", "\n", sep = ""))
    start = proc.time()
    predsE = filter('exact', XY)
    d = as.numeric(proc.time() - start)
    cat(paste("Exact filtering took ", d[3], "s\n", sep = ""))

    cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
    start = proc.time()
    predsLR  = filter('low.rank', XY)
    d = as.numeric(proc.time() - start)
    cat(paste("Low-rank filtering took ", d[3], "s\n", sep = ""))

    RRMSPE = calculateRRMSPE(predsMRA, predsLR, predsE, XY$x)
    LogSc = calculateLSs(predsMRA, predsLR, predsE, XY$x)
    write.csv(RRMSPE, file = paste(resultsDir, "/", data.model, "/RRMSPE.", iter, sep = ""))
    write.csv(LogSc, file = paste(resultsDir, "/", data.model, "/LogSc.", iter, sep = ""))

    if (iter == 1) {
      plotResults(XY, predsE, predsMRA, predsLR, resultsDir)
    }
    list(RRMSPE, LogSc)
}


avgRRMSPE = Reduce("+", lapply(scores, `[[`, 1))/length(scores)
avgdLS = Reduce("+", lapply(scores, `[[`, 2))/length(scores)
resultsAsString = list("===== avg. RRMSPE: ====\n")
resultsAsString = c(resultsAsString, capture.output(print(avgRRMSPE)))
resultsAsString = c(resultsAsString, "\n==== avg. dLS ====\n")
resultsAsString = c(resultsAsString, capture.output(print(avgdLS)))
resultsAsString = c(resultsAsString, "\n")
resultsAsString = paste(resultsAsString, collapse="\n")
cat(resultsAsString)
output = paste(c(AllParamsAsString, resultsAsString), sep="\n")
writeLines(output, paste(resultsDir, "/logs/", hash, sep=""))
