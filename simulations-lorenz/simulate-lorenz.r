setwd("~/HVLF")
rm(list = ls())
source("aux-functions.r")
source("scores.r")
resultsDir = "simulations-lorenz"
library(VEnKF)
library(rootSolve)
library(doParallel)
registerDoParallel(cores = 6)



filter = function(approx.name, XY){
  
  approx = approximations[[approx.name]]
  
  preds = list()
  L.tt = getL00(approx, Sig0, locs)
  mu.tt = x0
  preds[[1]] = list(state = x0, L = L.tt)
  
  for (t in 2:Tmax) {
    
    #cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    Et = Matrix::Matrix(gradient(evolFun, mu.tt, centered = TRUE))
    
    forecast = evolFun(mu.tt)
    Fmat = Et %*% L.tt
    covmodel = GPvecchia::getMatCov(approx, Fmat %*% Matrix::t(Fmat) + sig2*Sig0)
    
    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                                   likelihood_model = data.model, covmodel = covmodel,
                                                   covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W)
  }
  return( preds )
}



center_operator <- function(x) {
  n = nrow(x)
  ones = rep(1, n)
  H = diag(n) - (1/n) * (ones %*% t(ones))
  H %*% x
}


getCovariance = function(N, Force, dt, K){
  fileName.all = paste("simulations-lorenz/Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep = "_")
  
  X = Matrix::Matrix(scan(fileName.all, quiet = TRUE), nrow = N) 
  X = center_operator(X)
  S = (X %*% Matrix::t(X)) / (ncol(X) - 1)
}



######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 960
m = 50
frac.obs = 0.1
Tmax = 20


## evolution function ##
Force = 10
K = 32
dt = 0.005
M = 1
evolFun = function(X) Lorenz04M2Sim(as.numeric(X), Force, K, dt, M, iter = 1, burn = 0, order = 4)
max.iter = 100


## covariance function
sig2 = 0.1; range = .15; smooth = 1.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)


## likelihood settings
me.var = 0.2;
data.model = "gamma"
lik.params = list(data.model = data.model, me.var = me.var, alpha = 2)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(n))
locs = matrix(grid.oneside, ncol = 1)

## set initial state
Sig0 = getCovariance(n, Force, dt, K)
x0 = Matrix::t(chol(Sig0)) %*% matrix(rnorm(n), ncol = 1)#getX0(n, Force, K, dt)

## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)
exact = GPvecchia::vecchia_specify(locs, nrow(locs) - 1, conditioning = 'firstm')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm')
approximations = list(mra = mra, low.rank = low.rank, exact = exact)


RRMSPE = list(); LogSc = list()
foreach( iter=1:max.iter) %dopar% {
#for (iter in 1:max.iter) {  
  
  XY = simulate.xy(x0, evolFun, Sig0, frac.obs, lik.params, Tmax)
  
  cat(paste("iteration: ", iter, ", exact", "\n", sep = ""))
  start = proc.time()
  predsE = filter('exact', XY)
  d = as.numeric(proc.time() - start)
  cat(paste("Exact filtering took ", d[3], "s\n", sep = ""))
  cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
  start = proc.time()
  predsMRA = filter('mra', XY)
  d = as.numeric(proc.time() - start)
  cat(paste("MRA filtering took ", d[3], "s\n", sep = ""))
  cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
  start = proc.time()
  predsLR  = filter('low.rank', XY)
  d = as.numeric(proc.time() - start)
  cat(paste("Low-rank filtering took ", d[3], "s\n", sep = ""))
  
  RRMSPE = calculateRRMSPE(predsMRA, predsLR, predsE, XY$x)
  LogSc = calculateLSs(predsMRA, predsLR, predsE, XY$x)
  write.csv(RRMSPE, file = paste(resultsDir, "/", data.model, "/RRMSPE.", iter, sep = ""))
  write.csv(LogSc, file = paste(resultsDir, "/", data.model, "/LogSc.", iter, sep = ""))
  
  #data = list(XY = XY, predsMRA = predsMRA, predsE = predsE, predsLR = predsLR)
  #save(data, file = paste(resultsDir, "/", data.model, "/sim.", iter, sep = ""))
}
