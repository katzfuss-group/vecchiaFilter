setwd("~/HVLF")
rm(list = ls())
source("aux-functions.r")
source("scores.r")
resultsDir = "simulations-lorenz"
library(VEnKF)
library(rootSolve)


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
    covmodel = GPvecchia::getMatCov(approx, Fmat %*% Matrix::t(Fmat) + Q)
    
    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                                   likelihood_model = data.model, covmodel = covmodel,
                                                   covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W)
  }
  return( preds )
}



######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 120
m = 20
frac.obs = 0.1
Tmax = 20


## evolution function ##
Force = 10
K = 32
dt = 0.005
M = 1
evolFun = function(X) Lorenz04M2Sim(as.numeric(X), Force, K, dt, M, iter = 1, burn = 1, order = 4)
max.iter = 1


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
x0 = matrix(rnorm(n), ncol = 1)#getX0(n, Force, K, dt)
Q = covfun(locs)
Sig0 = (1/sig2)*covfun(locs)


XY = simulate.xy(x0, evolFun, Q, frac.obs, lik.params, Tmax)

## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)
exact = GPvecchia::vecchia_specify(locs, nrow(locs) - 1, conditioning = 'firstm')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm')
approximations = list(mra = mra, low.rank = low.rank, exact = exact)


RRMSPE = list(); LogSc = list()
#foreach( iter=1:max.iter) %dopar% {
for (iter in 1:max.iter) {  
  
  XY = simulate.xy(x0, evolFun, Q, frac.obs, lik.params, Tmax)
  
  cat(paste("iteration: ", iter, ", exact", "\n", sep = ""))
  predsE = filter('exact', XY)
  cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
  predsMRA = filter('mra', XY)
  cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
  predsLR  = filter('low.rank', XY)
  
  RRMSPE = calculateRRMSPE(predsMRA, predsLR, predsE, XY$x)
  LogSc = calculateLSs(predsMRA, predsLR, predsE, XY$x)
  write.csv(RRMSPE, file = paste(resultsDir, "/", data.model, "/RRMSPE.", iter, sep = ""))
  write.csv(LogSc, file = paste(resultsDir, "/", data.model, "/LogSc.", iter, sep = ""))
  
  data = list(XY = XY, predsMRA = predsMRA, predsE = predsE, predsLR = predsLR)
  save(data, file = paste(resultsDir, "/", data.model, "/sim.", iter, sep = ""))
  
  print(RRMSPE)
}


########## plot results ########## 
#for (t in 1:Tmax) {
#  zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm = TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm = TRUE)))))
#  nna.obs = which(!is.na(XY$y[[t]]))
#  plot( locs, XY$x[[t]], type = "l", xlim = c(0, 1), ylim = zrange, col = "red", main = paste("t =", t))
#  points( locs[nna.obs,], XY$y[[t]][nna.obs], pch = 16, col = "black")
#  lines( locs, predsVL[[t]]$state, type = "l", col = "blue")
#}