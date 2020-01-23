setwd("~/HVLF")
rm(list=ls())
source("aux-functions.r")
library(VEnKF)
library(rootSolve)


######### set parameters #########
set.seed(1988)
spatial.dim=2
n=120
m=2
frac.obs = 0.1
Tmax = 10


## evolution function ##
Force = 10
K = 32
dt = 0.005
M = 40
evolFun = function(X) Lorenz04M2Sim(as.numeric(X), Force, K, dt, M, iter=1, burn=1, order=4)


## covariance function
sig2=0.1; range=.15; smooth=1.5; 
covparms =c(sig2,range,smooth)
covfun <- function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)


## likelihood settings
me.var=0.1;
data.model = "gauss"
lik.params = list(data.model = data.model, me.var=me.var)


## generate grid of pred.locs
grid.oneside=seq(0,1,length=round(n))
locs=matrix(grid.oneside, ncol=1)

## set initial state
x0 = matrix(rnorm(n), ncol=1)#getX0(n, Force, K, dt)
Q = covfun(locs)

XY = simulate.xy(x0, evolFun, Q, frac.obs, lik.params, Tmax)

## define Vecchia approximation
vecchia.approx=GPvecchia::vecchia_specify(locs, m, conditioning='mra')#, mra.options=list(r=c(32)))
#vecchia.exact = GPvecchia::vecchia_specify(locs, conditioning='mra', verbose=TRUE, mra.options = list(r=c(n, 0)))#





########## filtering ##########
predsVL = list()
L.00.vl = getL00(vecchia.approx, covfun, locs)
predsVL[[1]] = list(state=x0, L=L.00.vl)

mu.tt=x0; L.tt=L.00.vl

for(t in 2:Tmax){
  cat(paste("filtering: t=", t, "\n", sep=""))
  obs.aux = as.numeric(XY$y[[t]])
  
  cat("\tcalculate gradient...")
  Et = Matrix::Matrix(gradient(evolFun, mu.tt, centered=TRUE))
  cat(" done.\n")
  Fmat = Et %*% L.tt
  covmodel = GPvecchia::getMatCov(vecchia.approx, Fmat %*% Matrix::t(Fmat) + Q)
  mu.tt1 = evolFun(mu.tt)
  
  preds.aux.vl = GPvecchia::calculate_posterior_VL( obs.aux, vecchia.approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
  
  L.tt = getLtt(vecchia.approx, preds.aux.vl)
  mu.tt = matrix(preds.aux.vl$mean, ncol=1)
  
  predsVL[[t]] = list(state=mu.tt, L=L.tt)
}


########## plot results ########## 
for(t in 1:Tmax){
  zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  nna.obs = which(!is.na(XY$y[[t]]))
  plot( locs, XY$x[[t]], type="l", xlim=c(0, 1), ylim=zrange, col="red", main=paste("t =", t))
  points( locs[nna.obs,], XY$y[[t]][nna.obs], pch=16, col="black")
  lines( locs, predsVL[[t]]$state, type="l", col="blue")
}