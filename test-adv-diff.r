setwd("~/HVLF")
rm(list=ls())
source("aux-functions.r")


######### set parameters #########
set.seed(1988)
spatial.dim=2
n=34**2
m=50
diffusion = 0.00004
advection = 0.01
frac.obs = 0.3
Tmax = 10


## covariance parameters
sig2=0.000001; range=.15; smooth=1.5; 
covparms =c(sig2,range,smooth)
covfun <- function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)

## likelihood settings
me.var=1e-8;
data.model = "gauss"
lik.params = list(data.model = data.model, me.var=me.var)




## generate grid of pred.locs
grid=seq(0,1,length=sqrt(n))
locs=as.matrix(expand.grid(grid,grid)) 

## set initial state
Q = covfun(locs)
x0 = matrix(rep(0, n), ncol=1); Sig0 = covfun(locs)
XY = simulate.xy(x0, function(x) evol(x, diff=diffusion, adv=advection), Q, frac.obs, lik.params, Tmax)

## define Vecchia approximation
vecchia.approx=GPvecchia::vecchia_specify(locs, m, conditioning='mra')#, mra.options=list(r=c(32)))
#vecchia.exact = GPvecchia::vecchia_specify(locs, conditioning='mra', verbose=TRUE, mra.options = list(r=c(n, 0)))#


########## filtering ##########
predsVL = list()
L.tt = getL00(vecchia.approx, covfun, locs)
mu.tt = x0
predsVL[[1]] = list(state=x0, L=L.tt)


for(t in 1:Tmax){

  cat(paste("filtering: t=", t, "\n", sep=""))
  obs.aux = as.numeric(XY$y[[t]])
  
  Fmat = evol(L.tt)
  covmodel = GPvecchia::getMatCov(vecchia.approx, Fmat %*% Matrix::t(Fmat) + Q)
  mu.tt1 = evol(mu.tt)
  
  preds.aux.vl = GPvecchia::calculate_posterior_VL( obs.aux, vecchia.approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
  
  L.tt = getLtt(vecchia.approx, preds.aux.vl)
  mu.tt = matrix(preds.aux.vl$mean, ncol=1)
  
  predsVL[[t]] = list(state=mu.tt, L=L.tt)

}


########## plot results ########## 
for(t in 1:Tmax){
  zrange = range(c(predsVL[[t]][["state"]], unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  #zrange = range(c(predsE[[t]][["state"]], predsVL[[t]][["state"]], unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  defpar = par(mfrow=c(2, 2), oma=c(0, 0, 2, 0))
  nna.obs = which(!is.na(XY$y[[t]]))
  fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
  fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
  #fields::quilt.plot( locs, predsVL[[t]]$state, zlim=zrange, nx=sqrt(n), ny=sqrt(n) )

  par(defpar)
}