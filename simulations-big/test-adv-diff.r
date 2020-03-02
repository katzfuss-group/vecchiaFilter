setwd("~/HVLF")
rm(list=ls())
source("aux-functions.r")


######### set parameters #########
set.seed(1996)
n=34**2
m=3
diffusion = 0.0000
advection = 0.01
frac.obs = 0.1
Tmax = 2


## covariance parameters
sig2=0.5; range=.15; smooth=0.5; 
covparms =c(sig2,range,smooth)
covfun <- function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)

## likelihood settings
me.var=0.1;
args = commandArgs(trailingOnly=TRUE)
if(length(args)==1 && !(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))){
    stop("One of the models has to be passed as argument")
} else if (length(args)==0) {
    data.model = "gamma"
} else {
    data.model = args[1]
}
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)




## generate grid of pred.locs
grid=seq(0,1,length=sqrt(n))
locs=as.matrix(expand.grid(grid,grid)) 

## set initial state
Q = covfun(locs)
#x0 = diffAdvVec2d(sqrt(n))
#fields::quilt.plot(locs, x0, nx=sqrt(n), ny=sqrt(n))
x0 = t(chol(Q)) %*% matrix(rnorm(n), ncol=1);
evolFun = function(x) evol(x, diff=diffusion, adv=advection)
XY = simulate.xy(x0, evolFun, Q, frac.obs, lik.params, Tmax)

## define Vecchia approximation
vecchia.approx=GPvecchia::vecchia_specify(locs, m, conditioning='mra')#, mra.options=list(r=c(32)))
#vecchia.exact = GPvecchia::vecchia_specify(locs, conditioning='mra', verbose=TRUE, mra.options = list(r=c(n, 0)))#


########## filtering ##########
predsVL = list()
L.tt = getL00(vecchia.approx, covfun, locs)
mu.tt = x0
predsVL[[1]] = list(state=x0, L=L.tt)


for(t in 2:Tmax){

  cat(paste("filtering: t=", t, "\n", sep=""))
  obs.aux = as.numeric(XY$y[[t]])

  E = evol(diag(n))
  Fmat = E %*% L.tt
  covmodel = GPvecchia::getMatCov(vecchia.approx, Fmat %*% Matrix::t(Fmat) + Q)
  mu.tt1 = E %*% mu.tt

  preds.aux.vl = GPvecchia::calculate_posterior_VL( obs.aux, vecchia.approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)

  L.tt = getLtt(vecchia.approx, preds.aux.vl)
  mu.tt = matrix(preds.aux.vl$mean, ncol=1)

  predsVL[[t]] = list(state=mu.tt, L=L.tt)

}


########## plot results ########## 
m = M = 0
for( t in 1:Tmax ){
  zrange = range(c(predsVL[[t]][["state"]], unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  #zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  m = min(m, zrange[1])
  M = max(M, zrange[2])
}

for(t in 1:Tmax){
  pdf(paste("filtering-", t, ".pdf", sep=""))
  defpar = par(mfrow=c(1, 3), oma=c(0, 0, 2, 0))
  nna.obs = which(!is.na(XY$y[[t]]))
  fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], zlim=c(m, M), nx=sqrt(n), ny=sqrt(n), main="obs" )
  fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim=c(m, M), nx=sqrt(n), ny=sqrt(n), main="truth" )
  fields::quilt.plot( locs, predsVL[[t]]$state, zlim=zrange, nx=sqrt(n), ny=sqrt(n), main="prediction" )
  par(defpar)
  dev.off()
}