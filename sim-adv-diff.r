## define the temporal evolution function
evol2 = function(state, adv=0, diff=0){
  # we assume that there is the same number of grid points
  # along each dimension

  N = dim(as.matrix(state))[1]
  Ny = Nx = sqrt(N)
  
  dx = dy = 1/Nx
  d = diff/(dx**2)
  
  c1 = 1 + 2*(d + d) - adv*(1/dx + 1/dy)
  c2 = - d + adv*(1/dy)
  c3 = - d
  
  diags = list(rep(c2, N-Nx), rep(c2, N-1), rep(c1, N), rep(c3, N-1), rep(c3,N-Nx) )
  E = Matrix::bandSparse(N, k=c(-Nx, -1, 0, 1, Nx), diag=diags)
  
  if(class(state)=='matrix' || class(state)=='dgCMatrix') return( E %*% state )
  else as.numeric(E %*% as.matrix(state))
}




setwd("~/HVLF")
source('aux-functions.r')
resultsDir="simulations-linear"
library(doParallel)
registerDoParallel(cores=5)
######### set parameters #########
set.seed(1988)
spatial.dim=2
n=34**2
m=50
frac.obs = 0.3
Tmax = 10
diffusion = 0.00004
advection = 0.01
evolFun = function(X) evol2(X, adv=advection, diff=diffusion)
max.iter = 5

## covariance parameters
sig2=0.5; range=.15; smooth=0.5; 
covparms =c(sig2,range,smooth)
covfun <- function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)

## likelihood settings
me.var=0.1;
data.model = "poisson"
lik.params = list(data.model = data.model, me.var=me.var)







## generate grid of pred.locs
grid.oneside=seq(0,1,length=round(sqrt(n)))
locs=as.matrix(expand.grid(grid.oneside,grid.oneside)) 

## set initial state
Q = covfun(locs)
Sig0 = (1/sig2)*covfun(locs)
x0 = t(chol(Sig0)) %*% matrix(rnorm(n), ncol=1); 


## define Vecchia approximation
mra=GPvecchia::vecchia_specify(locs, m, conditioning='mra')#, verbose=TRUE)#, mra.options=list(r=c(32)))
low.rank=GPvecchia::vecchia_specify(locs, m, conditioning='firstm')#, verbose=TRUE)#, mra.options=list(r=c(32)))
exact = GPvecchia::vecchia_specify(locs, conditioning='mra', verbose=TRUE, mra.options = list(r=c(n, 0)))#
vecchia = GPvecchia::vecchia_specify(locs, m, verbose=TRUE)

data=list()

foreach( iter=1:max.iter) %dopar% {
#for( iter in 1:max.iter) {  

  XY = simulate.xy(x0, evolFun, Q, frac.obs, lik.params, Tmax)
  
  ########## filtering ##########
  print(paste("iteration: ", iter, ", MRA", sep=""))
  predsMRA = list()
  L.tt = getL00(mra, Sig0, locs)
  mu.tt = x0
  predsMRA[[1]] = list(state=x0, L=L.tt)
  
  for(t in 2:Tmax){
    #cat(paste("filtering: t=", t, "\n", sep=""))
    forecast = evolFun(mu.tt)
    Fmat = evolFun(L.tt) 
    covmodel = GPvecchia::getMatCov(mra, Fmat %*% Matrix::t(Fmat) + Q)
    obs.aux = as.numeric(XY$y[[t]])
    
    preds.aux.vl = GPvecchia::calculate_posterior_VL( obs.aux, mra, prior_mean=forecast, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(mra, preds.aux.vl)
    mu.tt = matrix(preds.aux.vl$mean, ncol=1)
    predsMRA[[t]] = list(state=mu.tt, L = L.tt)
    
  }
  
  print(paste("iteration: ", iter, ", low-rank", sep=""))
  predsLR = list()
  L.tt = getL00(low.rank, Sig0, locs)
  mu.tt = x0
  predsLR[[1]] = list(state=x0, L = L.tt)

  for(t in 2:Tmax){
    #cat(paste("filtering: t=", t, "\n", sep=""))
    forecast = evolFun(mu.tt)
    Fmat = evolFun(L.tt)
    covmodel = GPvecchia::getMatCov(low.rank, Fmat %*% Matrix::t(Fmat) + Q)
    obs.aux = as.numeric(XY$y[[t]])

    preds.aux.lr = GPvecchia::calculate_posterior_VL( obs.aux, low.rank, prior_mean=forecast, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(low.rank, preds.aux.lr)
    mu.tt = matrix(preds.aux.lr$mean, ncol=1)
    predsLR[[t]] = list(state=mu.tt, L = L.tt)
  }
  

  print(paste("iteration: ", iter, ", exact", sep=""))
  predsE = list()
  L.tt = getL00(exact, Sig0, locs)
  mu.tt = x0
  predsE[[1]] = list(state=x0, L = L.tt)
  
  for(t in 2:Tmax){
    #cat(paste("filtering: t=", t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    
    forecast = evolFun(mu.tt)
    Fmat = evolFun(L.tt) 
    covmodel = GPvecchia::getMatCov(exact, Fmat %*% Matrix::t(Fmat) + Q)
    
    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, exact, prior_mean = forecast, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(exact, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol=1)
    predsE[[t]] = list(state=mu.tt, L = L.tt)
  }   
  
  
  data = list(XY=XY, predsMRA=predsMRA, predsE=predsE, predsLR=predsLR)
  save(data, file=paste(resultsDir, "/", data.model, ".", iter, sep=""))
  
}
#save(data, file=paste(data.model, "data", sep="."))
