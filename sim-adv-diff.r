## define the temporal evolution function
evolAdvDiff = function(state, adv=0, diff=0){
  
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



filter = function(approx.name, XY){
  
  approx = approximations[[approx.name]]
  
  preds = list()
  L.tt = getL00(approx, Sig0, locs)
  mu.tt = x0
  preds[[1]] = list(state=x0, L = L.tt)
  
  for(t in 2:Tmax){
    obs.aux = as.numeric(XY$y[[t]])
    
    forecast = evolFun(mu.tt)
    Fmat = evolFun(L.tt) 
    covmodel = GPvecchia::getMatCov(approx, Fmat %*% Matrix::t(Fmat) + Q)
    
    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = forecast, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol=1)
    preds[[t]] = list(state=mu.tt, W = preds.aux$W)
  }
  return( preds )
}



setwd("~/HVLF")
source('aux-functions.r')
source('score-adv-diff.r')
resultsDir="simulations-linear"
library(doParallel)
registerDoParallel(cores=5)

######### set parameters #########
set.seed(1988)
spatial.dim=2
n=20**2
m=50
frac.obs = 0.3
Tmax = 3
diffusion = 0.00004
advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv=advection, diff=diffusion)
max.iter = 2

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
save(locs, file=paste(resultsDir, "/locs", sep=""))


## set initial state
Q = covfun(locs)
Sig0 = (1/sig2)*covfun(locs)
x0 = t(chol(Sig0)) %*% matrix(rnorm(n), ncol=1); 


## define Vecchia approximation
mra=GPvecchia::vecchia_specify(locs, m, conditioning='mra')
low.rank=GPvecchia::vecchia_specify(locs, m, conditioning='firstm')
exact = GPvecchia::vecchia_specify(locs, n-1, conditioning='firstm')
approximations = list(mra=mra, low.rank = low.rank, exact=exact)


RRMSPE = list(); LogSc = list()
#foreach( iter=1:max.iter) %dopar% {
for( iter in 1:max.iter) {  

  XY = simulate.xy(x0, evolFun, Q, frac.obs, lik.params, Tmax)
  
  print(paste("iteration: ", iter, ", MRA", sep=""))
  predsMRA = filter('mra', XY)
  print(paste("iteration: ", iter, ", LR", sep=""))
  predsLR  = filter('low.rank', XY)
  print(paste("iteration: ", iter, ", exact", sep=""))
  predsE = filter('exact', XY)
 
  RRMSPE[[iter]] = calculateRRMSPE(predsMRA, predsLR, predsE, XY$x)
  LogSc[[iter]] = calculateLSs(predsMRA, predsLR, predsE, XY$x)
  
  data = list(RRMSPE=RRMSPE, LogScore=LogSc)
  save(data, file=paste(resultsDir, "/", data.model, ".scores.", iter, sep=""))
  
}
print(data)