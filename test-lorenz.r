########## Auxiliary functions ########## 



### This function is used mainly for testing.
### It takes the entire covariance matrix and creates
### a matrix of covariances
getL00 = function(vecchia.approx, covfun, locs){
  Sig.sel = GPvecchia::getMatCov(vecchia.approx, covfun(locs))
  inds = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(vecchia.approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(Sig.sel)))
  vals = GPvecchia::ic0(ptrs, inds, cov.vals)
  Laux = Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  ro = order(vecchia.approx$ord)
  return(Laux)
}


getLtt = function(vecchia.approx, preds){
  orig.order=order(vecchia.approx$ord)
  V = preds$V
  L.tt = (Matrix::solve(Matrix::t(V), sparse=TRUE)[seq(n, 1), ])[orig.order,]
  return(L.tt)
}



######### simulate and plot the data #########
## define the temporal evolution function
evol = function(state, f=0.1){
  N = dim(as.matrix(state))[1]
  
  diags = list(rep(f, 1), rep(-f, N-2), rep(-1, N), rep(f, N-1), rep(-f, 2))
  E = Matrix::bandSparse(N, k=c(-N+1, -2, 0, 1, N-2), diagonals = diags)
  
  if(dim(state)[2]>1) return( state + E %*% (state-f) )
  else as.numeric(state + E %*% as.matrix(state-f))
}


## simulate y given x 
simulate.y = function(x, frac.obs, lik.params){
  
  n.obs = round(n*frac.obs)
  obs.inds = sample(1:n, n.obs, replace = FALSE)
  data.model = lik.params["data.model"]
  # simulate data
  if(data.model=='poisson'){
    y.obs = rpois(n.obs, exp(x[obs.inds]))
  } else if(data.model=='logistic'){
    y.obs = rbinom(n.obs,1,prob = exp(x[obs.inds])/(1+exp(x[obs.inds])))
  } else if(data.model=='gamma'){
    #default_lh_params = list("alpha"=2, "sigma"=sqrt(.1), "beta"=.9, "phi"=1.5)
    #z = rgamma(n.obs, shape = default_lh_params$alpha, rate = default_lh_params$alpha*exp(-y[obs.inds]))
    y.obs = rgamma(n.obs, shape = params["alpha"], rate = lik.params[["alpha"]]*exp(-x[obs.inds]))
  } else if(data.model=='gauss'){
    y.obs = rnorm(n.obs, mean=x[obs.inds], sd=sqrt(lik.params[["me.var"]]))
  } else {
    print('Error: Distribution not implemented yet.')
  }
  y = rep(NA, n)
  y[obs.inds] = y.obs
  return(y)
}




## simulate x
simulate.xy = function(x0, E, Q, frac.obs, lik.params, Tmax, seed=NULL){
  
  if(!is.null(seed)) set.seed(seed)
  n=nrow(x0);
  x = list(); y = list()
  
  x[[1]] = x0
  y[[1]] = simulate.y(x0, frac.obs, lik.params)
  
  Qc = chol(Q)
  for(t in 2:Tmax){
    x[[t]] = E(x[[t-1]]) + t(Qc) %*% matrix(rnorm(n), ncol=1)
    y[[t]] = simulate.y(x[[t]], frac.obs, lik.params)
  }
  
  return(list(x=x, y=y))
}




######### set parameters #########
set.seed(1988)
spatial.dim=2
n=50
m=2
forcing=0.2
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
grid.oneside=seq(0,1,length=round(n))
locs=matrix(grid.oneside, ncol=1)

## set initial state
Q = covfun(locs)
x0 = matrix(rep(0, n), ncol=1); Sig0 = covfun(locs)
XY = simulate.xy(x0, function(x) evol(x, f=forcing), Q, frac.obs, lik.params, Tmax, forcing)

## define Vecchia approximation
vecchia.approx=GPvecchia::vecchia_specify(locs, m, conditioning='mra')#, mra.options=list(r=c(32)))
#vecchia.exact = GPvecchia::vecchia_specify(locs, conditioning='mra', verbose=TRUE, mra.options = list(r=c(n, 0)))#


########## filtering ##########
predsVL = list()
L.00.vl = getL00(vecchia.approx, covfun, locs)
predsVL[[1]] = list(state=x0, L=L.00.vl)

forecast=x0; Fmat=L.00.vl
for(t in 1:Tmax){
  cat(paste("filtering: t=", t, "\n", sep=""))
  obs.aux = as.numeric(XY$y[[t]])# - forecast)
  covmodel = GPvecchia::getMatCov(vecchia.approx, Fmat %*% Matrix::t(Fmat) + Q)
  preds.aux.vl = GPvecchia::calculate_posterior_VL( obs.aux, vecchia.approx, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
  L.tt = getLtt(vecchia.approx, preds.aux.vl)
  mu.tt = matrix(preds.aux.vl$mean, ncol=1)
  predsVL[[t]] = list(state=mu.tt, L=L.tt)
  if(t<Tmax){
    forecast = evol(mu.tt, f=forcing)
    Fmat = evol(L.tt, f=forcing)
  }
}


########## plot results ########## 
for(t in 1:Tmax){
  zrange = range(c(predsVL[[t]][["state"]], unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  #zrange = range(c(predsE[[t]][["state"]], predsVL[[t]][["state"]], unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  #defpar = par(mfrow=c(2, 2), oma=c(0, 0, 2, 0))
  nna.obs = which(!is.na(XY$y[[t]]))
  plot( locs[nna.obs,], XY$y[[t]][nna.obs], pch=16, ylim=c(0.195, 0.205), xlim=c(0, 1), main=paste("t =", t))#,zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
  lines( locs, as.numeric(XY$x[[t]]), type="l", col="red")#, zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
  #lines( locs, predsVL[[t]]$state, type="l", col="blue")#, zlim=zrange, nx=sqrt(n), ny=sqrt(n) )

  par(defpar)
}