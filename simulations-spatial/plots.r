library(GPvecchia)
rm(list=ls())

set.seed(1988)
nx = 7; ny = 7
n=nx*ny
m=10

# covariance function
sig2=1; range=.15; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*GPvecchia::MaternFun(fields::rdist(locs),covparms)
frac.obs = 0.5


# generate grid
if(ny==1){
  locs=matrix(runif(n),ncol=1)
} else {
  grid.x=seq(0,1,length=nx)
  grid.y=seq(0,1,length=ny)
  locs=as.matrix(expand.grid(grid.x,grid.y)) # grid of pred.locs
}

# simulate latent process
Om0 <- covfun(locs)
y=as.numeric(t(chol(Om0))%*%rnorm(n))


RMSPE = function(truth, prediction){
  n = length(truth$mean)
  sqrt(sum((truth$mean - prediction$mean)**2)/n)
}


# calculate MSE of VL vs just Laplace
compare.preds = function(data.model){
  print(data.model)
  
  # simulate data
  if(data.model=='poisson'){
    z = rpois(n, exp(y))
  } else if(data.model=='logistic'){
    z = rbinom(n,1,prob = exp(y)/(1+exp(y)))
  } else if(data.model=='gamma'){
    default_lh_params = list("alpha"=2, "sigma"=sqrt(.1), "beta"=.9, "phi"=1.5)
    z = rgamma(n, shape = default_lh_params$alpha, rate = default_lh_params$alpha*exp(-y))
  } else if(data.model=='gauss'){
    me.std=0.5
    z = rnorm(n, mean=y, sd=me.std)
  } else {
    print('Error: Distribution not implemented yet.')
  }
  obs.inds = sample(1:n, round(n*frac.obs), replace = FALSE)
  z[-obs.inds] = NA
  
  
  ######  specify Vecchia approximation   #######
  vecchia.approx = GPvecchia::vecchia_specify(locs, m, 'maxmin', conditioning='mra')
  posterior.vecchia = GPvecchia::calculate_posterior_VL(z, vecchia.approx, likelihood_model=data.model, covparms = covparms)
  
  exact = GPvecchia::vecchia_specify(locs, n-1)
  posterior.exact = GPvecchia::calculate_posterior_VL(z, exact, likelihood_model=data.model, covparms = covparms)
  
  low.rank = GPvecchia::vecchia_specify(locs, m, 'maxmin', conditioning='firstm')
  posterior.lr = GPvecchia::calculate_posterior_VL(z, low.rank, likelihood_model=data.model, covparms = covparms)
  
  output = list(vecchia = RMSPE(y, posterior.vecchia),
                lowrank = RMSPE(y, posterior.lr),
                exact   = RMSPE(y, posterior.exact))
  
  return(output)
}

rmspe = compare.preds('poisson')