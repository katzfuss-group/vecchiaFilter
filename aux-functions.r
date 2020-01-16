########## Auxiliary functions ########## 


RMSPE = function(pred, truth) sqrt(mean((pred - truth**2)))



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
evol = function(state, adv=0, diff=0){
  Nx = dim(as.matrix(state))[1]
  
  c1 = 1 - 2*diff*(Nx**2)
  c2 = 0.5*adv*Nx + diff*(Nx**2)
  c3 = -0.5*adv*Nx + diff*(Nx**2)
  
  diags = list(rep(c2, 1), rep(c3, Nx-1), rep(c1, Nx), rep(c2, Nx-1), rep(c3,1) )
  E = Matrix::bandSparse(Nx, k=c(-(Nx-1), -1, 0, 1, Nx-1), diag=diags)
  
  if(dim(state)[2]>1) return( E %*% state )
  else as.numeric(E %*% as.matrix(state))
}



getX0 = function(N, Force, K, dt, dir = '~/HVLF/models/'){
  
  fileName = paste("init_Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep="_")
  filePath = paste(dir, fileName, sep="")
  
  generateInit = !file.exists(filePath)

  if( generateInit ){
    message("Initial state for these parameters has not yet been generated. Generating now.")
    X0 = rnorm(N)  
    X1 = Lorenz04M2Sim(X0, Force, K, dt, M, iter=100, burn=0, order=1)
    x0 = X1[,100]
    write(x0, file=filePath)  
  } else {
    x0 = scan(filePath, quiet = TRUE)  
  }
  if(class(x0)=='numeric'){
    x0 = matrix(x0, ncol=1)
  }
  return(x0)

}



lorenz98 = function(state, f){
  N = dim(as.matrix(state))[1]
  
  diags = list(rep(f, 1), rep(-f, N-2), rep(-1, N), rep(f, N-1), rep(-f, 2))
  EL = Matrix::bandSparse(N, k=c(-N-1), -2, 0, 1, N-2)
  
  if(dim(state)[2]>1) return( E %*% state )
  else as.numeric(E %*% as.matrix(state))
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
    y.obs = rgamma(n.obs, shape = lik.params[["alpha"]], rate = lik.params[["alpha"]]*exp(-x[obs.inds]))
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
  
  if(Tmax>1){
    Qc = chol(Q)
    
    for(t in 2:Tmax){
      x[[t]] = E(x[[t-1]]) + t(Qc) %*% matrix(rnorm(n), ncol=1)
      y[[t]] = simulate.y(x[[t]], frac.obs, lik.params)
    } 
  }
  
  return(list(x=x, y=y))
}