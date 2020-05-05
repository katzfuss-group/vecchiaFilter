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
  n = nrow(vecchia.approx$locsord)
  orig.order=order(vecchia.approx$ord)
  V = preds$V
  L.tt = (Matrix::solve(Matrix::t(V), sparse=TRUE)[seq(n, 1), ])[orig.order,]
  return(L.tt)
}



getConfInt = function(preds, alpha){
  
  D = Matrix::diag(Matrix::solve(preds$W))
  mu = preds$state
  z = qnorm(1 - alpha/2)
  ub = mu + z*sqrt(D)
  lb = mu - z*sqrt(D)
  
  return(list(ub = ub, lb = lb))
}


getConfInt = function(preds, alpha){
  
  D = Matrix::diag(Matrix::solve(preds$W))
  mu = preds$state
  z = qnorm(1 - alpha/2)
  #browser()
  ub = mu + z*sqrt(D)
  lb = mu - z*sqrt(D)
  
  return(list(ub = ub, lb = lb))
}





######### simulate and plot the data #########
## define the temporal evolution function
evolAdvDiff = function(state, adv=0, diff=0){
  # we assume that there is the same number of grid points
  # along each dimension
  
  N = dim(state)[1]
  Ny = Nx = sqrt(N)
  
  dx = dy = 1/Nx
  d = diff/(dx**2)
  
  c1 = 1 + 2*(d + d) - adv*(1/dx + 1/dy)
  c2 = - d + adv*(1/dy)
  c3 = - d
  
  diags = list(rep(c2, N-Nx), rep(c2, N-1), rep(c1, N), rep(c3, N-1), rep(c3,N-Nx) )
  E = Matrix::bandSparse(N, k=c(-Nx, -1, 0, 1, Nx), diag=diags)
  
  if (class(state) == 'matrix' || methods::is(state, 'sparseMatrix')){
    return( E %*% state )
  } else {
    return( as.numeric(E %*% as.matrix(state)) )
  }
}

  


getX0 = function(N, Force, K, dt, dir = '~/HVLF/simulations-lorenz/'){
  
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




diffAdvVec2d = function(nx, ny=nx, height=1, rnge=4){
  v = matrix(rep(0, nx*ny), ncol=ny)
  if((nx %% 2)==0) mid_x = nx/2 else mid_x = nx/2+1
  if((ny %% 2)==0) mid_y = ny/2 else mid_y = ny/2+1
  v[round(mid_x-rnge/2):round(mid_x+rnge/2), round(mid_y-rnge/2):round(mid_y+rnge/2)] = height
  
  return(matrix(as.numeric(v), ncol=1))
}




## simulate y given x 
simulate.y = function(x, frac.obs, lik.params){
  
  n = nrow(x)
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

  } else if(data.model == 'gauss'){
    y.obs = rnorm(n.obs, mean = x[obs.inds], sd=lik.params[["sigma"]])


  } else {
    print('Error: Distribution not implemented yet.')
  }
  y = rep(NA, n)
  y[obs.inds] = y.obs
  return(y)
}



## simulate x
simulate.xy = function(x0, E, Q, frac.obs, lik.params, Tmax, seed=NULL, sig2=1, smooth = 0.5, range = 1, locs = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  n = nrow(x0);
  x = list(); y = list()
  x[[1]] = x0
  y[[1]] = simulate.y(x0, frac.obs, lik.params)
  
  if (Tmax > 1) { 
    
    if (!is.null(Q) && any(Q)) {
      Qc = Matrix::chol(Q)
    } 
    
    for (t in 2:Tmax) {
      if (sig2 > 0 || (!is.null(Q) && sum(abs(Q))>0)) {
        if (!is.null(Q)) {
          w =  t(Qc) %*% matrix(rnorm(n), ncol = 1)
        } else {
          w = matrix(sig2*RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range),
                                            x = locs[,1], y = locs[,2], spConform = FALSE), ncol=1)
        } 
      } else {
        w = matrix(rep(0, n), ncol=1)
      }
      
      x[[t]] = E(x[[t - 1]]) + w
      y[[t]] = simulate.y(x[[t]], frac.obs, lik.params)
    } 
  }
  
  return(list(x = x, y = y))

}
