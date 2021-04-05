########## Auxiliary functions ########## 



### This function is used mainly for testing.
### It takes the entire covariance matrix and creates
### a matrix of covariances
#getL00 = function(vecchia.approx, covfun, locs){
getL00 = function(vecchia.approx, covfun){
  Sig.sel = GPvecchia::getMatCov(vecchia.approx, covfun(vecchia.approx$locs))
  inds = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(vecchia.approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(Sig.sel)))
  vals = GPvecchia::ic0(ptrs, inds, cov.vals)
  Laux = Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  ro = order(vecchia.approx$ord)
  return(Laux[ro, ro])
  #return(L)
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





######### simulate and plot the data #########
## define the temporal evolution function
evolAdvDiff = function(state, adv=0, diff=0){
  # we assume that there is the same number of grid points
  # along each dimension
    if(is.null(dim(state))) {
        N = length(state)
    } else {
        N = dim(state)[1]
    }

    Ny = Nx = sqrt(N)
  
    dx = dy = 1/Nx
    D = diff/(dx**2)


    c1 = 1 + 2*(D + D) - adv*(1/dx + 1/dy)
    c2 = - D + adv*(1/dy)
    c3 = - D

    ## Let's consider the following layout on the grid
    ## . . . . .
    ## . . c . .
    ## . a x b .
    ## . . d . .
    ## . . . . .
    ## where the values "x" are given by d3, "d" - by d1,
    ## "a" - by d2, "b" - by d4 and "c" by d5.
    d = rep(c2, N-Nx)
    a = rep(c2, N-1)
    x = rep(c1, N)
    b = rep(c3, N-1)
    c = rep(c3,N-Nx)

    # Now notice that for k = 1:n, in each kn-th row there should
    # be no b and for every kn+1-st row there should be no a
    inds.b = Nx*(1:(Ny-1))
    inds.a = Nx*(1:(Ny-1))
    b[inds.b] = 0
    a[inds.a] = 0
    
    diags = list(d, a, x, b, c)
    E = Matrix::bandSparse(N, k=c(-Nx, -1, 0, 1, Nx), diag=diags)

    #browser()
    
    if (class(state) == 'matrix' || methods::is(state, 'sparseMatrix')){
        return( E %*% state )
    } else {
        return( as.numeric(E %*% as.matrix(state)) )
    }
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
simulate.xy = function(init_mean, sig02, range0, nu0, E, Q, frac.obs, lik.params, Tmax, seed=NULL, sig2=1, smooth = 0.5, range = 1, locs = NULL){

    if (!is.null(seed)) set.seed(seed)
    n = nrow(init_mean);
    x = list(); y = list()

    x0 = matrix(RandomFields::RFsimulate(model = RandomFields::RMwhittle(nu = nu0, scale = range0, var = sig02),
                                         x = locs[,1], y = locs[,2], spConform = FALSE), ncol=1) + init_mean
        
    
    if (!is.null(Q) && any(Q)) {
        Qc = Matrix::chol(Q)
    } 
    
    for (t in 1:Tmax) {
        cat(sprintf("simulating at time t=%d\n", t))
        if (sig2 > 0 || (!is.null(Q) && sum(abs(Q))>0)) {
            if (!is.null(Q)) {
                w =  t(Qc) %*% matrix(rnorm(n), ncol = 1)
            } else {
                w = matrix(RandomFields::RFsimulate(model = RandomFields::RMwhittle(nu = smooth, scale = range, var = sig2),
                                                    x = locs[,1], y = locs[,2], spConform = FALSE), ncol=1)
            } 
        } else {
            w = matrix(rep(0, n), ncol=1)
        }

        if (t==1) {
            x[[t]] = E(x0) + w
        } else {
            x[[t]] = E(x[[t - 1]]) + w
        }
        y[[t]] = simulate.y(x[[t]], frac.obs, lik.params)
    } 
    
    return(list(x = x, y = y, x0 = x0))
    
}
