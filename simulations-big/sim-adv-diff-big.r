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

  if (class(state) == 'matrix' || methods::is(state, 'sparseMatrix')){
    return( E %*% state )
  } else {
    return( as.numeric(E %*% as.matrix(state)) )
  }
    
  #if(class(state)=='matrix' || class(state)=='dgCMatrix') return( E %*% state )
  #else as.numeric(E %*% as.matrix(state))
}



filter.old = function(approx.name, XY){
  
  approx = approximations[[approx.name]]
  
  preds = list()
  L.tt = getL00(approx, Sig0, locs)
  mu.tt = x0
  preds[[1]] = list(state=x0, L = L.tt)

  for(t in 2:Tmax){

      #cat(paste("\tfiltering: t=",t, "\n", sep=""))
      obs.aux = as.numeric(XY$y[[t]])
      
      forecast = evolFun(mu.tt)
      Fmat = evolFun(L.tt)
      M = as.matrix(Fmat %*% Matrix::t(Fmat) + Q)
      covmodel = GPvecchia::getMatCov(approx, M)
      preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                                  likelihood_model = data.model, covmodel = covmodel,
                                                  covparms = covparms, likparms = lik.params, return_all = TRUE)
      L.tt = getLtt(approx, preds.aux)
      mu.tt = matrix(preds.aux$mean, ncol=1)
      preds[[t]] = list(state=mu.tt, W = preds.aux$W)
  }
  return( preds )
}



########## filtering ##########
filter = function(approx.name, XY){
  
  approx = approximations[[approx.name]]
  preds = list()
  cat(paste("filtering: t=1\n", sep = ""))
  cat("\tCalculate covariance elements from function: ")
  t0 = proc.time()
  covmodel = getMatCov(approx, covfun.d)
  t1 = proc.time()
  cat(paste((t1 - t0)[3], "\n"))
  
  mu.tt1 = rep(0, n)
  obs.aux = as.numeric(XY$y[[1]])
  
  cat("\tcalculate posterior: ")
  t0 = proc.time()
  preds.aux.vl = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
  t1 = proc.time()
  cat(paste((t1 - t0)[3], "\n"))
  
  cat("\tstore results\n")
  cat("\t\textract L.tt matrix\n")
  L.tt = getLtt(approx, preds.aux.vl)
  cat("\t\textract mean vector\n")
  mu.tt = matrix(preds.aux.vl$mean, ncol = 1)
  preds[[1]] = list(state = mu.tt, L = L.tt)
  
  if ( Tmax > 1 ) {
    
    for (t in 2:Tmax) {
      
      cat(paste("filtering: t=", t, "\n", sep = ""))
      obs.aux = as.numeric(XY$y[[t]])
      
      cat("\tevolve the L.tt matrix:\n")
      t0 = proc.time()
      cat("\t\tbuild the evolution matrix\n")
      E = evolFun(Matrix::Diagonal(n))
      cat("\t\tmultpily E by L\n")
      Fmat = E %*% L.tt
      t1 = proc.time()
      #cat(paste((t1 - t0)[3], "\n"))
      
      
      cat("\tCalculate covariance elements from factor: ")
      t0 = proc.time()
      
      M1 = getMatCov(approx, Matrix::t(Fmat), factor = TRUE)
      t1 = proc.time()
      cat(paste((t1 - t0)[3], "\n"))
      cat("\t... from function: ")
      t0 = proc.time()
      M2 = getMatCov(approx, covfun.d)
      t1 = proc.time()
      cat(paste((t1 - t0)[3], "\n"))
      covmodel = M1 + M2
      
      mu.tt1 = E %*% mu.tt
      
      cat("\tcalculate posterior: ")
      t0 = proc.time()
      preds.aux.vl = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
      t1 = proc.time()
      cat(paste((t1 - t0)[3], "\n"))
      
      cat("\tstore results\n")
      cat("\t\textract L.tt matrix\n")
      L.tt = getLtt(approx, preds.aux.vl)
      cat("\t\textract mean vector\n")
      mu.tt = matrix(preds.aux.vl$mean, ncol = 1)

      preds[[t]] = list(state = mu.tt, L = L.tt)
      
    }
    
  }
  return( preds )
}






setwd("~/HVLF")
source('aux-functions.r')
source('getMatCov.r')
Rcpp::sourceCpp('src/getMatCovFromFactor.cpp')
source('scores.r')
resultsDir = "simulations-linear"
library(doParallel)
library(Matrix)
registerDoParallel(cores = 5)

######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 100**2
m = 50
frac.obs = 0.1
Tmax = 20
diffusion = 0.0000001
advection = 0.001
#diffusion = 0.00004
#advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
max.iter = 1

## covariance parameters
sig2 = 1.0; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.25;
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  if (!(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))) {
    stop("One of the models has to be passed as argument")
  } else {
    data.model = args[1]
  }
} else {
  data.model = "gauss"  
}
lik.params = list(data.model = data.model, me.var = me.var)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 
save(locs, file = paste(resultsDir, "/locs", sep = ""))


## set initial state
Q = covfun(locs)
Sig0 = (1/sig2)*covfun(locs)
x0 = t(chol(Sig0)) %*% Matrix::Matrix(rnorm(n), ncol = 1); 


## define Vecchia approximation
#exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm', verbose = TRUE)
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
approximations = list(mra = mra, low.rank = low.rank)#, exact = exact)


RRMSPE = list(); LogSc = list()
#foreach( iter=1:max.iter) %dopar% {
for (iter in 1:max.iter) {  

    XY = simulate.xy(x0, evolFun, Q, frac.obs, lik.params, Tmax)

    #cat(paste("iteration: ", iter, ", exact", "\n", sep = ""))
    #predsE = filter('exact', XY)
    cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
    predsMRA = filter('mra', XY)
    cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
    predsLR  = filter('low.rank', XY)
    
    RMSPE = calculateRMSPE(predsMRA, predsLR, XY$x)
  
    write.csv(RRMSPE, file = paste(resultsDir, "/", data.model, "/RRMSPE.", iter, sep = ""))

    #data = list(XY = XY, predsMRA = predsMRA, predsE = predsE, predsLR = predsLR)
    #save(data, file = paste(resultsDir, "/", data.model, "/sim.", iter, sep = ""))
    
    print(RMSPE)
}
