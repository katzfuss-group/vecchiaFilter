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
resultsDir = "simulations-big"
library(doParallel)
library(Matrix)
registerDoParallel(cores = 25)

######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 300**2
m = 50
frac.obs = 0.1
Tmax = 20
diffusion = 0.0000001
advection = 0.001
#diffusion = 0.00004
#advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
max.iter = 10

## covariance parameters
sig2 = 0.2; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.05;
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
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 
save(locs, file = paste(resultsDir, "/locs", sep = ""))


## set initial state
x0 = matrix(0.5*RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range),
                                          x = locs[,1], y = locs[,2], spConform = FALSE), ncol = 1)


## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
approximations = list(mra = mra, low.rank = low.rank)


RMSPE = list(); LogSc = list()
foreach( iter=1:max.iter) %dopar% {
#for (iter in 1:max.iter) {  

    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
    
    cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
    predsMRA = filter('mra', XY)
    cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
    predsLR  = filter('low.rank', XY)
    
    RMSPE = calculateRMSPE(predsMRA, predsLR, XY$x)

    write.csv(RMSPE, file = paste(resultsDir, "/", data.model, "/RMSPE.", iter, sep = ""))
    
    print(RMSPE)

    if(iter==1){
        m = M = 0
        for (t in 1:Tmax) {
            zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE)))))
            m = min(m, zrange[1])
            M = max(M, zrange[2])
        }
        for (t in 1:Tmax) {
            
            if(t<10){
                number = paste("0", t, sep="")  
            } else {
                number = t
            }
            pdf(paste(resultsDir, "/", data.model, "/", number, ".pdf",sep=""), width=8, height=8)
            defpar = par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
            nna.obs = which(!is.na(XY$y[[t]]))
            fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], nx = sqrt(n), ny = sqrt(n), main = "obs" )
            fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim = c(m, M), nx = sqrt(n), ny = sqrt(n), main = "truth" )
            fields::quilt.plot( locs, predsMRA[[t]]$state, zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = "prediction MRA" )
            fields::quilt.plot( locs, predsLR[[t]]$state, zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = "prediction LR" )
            par(defpar)
            dev.off()
        }
    }

    
}


