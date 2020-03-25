rm(list = ls())
source("~/HVLF/aux-functions.r")
source("~/HVLF/getMatCov.r")
source("~/HVLF/scores.r")
resultsDir = "~/HVLF/simulations-big"
library(doParallel)
registerDoParallel(cores=25)


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
      M1 = GPvecchia::getMatCov(approx, Matrix::t(Fmat), factor = TRUE)
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



######### set parameters #########
set.seed(1996)
n = 50**2
m = 10
diffusion = 0.0000001
advection = 0.001
#diffusion = 0.00004
#advection = 0.01
frac.obs = 0.1
Tmax = 2
max.iter = 2

## covariance parameters
sig2 = 0.5; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)
covfun <- function(locs) covfun.d(fields::rdist(locs))


## likelihood settings
me.var = 1e-8;
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

lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha = 2)




## generate grid of pred.locs
grid = seq(0,1,length = sqrt(n))
locs = as.matrix(expand.grid(grid,grid)) 

## set initial state
#Q = covfun(locs)
x0 = matrix(sig2*RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range),
                                   x = locs[,1], y = locs[,2], spConform = FALSE), ncol = 1)
evolFun = function(x) evolAdvDiff(x, diff = diffusion, adv = advection)


## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm')
approximations = list(mra = mra, low.rank = low.rank)#, exact = exact)


RRMSPE = list(); LogSc = list()

#foreach( iter=1:max.iter) %dopar% {
for (iter in 1:max.iter) {  
  
    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
    
    cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
    predsMRA = filter('mra', XY)
    cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
    predsLR  = filter('low.rank', XY)
   
    RMSPE = calculateRMSPE(predsMRA, predsLR, XY$x)
    
    write.csv(RMSPE, file = paste(resultsDir, "/", data.model, "/RMSPE.", iter, sep=""))
    
    m = M = 0
    for (t in 1:Tmax) {
                                        
        zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
        m = min(m, zrange[1])
        M = max(M, zrange[2])
    }
    for (t in 1:Tmax) {
        #pdf(paste("~/HVLF/simulations-big/", data.model, "/", t, ".pdf", sep=""))
        defpar = par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
        nna.obs = which(!is.na(XY$y[[t]]))
        fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], zlim = c(m, M), nx = sqrt(n), ny = sqrt(n), main = "obs" )
        fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim = c(m, M), nx = sqrt(n), ny = sqrt(n), main = "truth" )
        fields::quilt.plot( locs, predsMRA[[t]]$state, zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = "prediction" )
        par(defpar)
        #dev.off()
    }
}
