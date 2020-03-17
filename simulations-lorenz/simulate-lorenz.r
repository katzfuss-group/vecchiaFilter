setwd("~/HVLF")
rm(list = ls())
source("aux-functions.r")
source("scores.r")
resultsDir = "simulations-lorenz"
#library(VEnKF)
library(rootSolve)
library(Matrix)
library(GPvecchia)



filter = function(approx.name, XY){
  
  approx = approximations[[approx.name]]
  
  MSigt = getMatCov(approx, Sigt)
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  covmodel = getMatCov(approx, Sig0)
  mu.tt1 = matrix(rep(0, nrow(XY$x[[1]])), ncol = 1)
  
  preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1,
                                      likelihood_model = data.model, covmodel = covmodel,
                                      covparms = covparms, likparms = lik.params, return_all = TRUE)
  
  L.tt  = getLtt(approx, preds.aux)
  mu.tt = matrix(preds.aux$mean, ncol = 1)
  preds[[1]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
  
  for (t in 2:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])
    Et = Matrix::Matrix(gradient(evolFun, mu.tt, centered = TRUE))
    
    forecast = evolFun(mu.tt)
    Fmat = Et %*% L.tt
    M = getMatCov(approx, t(Fmat), factor=TRUE)
    
    covmodel = M + MSigt
    
    preds.aux = calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                        likelihood_model = data.model, covmodel = covmodel,
                                        covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}


center_operator = function(x) {
  n = nrow(x)
  ones = rep(1, n)
  H = diag(n) - (1/n) * (ones %*% t(ones))
  H %*% x
}


getCovariance = function(N, Force, dt, K){
  fileName.all = paste("simulations-lorenz/Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep = "_")
  
  X = Matrix::Matrix(scan(fileName.all, quiet = TRUE), nrow = N) 
  X = center_operator(X)
  S = matrix((X %*% Matrix::t(X)) / (ncol(X) - 1), ncol=N)
  return(S)
}



######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 960
m = 50
frac.obs = 0.01
Tmax = 1


## evolution function ##
Force = 10
K = 32
dt = 0.005
M = 1
b = 0.33
evolFun = function(X) (1/b)*VEnKF::Lorenz04M2Sim(b*as.numeric(X), Force, K, dt, M, iter = 1, burn = 0, order = 4)
max.iter = 1


## covariance function
sig2 = 1e-6
range = .15; smooth = 1.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)


## likelihood settings
me.var = 1;
data.model = "logistic"
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha = 6)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(n))
locs = matrix(grid.oneside, ncol = 1)

## set initial state
Sig0 = (b**2)*getCovariance(n, Force, dt, K) + diag(1e-10, n)
x0 = b*getX0(n, Force, K, dt)
Sigt = sig2*Sig0



## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose = TRUE)
#exact = GPvecchia::vecchia_specify(locs, nrow(locs) - 1, conditioning = 'firstm')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm')
approximations = list(mra = mra, low.rank = low.rank)#, exact = exact)


RRMSPE = list(); LogSc = list()
#foreach( iter=1:max.iter) %dopar% {
for (iter in 1:max.iter) {

  XY = simulate.xy(x0, evolFun, Sigt, frac.obs, lik.params, Tmax)
  
  #cat(paste("iteration: ", iter, ", exact", "\n", sep = ""))
  #start = proc.time()
  #predsE = filter('exact', XY)
  #d = as.numeric(proc.time() - start)
  #cat(paste("Exact filtering took ", d[3], "s\n", sep = ""))
  cat(paste("iteration: ", iter, ", MRA", "\n", sep = ""))
  start = proc.time()
  predsMRA = filter('mra', XY)
  d = as.numeric(proc.time() - start)
  cat(paste("MRA filtering took ", d[3], "s\n", sep = ""))
  cat(paste("iteration: ", iter, ", LR", "\n", sep = ""))
  start = proc.time()
  predsLR  = filter('low.rank', XY)
  d = as.numeric(proc.time() - start)
  cat(paste("Low-rank filtering took ", d[3], "s\n", sep = ""))
  # 
  # RRMSPE = calculateRRMSPE(predsMRA, predsLR, predsE, XY$x)
  # LogSc = calculateLSs(predsMRA, predsLR, predsE, XY$x)
  # write.csv(RRMSPE, file = paste(resultsDir, "/", data.model, "/RRMSPE.", iter, sep = ""))
  # write.csv(LogSc, file = paste(resultsDir, "/", data.model, "/LogSc.", iter, sep = ""))
  
  #data = list(XY = XY, predsMRA = predsMRA, predsE = predsE, predsLR = predsLR)
  #save(data, file = paste(resultsDir, "/", data.model, "/sim.", iter, sep = ""))
}


######### plot results ##########
for (t in 1:Tmax) {
  zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm = TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm = TRUE)))))
  #zrange = range( unlist(lapply(XY$x, function(t) range(t, na.rm = TRUE))) )
  nna.obs = which(!is.na(XY$y[[t]]))
  
  plot( NULL, type = "l", xlim = c(0, 1), ylim = zrange, col = "red", main = paste("t =", t), xlab = "locations", ylab = "")
  ci.mra = getConfInt(predsMRA[[t]], 0.05)
  polygon(c(rev(locs), locs), c(rev(ci.mra$ub), ci.mra$lb), col = 'grey80', border = NA)
  lines(locs, XY$x[[t]], col = "red")
  points( locs[nna.obs,], XY$y[[t]][nna.obs], pch = 16, col = "black")
  
  #lines( locs, predsE[[t]]$state, type = "l", col = "black", lty = "dashed")
  lines( locs, predsMRA[[t]]$state, type = "l", col = "#500000")
  #lines( locs, predsLR[[t]]$state, type = "l", col = "black", lty = "dotted")
}
