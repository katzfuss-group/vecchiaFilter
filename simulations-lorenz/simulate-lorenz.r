setwd("~/HVLF")
rm(list = ls())
source("aux-functions.r")
source("scores.r")
resultsDir = "simulations-lorenz"
library(VEnKF)
library(rootSolve)
library(doParallel)
registerDoParallel(cores = 6)



filter = function(approx.name, XY){
  
  approx = approximations[[approx.name]]
  
  MSigt = GPvecchia::getMatCov(approx, as.matrix(Sigt))
  
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])
  covmodel = GPvecchia::getMatCov(approx, as.matrix(Sig0))
  mu.tt1 = mu
  
  preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1,
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
    covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat %*% Matrix::t(Fmat) + sig2*Sig0))
    
    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                                   likelihood_model = data.model, covmodel = covmodel,
                                                   covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W)
  }
  return( preds )
}



center_operator = function(x) {
  n = nrow(x)
  ones = rep(1, n)
  H = diag(n) - (1/n) * (ones %*% t(ones))
  H %*% x
}



getLRMuCovariance = function(N, Force, dt, K){
  fileName.all = paste("simulations-lorenz/Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep = "_")
  X = Matrix::Matrix(scan(fileName.all, quiet = TRUE), nrow = N) 
  Xbar = matrix(rowMeans(as.matrix(X)), ncol=1)
  X = center_operator(X)
  S = matrix((X %*% Matrix::t(X)) / (ncol(X) - 1), ncol=N)
  return(list(mu = Xbar, Sigma = S))
}




######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 960
m = 50
frac.obs = 0.1
Tmax = 10


## evolution function ##
Force = 10
K = 32
dt = 0.005
M = 5
b = 0.2
evolFun = function(X) b*Lorenz04M2Sim(as.numeric(X)/b, Force, K, dt, M, iter = 1, burn = 0, order = 4)
max.iter = 1



## covariance function
sig2 = 0.1; range = .15; smooth = 1.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)


## likelihood settings
me.var = 0.2;
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
grid.oneside = seq(0,1,length = round(n))
locs = matrix(grid.oneside, ncol = 1)


## set initial state
cat("Loading the moments of the long-run Lorenz\n")
moments = getLRMuCovariance(n, Force, dt, K)
Sig0 = (b**2)*moments[["Sigma"]] + diag(1e-10, n)
mu = b*moments[["mu"]]
x0 = b*getX0(n, Force, K, dt)
Sigt = sig2*Sig0


## define Vecchia approximation
cat("Calculating the approximations\n")
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', ordering = 'maxmin')
exact = GPvecchia::vecchia_specify(locs, nrow(locs) - 1, ordering = 'maxmin', conditioning = 'firstm')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, ordering = 'maxmin', conditioning = 'firstm')
approximations = list(mra = mra, low.rank = low.rank, exact = exact)


RRMSPE = list(); LogSc = list()
foreach( iter=1:max.iter) %dopar% {
#for (iter in 1:max.iter) {  

    cat("Simulating data")
    XY = simulate.xy(x0, evolFun, Sigt, frac.obs, lik.params, Tmax)
   
    cat(paste("iteration: ", iter, ", exact", "\n", sep = ""))
    start = proc.time()
    predsE = filter('exact', XY)
    d = as.numeric(proc.time() - start)
    cat(paste("Exact filtering took ", d[3], "s\n", sep = ""))

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
    
    RRMSPE = calculateRRMSPE(predsMRA, predsLR, predsE, XY$x)
    LogSc = calculateLSs(predsMRA, predsLR, predsE, XY$x)
    write.csv(RRMSPE, file = paste(resultsDir, "/", data.model, "/RRMSPE.", iter, sep = ""))
    write.csv(LogSc, file = paste(resultsDir, "/", data.model, "/LogSc.", iter, sep = ""))

    print(RRMSPE)
    print(LogSc)

    ## plot results for the first iteration
    if( iter==1 ){
        for (t in 1:Tmax) {
            if(t<10){
                number = paste("0", t, sep="")  
            } else {
                number = t
            }
            pdf(paste(resultsDir, "/", data.model, "/", number, ".pdf",sep=""), width=8, height=3.5)
            zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm = TRUE)))))
            if( data.model!='poisson' ){
                range2 = range(c(unlist(lapply(XY$y, function(t) range(t, na.rm = TRUE)))))
                zrange = c(zrange, range2)
            }
            nna.obs = which(!is.na(XY$y[[t]]))
            
            plot( NULL, type = "l", xlim = c(0, 1), ylim = zrange, col = "red", main = paste("t =", t), xlab = "", ylab = "")
            ci.mra = getConfInt(predsMRA[[t]], 0.05)
            polygon(c(rev(locs), locs), c(rev(ci.mra$ub), ci.mra$lb), col = 'grey80', border = NA)
            lines(locs, XY$x[[t]], col = "red")
            if( data.model!='poisson' ){
                points( locs[nna.obs,], XY$y[[t]][nna.obs], pch = 16, col = "black")
            }
            
            lines( locs, predsE[[t]]$state, type = "l", col = "black", lty = "dashed")
            lines( locs, predsMRA[[t]]$state, type = "l", col = "#500000")
            lines( locs, predsLR[[t]]$state, type = "l", col = "black", lty = "dotted")
            dev.off()
        }
    }
}
