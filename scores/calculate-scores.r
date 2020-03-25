rm(list=ls())
setwd("/home/marcin/HVLF/")
source("aux-functions.r")


######### set parameters #########
# args = commandArgs(trailingOnly=TRUE)
# if(length(args)!=1 || !(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))){
#     stop("One of the models has to be passed as argument")
# } else {
#     model = args[1]
# }

model="poisson"

set.seed(1988)
spatial.dim=2
n=34**2
ms=c(40)#c(1, 2, 5, 10, 20, 40, 60)
ms = c(60)
frac.obs = 0.3
Tmax = 2
diffusion = 0.00004
advection = 0.01
max.iter = 1

## covariance parameters
sig2=0.5; range=.15; smooth=1.5; 
covparms = c(sig2,range,smooth)
covfun <- function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)

## generate grid of pred.locs
grid.oneside=seq(0,1,length=round(sqrt(n)))
locs=as.matrix(expand.grid(grid.oneside,grid.oneside)) 

## Initialize approximation objects
message("=== Initialize approximations ===")
vecchia.exact = GPvecchia::vecchia_specify(locs, n-1, conditioning='firstm')#, mra.options = list(r=c(n, 0)))
vecchia.approxs = list(); low.ranks = list()
for(m in ms){
    vecchia.approxs[[m]] = GPvecchia::vecchia_specify(locs, m, conditioning='mra')
    low.ranks[[m]] = GPvecchia::vecchia_specify(locs, m, conditioning='firstm')
}





## likelihood settings
me.var=1e-4;

scores = data.frame(matrix(rep(0, length(ms)*2), nrow=2))
colnames(scores) = paste("m", ms, sep="")
row.names(scores) = c("lr", "vl")

message(paste("=== likelihood family: ", model, "==="))
lik.params = list(data.model = model, me.var=me.var, alpha=2)
  
## set initial state
Q = covfun(locs)
x0 = matrix(rep(0, n), ncol=1); Sig0 = covfun(locs)
  
## only Laplace, no Vecchia approximation

for( iter in 1:max.iter ){
  message(paste("Iteration ", iter))

  message("Simulating data")
  XY = simulate.xy(x0, evol, Q, frac.obs, lik.params, Tmax)
    
  cat("Calculating exact predictions")
  
  predsE = list()
  L.00.e = getL00(vecchia.exact, covfun, locs)
  predsE[[1]] = list(state=x0, L=L.00.e)
  forecast=x0; Fmat=L.00.e
    
  SSE.E = 0
    
  for(t in 1:Tmax){
    cat(".")
    obs.aux = as.numeric(XY$y[[t]])
    covmodel = GPvecchia::getMatCov(vecchia.exact, Fmat %*% Matrix::t(Fmat) + Q)
    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, vecchia.exact, likelihood_model = lik.params[["data.model"]],
                                                     covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
    L.tt = getLtt(vecchia.exact, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol=1)
    predsE[[t]] = list(state=mu.tt, L=L.tt)
    if(t<Tmax){
      forecast = evol(mu.tt)
      Fmat = evol(L.tt)
    }
    
    SSE.E = SSE.E + sum((predsE[[t]]$state - XY$x[[t]])**2)
  }
  cat("\n")
  RMSPE.E = sqrt(SSE.E/(n*Tmax))
  
  
  ## Laplace + Vecchia approximation

  for( m in ms ){
    message(paste("Calculating scores for VL m = ", m))

    RMSPEs.vl = c()  
    vecchia.approx = vecchia.approxs[[m]]
    L.00.vl = getL00(vecchia.approx, covfun, locs)

    predsVL = list()
    predsVL[[1]] = list(state=x0, L=L.00.vl)
    SSE.VL = 0
    
    forecast=x0; Fmat=L.00.vl
      
    for(t in 1:Tmax){
      obs.aux = as.numeric(XY$y[[t]])
      covmodel = GPvecchia::getMatCov(vecchia.approx, Fmat %*% Matrix::t(Fmat) + Q)
      preds.aux.vl = GPvecchia::calculate_posterior_VL( obs.aux, vecchia.approx, likelihood_model = lik.params[["data.model"]],
                                                          covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
      L.tt = getLtt(vecchia.approx, preds.aux.vl)
      mu.tt = matrix(preds.aux.vl$mean, ncol=1)
      predsVL[[t]] = list(state=mu.tt, L=L.tt)
      if(t<Tmax){
        forecast = evol(mu.tt)
        Fmat = evol(L.tt) 
      }
      SSE.VL = SSE.VL + sum((predsVL[[t]]$state - XY$x[[t]])**2)
    }
    #RMSPEs.vl = c(RMSPEs.vl, sqrt(SSE.VL/(n*Tmax))/RMSPE.E)
    RMSPE.vl = sqrt(SSE.VL/n)/RMSPE.E

    
    message(paste("Calculating scores for LR m = ", m))

    RMSPEs.lr = c()  
    vecchia.approx = low.ranks[[m]]
    L.00.lr = getL00(vecchia.approx, covfun, locs)

    predsLR = list()
    predsLR[[1]] = list(state=x0, L=L.00.lr)
    SSE.LR = 0

    forecast = x0; Fmat = L.00.lr

    for( t in 1:Tmax) {
        obs.aux = as.numeric(XY$y[[t]])
        covmodel = GPvecchia::getMatCov(vecchia.approx, Fmat %*% Matrix::t(Fmat) + Q)
        preds.aux.lr = GPvecchia::calculate_posterior_VL( obs.aux, vecchia.approx, likelihood_model = lik.params[["data.model"]],
                                                         covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
        L.tt = getLtt(vecchia.approx, preds.aux.lr)
        mu.tt = matrix(preds.aux.lr$mean, ncol=1)
        predsLR[[t]] = list(state=mu.tt, L=L.tt)
        if(t<Tmax){
            forecast = evol(mu.tt)
            Fmat = evol(L.tt)
        }
        SSE.LR = SSE.LR + sum((predsLR[[t]]$state - XY$x[[t]])**2)
    }
    #RMSPEs.lr = c(RMSPEs.lr, sqrt(SSE.LR/(n*Tmax))/RMSPE.E)
    RMSPE.lr = sqrt(SSE.LR/n)/RMSPE.E
    
    RMSPEs = c(RMSPE.lr, RMSPE.vl)
    colname = paste("m", m, sep="")
    scores[[colname]] = ((iter-1)/iter)*scores[[colname]] + RMSPEs/iter
    write.csv(scores,paste('scores-', model, '.csv', sep=""))
  }
}

print(scores)
