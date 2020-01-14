rm(list=ls())
setwd("/home/marcin/HVLF/")
source("aux-functions.r")


# ######### set parameters #########
# args = commandArgs(trailingOnly=TRUE)
# if(length(args)!=1 || !(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))){
#     stop("One of the models has to be passed as argument")
# } else {
#     model = args[1]
# }

model="gauss"

set.seed(1988)
spatial.dim=2
n=34**2
ms=c(1, 2, 5, 10, 20, 40, 60)
frac.obs = 0.3

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
vecchia.exact = GPvecchia::vecchia_specify(locs, n-1, conditioning='firstm')
vecchia.mra = list(); low.ranks = list(); vecchia = list()
for(m in ms){
    vecchia[[m]] = GPvecchia::vecchia_specify(locs, m)
    vecchia.mra[[m]] = GPvecchia::vecchia_specify(locs, m, conditioning='mra')
    low.ranks[[m]] = GPvecchia::vecchia_specify(locs, m, conditioning='firstm')
}





## likelihood settings
me.var=1e-4;

scores = data.frame(matrix(rep(0, length(ms)*3), nrow=3))
colnames(scores) = paste("m", ms, sep="")
row.names(scores) = c("lr", "vl", "vl-mra")

message(paste("=== likelihood family: ", model, "==="))
lik.params = list(data.model = model, me.var=me.var, alpha=2)
  
## set initial state
Q = covfun(locs)
Qc = chol(Q)
x0 = Qc %*% matrix(rnorm(n), ncol=1); Sig0 = covfun(locs)



getSSE = function(vecchia.approx, obs, truth){
  preds = GPvecchia::calculate_posterior_VL( obs, vecchia.approx, likelihood_model = lik.params[["data.model"]],
                                                    covmodel = 'matern', covparms = covparms, likparms = lik.params, return_all = TRUE)
  mu = matrix(preds$mean, ncol=1)
  SSE = sum((mu - truth)**2)
  return(SSE)
}





## only Laplace, no Vecchia approximation

for( iter in 1:max.iter ){
  message(paste("Iteration ", iter))

  message("Simulating data")
  XY = simulate.xy(x0, evol, Q, frac.obs, lik.params, 1)
  obs = as.numeric(XY$y[[1]])
  truth = XY$x[[1]] 
  
      
  message("Calculating exact predictions")
  SSE.E = getSSE(vecchia.exact, obs, truth)
  RMSPE.E = sqrt(SSE.E/n)
  
  
  

  for( m in ms ){
    message(paste("Calculating scores for VL m = ", m))
    SSE.VL = getSSE(vecchia[[m]], obs, truth)
  
    
    message(paste("Calculating scores for LR m = ", m))
    SSE.LR = getSSE(low.ranks[[m]], obs, truth)
  
    
    message(paste("Calculating scores for VL-MRA m = ", m))
    SSE.MRA = getSSE(vecchia.mra[[m]], obs, truth)
    
    
    RMSPEs = sqrt(c(SSE.LR, SSE.VL, SSE.MRA)/n)/RMSPE.E
    colname = paste("m", m, sep="")
    scores[[colname]] = ((iter-1)/iter)*scores[[colname]] + RMSPEs/iter
    write.csv(scores,paste('scores/scores-', model, '.csv', sep=""))
  }
}

print(scores)
