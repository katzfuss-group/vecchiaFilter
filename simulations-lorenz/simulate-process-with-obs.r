setwd("~/HVLF")
rm(list = ls())
source("aux-functions.r")
source("simulations-lorenz/aux-functions-Lorenz.r")
source("scores.r")
source("radial_plot.r")
resultsDir = "simulations-lorenz"
library(VEnKF)
library(Matrix)

## library(foreach)
## library(iterators)
## library(parallel)
## library(doParallel)
## registerDoParallel(cores=5)



AllParamsAsString = list()



######### set parameters #########
set.seed(1988)
n = 960
frac.obs = 0.1
Tmax = 20

## evolution function ##
Force = 10
K = 32
dt = 0.005
M = 1
Mplot = 5
b = 0.1
evolFun = function(X) b*Lorenz04M2Sim(as.numeric(X)/b, Force, K, dt, M, iter = 1, burn = 0)



## covariance function
sig2 = 0.01; range = .15; smooth = 0.5;
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)

data.model = 'poisson'
lik.params = list(data.model = data.model)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(n))
locs = matrix(grid.oneside, ncol = 1)


## set initial state
cat("Loading the moments of the long-run Lorenz\n")
moments = getLRMuCovariance(n, Force, dt, K)
Sig0 = (b**2)*moments[["Sigma"]] + diag(1e-10, n)
mu = b*moments[["mu"]]
#x0 = b*getX0(n, Force, K, dt)
x0 = t(chol(Sig0)) %*% matrix(rnorm(n), ncol=1) + mu
#Sigt = sig2*Sig0
Sigt = covfun(locs)

m = 50
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', ordering = 'maxmin')
approximations = list(mra = mra)#, low.rank = low.rank, exact = exact)

cat("Simulating data\n")
XY = simulate.xy(x0, evolFun, Sigt, frac.obs, lik.params, Tmax)

for( t in 1:Tmax ){
    if( t %% Mplot != 1){
        XY$y[[t]] = rep(NA, n)
    }
}

preds = filter('mra', XY)
predsState = lapply( X = preds, FUN = function(el) el$state )


rangeX = 2*range(unlist(XY$x))
#dfY = data.frame(matrix(unlist(XY$y), nrow=n, byrow=F))

angles = c( locs, locs[ 1 ] ) * 360


for( t in 1:Tmax ){

    scale = 1#switch( t %% Mplot + 1, 1, 0, 0.33, 0.66, 1 )
    
    ind = ceiling( t/Mplot ) * Mplot
    y = 5 * scale * XY$y[[ ind ]]
    y[ is.na(y) ] = 0
    y = c(y, y[1])

    x = c( XY$x[[ t ]], XY$x[[ t ]][ 1 ] )
    filtered = c( predsState[[ t ]], predsState[[ t ]][ 1 ] )
        
    
    df = data.frame( latent = x, obs = y, angles = angles, filtered = filtered )

    if( t < 10 ){  
        filename = paste("frames/0", t, ".jpeg", sep="")
    } else {
        filename = paste("frames/", t, ".jpeg", sep="")
    }
    graph_radial(df, rangeX, filename)

}

