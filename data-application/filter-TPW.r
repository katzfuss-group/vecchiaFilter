setwd("~/vecchiaFilter")
#source("plotting.r")
suppressPackageStartupMessages(library(tidyverse))
source("data-application/particle-filtering.r")
source('aux-functions.r')
source("data-application/data/process-data.r")


######### set parameters #########
set.seed(1996)
n.test = 34**2
m = 50
Np = 50

## starting values
testing = FALSE
testing.locs = TRUE



## temporal model settings
Tmax = 10
#diffusion = 0.000001
#advection = 0.01
diffusion = 0#0.000001
advection = 0#0.002
c = 0.8



## model error covariance
if( testing ) {
    sig2 = 0; range = 7; smooth = 1.5
} else {
    sig_02 = 25.0
    sig2 = 5.0; range = 7; smooth = 1.5
    prior.covparms = c(sig_02, range, smooth)
}

## likelihood settings
#data.model = "gauss"
#me.var = 0.25;
#lik.params = list(data.model = data.model, sigma = sqrt(me.var))#, alpha=a)
data.model = "gamma"
a = 2
lik.params = list(data.model = data.model, alpha=a)#sigma = sqrt(me.var), alpha=a)
frac.obs = 1.0


## temporal evolution
evolFun = function(X) c*evolAdvDiff(X, adv = advection, diff = diffusion)


## covariance
#covparms = c(sig2,range,smooth)
#covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
#covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## generate locations
if( testing.locs ){
    grid.oneside = seq(0,1,length = round(sqrt(n.test)))
    locs = as.matrix(expand.grid(grid.oneside,grid.oneside))
} else {
    locs = read_day(1) %>%
        select( x, y ) %>%
        dplyr::filter( x < -21.76952 ) %>%
        dplyr::filter( x < quantile(x, 0.05), y < quantile(y, 0.05) ) %>%
        data.matrix()
    n.test = nrow(data.matrix)
}
nx = length(unique(locs[,1]))
ny = length(unique(locs[,2]))
if( nx!= ny ){
    stop("Wrong dimensions")
}



## set initial state
if( testing ){
    indsx = which( locs[,1] > quantile(locs[,1], 0.45) & locs[,1] < quantile(locs[,1], 0.55) )
    indsy = which( locs[,2] > quantile(locs[,2], 0.45) & locs[,2] < quantile(locs[,2], 0.55) )
    x0 = rep(0, nrow(locs))
    x0[intersect(indsx,indsy)] = 1
} else {
    cat(sprintf("simulating at time t=%d\n", 1))
    x0 = RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range, var = sig_02),
                                  x = locs[,1], y = locs[,2], spConform = FALSE)
}
x0 = matrix(x0, ncol=1)

## simulate data
if( testing.locs ) {
    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
} else {
    print("load data")   
}

if( testing.locs) {
    ## Plot the simulated data
    #pdf("simulated-plots.pdf")
    #oldpar = par(mfrow=c(3, 8), oma=c(1, 1, 0, 0) + 1, mar=c(0, 0, 1, 1) + 1)
    #k = ceiling(Tmax / prod(par("mfrow")))
    oldpar = par(mfrow=c(2, 4), oma=c(1, 1, 0, 0) + 1, mar=c(0, 0, 1, 1) + 1)
    k = ceiling( Tmax / par("mfrow")[2] )
    zlim = range(sapply(XY$x, range))
    for( t in 1:Tmax ){
        if( t %% k == 1 ) {
            fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim = zlim, nx = nx, ny = ny, main = sprintf("true field at t=%d",t), axes=FALSE )
        }
    }
    for( t in 1:Tmax ){
        if(t %% k == 1) {
            inds.obs = which( !is.na(XY$y[[t]]) )
            fields::quilt.plot( locs[inds.obs,], as.numeric(XY$y[[t]])[inds.obs], nx = nx, ny = ny, main = sprintf("observations at t=%d",t) )
        }
    }
    par(oldpar)
    #dev.off()
}


## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
predsMRA = filter(mra, XY, Np, lik.params, prior.covparms, saveUQ="L")

par(mfrow=c(2, 2))
#pdf(sprintf("%s-particles-over-time.pdf", data.model), width=6, height=4)
for( par.name in colnames(predsMRA$particles[[1]]) ) {

    time = 1:Tmax
    
    uq = sapply(predsMRA$particles, function( t ) as.numeric(quantile(t[ , par.name ], 0.9)))
    means = sapply(predsMRA$particles, function( t ) as.numeric(mean(t[ , par.name ])))
    lq = sapply(predsMRA$particles, function( t ) as.numeric(quantile(t[ , par.name ], 0.1)))

    if( all(uq==lq) ){
        cat(sprintf("%s was not sampled\n", par.name))
        next()
    }
    
    ylim = range(c(uq, lq))
    ylim[2] = 1.1*ylim[2]
    plot(time, uq, type="l", lty=2, ylim=ylim, ylab=par.name, main=par.name)
    lines(time, means, type="l")
    lines(time, lq, type="l", lty=2)    
}
#dev.off()
