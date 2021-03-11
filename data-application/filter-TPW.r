setwd("~/vecchiaFilter")
#source("plotting.r")
suppressPackageStartupMessages(library(tidyverse))
source("data-application/particle-filtering.r")
source('aux-functions.r')
source("data-application/data/process-data.r")


######### set parameters #########
set.seed(1996)
#n.test = 34**2
m = 50
Np = 36

## starting values
testing = FALSE
testing.locs = FALSE

X_UPPER_LIMIT = -21.76952
LOC_QUANTILE = 0.5

## temporal model settings
Tmax = 2


## model error covariance
sig_02 = 25.0
sig2 = 5.0; range = 5; smooth = 1.5
covparms = c(sig_02, range, smooth)

## likelihood settings
data.model = "gamma"
a = 0.12
lik.params = list(data.model = data.model, alpha=a)
frac.obs = 1.0


## temporal evolution
#c = 0.8
#evolFun = function(X) c*X#evolAdvDiff(X, adv = advection, diff = diffusion)


TPW = readr::read_csv("TPW_10k.csv")
#TPW = readr::read_csv("data-application/data/TPW.csv") %>%
#    dplyr::filter( x < X_UPPER_LIMIT ) %>%
#    dplyr::filter( x < quantile(x, LOC_QUANTILE), y < quantile(y, LOC_QUANTILE) ) %>%
#    sample_n(size = 10000)

locs = TPW %>% select(x, y) %>% data.matrix()
Y = lapply(split.default(TPW, colnames(TPW))[as.character(1:Tmax)], function(t) remove.mean( t %>% pull()) )
nx = length(unique(locs[,1]))
ny = length(unique(locs[,2]))


## set initial state
#cat(sprintf("simulating at time t=%d\n", 1))
#x0 = RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range, var = sig_02),
#                              x = locs[,1], y = locs[,2], spConform = FALSE)
#x0 = matrix(x0, ncol=1)


## ## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
predsMRA = filter(mra, Y, Np, lik.params, covparms, saveUQ="L")

save(predsMRA, file="predsMRA2")

## oldpar = par(mfrow=c(2, 3))
## for( par.name in colnames(predsMRA$particles[[1]]) ) {
    
##     time = 1:Tmax
    
##     uq = sapply(predsMRA$particles, function( t ) as.numeric(quantile(t[ , par.name ], 0.9)))
##     means = sapply(predsMRA$particles, function( t ) as.numeric(mean(t[ , par.name ])))
##     lq = sapply(predsMRA$particles, function( t ) as.numeric(quantile(t[ , par.name ], 0.1)))

##     if( all(uq==lq) ){
##        cat(sprintf("%s was not sampled\n", par.name))
##        next()
##     }
    
##     ylim = range(c(uq, lq))
##     ylim[2] = 1.1*ylim[2]
##     plot(time, uq, type="l", lty=2, ylim=ylim, ylab=par.name, main=par.name)
##     lines(time, means, type="l")
##     lines(time, lq, type="l", lty=2)    
## }
## par(oldpar)


