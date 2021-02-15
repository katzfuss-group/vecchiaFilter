setwd("~/vecchiaFilter")
#source("plotting.r")
source("data-application/particle-filtering.r")
source('aux-functions.r')



######### set parameters #########
set.seed(1988)
n = 34**2
m = 50
frac.obs = 1.0
Tmax = 2
diffusion = 0.00004
advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
Np = 2


## covariance parameters
sig2 = 1.0; range = .15; smooth = 0.5;
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.25;
alpha = 2
data.model = "gamma"  
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=alpha)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## set initial state
x0 = matrix(RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range, var = sig2),
                                          x = locs[,1], y = locs[,2], spConform = FALSE), ncol = 1)


## define Vecchia approximation
#exact = GPvecchia::vecchia_specify(locs, n - 1, conditioning = 'firstm', verbose = TRUE)
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
#low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, conditioning = 'firstm', verbose = TRUE)
#approximations = list(mra = mra, low.rank = low.rank, exact = exact)
approximations = list(mra = mra)


XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)
## predsLR = filter('low.rank', XY)
## predsE = filter('exact', XY)
cat("=== Running the new filter ===\n")
predsMRAnew = filter('mra', XY, saveUQ="L", old=FALSE)
cat("=== Running the old filter ===\n")
predsMRAold = filter('mra', XY, saveUQ="L", old=TRUE)
#save(predsMRA, file=sprintf("%s.particles.data", data.model))


uq.old = apply(predsMRAold$particles, 2, `quantile`, 0.9)
means.old = colMeans(predsMRAold$particles)
lq.old = apply(predsMRAold$particles, 2, `quantile`, 0.1)

uq = sapply(predsMRAnew$particles, function( t ) as.numeric(quantile(t[,"range"], 0.9)))
means = sapply(predsMRAnew$particles, function( t ) as.numeric(mean(t[,"range"])))
lq = sapply(predsMRAnew$particles, function( t ) as.numeric(quantile(t[,"range"], 0.1)))


#pdf(sprintf("%s-particles-over-time.pdf", data.model), width=6, height=4)
time = 1:Tmax
ylim = range(c(uq.old, lq.old, uq, lq, range))
ylim[2] = 1.1*ylim[2]
plot(time, uq, type="l", lty=2, ylim=ylim, ylab="range")
lines(time, means, type="l")
lines(time, lq, type="l", lty=2)

lines( time, uq.old, type="l", lty=2, col="blue" )
lines( time, lq.old, type="l", lty=2, col="blue" )
lines( time, means.old, type="l", col="blue")

legend(1, 95, legend=c("true value", "old method", "new method"),
       col=c("red", "blue", "black"), cex=0.8)
abline(h = range, col="red")
#dev.off()
