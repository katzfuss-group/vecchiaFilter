setwd("~/vecchiaFilter")
source("plotting.r")
source("particle-filtering.r")
source('aux-functions.r')


######### set parameters #########
set.seed(1988)
n = 34**2
m = 50
frac.obs = 0.3
Tmax = 20
diffusion = 0.00004
advection = 0.01
evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
Nparticles = 1000
prop.sd = 0.1
prior.sd = 0.5


## covariance parameters
sig2 = 1.0; range = .15; smooth = 0.5;
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.25;
data.model = "gauss"  
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)


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
predsMRA = filter('mra', XY)
#save(predsMRA, file="particles.data")

uq = apply(predsMRA$particles, 2, `quantile`, 0.9)
means = colMeans(predsMRA$particles)
lq = apply(predsMRA$particles, 2, `quantile`, 0.1)


pdf("particles-over-time.pdf")
time = 1:Tmax
ylim = range(c(uq, lq, range))
plot(time, uq, type="l", lty=2, ylim=ylim, ylab="range")
lines(time, means, type="l")
lines(time, lq, type="l", lty=2)
abline(h = range, col="red")
dev.off()
