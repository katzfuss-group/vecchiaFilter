rm(list=ls())

## set parameters ##
set.seed(1988)
spatial.dim=2
n=15**2
frac.obs = 0.95

## covariance parameters ##
beta=2; sig2=0.2; range=.2; smooth=1.5; me=1e-6
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*GPvecchia::MaternFun(fields::rdist(locs),covparms)
nuggets=rep(me,n)


## grid of pred.locs ##
grid.oneside=seq(0,1,length=round(sqrt(n)))
locs=as.matrix(expand.grid(grid.oneside,grid.oneside))

# ## simulate data ##
Q = covfun(locs)
x = t(chol(Q)) %*% rnorm(n)
inds.obs = sample(1:n, round(frac.obs*n))
y = rep(NA, n)
eps = me*rnorm(n)
y[inds.obs] = x[inds.obs] + eps[inds.obs]
y.full = x+eps

cov.xy = covfun(locs)
nugg = nuggets; nugg[-inds.obs] = 1e6

var.y = covfun(locs[inds.obs,]) + nugg[inds.obs]
preds.exact = cov.xy[,inds.obs] %*% solve(var.y) %*% y[inds.obs]


## Vecchia prediction ##
m=20
vecchia.approx=GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', verbose=FALSE)
preds = GPvecchia::vecchia_prediction(y,vecchia.approx,covparms,nuggets)
preds.full = GPvecchia::vecchia_prediction(y.full,vecchia.approx,covparms,nuggets)




## plot results ##
zrange = range(c(x, Filter(function(t) !is.na(t), y), preds$mu.obs, preds.exact))
defpar = par(mfrow=c(2,2), oma=c(0, 0, 2, 0))
nna.obs = which(!is.na(y))
fields::quilt.plot(locs[nna.obs,], y[nna.obs], main="data", zlim=zrange, nx=sqrt(n), ny=sqrt(n))
fields::quilt.plot(locs, x, main="state", zlim=zrange, nx=sqrt(n), ny=sqrt(n))
fields::quilt.plot(locs, preds.full$mu.obs, main="GPvecchia pred.", zlim=zrange, nx=sqrt(n), ny=sqrt(n))
fields::quilt.plot(locs, as.numeric(preds.exact), main="exact pred.", zlim=zrange, nx=sqrt(n), ny=sqrt(n))
par(defpar)
