source('~/vecchiaFilter/aux-functions.r')
source('~/vecchiaFilter/getMatCov.r')
Rcpp::sourceCpp('~/vecchiaFilter/src/getMatCovFromFactor.cpp')


######### set parameters #########
set.seed(1988)
n = 34**2
m = 100
frac.obs = 1
Tmax = 1
Nparticles = 20



## covariance parameters
sig2 = 1.0; range = .15; smooth = 0.5;
covparms = c(sig2,range,smooth)
covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.25;
data.model = "logistic"  
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## set initial state
x0 = matrix(RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range, var = sig2),
                                          x = locs[,1], y = locs[,2], spConform = FALSE), ncol = 1)
XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)


## define Vecchia approximation
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra')
   
particles = seq(1e-4, 1.0, length.out=Nparticles)
loglik = rep(0, Nparticles)
    
yt = as.numeric(XY$y[[1]])

for( l in 1:Nparticles ) {
    print(l)
    covparms[2] = particles[l]    
    covfun.d = function(D) GPvecchia::MaternFun(D, covparms)        
    covmodel = getMatCov(mra, covfun.d)
    loglik[l] = GPvecchia::vecchia_laplace_likelihood( yt, mra, likelihood_model = data.model, covmodel = covmodel, likparms = lik.params)
}

defpar = par(mfrow=c(1, 3))
fields::quilt.plot(locs[,1], locs[,2], XY$y[[1]], nx=sqrt(n), ny=sqrt(n))
fields::quilt.plot(locs[,1], locs[,2], XY$x[[1]], nx=sqrt(n), ny=sqrt(n))
plot(particles, loglik, type="l", main=data.model)
par(defpar)
