## settings -----------------------------------
SEED = 1996
set.seed(SEED)
N = 34 ** 2
COND_SET_SIZE = 5
FRAC_OBS = 1.0
TMAX = 20
NSTEPS = 1
N_PARTS = 1


## evolution function
C = 0.99
DIFFUSION = 0.000002
ADVECTION = 0.01
#evolFun = function(X) C * X#evolAdvDiff(X, adv = ADVECTION, diff = DIFFUSION)
evolFun = function(X) {
    for (s in 1:NSTEPS) {
        #X = evolAdvDiff(X, adv = ADVECTION, diff = DIFFUSION)
        X = evolDiff(X, diff = DIFFUSION, rho = 1)
    }
    return(X)
}
#evolFun = function(X) evolDiff(X, adv = ADVECTION, diff = DIFFUSION)


## covariance parameters
#SIG_02 = 0.1
#SIG2 = 0.2
SIG_02 = 1e-8
SIG2 = 1e-8
RANGE = 0.1
SMOOTH = 1.5
MEAN_COEFS = matrix(c(rep(0, TMAX), rep(1, TMAX)), ncol = 2, byrow = FALSE)


## likelihood
## DATA_MODEL = "gamma"
ALPHA = 2
ALPHAS = rep(2, TMAX)
## lik.params = list(data.model = DATA_MODEL, alpha = ALPHA)
DATA_MODEL = "gauss"
ME_VAR = 1e-12


## other
NCORES = 10



## sampling
prior    = list(a     = list(mean = log(ALPHA),  sd = 0),
                c     = list(mean = log(C),      sd = 0),
                sig2  = list(mean = log(SIG2),   sd = 0),
                range = list(mean = log(RANGE),  sd = 0),
                nu    = list(mean = log(SMOOTH), sd = 0))


PROP     = list(a     = list(sd = 0),
                c     = list(sd = 0),
                sig2  = list(sd = 0),
                range = list(sd = 0),
                nu    = list(sd = 0))

