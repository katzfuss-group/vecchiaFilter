## settings -----------------------------------
#set.seed(1988)
N = 34 ** 2
COND_SET_SIZE = 50
FRAC_OBS = 1.0
TMAX = 5
N_PARTS = 400


## evolution function
C = 0.8
DIFFUSION = 0#.00004
ADVECTION = 0#.01
evolFun = function(X) C * X#evolAdvDiff(X, adv = ADVECTION, diff = DIFFUSION)



## covariance parameters
SIG_02 = 0.5/C
SIG2 = 1 - C*SIG_02
RANGE = 0.15
SMOOTH = 1.5




## likelihood 
## DATA_MODEL = "gamma"
ALPHA = 2
## lik.params = list(data.model = DATA_MODEL, alpha = ALPHA)
DATA_MODEL = "gauss"
ME_VAR = 1e-12

## other
NCORES = 16



## sampling
prior    = list(a     = list(mean = log(ALPHA),  sd = 0),
                c     = list(mean = C,           sd = 0.2),
                sig2  = list(mean = log(SIG2),   sd = 0.2),
                range = list(mean = log(RANGE),  sd = 0.1),
                nu    = list(mean = log(SMOOTH), sd = 0))


PROP     = list(a     = list(sd = 0),
                c     = list(sd = 0.1),
                sig2  = list(sd = 0.02),
                range = list(sd = 0.05),
                nu    = list(sd = 0))

