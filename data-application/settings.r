## settings -----------------------------------
SEED = 1996
set.seed(SEED)
N = 34 ** 2
COND_SET_SIZE = 50
FRAC_OBS = 0.3
TMAX = 10
N_PARTS = 1
SUB_SAMPLE = 1.0


X_UPPER_LIMIT = -21.76952
LOC_QUANTILE = 0.1


## evolution function
C = 0.5
DIFFUSION = 0#.00004
ADVECTION = 0#.01
evolFun = function(X) C * X#evolAdvDiff(X, adv = ADVECTION, diff = DIFFUSION)


## covariance parameters
SIG_02 = 0.2#0.045381#.5/C
SIG2 = C * SIG_02#(1 - C) * SIG_02
RANGE = 0.04
## SIG_02 = 0.5/C
## SIG2 = 1 - C*SIG_02
## RANGE = 0.15
SMOOTH = 2.5
MEAN = 6.5
MEAN_COEFS = matrix(c(-1.652611, 0.084358, -7.401559, 0.285934, -7.714005, 0.315610,
                      -1.358753, 0.069904, -5.086568, 0.218950, 1.059064, -0.024937,
                      -19.881873, 0.738750, 2.987584, -0.131574, -12.743291, 0.419174,
                      -0.230070, 0.026262), ncol = 2, byrow = TRUE)
## MEAN_COEFS = matrix(c(-7.401559,  0.285934, -7.714005,  0.315610,
##                       -1.358753,  0.069904, -5.086568,  0.218950,   1.059064, -0.024937,
##                       -19.881873, 0.738750,  2.987584, -0.131574, -12.743291,  0.419174,
##                       -0.230070,  0.026262), ncol = 2, byrow = TRUE)


## likelihood
## DATA_MODEL = "gamma"
## ALPHA = 2
#ALPHAS = c(19.390610, 3.054154, 4.595637, 2.654783, 6.494913, 2.575130, 3.318540, 0.244522, 4.863106)
ALPHAS = c(3.516784, 19.390610, 3.054154, 4.595637, 2.654783, 6.494913, 2.575130, 3.318540, 0.244522, 4.863106)


## lik.params = list(data.model = DATA_MODEL, alpha = ALPHA)
DATA_MODEL = "gamma"
ME_VAR = 1e-12


## other
NCORES = 10


## sampling
prior    = list(a     = list(mean = log(3),      sd = 0),
                c     = list(mean = log(C),      sd = 0),
                sig2  = list(mean = log(SIG2),   sd = 0),
                range = list(mean = log(RANGE),  sd = 0),
                nu    = list(mean = log(SMOOTH), sd = 0))


PROP     = list(a     = list(sd = 0),
                c     = list(sd = 0),
                sig2  = list(sd = 0),
                range = list(sd = 0),
                nu    = list(sd = 0))
