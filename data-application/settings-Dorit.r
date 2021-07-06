## settings -----------------------------------
SEED = 1996
set.seed(SEED)
N = 34 ** 2
COND_SET_SIZE = 50
FRAC_OBS = 1.0
TMAX = 5
NSTEPS = 1
N_PARTS = 1
SUB_SAMPLE = 1.0

X_UPPER_LIMIT = Inf
LOC_QUANTILE = 1.0

#X_UPPER_LIMIT = -21.76952
#LOC_QUANTILE = 0.1

FILENAME = "~/vecchiaFilter/data-application/data/katzfuss-dorit-data/MIRS.csv"

## evolution function
C = 0.5
DIFFUSION = 0.000002
ADVECTION = 0.001
#evolFun = function(X) C * X#evolAdvDiff(X, adv = ADVECTION, diff = DIFFUSION)
evolFun = function(X) {
    for (s in 1:NSTEPS) {
        #X = evolAdvDiff(X, adv = ADVECTION, diff = DIFFUSION)
        X = evolDiff(X, diff = DIFFUSION)
    }
    return(C * X)
}
#evolFun = function(X) evolDiff(X, adv = ADVECTION, diff = DIFFUSION)

#MRA_OPTIONS = list(r = c(100, 50, 50, 10, 1), J = c(2, 4, 8, 256), M = 4)
MRA_OPTIONS = NULL#list(r = c(10, 10, 5, 5, 5, 5, 5, 5), J = c(2, 2, 2, 2, 4, 4, 4, 4, 4, 4), M = 9)


params = readr::read_csv("~/vecchiaFilter/data-application/static-estimation/params_1_11.csv", col_types = readr::cols())
`%>%` = dplyr::`%>%`
## covariance parameters
#SIG_02 = 0.164#0.045381#.5/C
SIG_02 = params %>% dplyr::select(sig2) %>% dplyr::summarize(mean(sig2)) %>% dplyr::pull()
SIG2 = (1 - C) * SIG_02
#RANGE = 0.0469
RANGE = params %>% dplyr::select(range) %>% dplyr::summarize(mean(range)) %>% dplyr::pull()
## SIG_02 = 0.5/C
## SIG2 = 1 - C*SIG_02
##RANGE = 0.15
SMOOTH = 2.5
#MEAN = 6.5


#MEAN_COEFS = matrix(c(-1.652611, 0.084358, -7.401559, 0.285934, -7.714005, 0.315610,
#                      -1.358753, 0.069904, -5.086568, 0.218950, 1.059064, -0.024937,
#                      -19.881873, 0.738750, 2.987584, -0.131574, -12.743291, 0.419174,
#                      -0.230070, 0.026262), ncol = 2, byrow = TRUE)

MEAN_COEFS = dplyr::select(params, beta0, beta1) %>% data.matrix()

## MEAN_COEFS = matrix(c(-7.401559,  0.285934, -7.714005,  0.315610,
##                       -1.358753,  0.069904, -5.086568,  0.218950,   1.059064, -0.024937,
##                       -19.881873, 0.738750,  2.987584, -0.131574, -12.743291,  0.419174,
##                       -0.230070,  0.026262), ncol = 2, byrow = TRUE)


## likelihood
## DATA_MODEL = "gamma"
## ALPHA = 2
#ALPHAS = c(19.390610, 3.054154, 4.595637, 2.654783, 6.494913, 2.575130, 3.318540, 0.244522, 4.863106)
#ALPHAS = c(3.516784, 19.390610, 3.054154, 4.595637, 2.654783, 6.494913, 2.575130, 3.318540, 0.244522, 4.863106)
#ALPHAS = c(3.5194290, 19.4212010, 3.0625700, 4.5951373, 2.6493544, 6.5099610, 2.5757933, 3.3188420, 0.2444594, 4.8677700)
#ALPHAS = rep(10o00)#mean(ALPHAS), TMAX)
ALPHAS = dplyr::select(params, a) %>% dplyr::pull()


## lik.params = list(data.model = DATA_MODEL, alpha = ALPHA)
DATA_MODEL = "gamma_mean"
#ME_VAR = 0.1


## other
#NCORES = 10


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
