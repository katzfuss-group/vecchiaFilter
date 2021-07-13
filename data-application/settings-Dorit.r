## imports
`%>%` = dplyr::`%>%`


## settings -----------------------------------
COND_SET_SIZE = 80 
TMAX = 11
NSTEPS = 1
FRAC_OBS = 1.0
FILENAME = "~/vecchiaFilter/data-application/data/katzfuss-dorit-data/MIRS.csv"

#MRA_OPTIONS = list(r = c(100, 50, 50, 10, 1), J = c(2, 4, 8, 256), M = 4)
MRA_OPTIONS = NULL#list(r = c(20, 10, 10, 5, 5, 5, 3), J = 4, M = 7)


## evolution function -----------------------
C = 0.99
DIFFUSION = 0.000002
evolFun = function(X) {
    X = C * X#evolDiff(X, diff = DIFFUSION, nsteps = NSTEPS, c = C)
    return(X)
}


## covariance parameters -------------------
params = readr::read_csv("/home/marcin/vecchiaFilter/data-application/static-estimation/params_MIRS_1_5.csv", col_types = readr::cols())
#params = readr::read_csv("~/vecchiaFilter/data-application/static-estimation/params_1_11.csv", col_types = readr::cols())

SIG_02 = params %>% dplyr::select(sig2) %>% dplyr::filter(sig2 > 0) %>% dplyr::summarize(mean(sig2)) %>% dplyr::pull()
SIG2 = (1 - C) * SIG_02
RANGE = params %>% dplyr::select(range) %>% dplyr::filter(range > 0) %>% dplyr::summarize(mean(range)) %>% dplyr::pull()
SMOOTH = 1.5


## likelihood parameters -------------------
MEAN_COEFS = matrix(c(rep(2.56, TMAX), rep(0, TMAX)), ncol = 2, byrow = FALSE)

ALPHAS = dplyr::select(params, a) %>% dplyr::pull()
DATA_MODEL = "gamma_mean"



## ## other
## #NCORES = 10


## ## sampling
## prior    = list(a     = list(mean = log(3),      sd = 0),
##                 c     = list(mean = log(C),      sd = 0),
##                 sig2  = list(mean = log(SIG2),   sd = 0),
##                 range = list(mean = log(RANGE),  sd = 0),
##                 nu    = list(mean = log(SMOOTH), sd = 0))


## PROP     = list(a     = list(sd = 0),
##                 c     = list(sd = 0),
##                 sig2  = list(sd = 0),
##                 range = list(sd = 0),
##                 nu    = list(sd = 0))
