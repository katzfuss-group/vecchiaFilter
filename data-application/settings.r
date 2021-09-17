## imports
`%>%` = dplyr::`%>%`


## settings -----------------------------------
COND_SET_SIZE = 15
TMAX = 2
NSTEPS = 1
FRAC_OBS = 0.99
FILENAME = "~/vecchiaFilter/data-application/data/katzfuss-dorit-data/MIRS.csv"


MRA_OPTIONS = NULL#list(r = c(25, 15, 15, 5, 5, 5, 3), J = 4, M = 7)


## evolution function -----------------------
C = 0.99
DIFFUSION = 0.000003
evolFun = function(X) {
    X = evolDiff(X, diff = DIFFUSION, nsteps = NSTEPS, rho = C)
    return(X)
}


## covariance parameters -------------------
params = readr::read_csv("/home/marcin/vecchiaFilter/data-application/static-estimation/MIRS-Gaussian-params-smooth-1_5.csv", col_types = readr::cols())

SIG_02 = params %>% dplyr::select(sig2) %>% dplyr::filter(sig2 > 0) %>% dplyr::summarize(mean(sig2)) %>% dplyr::pull()
SIG2 = (1 - C) * SIG_02
RANGE = params %>% dplyr::select(range) %>% dplyr::filter(range > 0) %>% dplyr::summarize(mean(range)) %>% dplyr::pull()
SMOOTH = 1.5


## likelihood parameters -------------------
MEAN_COEFS = dplyr::select(params, mu) %>% data.matrix()
MEAN_COEFS[10,1] = sum(MEAN_COEFS[-10,1])/10

ME_VAR = 4.5 ** 2
DATA_MODEL = "gauss"
lik.params = list(data.model = DATA_MODEL, sigma = sqrt(ME_VAR))
