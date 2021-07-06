library(Matrix)
suppressPackageStartupMessages(library(tidyverse))
source("static-estimation.r")
source("../data/process-data.r")


ndays = 15
m = 50
n_sample = 100000
SMOOTH = 2.5
## these are the values fot the data from the Dorit paper
LOC_QUANTILE = 1
X_UPPER_LIMIT = Inf
## these are the values for the NASA modis data
#LOC_QUANTILE = 0.1
#X_UPPER_LIMIT = -21.76952 

maxiter = 20
set.seed(1996)
niter.mean = 40
params.names = c("a", "sig2", "range", "nu", "beta0", "beta1")
params = matrix(rep(0, ndays*length(params.names)), ncol=length(params.names))
colnames(params) = params.names



for( day.no in 1:ndays ){

    cat(sprintf("============== day %d =====================\n", day.no))
    TPW = read_day_MIRS( day.no ) %>%
    #TPW = read_day( day.no ) %>%
        dplyr::filter( x < X_UPPER_LIMIT ) #%>%
    #    dplyr::filter( x > quantile(x, 1-LOC_QUANTILE), y > quantile(y, 1-LOC_QUANTILE) ) #%>%
    #TPW = readr::read_csv("../../TPW_10k.csv")

    n = nrow(TPW)
    
    z = matrix( dplyr::select(TPW, value) %>% filter( !is.na(value) ) %>% pull(value), ncol=1 )
    x = TPW %>% filter( !is.na(value) ) %>% pull(x)
    y = TPW %>% filter( !is.na(value) ) %>% pull(y)
    locs = matrix(c(x, y), ncol=2, byrow=FALSE)
    #z = TPW %>% select(`1`) %>% data.matrix()
    #locs = TPW %>% select(x, y) %>% data.matrix()

    # test on subset
    sub_idx = sample(length(z), min(nrow(z), min(n, n_sample)), replace = FALSE)
    #pred_idx = sample(setdiff(1:length(z), sub_idx), n_sample, replace = FALSE)
    locs = locs[sub_idx,]
    obs = z[sub_idx]   

    #### trend estimation ####
    cat("Estimating mean\n")
    X = locs
    X[,1]= 1


    ## # Glm: IRLS method from taking deriv of llh
    beta = c(1,0)
    for(i in c(1:niter.mean)){
        XB = X %*% beta
        W  = -Matrix::Diagonal(x = as.numeric(exp(-XB)*obs) )
        A  = exp(-XB)*obs-1
        U  = W %*% XB - A
        beta = solve( t(X) %*% W %*% X , t(X) %*% U)
    }
    cat(sprintf("Mean coefficients: %f, %f\n", beta[1], beta[2]))
    XB = as.numeric(X %*% beta)
    
    ## Do parameter estimation for multiple m values
    cat(sprintf("Estimating parameters for m = %d\n", m))

    cat("Step 1, generating vecchia approximations\n")
    vecchia.approx = GPvecchia::vecchia_specify(locs, conditioning = 'mra', m = m)#, cond.yz = "zy")
    #vecchia.approx.IW = GPvecchia::vecchia_specify(sub_locs, m = m, conditioning = 'mra')

    ## Iterative method:  estimate a, then covparms, then a again
    cat("Step 2, optimizing parameters\n")
    a = a_prev = (mean(obs)/sd(obs))**2
    #covparms_prev = c(0.055376, 0.048875)#c(2, 0.1)
    covparms_prev = c(0.05, 0.1)#c(2, 0.1)

    cat("\tInitial parameter values:\n")
    cat(sprintf("\t\ta: %f\n", a_prev))
    #cat(sprintf("\t\tcovparms: %f, %f, %f\n", covparms_prev[1], covparms_prev[2], covparms_prev[3]))
    cat(sprintf("\t\tcovparms: %f, %f\n", covparms_prev[1], covparms_prev[2]))
    
    t_start = Sys.time()
    iter_count = 1

    for (i in 1:maxiter) {
        cat(sprintf("+++++ Iteration %d +++++\n", i))
        a = update_a(a_prev, c(covparms_prev, SMOOTH), vecchia.approx, NULL, XB, obs)
        covparms = fit_covparms(a, covparms_prev, vecchia.approx, NULL, XB, obs)
        cat(sprintf("Found shape parameter a = %f\n", a))
        cat(sprintf("Found covariance parameters (sig, range, smooth) = (%f, %f, %f)\n",covparms[1], covparms[2], covparms[3]))
        
        if (abs( a -a_prev) < 1e-3 &
           abs(covparms_prev[1] - covparms[1]) < 1e-2 &
           abs(covparms_prev[2] - covparms[2])/covparms_prev[2] < 1e-2 #&
           #abs(covparms_prev[3] - covparms[3]) < 1e-2
           ) {
               cat("Convergence criteria met (fitting all parameters)\n")
               iter_count = i
               break
           }
           a_prev = a
           covparms_prev = covparms
        }
    t_end = Sys.time()
    time_dur = as.double(difftime(t_end, t_start, units = "mins"))
    cat(sprintf("Estimation took %.2f minutes\n", time_dur))
    cat(sprintf("Parameters:\n"))
    cat(sprintf("a = %f, sig = %f, range = %f, smoothness = %f\n", a, covparms[1], covparms[2], SMOOTH))
    cat(sprintf("m = %d, time = %.2f, Niter = %d\n", m, time_dur, iter_count))

    params[ day.no, ] = c(a, covparms, SMOOTH, as.numeric(beta))
    
}

params.file = sprintf("params_1_%d.csv", ndays)
readr::write_csv(as.data.frame(params), file = params.file)
