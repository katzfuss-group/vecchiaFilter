## requires: settings.r, particle-sampling.r and doParallel

## Imports -----------------------------------------
#source('/home/marcin/vecchiaFilter/data-application/particle-sampling.r')

 

saveResults = function(preds.aux, L.tt, saveUQ) {

    results = list(state = matrix(preds.aux$mean, ncol = 1))
    if(saveUQ == "W"){
        results[["W"]] = preds.aux$W
    } else if(saveUQ == "L") {
        results[["L"]] = L.tt
    }
    return( results )

}




forecastStep = function(appr, prior_mean, prior_covmodel, Qcovparms, evolFun) {
    cat(sprintf("\tForecast step\n"))
    n = nrow(appr$locs)
    covfun.d = function(D) GPvecchia::MaternFun(D, Qcovparms)
    cat(sprintf("\t\tCalculating Q matcov\n"))
    Qcovmodel = getMatCov(appr, covfun.d)
    
    E = evolFun(Matrix::Diagonal(n))

    mu.tt1 = as.numeric(E %*% prior_mean)
    Fmat = E %*% prior_covmodel
    cat(sprintf("\t\tCalculating forecast matcov\n"))
    priorMatCov = getMatCov(appr, Matrix::t(Fmat), factor = TRUE)
    
    covmodel = Qcovmodel + priorMatCov
    
    return(list(mean = mu.tt1, covmodel = covmodel))
}



    
updateStep = function(y, appr, forecast_mean, forecast_covmodel, lik.params, saveUQ) {
    cat(sprintf("\tUpdate step\n"))
    loglik = NA    
    ## We need this because for some parameter values the posterior cannot be reasonably evaluated    
    tryCatch({
        cat(sprintf("\t\tCalculating posterior\n"))
        preds.aux = GPvecchia::calculate_posterior_VL( y, appr, prior_mean = forecast_mean, likelihood_model = lik.params[["data.model"]], covmodel = forecast_covmodel, likparms = lik.params, return_all = TRUE, max.iter=1000)
        if (!preds.aux$cnvgd) {
            warning("Posterior estimation did not converge")
        }
        L.tt = getLtt(appr, preds.aux)
        U = Matrix::t(L.tt)
        vars = as.numeric(sapply(split(U@x, cut(1:length(U@x), U@p, labels=FALSE)), function(v) sum(v**2)))
        mu.tt = matrix(preds.aux$mean, ncol = 1)
        preds.tt = saveResults(preds.aux, L.tt, saveUQ)
        msg = NULL
    },
    error = function(c) {
        msg = conditionMessage(c)
        browser()
        logweight <<- -Inf
        loglik <<- -Inf
        cat(sprintf("%s\n", msg))
    })

    return(list(preds = preds.tt, vars = vars, msg = msg))
}




# filtering ---------------------

# We now extend this function to the case of multiple parameters.
# Here is a list of all parameters: likelihood params, c (multiplicative constant
# of the temporal evolution), covparms (constants in the matern model); that's the
# total of k=6.
filter = function(appr, Y, lik.params, prior_covparms, prior_mean = NULL, saveUQ = "L"){

    if (is.null(appr)) {
        stop("No vecchia approximation provided")
    }
    
    n = nrow(appr$locsord)
    if (is.null(prior_mean)) {
        prior_mean = matrix(rep(0, n), ncol=1)
    }
    Tmax = length(Y)
    
    particles.all = list()
    preds.all = list()
    
    preds = list()
    logliks = list()
    indices = list()
    vars.all = list()
    
    # this can be modified to accomodate sampling particles for unknown
    # prior covariance models
    covfun.d = function(D) GPvecchia::MaternFun(D, prior_covparms)
    prior_covmodel = getL00(appr, covfun.d)
    Qcovparms = c(SIG2, RANGE, SMOOTH)
    
    for (t in 1:Tmax) {

        cat(sprintf("+++ Filtering for t = %d +++\n", t))
        forecasted = forecastStep(appr, prior_mean, prior_covmodel, Qcovparms, evolFun)

        if (lik.params[["data.model"]]=="gamma") {
            lik.params[["alpha"]] = ALPHAS[t]
            lik.params[["mu"]] = MEAN_COEFS[t, 1] + MEAN_COEFS[t, 2] * mean(appr$locsord[, 2])
        }
       
        updated = updateStep(Y[[t]], appr, forecasted$mean, forecasted$covmodel, lik.params, saveUQ)

        prior_mean = updated$preds$state
        prior_covmodel = updated$preds$L
        

        preds[[t]] = list()
        preds[[t]][[1]] = list(state = prior_mean)
        vars.all[[t]] = list(updated$vars)
    }

    return(list(preds = preds, vars = vars.all))
}
