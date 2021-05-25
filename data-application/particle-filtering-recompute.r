## requires: settings.r, particle-sampling.r and doParallel

## Imports -----------------------------------------
source('~/vecchiaFilter/data-application/particle-sampling.r')


# Systematic resampling scheme; this is generally better than multinomial
# resampling. which is what you get if you use
# sample(particles, Np, prob=weights, replace=TRUE)
# See: "A tutorial on particle filtering and smoothing"
resample = function(weights) {

    Rcpp::cppFunction('int min_index(NumericVector v, double a){
                            NumericVector::iterator low=std::lower_bound (v.begin(), v.end(), a);
                            return (low - v.begin());
                       }')
    
    N = length(weights)
    set.seed(round(100*weights[1]))
    u = runif(1) / N
    us = c(u, (1:(N - 1))/N + u)

    cumweights = cumsum(weights)
    assignments = sapply(us, function(t) min_index(cumweights, t))

    tab = table(assignments)
    counts = rep(0, N)
    inds = as.numeric(names(tab)) + 1
    counts[inds] = as.numeric(tab)
    resampled.indices = rep(1:N, counts)
    return(resampled.indices)
        
}
    


saveResults = function(preds.aux, L.tt, saveUQ) {

    results = list(state = matrix(preds.aux$mean, ncol = 1))
    if(saveUQ == "W"){
        results[["W"]] = preds.aux$W
    } else if(saveUQ == "L") {
        results[["L"]] = L.tt
    }
    return( results )

}



getWeight = function(p, samp.dist, l) {

    logweight = log.dist.eval("a", p, prior, l)      - log.dist.eval("a", p, samp.dist, l)
    logweight = log.dist.eval("c", p, prior, l)      - log.dist.eval("c", p, samp.dist, l)      + logweight
    logweight = log.dist.eval("sig2", p, prior, l)   - log.dist.eval("sig2", p, samp.dist, l)   + logweight
    logweight = log.dist.eval("range", p, prior, l)  - log.dist.eval("range", p, samp.dist, l)  + logweight
    logweight = log.dist.eval("nu", p, prior, l)     - log.dist.eval("nu", p, samp.dist, l)     + logweight

    if (length(logweight) > 1) {
        stop("something went wrong with calculating the weight in the getWeight() funciton")
    }

    return(logweight)
     
}


forecastStep = function(appr, prior_mean, prior_covmodel, Qcovparms, evolFun) {

    n = nrow(appr$locs)
    covfun.d = function(D) GPvecchia::MaternFun(D, Qcovparms)
    Qcovmodel = GPvecchia::getMatCov(appr, covfun.d)
    
    E = evolFun(Matrix::Diagonal(n))
        
    mu.tt1 = as.numeric(E %*% prior_mean)
    Fmat = E %*% prior_covmodel
    priorMatCov = GPvecchia::getMatCov(appr, Matrix::t(Fmat), factor = TRUE)
    
    covmodel = Qcovmodel + priorMatCov
    
    return(list(mean = mu.tt1, covmodel = covmodel))
}



    
updateStep = function(y, appr, forecast_mean, forecast_covmodel, lik.params, saveUQ) {

    loglik = NA
    preds.tt = list()
    ## We need this because for some parameter values the posterior cannot be reasonably evaluated    
    tryCatch({
        preds.aux = GPvecchia::calculate_posterior_VL(y, appr, prior_mean = forecast_mean, likelihood_model = lik.params[["data.model"]], covmodel = forecast_covmodel, likparms = lik.params, return_all = TRUE, max.iter=100)
        if (!preds.aux$cnvgd) {
            stop("Posterior estimation did not converge")
        }
        set.seed(2)
        loglik = GPvecchia::vecchia_laplace_likelihood_from_posterior( y, preds.aux, appr, prior_mean = forecast_mean, covmodel = forecast_covmodel)
        cat(sprintf("\t\titer: %d, loglik = %f\n", preds.aux$iter, loglik))

        L.tt = getLtt(appr, preds.aux)
        mu.tt = matrix(preds.aux$mean, ncol = 1)
        preds.tt = saveResults(preds.aux, L.tt, saveUQ)
        #logweight = loglik + getWeight(p, sampling.d, l)   
    },
    error = function(c) {
        errMsg = conditionMessage(c)
        if (errMsg == "Derivative of the loglikehood is infinite. Try different parameter values") {
            logweight <<- -Inf
            loglik <<- -Inf
        } else {
            loglik <<- -Inf
            #stop(errMsg)
        }
    })

    return(list(preds = preds.tt, loglik = loglik))
}



getResultsForParticle = function(p, l, indices, preds, appr, obs, saveUQ, t, sampling.d, prior_mean=NULL, prior_covmodel=NULL) {

    Qcovparms = p[c("sig2", "range", "nu")]
    if (t > 1) {
        previous_ind = indices[[t-1]][l]
        prior_mean = preds[[previous_ind]]$state
        prior_covmodel = preds[[previous_ind]]$L
    }
    
    evolFun = function(X) as.numeric(p["c"]) * X
    forecasted = forecastStep(appr, prior_mean, prior_covmodel, Qcovparms, evolFun)
    
    lik.params[["alpha"]] = as.numeric(p["a"])
    #if (t==2 & l==2){
    #    browser()
    #}
    updated = updateStep(obs, appr, forecasted$mean, forecasted$covmodel, lik.params, saveUQ)

    # We do not include the weight from the previous step because the particles
    # were resampled so weights are the same.
    logweight = updated$loglik + getWeight(p, sampling.d, l)

    return(list(logweight = logweight, loglik = updated$loglik, preds = updated$preds))
}




# filtering ---------------------

# We now extend this function to the case of multiple parameters.
# Here is a list of all parameters: likelihood params, c (multiplicative constant
# of the temporal evolution), covparms (constants in the matern model); that's the
# total of k=6.
filter_lean = function(appr, Y, Np, lik.params, prior_covparms, prior_mean = NULL, saveUQ = "L", seed = NULL){
    
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
    sampling.d = prior
    
    preds.all = list()
    logliks = list()
    indices = list()
    
    # this can be modified to accomodate sampling particles for unknown
    # prior covariance models
    covfun.d = function(D) GPvecchia::MaternFun(D, prior_covparms)
    prior_covmodel = getL00(appr, covfun.d)

    
    for (t in 1:Tmax) {

        cat(sprintf("+++ Filtering for t = %d +++\n", t))
        if (!is.null(seed)) {
            set.seed(seed + t)
        }
        particles = sample.particles( Np, sampling.d )

        #results = list();
        #for (l in 1:Np) {
        results = foreach(l = 1:Np) %dopar% {
            p = particles[l, ]
            cat(sprintf("\tWorking on particle %d, %f\n", l, p["c"]))
            res = getResultsForParticle(p, l, indices, preds, appr, Y[[t]], saveUQ, t, sampling.d, prior_mean, prior_covmodel)
            list(res$logweight, res$loglik)
            #results[[l]] = list(res$logweight, res$loglik)
        }
        if (any(sapply(results, is.null))) {
            cat(sprintf("The following particles led to errors:\n"))
            print(which(is.null(results)))
            warning("Parallelization error. Computations for certain particles failed")
            browser()
        }

        logweights = sapply(results, `[[`, 1)
        logliks[[t]] = sapply(results, `[[`, 2)

        # we need this because for some particles the likelihood is so small or has not been computed
        infinite.weights = which(!is.finite(logweights))
        logweights[infinite.weights] = 10*min(logweights[is.finite(logweights)])

        logweights = logweights - max(logweights, na.rm=TRUE)
        weights = exp(logweights) / sum(exp(logweights))
        resampled.indices = resample(weights)
        
        no.unique = length(unique(resampled.indices))
        cat(sprintf("\tNo. of unique particles after resampling: %d\n", no.unique))

        if( no.unique == 1 ){
            warning( "Only one unique particles is left." )
            #return(list(particles = particles, logliks = logliks, preds = preds, resampled.indices = indices))
        } #else {

        cat("\tRecomputing moments for remaining particles\n")
        preds.new = list()
        IndsOfParticlesToRecompute = unique(resampled.indices)
        preds.new = foreach(l=1:Np) %dopar% {
            res = list()
            if (l %in% IndsOfParticlesToRecompute) {
                res = getResultsForParticle(particles[l, ], l, indices, preds, appr, Y[[t]], saveUQ, t, sampling.d, prior_mean, prior_covmodel)
                res$preds
            } else {
                res[["preds"]] = list(NULL)
            }
            res$preds
        }

        preds = preds.new
        indices[[t]] = resampled.indices
        preds.all[[t]] = preds
        particles.all[[t]] = particles
        # these two might ultimately be swapped but this is easier for testing
        particles = particles[resampled.indices, ]
        
        sampling.d = update.sampling(particles, PROP)

        #}
    }

    return(list(particles = particles.all, logliks = logliks, preds = preds.all, resampled.indices = indices))
}
      
