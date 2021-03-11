source('~/vecchiaFilter/aux-functions.r')
source('~/vecchiaFilter/getMatCov.r')
source('~/vecchiaFilter/data-application/particle-sampling.r')
Rcpp::sourceCpp('~/vecchiaFilter/src/getMatCovFromFactor.cpp')
library(doParallel)
registerDoParallel(cores = 6)


## implements systematic resampling scheme; this is generally better than multinomial
## resampling. which is what you get if you use
## sample(particles, Np, prob=weights, replace=TRUE)
## See: "A tutorial on particle filtering and smoothing"
resample = function(weights){

    Rcpp::cppFunction('int min_index(NumericVector v, double a){
                            NumericVector::iterator low=std::lower_bound (v.begin(), v.end(), a);
                            return (low - v.begin());
                       }')

    N = length(weights)
    u = runif(1)/N
    us = c(u, (1:(N-1))/N + u)

    cumweights = cumsum(weights)
    assignments = sapply(us, function(t) min_index(cumweights, t))

    tab = table(assignments)
    counts = rep(0, N)
    inds = as.numeric(names(tab))+1
    counts[inds] = as.numeric(tab)
    resampled.indices = rep(1:N, counts)
    return( resampled.indices )
        
}
    


saveResults = function(preds.aux, L.tt, saveUQ){

    results = list(state = matrix(preds.aux$mean, ncol = 1))
    if(saveUQ == "W"){
        results[["W"]] = preds.aux$W
    } else if(saveUQ == "L") {
        results[["L"]] = L.tt
    }
    return( results )
}



getWeight = function( p, samp.dist, l ){

    logweight = log.dist.eval( "a", p, prior, l )      - log.dist.eval( "a", p, samp.dist, l )
    logweight = log.dist.eval( "c", p, prior, l )      - log.dist.eval( "c", p, samp.dist, l )      + logweight
    logweight = log.dist.eval( "sig2", p, prior, l )   - log.dist.eval( "sig2", p, samp.dist, l )   + logweight
    logweight = log.dist.eval( "range", p, prior, l )  - log.dist.eval( "range", p, samp.dist, l )  + logweight
    logweight = log.dist.eval( "nu", p, prior, l )     - log.dist.eval( "nu", p, samp.dist, l )     + logweight

    if( length( logweight ) > 1){
        stop("something went wrong with calculating the weight in the getWeight() funciton")
    }

    return( logweight )
     
}



########## filtering ##########
#### new code
## We now extend this function to the case of multiple parameters.
## Here is a list of all parameters: likelihood params, c (multiplicative constant
## of the temporal evolution), covparms (constants in the matern model); that's the
## total of k=6.
filter = function(appr, Y, Np, lik.params, covparms, saveUQ="L"){

    if( is.null(approx) ){
        stop("Approximation with this name does not exist")
    }
    n = nrow(appr$locsord)
    Tmax = length(Y)
    
    particles.all = list()
    preds.all = list()
    sampling.d = prior
    
    preds = list()
    logliks = list()
    
    for (t in 1:Tmax) {
        cat(sprintf("+++ Filtering for t=%d +++\n", t))

        particles = sample.particles( Np, sampling.d )

        yt = as.numeric(Y[[t]])

        #results = list()
        #for( l in 1:Np ){
        results = foreach( l = 1:Np ) %dopar% {

            p = particles[l,]

            lik.params[[ "alpha" ]] = as.numeric(p["a"])
            
            if( t>1 ) {
                covparms = p[c("sig2", "range", "nu")]
            }

            covfun.d = function(D) GPvecchia::MaternFun(D, covparms)
            covmodel = getMatCov(appr, covfun.d)
            mu.tt1 = rep(0, n)

            if( t>1 ) {
                #evolFun = function(X) as.numeric(p["c"]) * evolAdvDiff(X, adv = advection, diff = diffusion)
                evolFun = function(X) as.numeric(p["c"]) * X
                E = evolFun(Matrix::Diagonal(n))
                
                mu.tt1 = as.numeric(E %*% preds[[l]]$state)
                Fmat = E %*% preds[[l]]$L

                M1 = getMatCov(appr, Matrix::t(Fmat), factor = TRUE)
                covmodel = covmodel + M1
                
            }                

            preds.aux = GPvecchia::calculate_posterior_VL( yt, appr, prior_mean = mu.tt1, likelihood_model = lik.params[["data.model"]], covmodel = covmodel, likparms = lik.params, return_all = TRUE, max.iter=100)   
            loglik = GPvecchia::vecchia_laplace_likelihood_from_posterior( yt, preds.aux, appr, prior_mean = mu.tt1, likelihood_model = lik.params[["data.model"]], covmodel = covmodel, likparms = lik.params, max.iter=100)
            
            L.tt = getLtt(appr, preds.aux)
            mu.tt = matrix(preds.aux$mean, ncol = 1)
            preds.tt = saveResults(preds.aux, L.tt, saveUQ)
            
            # we do not include the weight from the previous step because the particles were resampled so weights are the same
            logweight = loglik + getWeight(p, sampling.d, l)   

            list(logweight, loglik, preds.tt)
            #results[[l]] = list(logweight, loglik, preds.tt)
        }


        if( any(sapply(results, is.null)) ){
            browser()
        }
        
        preds = lapply(results, `[[`, 3)       

        logliks[[t]] = sapply(results, `[[`, 2, simplify=TRUE)        
        logweights = sapply(results, `[[`, 1)

        infinite.weights = which(!is.finite(logweights))
        logweights[ infinite.weights ] = 10*min(logweights)
        logweights = logweights - mean(logweights, na.rm=TRUE)
        weights = exp(logweights) / sum(exp(logweights))
        #cat(sprintf("Weights at time %d:\n", t))
        #print(weights)
        resampled.indices = resample(weights)

        no.unique = length(unique(resampled.indices))
        cat(sprintf("\tNo. of unique particles after resampling: %d\n", no.unique))
        if( no.unique < 1 ){
            logliks = lapply(results, `[[`, 2)
            stop( "Less than 1 unique particles are left" )
        } else {
            preds = preds[ resampled.indices ]
            preds.all[[ t ]] = preds[[ floor(Np/2)+1 ]]

            particles.all[[ t ]] = particles            
            ## these two should be swapped and this is only for exploring the likelihood space
            particles = particles[ resampled.indices, ]

            sampling.d = update.sampling( particles, sampling.d )
        }

    }

    # for parallel for we have to transpose the list
    #preds = purrr::transpose(preds)
    #return( list(particles = particles.all) )
    #return( list(preds = preds, particles = particles.all, logliks = logliks) )
    return( list(particles = particles.all, logliks = logliks, preds = preds.all) )
}
