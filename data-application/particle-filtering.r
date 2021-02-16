source('~/vecchiaFilter/aux-functions.r')
source('~/vecchiaFilter/getMatCov.r')
source('~/vecchiaFilter/data-application/particle-sampling.r')
Rcpp::sourceCpp('~/vecchiaFilter/src/getMatCovFromFactor.cpp')
library(doParallel)
registerDoParallel(cores = 1)


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

    #logweight = log.dist.eval( "a", p, prior ) -         log.dist.eval( "a", p, samp.dist )
    #logweight = log.dist.eval( "c", p, prior ) -         log.dist.eval( "c", p, samp.dist ) + logweight
    #logweight = log.dist.eval( "sig2", p, prior ) -   log.dist.eval( "sig2", p, samp.dist ) + logweight
    logweight = log.dist.eval( "range", p, prior, l ) - log.dist.eval( "range", p, samp.dist, l ) #+ logweight
    #logweight = log.dist.eval( "nu", p, prior ) -       log.dist.eval( "nu", p, samp.dist ) + logweight

    if( length( logweight ) > 1 ){
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
filter = function(approx.name, XY, saveUQ="L", old=TRUE){

    approx = approximations[[approx.name]]
    if( is.null(approx) ){
        stop("Approximation with this name does not exist")
    }

    if( old ){
        #### old code
        particles.all = matrix(rep(range, Np*Tmax), ncol=Tmax)
        particles = rep(range, Np)
        sample.sd = prior[["range"]][["sd"]]
    } else {
        #### new code
        ## At each time point we now want to have a Np x k matrix
        ## where each column represent a draw from the joint parameter
        ## proposal distribution. We also need a list particles.all
        ## which holds all the draws through time
        particles.all = list()
        sampling.d = prior
    }

    
    
    preds = list()

    for (t in 1:Tmax) {
        cat(sprintf("+++ Filtering for t=%d +++\n", t))

        if( old ){
            prop.mean = log(particles)
            particles = exp(rnorm(Np, mean=prop.mean, sd=sample.sd))
            #cat(sprintf("Range particles sampled at time %d:\n", t))
            #print(particles)
        } else {
            particles = sample.particles( Np, sampling.d )
            #cat(sprintf("Range particles sampled at time %d:\n", t))
            #print(particles[,"range"])
        }
        
        yt = as.numeric(XY$y[[t]])
        
        #for( l in 1:Np ){
        results = foreach( l = 1:Np ) %dopar% {
            
            if( old ){
                covparms[2] = particles[l]
            } else {
                p = particles[l,]
                covparms = p[c("sig2", "range", "nu")]
                lik.params[[ "alpha" ]] = p["a"]                
            }

            covfun.d = function(D) GPvecchia::MaternFun(D, covparms)
            covmodel = getMatCov(approx, covfun.d)
            mu.tt1 = rep(0, n)

            if( t>1 ) {

                if( old ) {
                    E = evolFun(Matrix::Diagonal(n))
                } else {
                    E = p["c"] * evolFun(Matrix::Diagonal(n))
                }
                
                mu.tt1 = as.numeric(E %*% preds[[l]]$state)
                Fmat = E %*% preds[[l]]$L

                M1 = getMatCov(approx, Matrix::t(Fmat), factor = TRUE)
                covmodel = covmodel + M1

            }
            if( old ){
                loglik = GPvecchia::vecchia_laplace_likelihood( yt, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, likparms = lik.params)
                preds.aux = GPvecchia::calculate_posterior_VL( yt, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, likparms = lik.params, return_all = TRUE)
            } else {
                loglik = GPvecchia::vecchia_laplace_likelihood( yt, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, likparms = lik.params)
                preds.aux = GPvecchia::calculate_posterior_VL( yt, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, likparms = lik.params, return_all = TRUE)
            }   
            L.tt = getLtt(approx, preds.aux)
            mu.tt = matrix(preds.aux$mean, ncol = 1)
            preds.tt = saveResults(preds.aux, L.tt, saveUQ)

            # we do not include the weight from the previous step because the particles were resampled so weights are the same
            if( old ){
                logweight = loglik +
                    dnorm(log(particles[l]), prior[["range"]][["mean"]], prior[["range"]][["sd"]], log = TRUE) -
                    dnorm(log(particles[l]), prop.mean[l], prop[["range"]][["sd"]], log = TRUE)
            } else {
                logweight = loglik + getWeight(p, sampling.d, l)
                if( length(logweight) > 1 ){
                    stop("incorrect weight length")
                }
            }
            
            list(logweight, loglik, preds.tt)
        }
        
        preds = lapply(results, `[[`, 3)       
        
        logweights = sapply(results, `[[`, 1)        
        logweights = logweights - mean(logweights)
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
                if( old ){
                    particles = particles[ resampled.indices ]
                    particles.all[ ,t ]  = particles
                } else {
                    particles = particles[ resampled.indices, ]
                    particles.all[[ t ]] = particles
                    sampling.d = update.sampling( particles, sampling.d )
                }
        }
        
    }

    # for parallel for we have to transpose the list
    #preds = purrr::transpose(preds)
    return( list(particles = particles.all) )
    #return( list(preds = preds, particles = particles.all, logliks = logliks) )
}
