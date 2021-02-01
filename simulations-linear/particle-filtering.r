source('~/vecchiaFilter/aux-functions.r')
source('~/vecchiaFilter/getMatCov.r')
Rcpp::sourceCpp('~/vecchiaFilter/src/getMatCovFromFactor.cpp')
library(doParallel)
registerDoParallel(cores = 15)


## implements systematic resampling scheme; this is generally better than multinomial
## resampling. which is what you get if you use
## sample(particles, Nparticles, prob=weights, replace=TRUE)
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




########## filtering ##########
filter = function(approx.name, XY, saveUQ="L"){
    
    approx = approximations[[approx.name]]
    if( is.null(approx) ){
        stop("Approximation with this name does not exist")
    }
    
    particles.all = matrix(rep(range, Nparticles*Tmax), ncol=Tmax)
    particles = exp(rnorm(Nparticles, mean=log(range), sd=prior.sd))           

    preds = list()

    for (t in 1:Tmax) {
        cat(sprintf("Filtering for t=%d\n", t))
        
        prop.mean = log(particles)

        particles = exp(rnorm(Nparticles, sd=prop.sd) + prop.mean)
        
        yt = as.numeric(XY$y[[t]])
        

        results = foreach( l = 1:Nparticles ) %dopar% {

            covparms[2] = particles[l]    
            covfun.d = function(D) GPvecchia::MaternFun(D, covparms)        
            covmodel = getMatCov(approx, covfun.d)
            mu.tt1 = rep(0, n)
                
            if( t>1 ) {
 
                E = evolFun(Matrix::Diagonal(n))            
                mu.tt1 = as.numeric(E %*% preds[[l]]$state)
                Fmat = E %*% preds[[l]]$L
                
                M1 = getMatCov(approx, Matrix::t(Fmat), factor = TRUE)
                covmodel = covmodel + M1
            }
    
            loglik = GPvecchia::vecchia_laplace_likelihood( yt, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, likparms = lik.params)
            preds.aux = GPvecchia::calculate_posterior_VL( yt, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, likparms = lik.params, return_all = TRUE)      
            
            L.tt = getLtt(approx, preds.aux)
            mu.tt = matrix(preds.aux$mean, ncol = 1)
            preds.tt = saveResults(preds.aux, L.tt, saveUQ)

            # we do not include the weight from the previous step because the particles were resampled so weights are the same
            logweight = loglik + dnorm(log(particles[l]), log(range), prior.sd, log = TRUE) - dnorm(log(particles[l]), prop.mean[l], prop.sd, log = TRUE)
            
            list(logweight, loglik, preds.tt)
        }
       
        preds = lapply(results, `[[`, 3)       
        
        logweights = sapply(results, `[[`, 1)        
        logweights = logweights - mean(logweights)
        weights = exp(logweights) / sum(exp(logweights))
        resampled.indices = resample(weights)
        
        #resampled.indices = sample(1:Nparticles, Nparticles, prob=weights, replace=TRUE)
               
        no.unique = length(unique(resampled.indices))
        cat(sprintf("\tNo. of unique particles after resampling: %d\n", no.unique))
        if( no.unique < 5 ){
            logliks = lapply(results, `[[`, 2)
            stop( "Less than 5 unique particles are left" )
        } else {
            preds = preds[ resampled.indices ]
            particles = particles[ resampled.indices ]
            particles.all[,t] = particles
        }

    }

    # for parallel for we have to transpose the list
    #preds = purrr::transpose(preds)
    return( list(particles = particles.all) )
    #return( list(preds = preds, particles = particles.all, logliks = logliks) )
}
