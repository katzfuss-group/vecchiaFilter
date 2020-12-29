source('~/vecchiaFilter/aux-functions.r')
source('~/vecchiaFilter/getMatCov.r')
Rcpp::sourceCpp('~/vecchiaFilter/src/getMatCovFromFactor.cpp')
library(doParallel)
registerDoParallel(cores = 15)



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

    logliks = matrix(rep(0, Nparticles*Tmax), ncol=Tmax)
    preds = list()

    for (t in 1:Tmax) {
        cat(sprintf("Filtering for t=%d\n", t))
        
        prop.mean = log(particles)
        pdf(sprintf("histograms-%d.pdf", t))
        par(mfrow=c(1, 2))
        hist(particles, main=sprintf("resampled particles from previous time"), xlim=c(0,0.6), ylim=c(0, 1), freq=TRUE)
        particles = exp(rnorm(Nparticles, sd=prop.sd) + prop.mean)
        hist(particles, main=sprintf("sampled particles at t=%d", t), xlim=c(0, 0.6), ylim=c(0, 1), freq=TRUE)
        dev.off()
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
       
        logliks[,t] = sapply(results, `[[`, 2)
        preds = lapply(results, `[[`, 3)
        
        logweights = sapply(results, `[[`, 1)        
        logweights = logweights - mean(logweights)
        weights = exp(logweights) / sum(exp(logweights))
        resampled.indices = sample(1:Nparticles, Nparticles, prob=weights, replace=TRUE)
        
        preds = preds[ resampled.indices ]
        particles = particles[ resampled.indices ]
        particles.all[,t] = particles
        logliks[,t] = logliks[ resampled.indices, t ]
        cat(sprintf("\tNo. of unique particles after resampling: %d\n", length(unique(particles.all[, t]))))
        if( length(unique(particles.all[, t])) < 5 ){
            stop( "Less than 5 unique particles are left" )
        }
    }

    # for parallel for we have to transpose the list
    #preds = purrr::transpose(preds)
    return( list(particles = particles.all, logliks = logliks) )
    #return( list(preds = preds, particles = particles.all, logliks = logliks) )
}
