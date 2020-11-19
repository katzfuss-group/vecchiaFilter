source('~/vecchiaFilter/aux-functions.r')
source('~/vecchiaFilter/getMatCov.r')
Rcpp::sourceCpp('~/vecchiaFilter/src/getMatCovFromFactor.cpp')
source('~/vecchiaFilter/scores.r')


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
filter = function(approx.name, XY, saveUQ=""){
   
  approx = approximations[['exact']]
  preds = list()
  weights = matrix(rep(1/Nparticles, Nparticles*Tmax), ncol=Tmax)
  particles = matrix(rep(0, Nparticles*Tmax), ncol=Tmax)

  for(l in 1:Nparticles){

      preds[[l]] = list()
      particles[l,1] = covparms[2] = exp(rnorm(1, log(range), 0.25))
      covfun.d = function(D) GPvecchia::MaternFun(D, covparms)
      
      covmodel = getMatCov(approx, covfun.d)
  
      mu.tt1 = rep(0, n)
      obs.aux = as.numeric(XY$y[[1]])

      browser()
      posterior.lik = vecchia_laplace_likelihood( obs.aux, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params)
      preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
      
      L.tt = getLtt(approx, preds.aux)
      mu.tt = matrix(preds.aux$mean, ncol = 1)
      preds[[l]][[1]] = saveResults(preds.aux, L.tt, saveUQ)

      print(posterior.lik)
      
      weights[l,1] = posterior.lik
      
  }

   print("filtering for t>2")
  
  if ( Tmax > 1 ) {
    
    for (t in 2:Tmax) {

      for(l in 1:Nparticles){

          particles[l,t] = covparms[2] = exp(rnorm(1, log(particles[l,t-1]), 0.5))

          obs.aux = as.numeric(XY$y[[t]])
          E = evolFun(Matrix::Diagonal(n))

          mu.tt1 = E %*% preds[[l]][[t-1]]$state
          Fmat = E %*% preds[[l]][[t-1]]$L
     
          M1 = getMatCov(approx, Matrix::t(Fmat), factor = TRUE)
          
          covfun.d = function(D) GPvecchia::MaternFun(D, covparms)
          M2 = getMatCov(approx, covfun.d)
          covmodel = M1 + M2
          

          
          posterior.lik = GPvecchia::vecchia_laplace_likelihood( obs.aux, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params)
          preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)      

          L.tt = getLtt(approx, preds.aux)
          mu.tt = matrix(preds.aux$mean, ncol = 1)
          preds[[l]][[t]] = saveResults(preds.aux, L.tt, saveUQ)
          
          weights[l,t] = weights[l,t-1] * posterior.lik * dnorm(log(particles[l,t]), log(range), 0.25) / dnorm(log(particles[l,t]), log(particles[l,t-1]), 0.5)
          
      }
    
    }

  }
  return( preds )
}
