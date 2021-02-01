source('getMatCov.r')


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

  approx = approximations[[approx.name]]
  preds = list()
    
  for (t in 1:Tmax) {
      
      cat(paste("filtering: t=", t, "\n", sep = ""))
      yt = as.numeric(XY$y[[t]])
     
      covmodel = getMatCov(approx, covfun.d)
      mu.tt1 = rep(0, n)
      
      if(t > 1) {
          Fmat = apply(L.tt, 2, evolFun)
          M = getMatCov(approx, Matrix::t(Fmat), factor = TRUE)
          covmodel = covmodel + M
          mu.tt1 = evolFun( mu.tt )
      }
      
      preds.aux = GPvecchia::calculate_posterior_VL( yt, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, likparms = lik.params, return_all = TRUE, covparms = covparms)
      
      L.tt = getLtt(approx, preds.aux)
      mu.tt = matrix(preds.aux$mean, ncol = 1)
        
      preds[[t]] = saveResults(preds.aux, L.tt, saveUQ)
      
  }
    
  return( preds )
}
