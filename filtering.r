setwd("~/vecchiaFilter")
source('aux-functions.r')
source('getMatCov.r')
Rcpp::sourceCpp('src/getMatCovFromFactor.cpp')
source('scores.r')



## define the temporal evolution function
evolAdvDiff = function(state, adv=0, diff=0){
  
  # we assume that there is the same number of grid points
  # along each dimension
  N = dim(as.matrix(state))[1]
  Ny = Nx = sqrt(N)
  
  dx = dy = 1/Nx
  d = diff/(dx**2)
  
  c1 = 1 + 2*(d + d) - adv*(1/dx + 1/dy)
  c2 = - d + adv*(1/dy)
  c3 = - d
  
  diags = list(rep(c2, N-Nx), rep(c2, N-1), rep(c1, N), rep(c3, N-1), rep(c3,N-Nx) )
  E = Matrix::bandSparse(N, k=c(-Nx, -1, 0, 1, Nx), diag=diags)

  if (class(state) == 'matrix' || methods::is(state, 'sparseMatrix')){
    return( E %*% state )
  } else {
    return( as.numeric(E %*% as.matrix(state)) )
  }
    
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
filter = function(approx.name, XY, saveUQ=""){
  
  approx = approximations[[approx.name]]
  preds = list()
  cat(paste("filtering: t=1\n", sep = ""))
  cat("\tCalculate covariance elements from function: ")
  t0 = proc.time()
  covmodel = getMatCov(approx, covfun.d)
  t1 = proc.time()
  cat(paste((t1 - t0)[3], "\n"))
  
  mu.tt1 = rep(0, n)
  obs.aux = as.numeric(XY$y[[1]])
  
  cat("\tcalculate posterior: ")
  t0 = proc.time()
  preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
  t1 = proc.time()
  cat(paste((t1 - t0)[3], "\n"))


  cat("\tstore results\n")
  cat("\t\textract L.tt matrix\n")
  L.tt = getLtt(approx, preds.aux)
  cat("\t\textract mean vector\n")
  mu.tt = matrix(preds.aux$mean, ncol = 1)
  preds[[1]] = saveResults(preds.aux, L.tt, saveUQ)

  
  if ( Tmax > 1 ) {
    
    for (t in 2:Tmax) {
      
      cat(paste("filtering: t=", t, "\n", sep = ""))
      obs.aux = as.numeric(XY$y[[t]])
      
      cat("\tevolve the L.tt matrix:\n")
      t0 = proc.time()
      cat("\t\tbuild the evolution matrix\n")
      E = evolFun(Matrix::Diagonal(n))
      cat("\t\tmultpily E by L\n")
      Fmat = E %*% L.tt
      t1 = proc.time()      
      
      cat("\tCalculate covariance elements from factor: ")
      t0 = proc.time()
      
      M1 = getMatCov(approx, Matrix::t(Fmat), factor = TRUE)
      t1 = proc.time()
      cat(paste((t1 - t0)[3], "\n"))
      cat("\t... from function: ")
      t0 = proc.time()
      M2 = getMatCov(approx, covfun.d)
      t1 = proc.time()
      cat(paste((t1 - t0)[3], "\n"))
      covmodel = M1 + M2
      
      mu.tt1 = E %*% mu.tt
      
      cat("\tcalculate posterior: ")
      t0 = proc.time()
      preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
      t1 = proc.time()
      cat(paste((t1 - t0)[3], "\n"))
      
      cat("\tstore results\n")
      cat("\t\textract L.tt matrix\n")
      L.tt = getLtt(approx, preds.aux)
      cat("\t\textract mean vector\n")
      mu.tt = matrix(preds.aux$mean, ncol = 1)
        
      preds[[t]] = saveResults(preds.aux, L.tt, saveUQ)
      
    }
    
  }
  return( preds )
}
