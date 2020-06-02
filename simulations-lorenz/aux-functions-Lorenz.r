getX0 = function(N, Force, K, dt, dir = '~/HVLF/simulations-lorenz/'){
  
  fileName = paste("init_Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep="_")
  filePath = paste(dir, fileName, sep="")
  
  generateInit = !file.exists(filePath)

  if( generateInit ){
    message("Initial state for these parameters has not yet been generated. Generating now.")
    X0 = rnorm(N)  
    X1 = Lorenz04M2Sim(X0, Force, K, dt, M, iter=100, burn=0, order=1)
    x0 = X1[,100]
    write(x0, file=filePath)  
  } else {
    x0 = scan(filePath, quiet = TRUE)  
  }
  if(class(x0)=='numeric'){
    x0 = matrix(x0, ncol=1)
  }
  return(x0)

}





center_operator = function(x) {
  n = nrow(x)
  ones = rep(1, n)
  H = diag(n) - (1/n) * (ones %*% t(ones))
  H %*% x
}



getLRMuCovariance = function(N, Force, dt, K){

  fileName.all = paste("~/HVLF/simulations-lorenz/Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep = "_")
  X = Matrix::Matrix(scan(fileName.all, quiet = TRUE), nrow = N)
  Xbar = matrix(rowMeans(as.matrix(X)), ncol=1)
  X = center_operator(X)
  S = matrix((X %*% Matrix::t(X)) / (ncol(X) - 1), ncol=N)

  return(list(mu = Xbar, Sigma = S))

}





filter = function(approx.name, XY){
  
  approx = approximations[[approx.name]]
    
  Tmax = length(XY$x)
  preds = list()
  obs.aux = as.numeric(XY$y[[1]])

  cat(paste("\tfiltering: t=1\n"))
  
  cat("\t\tcalculating forecast moments\n")
  covmodel = GPvecchia::getMatCov(approx, as.matrix(Sig0))
  mu.tt1 = mu
  
  cat("\t\tcalculating posterior\n")
  preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = mu.tt1,
                                      likelihood_model = data.model, covmodel = covmodel,
                                      covparms = NULL, likparms = lik.params, return_all = TRUE)
  
  cat("\t\tsaving the moments\n")
  L.tt  = getLtt(approx, preds.aux)
  mu.tt = matrix(preds.aux$mean, ncol = 1)
  preds[[1]] = list(state = mu.tt, W = preds.aux$W)#, V = preds.aux$V)
  
  if (Tmax == 1) { 
    return( preds )
  } 
    
  for (t in 2:Tmax) {
    
    cat(paste("\tfiltering: t=",t, "\n", sep=""))
    obs.aux = as.numeric(XY$y[[t]])

    cat("\t\tcalculating gradient...\n")
    Et = Matrix::Matrix(exactGradient(mu.tt, K, M, dt, Force))

    cat("\t\tcalculating forecast moments\n")
    forecast = evolFun(mu.tt)
    Fmat = Et %*% L.tt
    
    covmodel = GPvecchia::getMatCov(approx, as.matrix(Fmat %*% Matrix::t(Fmat) + Sigt))
    cat("\t\tcalculating posterior\n")
    preds.aux = GPvecchia::calculate_posterior_VL( obs.aux, approx, prior_mean = forecast,
                                                   likelihood_model = data.model, covmodel = covmodel,
                                                   covparms = covparms, likparms = lik.params, return_all = TRUE)
    cat("\t\tsaving the moments\n")

    L.tt = getLtt(approx, preds.aux)
    mu.tt = matrix(preds.aux$mean, ncol = 1)
    preds[[t]] = list(state = mu.tt, W = preds.aux$W, V = preds.aux$V)
  }
  return( preds )
}


plotResults = function(XY, predsE, predsMRA, predsLR, resultsDir = NULL){
  Tmax = length(XY$x)
  for (t in 1:Tmax) {
    if(t<10){
      number = paste("0", t, sep="")  
    } else {
      number = t
    }
    
    if(!is.null(resultsDir)){
      pdf(paste(resultsDir, "/", data.model, "/", number, ".pdf",sep=""), width=8, height=3.5)  
    }
    
    zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm = TRUE)))))
    
    if( data.model!='poisson' ){
      range2 = range(c(unlist(lapply(XY$y, function(t) range(t, na.rm = TRUE)))), na.rm=TRUE)
      zrange = range(c(zrange, range2))
    }
    nna.obs = which(!is.na(XY$y[[t]]))
    
    plot( NULL, type = "l", xlim = c(0, 1), ylim = zrange, col = "red", main = paste("t =", t), xlab = "", ylab = "")
    
    ci.mra = getConfInt(predsMRA[[t]], 0.05)
    polygon(c(rev(locs), locs), c(rev(ci.mra$ub), ci.mra$lb), col = 'grey80', border = NA)
    lines(locs, XY$x[[t]], col = "red")
    if( data.model!='poisson' ){
      points( locs[nna.obs,], XY$y[[t]][nna.obs], pch = 16, col = "black")
    }
    
    lines( locs, predsE[[t]]$state, type = "l", col = "black", lty = "dashed")
    lines( locs, predsMRA[[t]]$state, type = "l", col = "#500000")
    lines( locs, predsLR[[t]]$state, type = "l", col = "black", lty = "dotted")
    if(!is.null(resultsDir)){
      dev.off()
    }
  }
}  
