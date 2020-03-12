RMSE = function(x, y){
  N = length(x)
  res = sqrt(sum((x-y)**2)/N)
  return(res)
}



calculateRMSPE = function(predsMRA, predsLR, x){

  Tmax = length(x)
  RMSPEs = c()
  
  for(t in 1:Tmax){
    MRA = RMSE(x[[t]], predsMRA[[t]]$state)
    LR = RMSE(x[[t]], predsLR[[t]]$state)
    RMSPEs.t = c(t, MRA, LR)
    RMSPEs = c(RMSPEs, RMSPEs.t)
  }
  RMSPEs = matrix( RMSPEs, ncol=3, byrow=TRUE)
  colnames(RMSPEs) = c("time", "MRA", "LR")
  return( RMSPEs )
    
}



calculateRRMSPE = function(predsMRA, predsLR, predsE, x){
  
  Tmax = length(x)
  RRMSPEs = c()
  
  for(t in 1:Tmax){
    MRA = RMSE(x[[t]], predsMRA[[t]]$state)
    LR = RMSE(x[[t]], predsLR[[t]]$state)
    E  = RMSE(x[[t]], predsE[[t]]$state)
    RRMSPEs.t = c(t, MRA/E, LR/E)
    RRMSPEs = c(RRMSPEs, RRMSPEs.t)
  }
  RRMSPEs = matrix( RRMSPEs, ncol=3, byrow=TRUE)
  colnames(RRMSPEs) = c("time", "MRA", "LR")
  return( RRMSPEs )
}




LogS = function(preds, y){
  
  Tmax = length(preds)
  LS = rep(0, Tmax-1)
  for( t in 1:Tmax ){
    LS[t] = dposterior(y[[t]], preds[[t]])
  }
  return(LS)
}


dposterior = function(y, pred){
  mu = pred$state
  W= pred$W
  # VL calculates V, Laplace calculates W                                                                                                                                                                                        
  if("V" %in% names(pred)){
    det_term = sum(log(Matrix::diag(pred$V)))
  }else{
    V=t(chol(Matrix::forceSymmetric(pred$W)))
    #V = t(chol((W + t(W))/2))                                                                                                                                                                                                   
    det_term = sum(log(diag(V)))
    #det_term =log(det(W))/2                                                                                                                                                                                                     
  }
  quad_term = -Matrix::t(y-mu) %*% W %*% (y-mu)/2
  pi_term = -length(y)/2*log(2*pi)
  # value is of class "dgeMatrix"                                                                                                                                                                                                
  return((quad_term+det_term+pi_term)[1,1])
}


calculateLSs = function(predsMRA, predsLR, predsE, x){
  LS.MRA = LogS(predsMRA, x)
  LS.LR = LogS(predsLR, x)
  LS.E = LogS(predsE, x)
  time = seq(1:length(LS.MRA))
  results = cbind(time, LS.E-LS.MRA, LS.E - LS.LR)
  colnames(results) = c("time", "MRA", "LR")
  return(results)
}