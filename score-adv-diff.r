RMSE = function(x, y){
  N = length(x)
  res = sqrt(sum((x-y)**2)/N)
  return(res)
}


calculateRRMSPE = function(predsMRA, predsLR, predsE, x){
  
  Tmax = length(x)
  RRMSPEs = c()
  
  for(t in 2:Tmax){
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
  for( t in 2:Tmax ){
    LS[t] = dposterior(y[[t]], preds[[t]])
  }
  return(LS)
}


dposterior = function(y, pred){
  
  mu = pred$state
  W= pred$W
  # VL calculates V, Laplace calculates W                                                                                                                                                                                        
  if("V" %in% names(pred)){
    det_term = sum(log(diag(pred$V)))
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
  results = cbind(time, LS.MRA, LS.LR)
  colnames(results) = c("time", "MRA", "LR")
  return(results)
}



plotScores = function(resultsDir, data.models){
  
  for( model in data.models ){
    
    RRMSPE = 0 
    LogS = 0
    iter = 1
    filePath = paste(resultsDir, "/", model, ".scores.", iter, sep="")
    
    while(file.exists(filePath)){
      load(filePath)
      RRMSPE = ((iter-1)/iter)*RRMSPE + data$RRMSPE[[iter]]/iter
      LogS = ((iter-1)/iter)*LogS + data$LogS[[iter]]/iter
      iter = iter+1
      filePath = paste(resultsDir, "/", model, ".", iter, sep="")
      
    }
    
    defpar = par(mfrow=c(2, 1))
    plot(RRMSPE[,1], RRMSPE[,2], type="l", lwd=2, col="#500000", main="RRMSPE", 
         ylab="RRMSPE", xlab="time", ylim=range(RRMSPE[,-1]))
    lines(RRMSPE[,1], RRMSPE[,3], col="black", lwd=2)
    legend("topright", c("MRA", "low-rank"), col=c("#500000", "black"))
    
    plot(LogS[,1], LogS[,2], type="l", lwd=2, col="#500000", main="Log Score", 
         ylab="Log score", xlab="time", ylim=range(LogS[,-1]))
    lines(LogS[,1], LogS[,3], col="black", lwd=2)
    legend("topright", c("MRA", "low-rank"), col=c("#500000", "black"))
    
    
    print(RRMSPE)
  } 
  
}




# plotSims = function(resultsDir, lik.model, iteration){
#   
#   dataPath = paste(resultsDir, "/", lik.model, ".", iteration, sep="")
#   locsPath = paste(resultsDir, "/locs", sep="")
#   load(dataPath)
#   load(locsPath)
#   Tmax = length(data$XY$x)
#   n = nrow(locs)
#   
#   
#   ########## plot results ########## 
#   for(t in 1:Tmax){
#     
#     X = as.numeric(data$XY$x[[t]])
#     MRA = data$predsMRA[[t]]$state
#     exact = data$predsE[[t]]$state
#     LR = data$predsLR[[t]]$state
#     
#     zrange = range(c(data$predsE[[t]][["state"]], data$predsVL[[t]][["state"]], unlist(lapply(data$XY$x, function(t) range(t, na.rm=TRUE)))))
#     
#     defpar <- par(mfrow = c(2,2),
#                   oma = c(5,4,2,0) + 0.1,
#                   mar = c(1,1,3,3) + 0.1)
#     
#     fields::quilt.plot( locs, X, main="truth", zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
#     fields::quilt.plot( locs, exact-MRA, main="MRA", nx=sqrt(n), ny=sqrt(n))#, zlim=zrange )
#     fields::quilt.plot( locs, exact, main="exact", nx=sqrt(n), ny=sqrt(n), zlim=zrange )
#     fields::quilt.plot( locs, exact-LR, main="low rank", nx=sqrt(n), ny=sqrt(n))#, zlim=zrange)
#     mtext(paste("t=", t, sep=""), outer = TRUE, cex = 1.5)
#     par(defpar)
#   }
#   
# }





setwd("~/HVLF")
resultsDir="simulations-linear"
data.models = c("poisson")
plotScores(resultsDir, data.models)
#calculateRRMSPE(resultsDir, data.models)