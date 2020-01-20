rm(list=ls())

RMSE = function(x, y){
  N = length(x)
  res = sqrt(sum((x-y)**2)/N)
  return(res)
}


updateRRMSPE = function(data, RRMSPE, iter){
  Tmax = length(data$XY$x)
  RRMSPEs = c()
  for(t in 2:Tmax){
    MRA = RMSE(data$XY$x[[t]], data$predsMRA[[t]]$state)
    LR = RMSE(data$XY$x[[t]], data$predsLR[[t]]$state)
    E  = RMSE(data$XY$x[[t]], data$predsE[[t]]$state)
    RRMSPEs.t = c(MRA/E, LR/E)
    RRMSPEs = rbind(RRMSPEs, RRMSPEs.t)
  }
  RRMSPE = RRMSPE*(iter-1)/iter + RRMSPEs/iter
  return( RRMSPE )
}


calculateRRMSPE = function(resultsDir, data.models){
  
  
  for( model in data.models ){
    
    RRMSPE = 0  
    iter = 1
    filePath = paste(resultsDir, "/", model, ".", iter, sep="")
    
    while(file.exists(filePath)){
      load(filePath)
      RRMSPE = updateRRMSPE(data, RRMSPE, iter)
      iter = iter+1
      filePath = paste(resultsDir, "/", model, ".", iter, sep="")
      
    }
    timeSeq = seq(from=2, to=(nrow(RRMSPE)+1))
    RRMSPE = cbind(RRMSPE, timeSeq)
    colnames(RRMSPE) = c("MRA", "low-rank", "time")
    row.names(RRMSPE) = c()
    print(RRMSPE)
  } 
}



plotSims = function(resultsDir, lik.model, iteration){
  
  dataPath = paste(resultsDir, "/", lik.model, ".", iteration, sep="")
  locsPath = paste(resultsDir, "/locs", sep="")
  load(dataPath)
  load(locsPath)
  Tmax = length(data$XY$x)
  n = nrow(locs)
  
  
  ########## plot results ########## 
  for(t in 1:Tmax){
    
    X = as.numeric(data$XY$x[[t]])
    MRA = data$predsMRA[[t]]$state
    exact = data$predsE[[t]]$state
    LR = data$predsLR[[t]]$state
    
    zrange = range(c(data$predsE[[t]][["state"]], data$predsVL[[t]][["state"]], unlist(lapply(data$XY$x, function(t) range(t, na.rm=TRUE)))))
    
    defpar <- par(mfrow = c(2,2),
                  oma = c(5,4,2,0) + 0.1,
                  mar = c(1,1,3,3) + 0.1)
    
    fields::quilt.plot( locs, X, main="truth", zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
    fields::quilt.plot( locs, exact-MRA, main="MRA", nx=sqrt(n), ny=sqrt(n))#, zlim=zrange )
    fields::quilt.plot( locs, exact, main="exact", nx=sqrt(n), ny=sqrt(n), zlim=zrange )
    fields::quilt.plot( locs, exact-LR, main="low rank", nx=sqrt(n), ny=sqrt(n))#, zlim=zrange)
    mtext(paste("t=", t, sep=""), outer = TRUE, cex = 1.5)
    par(defpar)
  }
  
}





setwd("~/HVLF")
resultsDir="simulations-linear"
data.models = c("poisson")
calculateRRMSPE(resultsDir, data.models)