rm(list=ls())
setwd("~/HVLF")
resultsDir="simulations-linear"

data.models = c("poisson")


updateRRMSPE = function(){
  
  Tmax = length(data$XY$x)
  N = nrow(data$XY$x[[1]])
  RRMSPEs = c()
  for(t in 2:Tmax){
    MRA = sum((data$XY$x[[t]] - data$predsMRA[[t]]$state)**2)/N
    LR = sum((data$XY$x[[t]] - data$predsLR[[t]]$state)**2)/N
    E  = sum((data$XY$x[[t]] - data$predsE[[t]]$state)**2)/N
    RRMSPEs.t = c(MRA/E, LR/E)
    RRMSPEs = rbind(RRMSPEs, RRMSPEs.t)
  }

  RRMSPE = RRMSPE*(iter-1)/iter + RRMSPEs/iter
  return( RRMSPE )
}



for( model in data.models ){

  RRMSPE = 0  
  iter = 1
  filePath = paste(resultsDir, "/", model, ".", iter, sep="")
  
  while(file.exists(filePath)){
    load(filePath)
    
    RRMSPE = updateRRMSPE()
    
    iter = iter+1
    filePath = paste(resultsDir, "/", model, ".", iter, sep="")
    
  }
  RRMSPE = cbind(RRMSPE, seq(1:nrow(RRMSPE)))
  colnames(RRMSPE) = c("MRA", "low-rank", "time")
  row.names(RRMSPE) = c()
  print(RRMSPE)
}