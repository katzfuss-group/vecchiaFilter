plotScores = function(resultsDir, data.models){
  
  for( model in data.models ){
    
    RRMSPE = 0 
    LogS = 0
    iter = 1
    fileRRMSPE = paste(resultsDir, "/", model, ".RRMSPE.", iter, sep="")
    fileLogSc = paste(resultsDir, "/", model, ".LogSc.", iter, sep="")
    
    while(file.exists(fileRRMSPE)){

      RRMSPEit = read.csv(fileRRMSPE, header=TRUE)
      LogSit = read.csv(fileLogSc, header=TRUE)

      RRMSPE = ((iter-1)/iter)*RRMSPE + RRMSPEit/iter
      LogS = ((iter-1)/iter)*LogS + LogSit/iter
      iter = iter+1
      
      fileRRMSPE = paste(resultsDir, "/", model, ".RRMSPE.", iter, sep="")
      fileLogSc = paste(resultsDir, "/", model, ".LogSc.", iter, sep="")
      
    }
    
    defpar = par(mfrow=c(2, 1))
    plot(RRMSPE[,"time"], RRMSPE[,"MRA"], type="l", lwd=2, col="#500000", main="RRMSPE", 
         ylab="RRMSPE", xlab="time", ylim=range(RRMSPE[,c("MRA", "LR")]))
    lines(RRMSPE[,"time"], RRMSPE[,"LR"], col="black", lwd=2)
    legend("topright", c("MRA", "low-rank"), col=c("#500000", "black"), lty=c(1, 1), lwd=c(2, 2))
    
    plot(LogS[,"time"], LogS[,"MRA"], type="l", lwd=2, col="#500000", main="Log Score", 
         ylab="Log score", xlab="time", ylim=range(LogS[,c("MRA", "LR")]))
    lines(LogS[,"time"], LogS[,"LR"], col="black", lwd=2)
    legend("topright", c("MRA", "low-rank"), col=c("#500000", "black"), lty=c(1, 1), lwd=c(2, 2))

  } 
  
  print(paste("Loaded ", iter, " files."))
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