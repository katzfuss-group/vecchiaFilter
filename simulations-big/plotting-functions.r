library(tcltk)
library(Matrix)
plotScores = function(resultsDir, data.models){
  
    for( model in data.models ){
    
        RRMSPE = 0 
        LogS = 0
        iter = 1
        fileRRMSPE = paste(resultsDir, "/", model, "/RRMSPE.", iter, sep="")
        fileLogSc = paste(resultsDir, "/", model, "/LogSc.", iter, sep="")

    
        while(file.exists(fileRRMSPE)){

            RRMSPEit = read.csv(fileRRMSPE, header=TRUE)[,-1]
            LogSit = read.csv(fileLogSc, header=TRUE)[,-1]

            if(model=="poisson" && max(abs(LogSit))>1e5){
                print(iter)
                M = max(abs(LogSit))
                print(LogSit[14:20,])
            }
            
            RRMSPE = ((iter-1)/iter)*RRMSPE + RRMSPEit/iter
            LogS = ((iter-1)/iter)*LogS - LogSit/iter
            iter = iter+1
        
            fileRRMSPE = paste(resultsDir, "/", model, "/RRMSPE.", iter, sep="")
            fileLogSc = paste(resultsDir, "/", model, "/LogSc.", iter, sep="")
            
        }
        
        LogS[,1] = -LogS[,1]

        X11()
        defpar = par(mfrow=c(2, 1))
        
        time = RRMSPE[,1]
        plot(time, RRMSPE[,"MRA"], type="l", lwd=2, col="#500000", main=paste(model, "RRMSPE", sep=" - "),
             ylab="RRMSPE", ylim=range(RRMSPE[,c("MRA", "LR")]), xaxt="n")
        lines(RRMSPE[,1], RRMSPE[,"LR"], col="black", lwd=2, lty=2)
        legend("topleft", c("MRA", "low-rank"), col=c("#500000", "black"), lty=c(1, 2), lwd=c(2, 2))
        
        plot(LogS[,1], LogS[,"MRA"], type="l", lwd=2, col="#500000", #main=paste(model, "Log Score", sep=" - "),
             ylab="Log score", xlab="time", ylim=range(LogS[,c("MRA", "LR")]))
        lines(LogS[,1], LogS[,"LR"], col="black", lwd=2, lty=2)
        legend("topleft", c("MRA", "low-rank"), col=c("#500000", "black"), lty=c(1, 2), lwd=c(2, 2))
        
        while (!is.null(dev.list())) Sys.sleep(1)  
    } 

}



plotSims = function(resultsDir, lik.model, iteration){

    dataPath = paste(resultsDir, "/", lik.model, "/sim.", iteration, sep="")
    locsPath = paste(resultsDir, "/locs", sep="")
    load(dataPath)
    load(locsPath)
    Tmax = length(data$XY$x)
    n = nrow(locs)

    zrange=c(0, 0)
    for(t in 1:Tmax){
        Xt = as.numeric(data$XY$x[[t]])
        pMRA = as.numeric(data$predsMRA[[t]]$state)
        pE = as.numeric(data$predsE[[t]]$state)
        pLR = as.numeric(data$predsLR[[t]]$state)
        zrange = range(c(range(pE), range(pLR), range(pMRA), range(Xt)), zrange)
    }
  
    ## plot results ##
    for(t in 1:Tmax){
        
        Xt = as.numeric(data$XY$x[[t]])
        pMRA = as.numeric(data$predsMRA[[t]]$state)
        pE = as.numeric(data$predsE[[t]]$state)
        pLR = as.numeric(data$predsLR[[t]]$state)

        x11()
        defpar <- par(mfrow = c(2,2),
                      oma = c(5,4,2,0) + 0.1,
                      mar = c(1,1,3,3) + 0.1)

        fields::quilt.plot( locs, Xt, main="truth", nx=sqrt(n), ny=sqrt(n), zlim=zrange )
        fields::quilt.plot( locs, pMRA, main="MRA", nx=sqrt(n), ny=sqrt(n), zlim=zrange )
        fields::quilt.plot( locs, pE, main="exact", nx=sqrt(n), ny=sqrt(n), zlim=zrange )
        fields::quilt.plot( locs, pLR, main="low rank", nx=sqrt(n), ny=sqrt(n), zlim=zrange)
        mtext(paste("t=", t, sep=""), outer = TRUE, cex = 1.5)
        par(defpar)
        while (!is.null(dev.list())) Sys.sleep(1)  
    }
    
}
