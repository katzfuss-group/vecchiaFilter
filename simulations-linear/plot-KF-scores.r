load("EnKF.scores")

Niter = length(scores)
Tmax = nrow( scores[[1]][[1]] )


LogSc  = matrix( rep(0, Tmax*2), byrow=FALSE, ncol=2)
RRMSPE = matrix( rep(0, Tmax*2), byrow=FALSE, ncol=2)

for( i in 1:Niter ){

    if(is.null(scores[[i]]))
        next
    RRMSPE = RRMSPE + scores[[i]][[1]][,-1]
    LogSc  = LogSc + scores[[i]][[2]][,-1]

}


RRMSPE = RRMSPE/( Niter - sum(sapply(scores, is.null)) )
LogSc = LogSc/( Niter - sum(sapply(scores, is.null)) )

time = 1:Tmax


pdf("EnKF-RRMSPE.pdf")
plot(time, RRMSPE[,2], type="l", col="black", lty=2, ylim = range(RRMSPE), ylab="RRMSPE")
lines(time, RRMSPE[,1], type="l", col="red")
dev.off()

pdf("EnKF-LogSc.pdf")
plot(time, LogSc[,2], type="l", col="black", lty=2, ylim = range(LogSc), ylab="LogSc")
lines(time, LogSc[,1], type="l", col="red")
dev.off()

