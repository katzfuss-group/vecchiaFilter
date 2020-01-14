rm(list=ls())


models = c("gauss", "poisson", "logistic", "gamma")
scores = list()
ms = c(1, 2, 5, 10, 20, 40, 60)
defpar = par(mfrow=c(1, 4), oma=c(0, 0, 0, 0), mar=c(4, 4, 2, 1))
m = M = 1.5
for( model in models) {
  filename = paste("scores/scores-", model, ".csv", sep="")
  scores = read.csv(filename)
  
  lr = as.numeric(scores[which(scores$X=='lr'),-1])
  vl.mra = as.numeric(scores[which(scores$X=='vl-mra'),-1])
  vl = as.numeric(scores[which(scores$X=='vl'), -1])
  
  m = min(m, min(c(lr, vl.mra, vl)))
  M = max(M, max(c(lr, vl.mra, vl)))
  
}




idx=0
for( model in models){
  idx = idx + 1
  filename = paste("scores/scores-", model, ".csv", sep="")
  scores = read.csv(filename)

  lr = as.numeric(scores[which(scores$X=='lr'),-1])
  vl.mra = as.numeric(scores[which(scores$X=='vl-mra'),-1])
  vl = as.numeric(scores[which(scores$X=='vl'), -1])

  if(idx==1){
    plot(ms, vl, type="o",, lty="dashed", ylab = "RRMSPE", main=model, ylim=c(m, M))
  }  else {
    plot(ms, vl, type="o",, lty="dashed", ylab="", main=model, ylim=c(m, M))
  }
  lines(ms, rep(1, length(ms)))
  lines(ms, vl.mra, type="o", col="#800000", lwd=2)
  lines(ms, lr, type="o", col="darkgrey", lwd=2)
}


par(defpar)