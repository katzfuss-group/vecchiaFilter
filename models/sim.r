rm(list=ls())
library(VEnKF)
library(rootSolve)
setwd("~/HVLF/models")

generateInit = F

N = 120
Force = 10
dt = 0.005
M = 40
K = 32
Tmax = 100
seed = 1988
fileName = paste("init_Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep="_")


#### generate solutions ####

if( generateInit ){
  set.seed(seed)
  X0 = rnorm(N)  
} else {
  X0 = scan(fileName, quiet = TRUE)  
}

X1 = Lorenz04M2Sim(X0, Force, K, dt, M, iter=Tmax, burn=0, order=1)

if( generateInit ){
  last = X1[,Tmax]
  write(last, file=fileName)  
} 

X4 = Lorenz04M2Sim(X0, Force, K, dt, M, iter=Tmax, burn=0, order=4)



#### set up functions and constants for the linear approximation ####

taylor_x0 = 0#mean(X0)
taylor_x0_vec = rep(taylor_x0, N)

f1 = function(X) DeltaLorenz04M2Cpp(X, Force, K, dt, 1, order=1)
Jf1 = Matrix::Matrix(gradient(f1, taylor_x0_vec, centered=TRUE))
E1 = Matrix::diag(N) + Jf1
const1 = -Jf1 %*% taylor_x0_vec + f1(taylor_x0_vec)


f4 = function(X) DeltaLorenz04M2Cpp(X, Force, K, dt, 1, order=4)
Jf4 = Matrix::Matrix(gradient(f4, taylor_x0_vec, centered=TRUE))
E4 = Matrix::diag(N) + Jf4
const4 = -Jf4 %*% taylor_x0_vec + f4(taylor_x0_vec)



#### simulate linear model ####

Xl1 = matrix(rep(NA, N*Tmax), ncol=Tmax)
Xl4 = matrix(rep(NA, N*Tmax), ncol=Tmax)
Xl1[,1] = X0
Xl4[,1] = X0
oldXl1 = X0
oldXl4 = X0
for(t in 2:(Tmax*M)){

  newXl1 = E1 %*% oldXl1 + const1
  oldXl1 = newXl1

  newXl4 = E4 %*% oldXl4 + const4
  oldXl4 = newXl4

  # record only every Mth observation
  if((t %% M) == 0){
    #print(t/M)
    Xl1[,round(t/M)] = as.numeric(newXl1)
    Xl4[,round(t/M)] = as.numeric(newXl4)
  }

}



#### plot results ####
for(t in 1:Tmax){
  m = min(min(X1), min(X4)); M = max(max(X1), max(X4))
  #m = min(min(X1), min(Xl1), min(X4), min(Xl4)); M = max(max(X1), max(Xl1), max(X4), max(Xl4))
  plot(X1[,t], type="l", main=paste("t=", t, sep=""), col="black", ylim=c(m, M))
  #lines(Xl1[,t], type="l", lty=2, col="black")
  
  lines(X4[,t], col="red", type="l")
  #lines(Xl4[,t], col="red", type="l", lty=2)
  
  legend("topright", c("order 1", "order 4", "linear order 1", "linear order 4"),
         col=c("black", "red", "black", "red"), lty=c(1, 1, 2, 2))
}