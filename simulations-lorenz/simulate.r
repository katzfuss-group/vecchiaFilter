rm(list=ls())
library(VEnKF)
library(rootSolve)
setwd("~/HVLF/simulations-lorenz")

generateInit = T
    
N = 960
Force = 10
dt = 0.005
M = 40
K = 32
Tmax = 80
seed = 1988
fileName.init = paste("init_Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep = "_")
fileName.all = paste("Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep = "_")





#### generate solutions ####

if (generateInit) {
  set.seed(seed)
  X0 = rnorm(N)  
} else {
  X0 = scan(fileName.init, quiet = TRUE)  
}
X = Lorenz04M2Sim(X0, Force, K, dt, M, iter = Tmax, burn = 0, order = 4)
#X1 = Lorenz04M2Sim(X0, Force, K, dt, M, iter=Tmax, burn=0, order=1)

browser()

if ( generateInit ) {
  Xlast = X[,Tmax]
  write(Xlast, file = fileName.init)  
} else {
  write(X, file = fileName.all)
}
