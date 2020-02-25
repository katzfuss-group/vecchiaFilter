rm(list = ls())
library(VEnKF)
library(Matrix)
setwd("~/HVLF/simulations-lorenz")

generateInit = T

N = 960
Force = 10
dt = 0.005
M = 40
K = 32
Tmax = 1
seed = 1988

fileName.all = paste("Lorenz04_N", N, "F", Force, "dt", dt, "K", K, sep = "_")
cov.file.name = paste("Lorenz04_covariance_N", N, "F", Force, "dt", dt, "K", K, sep = "_")

center_operator <- function(x) {
  n = nrow(x)
  ones = rep(1, n)
  H = diag(n) - (1/n) * (ones %*% t(ones))
  H %*% x
}



#### generate solutions ####

X = Matrix(scan(fileName.all, quiet = TRUE), nrow = N) 
X = center_operator(X)
S = (X %*% t(X)) / (ncol(X) - 1)
write(as.matrix(S), file = cov.file.name)