N = 960
Force = 10
dt = 0.005
M = 5
K = 31
Tmax = 1

X0 = matrix(rep(1,N))

evolFun = function(X) Lorenz04M2SimCpp(X, Force, K, dt, M, iter = 1, burn = 0, newAlgo=FALSE)


t0 = proc.time()
G1 = rootSolve::gradient(evolFun, X0)
d1 = proc.time() - t0


t1 = proc.time()
G2 = exactGradient(X0, K, M, dt, Force)  
d2 = proc.time() - t1


print(paste("Old gradient calculated in ", d1[3]))
print(paste("New gradient calculated in ", d2[3]))
