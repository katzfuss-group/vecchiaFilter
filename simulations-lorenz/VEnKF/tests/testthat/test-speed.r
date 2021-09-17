N = 960
Force = 10
dt = 0.005
M = 1
K = 3
Tmax = 500

X0 = rnorm(N)

t0 = proc.time()
X = Lorenz04M2SimCpp(X0, Force, K, dt, M, iter = Tmax, burn = 0, newAlgo=TRUE)
t1 = proc.time() - t0
print(paste("Simulation took", t1[3]))
print(dim(X))


t0 = proc.time()
X = Lorenz04M2SimCpp(X0, Force, K, dt, M, iter = Tmax, burn = 0, newAlgo=FALSE)
t1 = proc.time() - t0
print(paste("Simulation took", t1[3]))
print(dim(X))
plot(X[,Tmax], type="l")