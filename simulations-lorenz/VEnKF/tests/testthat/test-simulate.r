N = 960
Force = 10
dt = 0.005
M = 1
K = 1
Tmax = 100

X0 = rep(1, N)
X = VEnKF::Lorenz04M2SimCpp(X0, Force, K, dt, M, iter = Tmax, burn = 0, newAlgo = FALSE)

test_that("Lorenz04M2Sim works", {
  expect_gt(sum(abs(X - X0)), 0)
})
