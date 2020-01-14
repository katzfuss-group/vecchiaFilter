library(numDeriv)

f <- function(x) c(sin(x), cos(x))
fp = function(x) c(cos(x), -sin(x))


x <- (0:1)*2*pi
J = jacobian(sin, x)
print(J)