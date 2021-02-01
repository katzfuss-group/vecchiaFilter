
N = 10

weights = rgamma(N, shape=2)
                                        #weights = rnorm(N, mean=N/2, sd=N/1000)
weights = weights/sum(weights)

Rcpp::cppFunction('int min_index(NumericVector v, double a){
                      NumericVector::iterator low=std::lower_bound (v.begin(), v.end(), a);
                      return (low - v.begin());
                  }')

u = runif(1)/N
us = c(u, (1:(N-1))/N + u)

cumweights = cumsum(weights)
assignments = sapply(us, function(t) min_index(cumweights, t))


    
tab = table(assignments)
counts = rep(0, N)
inds = as.numeric(names(tab))+1
counts[inds] = as.numeric(tab)
resampled.indices = rep(1:N, counts)
print(resampled.indices)



