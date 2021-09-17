## Vecchia Update Source File


###########################################################
##  This file contains only the Vecchia Update function  ##
###########################################################



#' EnKF update function with Vecchia covariance approximation
#' 
#' Given an observation vector and an ensemble of state vectors, this function returns 
#' an updated ensemble of state vectors and the number of nearest neighbors used to
#' determine the update.
#' 
#'
#' @param x.mult.prior Ensemble of state data vectors at time t-1 to update
#' @param y.mult.i Ensemble of observation data at time t
#' @param H.mat The observation matrix relating $x_t$ to $y_t$
#' @param tau.mat The observation noise matrix
#' @param S The matrix of spatial locations in Max-min ordering
#' @param eps The threshold for selecting m based upon our hyperpriors. The default value 
#' is 0.001, and we do not recommend a value larger than 0.5 or smaller than 0.000001.
#' @param NN The maximum size wanted for creating a Nearest Neighbors array. This will 
#' put an upper limit on the number of possible neighbors the method will consider. 
#' The default value is 30.
#' @param m The number of neighbors to consider in the regression models for calculating
#' the modified Cholesky decomposition of the Covariance matrix. If left \code{NULL},
#' the method will select m as described above under the threshold \code{eps}. The default
#' value is set to \code{NULL}.
#'
#' @return A list of returned objects
#' \itemize{
#' \item{update}{The updated ensemble of state vectors}
#' \item{m}{The number of nearest neighbors used to determine the update}
#' }
#' @export
#'
#' @examples
#' data(x.mult.prior)
#' data(y.mult.i)
#' data(H.mat)
#' data(tau.mat)
#' data(S2d.ord)
#' 
#' vec.update(x.mult.prior, y.mult.i, H.mat, tau.mat, S2d.ord)
#' # Choose the m you want
#' vec.update(x.mult.prior, y.mult.i, H.mat, tau.mat, S2d.ord, m = 20)
#' 
vec.update = function(x.mult.prior, y.mult.i, H.mat, tau.mat, S, eps = 1e-03, NN = 30, m = NULL){
  # Vecchia Prior EnKF update
  ## set up priors from Brian's code
  ## Get the number of rows in our x.mult.prior
  n = nrow(x.mult.prior)
  x.mult.prior = as.matrix(x.mult.prior)
  
  ## Make sure eps is within a certain range
  if(eps < 1e-06 | eps > 5e-01)
  {
    stop("eps is not an acceptable value. Might cause numerical issues")
  }
  
  ## Have to generate our Nearest Neighbors array
  
  NNarray = find_ordered_nn(S, min(NN, m))
  NNarray = as.matrix(NNarray[, -1])
  
  
  # Get maximum a posteriori values for theta
  thets = optim(c(1,0,-1), loglikeli, datum = t(x.mult.prior), NNarray = NNarray, 
                eps = eps, m = m, method="L-BFGS-B", lower=-2, upper=3)$par
  
  # Use M.a.p. theta values to generate prior values for alpha, beta, and gamma
  thetps = thetas_to_priors(thets, n^2, eps = eps, m)
  alpha = thetps[[1]]
  beta = thetps[[2]]
  gamma = thetps[[3]]
  
  # Make sure m remains under the user allowed limit
  m = min(ncol(gamma), ncol(NNarray))
  
  ## Generate posterior values from sourced Brian's code
  # This sourced function generates the posterior parameters
  posts = get_posts(t(x.mult.prior), alpha, beta, gamma, NNarray)
  
  # This function use posterior parameters to sample from posteriors
  # Returns an upper triangular matrix
  vec.chol.mat = as.matrix(samp_posts(posts, NNarray, m))
  
  
  ## Calculate update using formula described in paper
  # Calculate new sigma and t(H.mat)%*%tau.mat.inv using crossprod because they 
  # need to be used twice in the code beneath
  tau.mat.inv=solve(tau.mat)
  new.sigma.vec = tcrossprod(vec.chol.mat)
  H.tau.mat = crossprod(H.mat, tau.mat.inv)
  
  # Calculations using update formula
  temp.vec = new.sigma.vec + H.tau.mat %*% H.mat
  post.vec.mat = chol(temp.vec)
  temp = new.sigma.vec %*% x.mult.prior + H.tau.mat %*% y.mult.i
  
  # Use solve against a vector instead of just a matrix to speed up calculation
  step.vec = solve(t(post.vec.mat), temp)
  x.mult.post.vec = solve(post.vec.mat, step.vec)
  
  # Return the updated sample
  return(list("update" = x.mult.post.vec, "m" = m))
}
