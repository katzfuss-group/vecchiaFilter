

# Base update function
g.mult=function(x.samp, y.samp, K){
  x.samp + K %*% (y.samp - H.mat %*% x.samp)
}



# Simple EnKF update function
#' EnKF Update with Covariance Matrix
#' 
#' This function is given an ensemble state matrix, an ensemble observation matrix,
#' and the estimated state covariance matrix, and it updates the ensemble matrix.
#'
#' @param x.mult.prior Ensemble of state data vectors at time t-1 to update
#' @param y.mult.i Ensemble of observation data at time t
#' @param sigma.mat The covariance matrix for the state ensemble
#' @param H.mat The observation matrix relating $x_t$ to $y_t$
#' @param tau.mat The observation noise matrix
#'
#' @return This function returns an updated state ensemble matrix
#' @export
#'
#' @examples
#' data(x.mult.prior)
#' data(y.mult.i)
#' data(H.mat)
#' data(tau.mat)
#' sig.mat = cov(t(x.mult.prior))
#' 
#' simple.update(x.mult.prior, y.mult.i, sigma.mat, H.mat, tau.mat)
#' 
simple.update=function(x.mult.prior, y.mult.i, sigma.mat, H.mat, tau.mat){
  # Compatibility Checks
  if(nrow(x.mult.prior) != nrow(y.mult.i) | ncol(x.mult.prior) != ncol(y.mult.i)){
    stop("Dimensions of x.mult.prior and y.mult.i are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(sigma.mat)){
    stop("Dimensions of x.mult.prior and sigma.mat are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(H.mat)){
    stop("Dimensions of x.mult.prior and H.mat are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(tau.mat)){
    stop("Dimensions of x.mult.prior and tau.mat are not compatible")
  }
  
  # Update calculations
  x.mult.prior = as.matrix(x.mult.prior)
  y.mult.i = as.matrix(y.mult.i)
  temp = H.mat %*% sigma.mat %*% t(H.mat) + tau.mat
  K = sigma.mat %*% t(H.mat) %*% solve(temp)
  x.mult.post = g.mult(x.mult.prior, y.mult.i, K)
  return(x.mult.post)
}



# EnKF update with cholesky decomps
#' EnKF Update with Cholesky Decomposition of Covariance Matrix
#' 
#' This function is given an ensemble state matrix, an ensemble observation matrix,
#' and the Cholesky decomposition of the estimated state covariance matrix, and 
#' it updates the ensemble matrix.
#'
#' @param x.mult.prior Ensemble of state data vectors at time t-1 to update
#' @param y.mult.i Ensemble of observation data at time t
#' @param chol.mat Cholesky decomposition of the state ensemble covariance matrix
#' @param H.mat The observation matrix relating $x_t$ to $y_t$
#' @param tau.mat The observation noise matrix
#'
#' @return This function returns the updated state ensemble matrix
#' @export
#'
#' @examples
#' data(x.mult.prior)
#' data(y.mult.i)
#' data(H.mat)
#' data(tau.mat)
#' sig.mat = chol(cov(t(x.mult.prior)))
#' 
#' chol.update(x.mult.prior, y.mult.i, sig.mat, H.mat, tau.mat)
#' 
chol.update=function(x.mult.prior,y.mult.i,chol.mat, H.mat, tau.mat){
  # Compatibility Checks
  if(nrow(x.mult.prior) != nrow(y.mult.i) | ncol(x.mult.prior) != ncol(y.mult.i)){
    stop("Dimensions of x.mult.prior and y.mult.i are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(chol.mat)){
    stop("Dimensions of x.mult.prior and chol.mat are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(H.mat)){
    stop("Dimensions of x.mult.prior and H.mat are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(tau.mat)){
    stop("Dimensions of x.mult.prior and tau.mat are not compatible")
  }
  
  # Update calculations
  new.sigma.chol = crossprod(chol.mat)
  
  temp.chol = new.sigma.chol + crossprod(H.mat, solve(tau.mat, H.mat))
  post.chol.mat = chol(temp.chol)
  
  temp.step = new.sigma.chol %*% x.mult.prior + crossprod(H.mat, solve(tau.mat, y.mult.i))
  step.chol = solve(t(post.chol.mat), temp.step)
  x.mult.post.chol = solve(post.chol.mat, step.chol)
  return(x.mult.post.chol)
}



# EnKF update with tapered cov matrix
#' EnKF Update with Tapered Covariance Matrix
#' 
#' This function is given an ensemble state matrix, an ensemble observation matrix,
#' and the estimated state covariance matrix, and it tapers the covariance matrix 
#' before updating the ensemble matrix.
#'
#' @param x.mult.prior Ensemble of state data vectors at time t-1 to update
#' @param y.mult.i Ensemble of observation data at time t
#' @param sigma.mat The estimated state ensemble covariance matrix
#' @param H.mat The observation matrix relating $x_t$ to $y_t$
#' @param tau.mat The observation noise matrix 
#' @param taper.radius The radius used to determine the tapered covariance matrix
#'
#' @return This function returns the updated state ensemble matrix 
#' @export
#'
#' @examples
#' 
#' data(x.mult.prior)
#' data(y.mult.i)
#' data(H.mat)
#' data(tau.mat)
#' sig.mat = cov(t(x.mult.prior))
#' 
#' tap.update(x.mult.prior, y.mult.i, sig.mat, H.mat, tau.mat, taper.radius = 0.3)
#' 
tap.update = function(x.mult.prior, y.mult.i, sigma.mat, H.mat, tau.mat, taper.radius){
  # Compatability Checks
  if(nrow(x.mult.prior) != nrow(y.mult.i) | ncol(x.mult.prior) != ncol(y.mult.i)){
    stop("Dimensions of x.mult.prior and y.mult.i are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(sigma.mat)){
    stop("Dimensions of x.mult.prior and sigma.mat are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(H.mat)){
    stop("Dimensions of x.mult.prior and H.mat are not compatible")
  }
  if(nrow(x.mult.prior) != nrow(tau.mat)){
    stop("Dimensions of x.mult.prior and tau.mat are not compatible")
  }
  
  # Build the tapered covariance matrix
  taper = Wendland(dist.mat, taper.radius, 1, 1)
  taper.sigma.mat = sigma.mat * taper
  
  # Use simple update function to calculate the update step
  simple.update(x.mult.prior, y.mult.i, taper.sigma.mat, H.mat, tau.mat)
}