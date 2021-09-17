# Rfile to export Lorenz functions from cpp without allowing access to those functions
# And with wrapping to check compatability

#' Single Lorenz simulation
#'
#' This function simulates one iteration through a Lorenz chaos model 
#' described as the second model in 
#' \href{https://journals.ametsoc.org/doi/pdf/10.1175/JAS3430.1}{Lorenz (2005)}
#'
#' @param X0 The vector of points you want to pass through a Lorenz chaos model.
#' @param F_Lor A Forcing constant for the chaos model. Recommended values are
#' @param K_Lor A constant that determines the wavelength within the model. A recommended 
#' value is the length of your data vector divided by 30.
#' @param dt The internal step size. The default value is set to 0.005 
#' @param M The number of internal steps. The default value is set to 40
#' @param order The order of the RK solution method to use.
#' The default value is 1
#'
#' @return The function returns a new vector that represents a single step in a Lorenz chaos
#' model as described in the paper above.
#' @export
#'
#' @examples
#' X0 = rnorm(960)
#' Lorenz04M2(X0, 10, 32, 0.005, 40)
#' 
Lorenz04M2 = function(X0, F_Lor, K_Lor, dt = 0.005, M = 40, order = 1){
  N_Lor = length(X0)
  
  if(!is.vector(X0)){ stop("X0 should be a vector")}
  
  Xout = X0 + DeltaLorenz04M2Cpp(X0, F_Lor, K_Lor, dt, M, order)
  return(Xout)
}



#' Lorenz Model Simulation
#' 
#' This function simulates a specified number of iterations through a Lorenz chaos model 
#' described as the second model in 
#' \href{https://journals.ametsoc.org/doi/pdf/10.1175/JAS3430.1}{Lorenz (2005)}
#' 
#'
#' @param X0 The vector of points you want to pass through a Lorenz chaos model.
#' @param F_Lor A Forcing constant for the chaos model. Recommended values are
#' @param K_Lor A constant that determines the wavelength within the model. A recommended 
#' value is the length of your data vector divided by 30.
#' @param dt The internal step size. The default value is set to 0.005 
#' @param M The number of internal steps. The default value is set to 40
#' @param iter The number of iterations of the Lorenz model you want to run. 
#' The default value is set to 500
#' @param burn The amount of burn-in you want to run based on your starting point.
#' The default value is set to 100
#' @param order The order of the RK solution method to use.
#' The default value is 1
#'
#' @return This function returns a matrix containg iter serial simulations of 
#' Lorenz chaos model updates.
#' @export
#'
#' @examples
#' X0 = rnorm(960)
#' Lorenz04M2Sim(X0, 10, 32, 0.005, 40)
#' 
Lorenz04M2Sim = function(X0, F_Lor, K_Lor, dt, M, iter = 500, burn = 100, order = 1){
  N_Lor = length(X0)
  
  if(!is.vector(X0)){ stop("X0 should be a vector")}
  
  Xout = Lorenz04M2SimCpp(X0, F_Lor, K_Lor, dt, M, iter, burn, order)
  return(Xout)
}

