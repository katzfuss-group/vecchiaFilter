#' Lorenz Model N = 480
#' 
#' 500 simulations after a burn-in of 200 from a Lorenz04M2 Model with N = 480, F_Lor = 10, K_Lor = 16,
#' dt = 0.005, and M = 40. Figure 4(b) gives example of what this model is supposed to
#' look similar to, up to 960 locations.
#'  
#' @format A matrix with 500 columns for the simulations and 480 rows for the locations
#'
#'
#' @source \href{https://journals.ametsoc.org/doi/pdf/10.1175/JAS3430.1}{Lorenz (2005)}
#'
"Lorenz480"

#' Example State Ensemble Matrix
#'
#' Ensemble state matrix for a 5x5 state-space model with ensemble size N=50 for
#' examples.
#'
#' @format 25x50 matrix
#'
"x.mult.prior"

#' Example Observation Ensemble Matrix
#'
#' Ensemble observation matrix for a 5x5 state-space model with ensemble size N=50
#' for examples.
#'
#' @format 25x50 matrix
#'
"y.mult.i"

#' Example Observation Matrix
#'
#' Observation matrix for a 5x5 state-space model for examples.
#'
#' @format 25x25 matrix
#'
"H.mat"

#' Example Observation Noise Matrix
#'
#' Observation noise matrix for a 5x5 state-space model for examples.
#'
#' @format 25x25 matrix
#'
"tau.mat"

#' Example Spatial Location Matrix
#'
#' Spatial locations for a 5x5 state-space model for examples ordered using the 
#' max-min technique.
#'
#' @format 25x2 matrix
#'
"S2d.ord"
