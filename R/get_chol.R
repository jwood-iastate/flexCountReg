#' Create Cholesky Matrix from Parameters in a Random Parameter Model
#'
#' This function creates a Cholesky matrix from parameters in a random parameter model.
#'
#' @param pars A numeric vector of parameters to include in the Cholesky matrix.
#' @param Nvars An integer specifying the number of random variables for the Cholesky matrix.
#'
#' @return A matrix with \code{Nvars} columns and rows containing the upper Cholesky matrix.
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
get_chol <- function(pars, Nvars) {
  sourceCpp("get_chol.cpp")
  get_chol_cpp(pars, Nvars)
}
