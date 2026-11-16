#' Generate a covariance matrix using a correlation matrix and vector of
#' standard deviations
#'
#' @param C A correlation matrix.
#' @param S A vector of standard deviations.
#'
#' @returns A covariance matrix
#' @examples
#' C <- matrix(c(1,-0.3,0.7,-0.3,1,-0.2,0.7,-0.2,1), 3, 3)
#' S <- c(0.5, 2, 1.25)
#' cor2cov(C,S)
#' @export
cor2cov <- function(C, S) {
  if (!is.matrix(C) || !is.numeric(C) || nrow(C) != ncol(C)) {
    stop("`C` must be a square numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(S) || !is.null(dim(S))) {
    stop("`S` must be a numeric vector.", call. = FALSE)
  }
  if (length(S) != nrow(C)) {
    stop("`S` must have the same length (", length(S), ") as the number of ",
         "rows/columns in `C` (", nrow(C), ").", call. = FALSE)
  }
  tol <- sqrt(.Machine$double.eps)
  if (max(abs(C - t(C))) > tol) {
    stop("`C` must be symmetric.", call. = FALSE)
  }
  if (any(abs(diag(C) - 1) > tol)) {
    stop("`C` must have a unit diagonal (it must be a correlation matrix).",
         call. = FALSE)
  }
  if (any(abs(C) > 1 + tol)) {
    stop("Correlations in `C` must be between -1 and 1.", call. = FALSE)
  }
  if (any(S < 0)) {
    stop("Standard deviations in `S` must be non-negative.", call. = FALSE)
  }
  sweep(sweep(C, 1, S, "*"), 2, S, "*")
}
