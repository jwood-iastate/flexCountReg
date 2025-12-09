#' Generate a covariance matrix using a correlation matrix and vector of
#' standard deviations
#'
#' @param C A correlation matrix.
#' @param S A vector of standard deviations.
#'
#' @return A covariance matrix
#' @examples
#' C <- matrix(c(1,-0.3,0.7,-0.3,1,-0.2,0.7,-0.2,1), 3, 3)
#' S <- c(0.5, 2, 1.25)
#' cor2cov(C,S)
#' @export
cor2cov <- function(C, S) {
  if (!is.matrix(C)) {
    warning("C must be a matrix") 
    return(NULL) # Exit to prevent crash
  }
  if (!is.vector(S)){
    warning("S must be a vector") 
    return(NULL) # Exit to prevent crash
  } 
  if (length(S) != nrow(C)) {
    msg <- paste0("S must have the same length (", length(S), 
                  ") as the number of rows/columns in C (", nrow(C), ")")
    warning(msg) 
    return(NULL) # Exit to prevent crash
  }
  sweep(sweep(C, 1, S, "*"), 2, S, "*")
}
