


#' Create a function that returns an object of class veryCoolMethod
#'
#' @param mydata Dataframe with at least two columns
#' @param ... Additional arguments
#'
#' @return A list
#' @export
#'
#' @examples
#' sample_object <- veryCoolMethod(cars)
#' sample_object
#' predict(sample_object)
veryCoolMethod <- function(mydata, ...) {

    model <- structure(
      list(x = mydata[, -1], y = mydata[, 1]),
      class = "veryCoolMethod")
    return(model)
}


#' Create a method for function print for class veryCoolMethod
#'
#' @param veryCoolMethodObject 
#'
#' @return Predicted valued for veryCoolMethod
#' @export
#'
#' @examples
#' sample_object <- veryCoolMethod(cars)
#' sample_object
#' predict(sample_object)
predict.veryCoolMethod = function(veryCoolMethodObject) {

  idx <- seq_along(veryCoolMethodObject$y)
  idx_mod <- (idx %% 24) + 1

  x <- c(
    " ", "n", "e", "v", "e", "r", " ", "g", "o", "n", "n", "a", " ", "g",
    "i", "v", "e", " ", "y", "o", "u", " ", "u", "p")

  return(x[idx_mod])
}

# 
# sample_object <- veryCoolMethod(cars)
# sample_object
# predict(sample_object)

