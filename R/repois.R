#' Random-Effects Poisson Panel Model
#'
#' Estimate a random-effects Poisson panel model for clustered event counts.
#' This function follows the gamma random-effects representation discussed by
#' Guo (1996), which yields the same likelihood form used for negative
#' multinomial clustered-count models.
#'
#' The function is intended for panel or clustered count data where repeated
#' observations share an unobserved cluster-specific effect.
#'
#' @param formula An \code{R} formula.
#' @param group_var The grouping variable(s) for the random effects (for
#'   example, an individual ID, corridor ID, or other panel identifier).
#' @param data A data frame containing all variables used in \code{formula}.
#' @param method A method used for maximum-likelihood optimization. The default
#'   is \code{"NM"}.
#' @param max.iters The maximum number of optimizer iterations. Default is
#'   \code{1000}.
#' @param print.level Integer controlling optimizer verbosity. Default is
#'   \code{0}.
#' @param bootstraps Optional integer specifying the number of bootstrap
#'   replications for standard errors. If \code{NULL}, no bootstrap is used.
#' @param offset Optional offset term provided as a string.
#'
#' @details
#' This is a thin user-facing wrapper around the existing random-effects count
#' regression engine used in \code{renb()}. The parameterization is presented as
#' a Poisson panel model with a cluster-level random effect, following the
#' derivation in Guo (1996). The fitted object is returned as class
#' \code{countreg}, so it remains compatible with existing methods such as
#' \code{summary()}, \code{predict()}, \code{regCompTable()}, and
#' \code{regCompTest()}.
#'
#' @return An object of class \code{countreg} with components including:
#' \describe{
#'   \item{\code{model}}{The fitted model object.}
#'   \item{\code{data}}{The data frame used to fit the model.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{formula}}{The formula used to fit the model.}
#' }
#'
#' @references
#' Guo, G. (1996). Negative multinomial regression models for clustered event
#' counts. \emph{Sociological Methodology}, 26, 113--132.
#'
#' @seealso \code{\link{renb}}, \code{\link{predict.flexCountReg}},
#' \code{\link{summary.flexCountReg}}, \code{\link{regCompTable}},
#' \code{\link{regCompTest}}
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)
#'
#' repois.mod <- repois(
#'   Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
#'   data = washington_roads,
#'   offset = "lnlength",
#'   group_var = "ID",
#'   method = "NM",
#'   max.iters = 1000
#' )
#'
#' summary(repois.mod)
#' }
#'
#' @export
repois <- function(formula,
                   group_var,
                   data,
                   method = "NM",
                   max.iters = 1000,
                   print.level = 0,
                   bootstraps = NULL,
                   offset = NULL) {
  
  fit <- renb(
    formula = formula,
    group_var = group_var,
    data = data,
    method = method,
    max.iters = max.iters,
    print.level = print.level,
    bootstraps = bootstraps,
    offset = offset
  )
  
  if (!is.null(fit$call) && length(fit$call) >= 1L) {
    fit$call[[1L]] <- as.name("repois")
  }
  
  if (!is.null(fit$model) && is.list(fit$model)) {
    fit$model$model_name <- "Random Effects Poisson Panel Model"
    fit$model$model_family <- "PoissonRE"
    fit$model$derivation <- "Guo (1996)"
    fit$model$likelihood_form <- "Negative multinomial equivalent"
  }
  
  class(fit) <- unique(c("countreg", class(fit)))
  fit
}