# Evaluates a fit of an underreporting model on a nearly flat likelihood
# ridge, muffling the warnings whose appearance depends on where the
# optimizer stops on the ridge (and so on the platform's floating point):
# weak identification, an indefinite Hessian and non-convergence. Any other
# warning is passed on.
fcr_quiet_ridge <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl(paste0("weakly identified|Hessian is not negative definite|",
                       "did not report convergence"), msg)) {
        invokeRestart("muffleWarning")
      }
    })
}
