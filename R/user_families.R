# User-defined count families: new_count_family(), count_param() and
# check_count_family(), and their registration in the family registry.
#
# A user family is an object of class "fcr_count_family". Estimators accept
# it as `family`; it is then registered under a key that includes a hash of
# its definition ("user:<name>:<hash>"), so two families with the same name
# never collide, and the fit stores the object (model$family_spec). When a
# saved fit is used in a new R session, the family is registered again from
# the fit the first time it is needed (.fcr_recover_family()).
#
# Internally every family has up to two distribution parameters, each the
# exponential of its linear predictor (see families.R). A user parameter
# with link g is mapped to that slot as par = g^{-1}(log(slot)), so the
# linear predictor of `dis_param_formula_1` (or `_2`) is g(par).

# Registered user families: key -> registry entry.
.fcr_user_families <- new.env(parent = emptyenv())

#' A distribution parameter of a user-defined count family
#'
#' @param link the link of the parameter's linear predictor: \code{"log"}
#'   (default; a positive parameter), \code{"identity"} (any real value) or
#'   \code{"logit"} (a parameter in (0, 1)).
#' @param start the starting value of the linear predictor (on the link
#'   scale).
#' @returns An object of class \code{"fcr_count_param"}.
#' @seealso \code{\link{new_count_family}}
#' @export
count_param <- function(link = c("log", "identity", "logit"), start = 0) {
  link <- match.arg(link)
  if (!is.numeric(start) || length(start) != 1L || !is.finite(start)) {
    stop("`start` must be a single finite number.", call. = FALSE)
  }
  structure(list(link = link, start = start), class = "fcr_count_param")
}

#' Define a count distribution for the regression functions
#'
#' Creates a count family from its log probability mass function, for use
#' as \code{family} in \code{\link{countreg}()}, \code{\link{countreg_rp}()}
#' and \code{\link{countreg_mix}()}. Only \code{logpmf} is required: the
#' mean, variance, CDF and random draws are computed from it when they are
#' not given.
#'
#' @param name a short name for the family (used in printed output).
#' @param logpmf \code{function(y, mu, par)} returning
#'   \eqn{\log P(Y = y)} for vectors \code{y} and \code{mu} and a named list
#'   \code{par} of parameter vectors of the same length. Invalid parameter
#'   values should give \code{-Inf} or \code{NaN} (errors are also caught).
#' @param params a named list of up to two \code{\link{count_param}()}
#'   objects. The first is modelled by \code{dis_param_formula_1}, the second
#'   by \code{dis_param_formula_2}. The coefficients are labelled
#'   \code{ln(name)}, \code{logit(name)} or \code{name}, by link.
#' @param mean_is_mu \code{TRUE} (default) if \eqn{\mu = \exp(X\beta)} is
#'   the mean of the distribution. Otherwise give \code{mean}, or the mean
#'   is computed by summing the PMF.
#' @param mean,variance optional \code{function(mu, par)}: the mean and
#'   variance. Without them they are computed by summing the PMF.
#' @param cdf optional \code{function(q, mu, par)}: \eqn{P(Y \le q)}. If it
#'   has arguments \code{lower.tail} and \code{log.p} they are used, which is
#'   more accurate in the tails. Without it the PMF is summed.
#' @param rand optional \code{function(n, mu, par)}: random draws. Without
#'   it draws are made by inverting the CDF.
#' @param valid optional \code{function(mu, par)}: \code{TRUE} where the
#'   parameters define a distribution. Without it every finite value is
#'   accepted and \code{logpmf} must signal invalid values.
#'
#' @details The family supports everything that needs only probabilities:
#'   fixed and random-parameter models, zero-inflated and hurdle models,
#'   finite mixtures, predictions of probabilities, quantile residuals,
#'   rootograms and simulation. Underreporting (binomial thinning) is not
#'   available, because a user family is not known to be closed under
#'   thinning.
#'
#'   The fit stores the family, so a saved model works in a new R session.
#'   Two families with the same \code{name} but different definitions are
#'   kept apart.
#'
#'   Use \code{\link{check_count_family}()} to test a new family before
#'   fitting.
#'
#' @returns An object of class \code{"fcr_count_family"}.
#' @examples
#' # Negative binomial with size theta (variance mu + mu^2 / theta)
#' nb_theta <- new_count_family(
#'   name = "NB-theta",
#'   params = list(theta = count_param("log")),
#'   logpmf = function(y, mu, par)
#'     dnbinom(y, size = par$theta, mu = mu, log = TRUE),
#'   variance = function(mu, par) mu + mu^2 / par$theta)
#' \donttest{
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
#'               family = nb_theta, method = "BFGS")
#' summary(m)
#' }
#' @export
new_count_family <- function(name, logpmf, params = list(),
                             mean_is_mu = TRUE, mean = NULL,
                             variance = NULL, cdf = NULL, rand = NULL,
                             valid = NULL) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be a single non-empty string.", call. = FALSE)
  }
  if (!is.function(logpmf)) {
    stop("`logpmf` must be a function(y, mu, par).", call. = FALSE)
  }
  if (!is.list(params) || length(params) > 2L ||
      (length(params) && (is.null(names(params)) ||
                          any(!nzchar(names(params))) ||
                          anyDuplicated(names(params))))) {
    stop("`params` must be a named list of at most two count_param() ",
         "objects (more parameters are not yet supported).", call. = FALSE)
  }
  if (!all(vapply(params, inherits, logical(1), "fcr_count_param"))) {
    stop("Each element of `params` must be made with count_param().",
         call. = FALSE)
  }
  hooks <- list(mean = mean, variance = variance, cdf = cdf, rand = rand,
                valid = valid)
  for (h in names(hooks)) {
    if (!is.null(hooks[[h]]) && !is.function(hooks[[h]])) {
      stop("`", h, "` must be a function or NULL.", call. = FALSE)
    }
  }
  if (!isTRUE(mean_is_mu) && !isFALSE(mean_is_mu)) {
    stop("`mean_is_mu` must be TRUE or FALSE.", call. = FALSE)
  }
  fam <- c(list(name = name, logpmf = logpmf, params = params,
                mean_is_mu = mean_is_mu), hooks)
  # The key identifies the definition, not only the name
  def <- lapply(fam, function(x) if (is.function(x)) deparse(x) else x)
  fam$key <- paste0("user:", name, ":",
                    substr(rlang::hash(def), 1, 12))
  class(fam) <- "fcr_count_family"
  fam
}

#' @export
print.fcr_count_family <- function(x, ...) {
  cat("Count family \"", x$name, "\"", sep = "")
  if (length(x$params)) {
    cat(" with parameter(s)",
        toString(paste0(names(x$params), " (",
                        vapply(x$params, `[[`, "", "link"), " link)")))
  }
  cat("\n  mean:", if (x$mean_is_mu) "mu" else if (!is.null(x$mean))
    "user function" else "summed from the PMF")
  cat("\n  supplied:", toString(c("logpmf", names(Filter(Negate(is.null),
    x[c("mean", "variance", "cdf", "rand", "valid")])))), "\n")
  invisible(x)
}

# ---- Mapping to the registry ------------------------------------------------

.fcr_param_label <- function(nm, link) {
  switch(link, log = paste0("ln(", nm, ")"), logit = paste0("logit(", nm, ")"),
         identity = nm)
}

.fcr_linkinv <- function(link) {
  switch(link, log = function(eta) exp(eta), logit = stats::plogis,
         identity = function(eta) eta)
}

# The user's named parameter list from the internal slots (alpha, sigma =
# exp of the linear predictors).
.fcr_user_par <- function(fam, alpha, sigma) {
  slots <- list(alpha, sigma)
  out <- list()
  for (i in seq_along(fam$params)) {
    eta <- log(slots[[i]])
    out[[names(fam$params)[i]]] <- .fcr_linkinv(fam$params[[i]]$link)(eta)
  }
  out
}

# log P(Y = y) of a user family, with invalid values (by `valid`, NaN, NA
# or an error) at -Inf; inputs recycled to a common length.
.fcr_user_logpmf <- function(fam, y, mu, par) {
  n <- max(length(y), length(mu), lengths(par))
  if (n == 0L) return(numeric(0))
  y <- rep_len(y, n)
  mu <- rep_len(mu, n)
  par <- lapply(par, rep_len, n)
  lp <- tryCatch(suppressWarnings(as.numeric(fam$logpmf(y, mu, par))),
                 error = function(e) rep(NA_real_, n))
  if (length(lp) != n) lp <- rep_len(lp, n)
  ok <- !is.na(mu) & mu > 0
  if (!is.null(fam$valid)) {
    v <- tryCatch(fam$valid(mu, par), error = function(e) rep(FALSE, n))
    ok <- ok & !is.na(v) & rep_len(v, n)
  }
  lp[!ok | is.na(lp)] <- -Inf
  lp
}

# Moments by summing the PMF to a tail mass below `tol` (each row).
.fcr_user_moments <- function(fam, mu, par, tol = 1e-12, max_y = 1e6) {
  n <- max(length(mu), lengths(par))
  mu <- rep_len(mu, n)
  par <- lapply(par, rep_len, n)
  s0 <- s1 <- s2 <- numeric(n)
  active <- which(!is.na(mu) & mu > 0)
  y <- 0
  while (length(active) && y <= max_y) {
    p <- exp(.fcr_user_logpmf(fam, y, mu[active],
                              lapply(par, `[`, active)))
    s0[active] <- s0[active] + p
    s1[active] <- s1[active] + y * p
    s2[active] <- s2[active] + y^2 * p
    done <- (1 - s0[active] < tol & y > mu[active]) | !is.finite(p)
    active <- active[!done]
    y <- y + 1
  }
  m <- s1 / s0
  v <- s2 / s0 - m^2
  bad <- is.na(mu) | !(mu > 0) | s0 < 0.5
  m[bad] <- NaN
  v[bad] <- NaN
  list(mean = m, var = pmax(v, 0))
}

# Random draws by inverting the CDF (summing the PMF).
.fcr_user_rand <- function(fam, n, mu, par, max_y = 1e6) {
  mu <- rep_len(mu, n)
  par <- lapply(par, rep_len, n)
  u <- stats::runif(n)
  out <- rep(NA_real_, n)
  cdf <- numeric(n)
  active <- which(!is.na(mu) & mu > 0)
  y <- 0
  while (length(active) && y <= max_y) {
    p <- exp(.fcr_user_logpmf(fam, y, mu[active],
                              lapply(par, `[`, active)))
    cdf[active] <- cdf[active] + p
    hit <- cdf[active] >= u[active]
    out[active[hit]] <- y
    active <- active[!hit]
    y <- y + 1
  }
  out
}

# The registry entry (see .fcr_family() in families.R) of a user family.
.fcr_user_entry <- function(fam) {
  pars <- function(alpha, sigma) .fcr_user_par(fam, alpha, sigma)
  labels <- lapply(seq_len(2), function(i) {
    if (i > length(fam$params)) return(NULL)
    .fcr_param_label(names(fam$params)[i], fam$params[[i]]$link)
  })
  starts <- vapply(fam$params, `[[`, 0, "start")
  cdf_hook <- fam$cdf
  cdf <- if (!is.null(cdf_hook)) {
    tail_args <- all(c("lower.tail", "log.p") %in% names(formals(cdf_hook)))
    function(q, predicted, alpha, sigma, lower.tail = TRUE) {
      par <- pars(alpha, sigma)
      if (tail_args) {
        return(cdf_hook(q, predicted, par, lower.tail = lower.tail,
                        log.p = TRUE))
      }
      p <- cdf_hook(q, predicted, par)
      if (lower.tail) log(p) else log1p(-p)
    }
  }
  .fcr_family(
    params = labels,
    dens = function(y, predicted, alpha, sigma, ..., log = FALSE) {
      lp <- .fcr_user_logpmf(fam, y, predicted, pars(alpha, sigma))
      if (log) lp else exp(lp)
    },
    variance = function(mu, alpha, sigma) {
      par <- pars(alpha, sigma)
      if (!is.null(fam$variance)) fam$variance(mu, par) else
        .fcr_user_moments(fam, mu, par)$var
    },
    rand = function(n, mu, alpha, sigma) {
      par <- pars(alpha, sigma)
      if (!is.null(fam$rand)) fam$rand(n, mu, par) else
        .fcr_user_rand(fam, n, mu, par)
    },
    mean = if (!fam$mean_is_mu) function(mu, alpha, sigma) {
      par <- pars(alpha, sigma)
      if (!is.null(fam$mean)) fam$mean(mu, par) else
        .fcr_user_moments(fam, mu, par)$mean
    },
    cdf = cdf,
    start = if (length(starts)) unname(starts) else 0,
    start_rp = if (length(starts)) unname(starts),
    valid = if (!is.null(fam$valid)) function(mu, alpha, sigma)
      fam$valid(mu, pars(alpha, sigma)))
}

# Register a user family (once per definition); returns its key.
.fcr_register_family <- function(fam) {
  if (is.null(.fcr_user_families[[fam$key]])) {
    entry <- .fcr_user_entry(fam)
    entry$user <- fam
    assign(fam$key, entry, envir = .fcr_user_families)
  }
  fam$key
}

# The family key for an estimator's `family` argument: a user family object
# (registered here), or a built-in family name in any case and spelling.
.fcr_resolve_family <- function(family) {
  if (inherits(family, "fcr_count_family")) {
    return(.fcr_register_family(family))
  }
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop("`family` must be a family name or a family made with ",
         "new_count_family().", call. = FALSE)
  }
  if (startsWith(family, "user:")) {
    .fcr_get_family(family)
    return(family)
  }
  key <- toupper(gsub("[^[:alnum:]]", "", family))
  if (!key %in% .fcr_families) {
    stop("Unknown `family` \"", family, "\". Must be one of: ",
         toString(.fcr_families), ", or a family made with ",
         "new_count_family().", call. = FALSE)
  }
  key
}

# The user family object of a fit or family (NULL for built-in families).
.fcr_family_spec <- function(x) {
  if (inherits(x, "fcr_count_family")) return(x)
  spec <- x$model$family_spec %||% x$family_spec
  if (inherits(spec, "fcr_count_family")) spec
}

# A user family that is not registered in this session (a fit saved in
# another one): find the fit or family object that holds it among the
# objects of the calling functions, and register it again.
.fcr_recover_family <- function(key) {
  if (!is.character(key) || length(key) != 1L || !startsWith(key, "user:")) {
    return(NULL)
  }
  holds <- function(obj) {
    if (is.null(obj) || !is.list(obj)) return(NULL)
    specs <- c(list(.fcr_family_spec(obj)),
               if (is.list(obj$model$family_spec)) obj$model$family_spec,
               if (is.list(obj$family_spec)) obj$family_spec)
    for (s in specs) {
      if (inherits(s, "fcr_count_family") && identical(s$key, key)) return(s)
    }
    NULL
  }
  for (fr in rev(sys.frames())) {
    nms <- ls(fr, all.names = TRUE)
    if (!length(nms)) next
    # Values already computed only: forcing a promise (an argument not yet
    # evaluated) or an active binding could run arbitrary code
    skip <- rlang::env_binding_are_lazy(fr, nms) |
      rlang::env_binding_are_active(fr, nms)
    for (nm in nms[!skip]) {
      obj <- tryCatch(get(nm, envir = fr, inherits = FALSE),
                      error = function(e) NULL)
      s <- tryCatch(holds(obj), error = function(e) NULL)
      if (!is.null(s)) {
        .fcr_register_family(s)
        return(.fcr_user_families[[key]])
      }
    }
  }
  NULL
}

# Display name of a family key.
.fcr_family_label <- function(family) {
  if (!is.character(family)) return(family)
  user <- !is.na(family) & startsWith(family, "user:")
  family[user] <- sub("^user:(.*):[^:]*$", "\\1", family[user])
  family
}

#' Check a user-defined count family
#'
#' Runs numerical checks of a family made with \code{\link{new_count_family}()}
#' over a grid of means and parameter values: that the probabilities sum to
#' one, that the mean is \eqn{\mu} (when \code{mean_is_mu = TRUE}) or matches
#' the \code{mean} function, that the \code{variance}, \code{cdf} and
#' \code{rand} functions (when given) agree with the PMF, that inputs are
#' recycled, and that missing and invalid values give \code{-Inf} rather than
#' errors.
#'
#' @param family a family from \code{new_count_family()}.
#' @param mu means to check.
#' @param par a named list with values of each parameter (on the parameter's
#'   own scale); every combination with \code{mu} is checked.
#' @param tol tolerance for the comparisons.
#' @param nsim draws per grid point to check \code{rand} (0 to skip).
#'
#' @returns A data frame with one row per check (\code{check}, \code{ok},
#'   \code{detail}), printed, and returned invisibly. A warning lists the
#'   checks that failed.
#' @examples
#' nb_theta <- new_count_family(
#'   name = "NB-theta",
#'   params = list(theta = count_param("log")),
#'   logpmf = function(y, mu, par)
#'     dnbinom(y, size = par$theta, mu = mu, log = TRUE),
#'   variance = function(mu, par) mu + mu^2 / par$theta)
#' check_count_family(nb_theta, mu = c(0.1, 2, 30),
#'                    par = list(theta = c(0.5, 5)))
#' @export
check_count_family <- function(family, mu = c(0.05, 1, 20), par = list(),
                               tol = 1e-6, nsim = 2000) {
  if (!inherits(family, "fcr_count_family")) {
    stop("`family` must be made with new_count_family().", call. = FALSE)
  }
  miss <- setdiff(names(family$params), names(par))
  if (length(miss)) {
    stop("Give values of the parameter(s) ", toString(miss), " in `par`.",
         call. = FALSE)
  }
  grid <- expand.grid(c(list(mu = mu), par[names(family$params)]),
                      KEEP.OUT.ATTRS = FALSE)
  gp <- as.list(grid[names(family$params)])
  res <- list()
  add <- function(check, ok, detail) {
    res[[length(res) + 1L]] <<- data.frame(check = check, ok = isTRUE(ok),
                                            detail = detail)
  }
  mom <- .fcr_user_moments(family, grid$mu, gp, tol = 1e-12)
  # PMF sums (tail mass at the cutoff used for the moments)
  tot <- vapply(seq_len(nrow(grid)), function(i) {
    hi <- max(100, ceiling(grid$mu[i] + 50 * sqrt(max(mom$var[i], 1))))
    sum(exp(.fcr_user_logpmf(family, 0:hi, grid$mu[i],
                             lapply(gp, `[`, i))))
  }, numeric(1))
  add("PMF sums to 1", all(abs(tot - 1) < tol),
      sprintf("largest |sum - 1| = %.2g", max(abs(tot - 1))))
  if (family$mean_is_mu) {
    e <- abs(mom$mean - grid$mu) / pmax(grid$mu, 1)
    add("mean equals mu", all(e < tol),
        sprintf("largest relative error = %.2g", max(e)))
  } else if (!is.null(family$mean)) {
    e <- abs(family$mean(grid$mu, gp) - mom$mean) / pmax(abs(mom$mean), 1)
    add("mean() matches the PMF", all(e < tol),
        sprintf("largest relative error = %.2g", max(e)))
  }
  if (!is.null(family$variance)) {
    e <- abs(family$variance(grid$mu, gp) - mom$var) / pmax(mom$var, 1)
    add("variance() matches the PMF", all(e < tol),
        sprintf("largest relative error = %.2g", max(e)))
  }
  if (!is.null(family$cdf)) {
    e <- vapply(seq_len(nrow(grid)), function(i) {
      q <- 0:max(10, ceiling(3 * grid$mu[i]))
      pi <- lapply(gp, `[`, i)
      direct <- family$cdf(q, grid$mu[i], lapply(pi, rep_len, length(q)))
      if (isTRUE(all(c("lower.tail", "log.p") %in%
                     names(formals(family$cdf))))) {
        direct <- exp(family$cdf(q, grid$mu[i],
                                 lapply(pi, rep_len, length(q)),
                                 lower.tail = TRUE, log.p = TRUE))
      }
      summed <- cumsum(exp(.fcr_user_logpmf(family, q, grid$mu[i], pi)))
      max(abs(direct - summed))
    }, numeric(1))
    add("cdf() matches the summed PMF", all(e < tol),
        sprintf("largest difference = %.2g", max(e)))
  }
  if (!is.null(family$rand) && nsim > 0) {
    z <- vapply(seq_len(nrow(grid)), function(i) {
      pi <- lapply(gp, `[`, i)
      x <- family$rand(nsim, grid$mu[i], lapply(pi, rep_len, nsim))
      se <- sqrt(max(mom$var[i], 1e-12) / nsim)
      abs(mean(x) - mom$mean[i]) / se
    }, numeric(1))
    add("rand() has the PMF's mean", all(z < 4.5),
        sprintf("largest |z| = %.2g (of %d draws per point)", max(z), nsim))
  }
  # Recycling: a scalar parameter with a vector of counts
  i <- 1
  lp1 <- .fcr_user_logpmf(family, 0:3, grid$mu[i], lapply(gp, `[`, i))
  lp2 <- vapply(0:3, function(y)
    .fcr_user_logpmf(family, y, grid$mu[i], lapply(gp, `[`, i)), 0)
  add("vectorized over y", isTRUE(all.equal(lp1, lp2)),
      if (isTRUE(all.equal(lp1, lp2))) "ok" else "elementwise values differ")
  raw <- tryCatch(family$logpmf(c(0, 1), c(1, NA), lapply(gp, function(v)
    rep_len(v[1], 2))), error = function(e) e)
  add("NA mean handled", !inherits(raw, "error"),
      if (inherits(raw, "error")) conditionMessage(raw) else "no error")
  bad <- tryCatch(suppressWarnings(
    family$logpmf(c(1, 1), c(-1, 1), lapply(gp, function(v) rep_len(v[1], 2)))),
    error = function(e) e)
  add("invalid mean gives no error", !inherits(bad, "error"),
      if (inherits(bad, "error")) conditionMessage(bad) else
        "no error (the row is given probability 0)")
  out <- do.call(rbind, res)
  print(out, row.names = FALSE)
  if (!all(out$ok)) {
    warning("Failed checks: ", toString(out$check[!out$ok]), call. = FALSE)
  }
  invisible(out)
}
