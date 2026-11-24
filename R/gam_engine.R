# Penalized likelihood estimation with smoothing-parameter selection, shared
# by the estimators that accept smooth terms.
#
# The objective is
#   l_p(theta) = l(theta) - 1/2 sum_m lambda_m theta' S_m theta,
# where each S_m is a smooth's penalty matrix embedded at that smooth's
# coefficients. For fixed lambda it is maximized by BFGS. The smoothing
# parameters are chosen by the generalized Fellner-Schall update (Wood and
# Fasiolo 2017),
#   lambda_m <- lambda_m [tr(S_lambda^- S_m) - tr(V S_m)] / (theta' S_m theta),
#   V = (H + S_lambda)^{-1},
# with H the negative Hessian of the unpenalized log-likelihood, which
# maximizes the Laplace-approximate restricted marginal likelihood. It needs
# only the Hessian at each fit.
#
# Inference is Bayesian: the covariance of the estimates is
# V_b = (H + S_lambda)^{-1}, and the effective degrees of freedom are the
# diagonal of F = V_b H, summed over each smooth's coefficients.

# Penalty terms on the full parameter vector: one per penalty matrix.
#   smooths  list of component smooth entries (see smooth_design.R), each
#            with `cols` relative to the component's design matrix
#   offsets  the position of each component's first coefficient minus 1
#   prefix   coefficient-name prefix of each component (for labels)
.fcr_pen_terms <- function(smooths, offsets, components) {
  out <- list()
  for (cmp in names(smooths)) {
    for (e in smooths[[cmp]]) {
      idx <- offsets[[cmp]] + e$cols
      for (m in seq_along(e$S)) {
        out[[length(out) + 1L]] <- list(
          idx = idx, S = e$S[[m]], smooth = e$label, component = cmp,
          label = components[[cmp]] %||% cmp, which = m)
      }
    }
  }
  out
}

# S_lambda = sum_m lambda_m S_m, as a P x P matrix.
.fcr_pen_matrix <- function(terms, lambda, P) {
  S <- matrix(0, P, P)
  for (m in seq_along(terms)) {
    i <- terms[[m]]$idx
    S[i, i] <- S[i, i] + lambda[m] * terms[[m]]$S
  }
  S
}

# Maximize the penalized log-likelihood for fixed lambda.
.fcr_pen_optimize <- function(fn, grad, start, terms, lambda, method,
                              iterlim) {
  P <- length(start)
  Sl <- .fcr_pen_matrix(terms, lambda, P)
  obj <- function(p) {
    v <- fn(p)
    if (any(is.na(v))) return(NA_real_)
    sum(v) - 0.5 * drop(crossprod(p, Sl %*% p))
  }
  gr <- if (!is.null(grad)) function(p) {
    g <- grad(p)
    g <- if (is.matrix(g)) colSums(g) else as.vector(g)
    g - drop(Sl %*% p)
  }
  if (method %in% c("BHHH", "NR", "BFGSR")) method <- "BFGS"
  fit <- maxLik::maxLik(obj, grad = gr, start = start, method = method,
                        finalHessian = FALSE,
                        control = list(iterlim = iterlim, reltol = 1e-10))
  names(fit$estimate) <- names(start)
  # Newton polish with the penalized Hessian: BFGS stops on a relative
  # change of the objective, which leaves directions of low curvature
  # (weakly penalized smooths) short of the optimum
  if (!is.null(gr)) {
    theta <- fit$estimate
    f0 <- obj(theta)
    for (it in 1:5) {
      g <- gr(theta)
      Hl <- .fcr_hessian_step(fn, theta, grad, 1e-4 * pmax(1, abs(theta)))
      if (is.null(Hl) || !all(is.finite(g))) break
      A <- -Hl + Sl
      step <- tryCatch(solve(A, g), error = function(e) NULL)
      if (is.null(step) || !all(is.finite(step))) break
      s <- 1
      repeat {
        cand <- theta + s * step
        f1 <- obj(cand)
        if (is.finite(f1) && f1 >= f0 - 1e-12 * abs(f0)) break
        s <- s / 2
        if (s < 1e-8) break
      }
      if (s < 1e-8) break
      theta <- cand
      done <- abs(f1 - f0) < 1e-12 * max(1, abs(f0)) &&
        max(abs(s * step)) < 1e-9
      f0 <- f1
      if (done) break
    }
    names(theta) <- names(start)
    fit$estimate <- theta
    fit$maximum <- f0
    fit$gradient <- gr(theta)
  }
  fit
}

# One generalized Fellner-Schall update of the smoothing parameters, at
# estimates theta with H the negative Hessian of the unpenalized
# log-likelihood.
.fcr_fs_step <- function(theta, H, terms, lambda) {
  P <- length(theta)
  Sl <- .fcr_pen_matrix(terms, lambda, P)
  V <- tryCatch(solve(H + Sl), error = function(e) MASS::ginv(H + Sl))
  # S_lambda is block diagonal by smooth (only the penalties of one
  # smooth, as in te(), share coefficients), so its generalized inverse
  # is taken block by block: one inverse of the whole matrix would treat
  # a weakly penalized block as zero next to a heavily penalized one
  Slg <- matrix(0, P, P)
  for (blk in unique(lapply(terms, `[[`, "idx"))) {
    Slg[blk, blk] <- MASS::ginv(Sl[blk, blk, drop = FALSE])
  }
  new <- lambda
  for (m in seq_along(terms)) {
    i <- terms[[m]]$idx
    Sm <- terms[[m]]$S
    num <- sum(Slg[i, i] * Sm) - sum(V[i, i] * Sm)
    den <- drop(crossprod(theta[i], Sm %*% theta[i]))
    new[m] <- lambda[m] * max(num, 1e-12) / max(den, 1e-12)
  }
  pmin(pmax(new, 1e-8), 1e10)
}

# Fit with smoothing-parameter selection (or fixed `lambda`).
# Returns the penalized fit and the quantities for inference.
.fcr_gam_fit <- function(fn, grad, start, terms, method = "BFGS",
                         iterlim = 1000, lambda = NULL, fs_maxit = 200,
                         fs_tol = 1e-4) {
  P <- length(start)
  nm <- names(start)
  fixed <- !is.null(lambda)
  if (!fixed) lambda <- rep(1, length(terms))
  if (length(lambda) != length(terms)) {
    stop("`lambda` must have one value per smoothing parameter (",
         length(terms), ").", call. = FALSE)
  }
  theta <- start
  hist <- NULL
  conv <- fixed
  it <- 0L
  repeat {
    it <- it + 1L
    fit <- .fcr_pen_optimize(fn, grad, theta, terms, lambda, method,
                             iterlim)
    theta <- fit$estimate
    if (fixed || it > fs_maxit) break
    h <- 1e-4 * pmax(1, abs(theta))
    Hl <- .fcr_hessian_step(fn, theta, grad, h)
    if (is.null(Hl)) break
    new <- .fcr_fs_step(theta, -Hl, terms, lambda)
    hist <- rbind(hist, new)
    done <- max(abs(log(new) - log(lambda))) < fs_tol * 10 ||
      max(abs(new - lambda) / lambda) < fs_tol
    lambda <- new
    if (done) {
      conv <- TRUE
      fit <- .fcr_pen_optimize(fn, grad, theta, terms, lambda, method,
                               iterlim)
      theta <- fit$estimate
      break
    }
  }
  # Final Hessian (adaptive) and Bayesian covariance
  Hf <- .fcr_final_hessian(fn, theta, grad)
  if (is.null(Hf)) {
    Hf <- .fcr_hessian_step(fn, theta, grad, 1e-4 * pmax(1, abs(theta)))
  }
  H <- -as.matrix(Hf)
  H <- (H + t(H)) / 2
  Sl <- .fcr_pen_matrix(terms, lambda, P)
  A <- H + Sl
  Vb <- tryCatch(solve(A), error = function(e) MASS::ginv(A))
  Vb <- (Vb + t(Vb)) / 2
  edf <- diag(Vb %*% H)
  names(edf) <- nm
  dimnames(Vb) <- list(nm, nm)
  # Negative-definite "Hessian" of the penalized objective, so that the
  # Hessian-based covariance of the fit is V_b
  hess <- -A
  dimnames(hess) <- list(nm, nm)
  list(fit = fit, theta = theta, lambda = lambda, H = H, Sl = Sl, Vb = Vb,
       edf = edf, hessian = hess, iterations = it, converged = conv,
       lambda_history = hist, terms = terms,
       accuracy = attr(Hf, "accuracy"))
}

# Per-smooth summaries: label, component, EDF, number of coefficients,
# smoothing parameter(s) and an approximate test that the smooth is zero.
# The test follows Wood (2013) in the space of the smooth's fitted values:
# with X_j its basis at the data rows (X_j = Q R), the statistic is
# f' V_f^{r-} f computed as (R b)' (R V R')^{r-} (R b), a rank-r pseudo-
# inverse with r the EDF rounded up, referred to a chi-squared with EDF
# degrees of freedom. It is approximate (smoothing-parameter uncertainty is
# ignored). `bases` holds X_j for each smooth (named "component\rlabel").
.fcr_smooth_table <- function(gam, bases = NULL) {
  terms <- gam$terms
  if (!length(terms)) return(NULL)
  key <- vapply(terms, function(t) paste(t$component, t$smooth, sep = "\r"),
                "")
  rows <- lapply(unique(key), function(k) {
    ts <- terms[key == k]
    idx <- ts[[1]]$idx
    edf <- sum(gam$edf[idx])
    b <- gam$theta[idx]
    V <- gam$Vb[idx, idx, drop = FALSE]
    Xj <- bases[[k]]
    if (!is.null(Xj)) {
      R <- qr.R(qr(Xj))
      b <- drop(R %*% b)
      V <- R %*% V %*% t(R)
    }
    r <- min(length(idx), max(1, ceiling(edf - 0.05)))
    e <- eigen((V + t(V)) / 2, symmetric = TRUE)
    keep <- seq_len(r)
    z <- crossprod(e$vectors[, keep, drop = FALSE], b)
    stat <- sum(z^2 / pmax(e$values[keep], .Machine$double.eps))
    data.frame(
      parameter = ts[[1]]$label, term = ts[[1]]$smooth,
      edf = edf, k = length(idx),
      lambda = paste(signif(gam$lambda[key == k], 3), collapse = ", "),
      chisq = stat, p.value = stats::pchisq(stat, df = max(edf, 1),
                                            lower.tail = FALSE),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Smooth summaries and plotting information after a penalized fit.
#   ents  list: component -> its smooth entries (see smooth_design.R)
#   mats  list: component -> its design matrix at the estimation rows
# Components named "het_mean_<k>" use the entries and matrix of "het_mean".
# Returns list(table, info), with `info` one element per smooth: term,
# parameter label, component, coefficient indices, mgcv smooth object.
.fcr_smooth_outputs <- function(gam, ents, mats) {
  bases <- list()
  info <- list()
  base_cmp <- function(cmp) if (startsWith(cmp, "het_mean_")) "het_mean" else
    cmp
  for (t in gam$terms) {
    key <- paste(t$component, t$smooth, sep = "\r")
    if (!is.null(bases[[key]])) next
    cmp <- t$base %||% base_cmp(t$component)
    es <- ents[[cmp]]
    e <- es[[which(vapply(es, `[[`, "", "label") == t$smooth)[1]]]
    bases[[key]] <- mats[[cmp]][, e$cols, drop = FALSE]
    info[[key]] <- list(term = t$smooth, parameter = t$label,
                        component = t$component, base = cmp, idx = t$idx,
                        sm = e$sm, vars = e$vars, by = e$by)
  }
  list(table = .fcr_smooth_table(gam, bases), info = unname(info))
}

# Does a prepare_model_data() result have smooth terms?
.fcr_prep_has_smooth <- function(prep) {
  length(prep$smooth$main) > 0 || any(lengths(prep$smooth$aux) > 0)
}

# A penalized fit with the fields of a maxLik fit, for the estimators whose
# covariance and summaries work from maxLik objects.
#   ents     component -> smooth entries; offsets, labels, mats as in
#            .fcr_pen_terms() and .fcr_smooth_outputs()
#   fn, grad the (weighted) log-likelihood contributions and their scores
# The returned fit has the penalized estimates, the unpenalized maximum,
# the Hessian of the penalized objective (so the Hessian-based covariance
# is the Bayesian V_b), the unpenalized scores (for robust covariances) and
# the smooth summaries.
.fcr_smooth_fit <- function(ents, fn, grad, start, method, max.iters, sp,
                            offsets, labels, mats) {
  pen <- .fcr_pen_terms(ents, offsets, labels)
  if (method %in% c("NM", "SANN", "BHHH", "NR", "BFGSR")) method <- "BFGS"
  gam <- .fcr_gam_fit(fn, grad, start, pen, method = method,
                      iterlim = max.iters, lambda = sp)
  fit <- gam$fit
  fit$estimate <- gam$theta
  fit$maximum <- sum(fn(gam$theta))
  fit$hessian <- gam$hessian
  fit$gradientObs <- if (!is.null(grad)) grad(gam$theta) else
    maxLik::numericGradient(fn, gam$theta)
  colnames(fit$gradientObs) <- names(start)
  fit$hessian_accuracy <- gam$accuracy
  so <- .fcr_smooth_outputs(gam, ents, mats)
  fit$smooth_table <- so$table
  fit$smooth_info <- so$info
  fit$smooth_lambda <- gam$lambda
  fit$smooth_converged <- gam$converged
  fit$smooth_iterations <- gam$iterations
  fit$smooth_coefs <- unique(unlist(lapply(gam$terms, `[[`, "idx")))
  fit$edf <- gam$edf
  fit$df_total <- sum(gam$edf)
  fit$code <- if (isTRUE(gam$converged) || !is.null(sp)) 0L else 1L
  fit$method <- method
  if (!isTRUE(gam$converged) && is.null(sp)) {
    warning("The smoothing parameters did not converge in ",
            gam$iterations - 1, " Fellner-Schall iterations.", call. = FALSE)
  }
  fit
}

# Bootstrap refits would re-select the smoothing parameters; not yet.
.fcr_smooth_no_boot <- function(stderr) {
  if (identical(stderr, "boot")) {
    stop("Bootstrap standard errors are not yet available with smooth ",
         "terms; use stderr = \"normal\" (Bayesian) or \"robust\".",
         call. = FALSE)
  }
}
