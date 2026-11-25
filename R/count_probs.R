# Count probabilities P(Y = k) and intervals for every estimator.
#
#   .fcr_count_set()      the counts requested by `counts` or `max_count`;
#   .fcr_prob_table()     the matrix returned by predict(type = "prob"), with
#                         the optional lower and upper tail columns;
#   .fcr_count_pmf()      P(Y = k) of any univariate model for the rows of a
#                         data frame (a list of matrices, one per outcome or
#                         category, for countreg_mv() and countreg_alloc());
#   .fcr_general_interval()  confidence and prediction intervals for the
#                         estimators other than countreg().

# The counts at which probabilities are evaluated: `counts` (distinct whole
# numbers, in increasing order) or 0:max_count.
.fcr_count_set <- function(counts = NULL, max_count = NULL, default = 20) {
  whole <- function(x) {
    is.numeric(x) && length(x) && !anyNA(x) && all(is.finite(x)) &&
      all(x >= 0) && all(x == round(x))
  }
  if (!is.null(counts)) {
    if (!is.null(max_count)) {
      stop("Give either `counts` or `max_count`, not both.", call. = FALSE)
    }
    if (!whole(counts)) {
      stop("`counts` must be non-negative whole numbers, such as 0:10 or ",
           "5:8.", call. = FALSE)
    }
    return(as.numeric(sort(unique(counts))))
  }
  if (is.null(max_count)) max_count <- default
  if (!whole(max_count) || length(max_count) != 1L) {
    stop("`max_count` must be a single non-negative whole number.",
         call. = FALSE)
  }
  as.numeric(seq.int(0, max_count))
}

# `pmf(k)` gives the n x length(k) matrix of P(Y = k). `upper` (optional)
# gives P(Y > q) directly where a family has an accurate upper tail.
.fcr_prob_table <- function(pmf, counts, overflow = FALSE, underflow = FALSE,
                            upper = NULL) {
  overflow <- isTRUE(overflow)
  underflow <- isTRUE(underflow)
  lo <- min(counts)
  hi <- max(counts)
  need <- if (overflow && is.null(upper)) union(seq.int(0, hi), counts) else
    if (underflow) union(seq_len(lo) - 1, counts) else counts
  need <- sort(need)
  A <- pmf(need)
  P <- A[, match(counts, need), drop = FALSE]
  colnames(P) <- format(counts, scientific = FALSE, trim = TRUE)
  if (underflow) {
    P <- cbind(rowSums(A[, need < lo, drop = FALSE]), P)
    colnames(P)[1L] <- paste0("<", lo)
  }
  if (overflow) {
    up <- if (!is.null(upper)) upper(hi) else pmax(0, 1 - rowSums(A))
    P <- cbind(P, up)
    colnames(P)[ncol(P)] <- paste0(">", hi)
  }
  rownames(P) <- NULL
  P
}

# Weighted average over nodes of the density of `family`: `E` is the n x Q
# matrix of log mean parameters, `w` the Q node weights.
.fcr_mix_pmf <- function(family, E, w, alpha, sigma, counts, dens_args) {
  n <- nrow(E)
  Q <- ncol(E)
  out <- matrix(NA_real_, n, length(counts))
  ok <- stats::complete.cases(E)
  if (!is.null(alpha)) ok <- ok & !is.na(alpha)
  if (!is.null(sigma)) ok <- ok & !is.na(sigma)
  if (!any(ok)) return(out)
  no <- sum(ok)
  dens <- get_probFunc(family)
  size <- max(1L, min(Q, floor(1e6 / no)))
  acc <- matrix(0, no, length(counts))
  for (st in seq.int(1L, Q, by = size)) {
    q <- st:min(Q, st + size - 1L)
    nq <- length(q)
    mu <- exp(as.vector(E[ok, q, drop = FALSE]))
    a <- if (!is.null(alpha)) rep(alpha[ok], nq)
    s <- if (!is.null(sigma)) rep(sigma[ok], nq)
    zero <- mu == 0     # an all-zero group: the count is 0
    for (j in seq_along(counts)) {
      p <- numeric(length(mu))
      p[zero] <- as.numeric(counts[j] == 0)
      if (any(!zero)) {
        p[!zero] <- do.call(dens, c(
          list(y = rep(counts[j], sum(!zero)), predicted = mu[!zero],
               alpha = a[!zero], sigma = s[!zero]), dens_args,
          list(log = FALSE)))
      }
      acc[, j] <- acc[, j] + as.vector(matrix(p, no, nq) %*% w[q])
    }
  }
  out[ok, ] <- acc
  out
}

.fcr_dens_args <- function(ndraws = 1500, engine = "quadrature") {
  h <- randtoolbox::halton(ndraws)
  list(haltons = h, normed_haltons = stats::qnorm(h), engine = engine)
}

# P(Y = k), one row per row of `data`, at the model's current estimates.
.fcr_count_pmf <- function(object, data, counts, ...) {
  if (inherits(object, "countreg_mix")) {
    r <- .fcr_fm_rows(object, data)
    return(.fcr_fm_prob_rows(r, exp(r$logpi), counts))
  }
  if (inherits(object, "countreg_mv")) {
    return(.fcr_mv_prob(object, data, counts, ...))
  }
  if (inherits(object, "countreg_alloc")) {
    return(.fcr_alloc_pmf(object, data, counts, ...))
  }
  if (inherits(object, "countreg_hier")) {
    return(.fcr_hier_pmf(object, data, counts, ...))
  }
  model <- object$model
  data <- as.data.frame(data)
  if (!is.null(model$cre)) data <- .fcr_cre_newdata(model$cre, data)
  type <- object$modelType %||% model$modelType %||% "countreg"
  switch(type,
    countreg.rp = .fcr_rp_pmf(object, data, counts, ...),
    countreg_fe = .fcr_fe_pmf(object, data, counts, ...),
    RENB = ,
    REPoisson = ,
    poisLindRE = .fcr_re_pmf(object, data, counts, type),
    NegMultinom = .fcr_nm_pmf(object, data, counts),
    .fcr_prob_matrix(object, data, counts = counts))
}

# Random-effects panel models: the marginal distribution of one count of a
# new unit, with the unit's effect integrated out in closed form.
.fcr_re_pmf <- function(object, data, counts, type) {
  mu <- as.vector(stats::predict(object, newdata = data))
  est <- object$model$estimate
  ok <- !is.na(mu)
  one <- switch(type,
    REPoisson = {
      a <- exp(est[[intersect(c("lnalpha", "ln(alpha)"), names(est))[1]]])
      function(k) stats::dnbinom(k, size = 1 / a, mu = mu[ok])
    },
    RENB = {
      a <- 1 + exp(est[["ln(a-1)"]])
      b <- exp(est[["ln(b)"]])
      lambda <- mu[ok] * (a - 1) / b
      function(k) exp(renb_ll(rep(k, length(lambda)), lambda, a, b,
                              seq_along(lambda)))
    },
    poisLindRE = {
      theta <- exp(est[["ln(theta)"]])
      function(k) dplind(rep(k, sum(ok)), mean = mu[ok], theta = theta)
    })
  out <- matrix(NA_real_, length(mu), length(counts))
  if (any(ok)) {
    out[ok, ] <- vapply(counts, one, numeric(sum(ok)))
  }
  out
}

# Negative multinomial: the frailty's shape belongs to the cluster (for NM1
# and NM-P it depends on the total mean of the cluster's rows in `data`;
# without the grouping column(s) every row is a cluster).
.fcr_nm_pmf <- function(object, data, counts) {
  model <- object$model
  fd <- forge_new_data(model$blueprint %||% .fcr_blueprint(object), data)
  n <- nrow(data)
  out <- matrix(NA_real_, n, length(counts))
  gv <- model$group_var
  has_g <- all(gv %in% names(data))
  ok <- stats::complete.cases(fd$X) & !is.na(fd$offset)
  if (has_g) ok <- ok & stats::complete.cases(data[, gv, drop = FALSE])
  for (m in fd$aux) if (!is.null(m)) ok <- ok & stats::complete.cases(m)
  if (!any(ok)) return(out)
  g <- if (!has_g) seq_len(sum(ok)) else
    as.integer(droplevels(interaction(data[ok, gv, drop = FALSE],
                                      drop = TRUE, sep = ":")))
  G <- max(g)
  one_col <- matrix(1, G, 1, dimnames = list(NULL, "(Intercept)"))
  clus <- function(M, nm) if (is.null(M)) one_col else
    suppressWarnings(.nm_cluster_design(M[ok, , drop = FALSE], g, G, nm))
  Zc <- clus(fd$aux$dis1, "dis_param_formula_1")
  Wc <- if (model$family == "NMP") clus(fd$aux$dis2, "dis_param_formula_2")
  prep <- list(y = rep(0, sum(ok)), X = fd$X[ok, , drop = FALSE],
               offset = fd$offset[ok], g = g, G = G, Zc = Zc, Wc = Wc,
               n_beta = ncol(fd$X), n_alpha = ncol(Zc),
               n_p = if (is.null(Wc)) 0L else ncol(Wc), Y = rep(0, G))
  cp <- .nm_components(as.numeric(model$estimate), prep, model$family)
  out[ok, ] <- vapply(counts, function(k)
    stats::dnbinom(k, size = cp$r[g], mu = cp$mu), numeric(sum(ok)))
  out
}

# Fixed effects: the family at the group's estimated effect.
.fcr_fe_pmf <- function(object, data, counts, newgroup = "none") {
  m <- object$model
  ok <- .fcr_fe_predict(object, data, newgroup, link = TRUE)
  b <- m$estimate
  pars <- get_params(m$family)
  n <- nrow(data)
  par <- function(nm) if (!is.null(nm)) rep(exp(b[[nm]]), n)
  .fcr_mix_pmf(m$family, matrix(ok, n, 1L), 1, par(pars[[1]]), par(pars[[2]]),
               counts, .fcr_dens_args())
}

# Random parameters: averaged over the draws of the random coefficients (a
# new unit). `hdraws`: uniform draws, one column per random parameter.
.fcr_rp_pmf <- function(object, data, counts, hdraws = NULL, method = NULL) {
  if (!is.null(method) && !identical(method, "Simulated")) {
    stop("Probabilities of countreg_rp() models are averages over the ",
         "distribution of the random parameters (method = \"Simulated\").",
         call. = FALSE)
  }
  P <- .fcr_rp_parts(object, data)
  w <- NULL
  if (is.null(hdraws)) {
    if (P$N_rand <= 4L &&
        (isTRUE(P$correlated) || all(P$rpardists %in% c("n", "ln")))) {
      # Gauss-Hermite nodes, as the uniforms of the draw functions
      nd <- .fcr_normal_nodes(P$N_rand)
      hdraws <- stats::pnorm(nd$Z)
      hdraws[] <- pmin(pmax(hdraws, 1e-300), 1 - 1e-15)
      w <- nd$w
    } else {
      hdraws <- .fcr_halton(max(2000, P$ndraws), P$N_rand, TRUE, 1L)
    }
  }
  xb <- generate_random_draws(
    hdraws = hdraws, random_coefs_means = P$random_coefs_means,
    rand_var_params = P$rand_var_params, rpardists = P$rpardists,
    rpar = P$rpar, X_rand = P$X_rand, het_mean_coefs = P$het_mean_coefs,
    X_het_mean = P$X_het_mean, het_var_coefs = P$het_var_coefs,
    X_het_var = P$X_het_var, correlated = isTRUE(P$correlated),
    het_form = P$het_form)$xb_rand_mat
  m <- object$model
  .fcr_mix_pmf(P$family, log(P$mu_fixed) + xb,
               w %||% rep(1 / ncol(xb), ncol(xb)), P$vec_1, P$vec_2, counts,
               .fcr_dens_args(m$mix_draws %||% m$ndraws,
                              m$engine %||% "halton"))
}

# Hierarchical models: given the effects conditioned on (`re.form`); other
# effects are zero, or integrated out with `marginal`. `cond_var` also
# integrates over the uncertainty of the estimated effects.
.fcr_hier_pmf <- function(object, data, counts, re.form = NULL,
                          allow.new.levels = TRUE, marginal = FALSE,
                          history = NULL, cond_var = FALSE) {
  r <- .fcr_hier_nodes(object, data, re.form, allow.new.levels, marginal,
                       history, cond_var)
  .fcr_mix_pmf(object$model$family, r$E, r$w, r$alpha, r$sigma, counts,
               .fcr_dens_args(object$model$ndraws %||% 1500))
}

# Linear predictors of the rows at nodes over the integrated effects.
.fcr_hier_nodes <- function(object, data, re.form, allow.new.levels,
                            marginal, history, cond_var) {
  m <- object$model
  data <- as.data.frame(data)
  eff <- if (!is.null(history)) .fcr_hier_history_effects(object, history)
  parts <- .fcr_hier_parts(object, data, eff)
  if (!is.null(history)) allow.new.levels <- TRUE
  incl <- .fcr_hier_reform(re.form, m$re_terms)
  n <- nrow(data)
  eta <- parts$eta0
  v <- numeric(n)                     # variance of integrated normal effects
  other <- NULL                       # n x R effects of other distributions
  for (i in seq_along(parts$re)) {
    r <- parts$re[[i]]
    integ <- if (incl[i]) r$new else rep(TRUE, n)
    if (incl[i]) {
      if (any(r$new) && !allow.new.levels) {
        stop("`newdata` has groups of (", m$re_terms[[i]]$label, ") not ",
             "seen in estimation; set `allow.new.levels = TRUE`.",
             call. = FALSE)
      }
      eta <- eta + r$contrib
      if (cond_var) v <- v + .fcr_hier_cond_var(m, i, r)
    }
    if (!marginal || !any(integ)) next
    if (all(r$dist == "normal")) {
      v <- v + integ * rowSums((r$X %*% r$Sigma) * r$X)
    } else {
      # Only single-term models have other distributions
      q <- ncol(r$X)
      R <- 512L
      U <- if (q == 1L) matrix((seq_len(R) - 0.5) / R) else
        .fcr_halton(R, q, TRUE, 1L)
      S <- lapply(seq_len(q), function(k)
        if (r$dist[k] == "gamma") matrix(U[, k], 1L) else
          matrix(.fcr_re_std(U[, k], r$dist[k]), 1L))
      B <- t(vapply(seq_len(R), function(j) as.vector(
        .fcr_re_draw_effects(S, matrix(U[, 1L], 1L), r$dist, r$L, j)),
        numeric(q)))
      if (q == 1L) B <- matrix(B, R, 1L)
      other <- integ * (r$X %*% t(B))
    }
  }
  if (!is.null(other)) {
    # Rows whose effect is estimated carry its uncertainty (v); spread it
    # over the nodes by normal quantiles in a scrambled order
    R <- ncol(other)
    z <- stats::qnorm((order(.fcr_seeded_runif(R, 1L)) - 0.5) / R)
    E <- eta + other + outer(sqrt(v), z)
    w <- rep(1 / R, R)
  } else if (any(v > 0, na.rm = TRUE)) {
    g <- .fcr_gh_rule(40L)
    E <- eta + outer(sqrt(v), g$z)
    w <- g$w
  } else {
    E <- matrix(eta, n, 1L)
    w <- 1
  }
  list(E = E, w = w, alpha = parts$alpha, sigma = parts$sigma, eta = eta,
       v = v)
}

# z' Var(b | y) z of a term's estimated effects (0 for groups not seen).
.fcr_hier_cond_var <- function(m, i, r) {
  cv <- m$condvar[[i]]
  n <- nrow(r$X)
  if (is.null(cv)) return(numeric(n))
  out <- numeric(n)
  q <- ncol(r$X)
  j <- r$j
  ok <- !is.na(j)
  for (a in seq_len(q)) for (b in seq_len(q)) {
    out[ok] <- out[ok] + r$X[ok, a] * r$X[ok, b] * cv[a, b, j[ok]]
  }
  pmax(out, 0)
}

# Multivariate models: one matrix per outcome (marginal of the others).
.fcr_mv_prob <- function(object, data, counts, random = "marginal") {
  data <- as.data.frame(data)
  E <- .fcr_mv_design(object, data)
  p <- object$coefficients
  out <- lapply(seq_along(object$outcomes), function(j)
    .fcr_mv_integrate(object, E, data, function(adds)
      .fcr_mv_pmf(p, E$eqs[[j]], adds[[j]], counts), random))
  stats::setNames(out, object$outcomes)
}

# Total-plus-allocation models. what = "total": P(total = k); "counts": a
# list with the probabilities of each category's count,
#   P(y_j = k) = sum_{c >= k} P(total = c) Binomial(k; c, pi_j),
# the total thinned by the category's probability (multinomial allocation).
.fcr_alloc_pmf <- function(object, data, counts, what = c("counts", "total"),
                           conditional = FALSE, outcome = NULL) {
  what <- match.arg(what)
  if (what == "counts" && object$likelihood == "fractional") {
    stop("The fractional quasi-likelihood does not define a distribution of ",
         "the category counts; use type = \"total_prob\" for the total.",
         call. = FALSE)
  }
  E <- .fcr_alloc_design(object, as.data.frame(data))
  p <- object$coefficients
  n <- nrow(E$X)
  cats <- object$categories
  fam <- object$family
  by_draw <- isTRUE(object$truncated)
  js <- if (is.null(outcome)) seq_along(cats) else match(outcome, cats)
  tot <- function(pc, cs) .fcr_mix_pmf(
    fam, matrix(log(pc$mu)), 1, pc$alpha, pc$sigma, cs, E$dens_args)
  # The pieces at one value of the random coefficients: a list of matrices
  cells <- function(add) {
    Er <- E
    idx <- seq_len(n)
    wv <- rep(1, n)
    if (!is.null(E$S$SH)) {
      Er <- E$Erep %||% .fcr_alloc_rep(E)
      idx <- Er$idx
      add <- lapply(add, function(a)
        if (is.matrix(a)) a[idx, , drop = FALSE] else a[idx])
      add <- .fcr_alloc_shared_add(add, .fcr_alloc_loadings(p, Er), Er$u)
      wv <- rep(E$GH$w, each = n)
    }
    pc <- .fcr_alloc_predict_core(p, Er, add)
    p0 <- tot(pc, 0)
    if (what == "total") {
      M <- list(tot(pc, counts))
    } else {
      # The total's distribution, until its tail is negligible
      hi <- max(counts)
      top <- max(hi, ceiling(3 * max(pc$total[is.finite(pc$total)], 1)) + 20)
      Pt <- tot(pc, 0:top)
      while (top < 20000 &&
             any(1 - rowSums(Pt) > 1e-8, na.rm = TRUE)) {
        more <- (top + 1):(2 * top)
        Pt <- cbind(Pt, tot(pc, more))
        top <- 2 * top
      }
      M <- lapply(js, function(j) {
        pj <- pc$prob[, j]
        vapply(counts, function(k) {
          cs <- k:top
          rowSums(Pt[, cs + 1, drop = FALSE] *
                    matrix(stats::dbinom(k, rep(cs, each = length(pj)), pj),
                           length(pj)))
        }, numeric(length(pj)))
      })
    }
    if (conditional) {
      # Given at least one event: remove the probability of a zero total
      M <- lapply(M, function(m) {
        m <- matrix(m, ncol = length(counts))
        z <- which(counts == 0)
        if (length(z)) m[, z] <- m[, z] - p0
        m
      })
    }
    M <- c(M, list(p0))
    # Average over the shared term first, as in the likelihood
    M <- lapply(M, function(m) {
      m <- matrix(m, nrow = length(idx))
      unname(rowsum(m * wv, idx, reorder = TRUE))
    })
    if (conditional && by_draw) {
      P0 <- M[[length(M)]][, 1L]
      M[-length(M)] <- lapply(M[-length(M)], function(m) m / (1 - P0))
    }
    M
  }
  if (is.null(object$rpar)) {
    M <- cells(NULL)
  } else {
    RS <- object$rpar$layout
    L <- .fcr_alloc_chol(p, RS)
    nd <- .fcr_normal_nodes(RS$K)
    M <- NULL
    for (r in seq_along(nd$w)) {
      B <- matrix(nd$Z[r, ] %*% t(L), n, RS$K, byrow = TRUE)
      mr <- lapply(cells(.fcr_alloc_add(B, E)), `*`, nd$w[r])
      M <- if (is.null(M)) mr else Map(`+`, M, mr)
    }
  }
  P0 <- M[[length(M)]][, 1L]
  M <- M[-length(M)]
  if (conditional && !by_draw) M <- lapply(M, function(m) m / (1 - P0))
  M <- lapply(M, function(m) if (conditional) pmax(m, 0) else m)
  if (what == "total") return(M[[1L]])
  stats::setNames(M, cats[js])
}

# ---- intervals ----------------------------------------------------------------

# Parameter vectors from the asymptotic normal distribution of the estimates.
# Parameters without a usable standard error (at a boundary) are held fixed.
.fcr_general_draws <- function(object, ndraws, seed, scrambled = TRUE) {
  theta <- .margEff_theta(object)
  V <- .margEff_vcov(object)
  if (is.null(V) || nrow(V) != length(theta)) {
    stop("The model has no covariance matrix for its estimates, so ",
         "intervals are not available.", call. = FALSE)
  }
  free <- is.finite(diag(V)) & diag(V) > 0
  # The log of a distribution parameter at its boundary (such as ln(alpha) at
  # the Poisson limit) has an enormous standard error, with no information:
  # the likelihood is flat, and draws of it would be degenerate
  flat <- grepl("(^|:)(ln\\([^)]*\\)|lnalpha)$", names(theta)) & free &
    diag(V) > 10^2
  free <- free & !flat
  if (!any(free) || !.fcr_is_pd(V[free, free, drop = FALSE])) {
    stop("The covariance matrix of the estimates is not positive definite, ",
         "so intervals are not available.", call. = FALSE)
  }
  if (any(!free)) {
    message("Held at their estimates, because they have no usable standard ",
            "error: ",
            toString(names(theta)[!free]), ".")
  }
  k <- sum(free)
  z <- stats::qnorm(.fcr_halton(ndraws, k + 1L, scrambled,
                                seed = if (is.null(seed)) 1L else seed))
  D <- matrix(theta, ndraws, length(theta), byrow = TRUE,
              dimnames = list(NULL, names(theta)))
  D[, free] <- D[, free] + z[, seq_len(k), drop = FALSE] %*%
    chol(V[free, free, drop = FALSE])
  list(draws = D, z = z[, k + 1L])
}

# Confidence intervals for the expected count and prediction intervals for a
# new count, for the estimators other than countreg(). Returns a data frame
# (fit, lwr, upr), or a list of them by outcome.
.fcr_general_interval <- function(object, data, interval, level,
                                  ndraws = NULL, seed = NULL,
                                  scrambled = TRUE, max_count = 1e4,
                                  uncertainty = "normal", outcome_type = NULL,
                                  ...) {
  .fcr_check_prob(level, "level")
  if (!identical(uncertainty, "normal")) {
    stop("uncertainty = \"", uncertainty, "\" is available for countreg() ",
         "models only.", call. = FALSE)
  }
  heavy <- inherits(object, c("countreg_mv", "countreg_alloc")) &&
    !is.null(object$rpar)
  if (is.null(ndraws)) {
    ndraws <- if (interval == "confidence") 500L else if (heavy) 50L else 200L
  }
  if (!is.numeric(ndraws) || length(ndraws) != 1L || is.na(ndraws) ||
      ndraws < 10) {
    stop("`ndraws` must be a number of at least 10.", call. = FALSE)
  }
  ndraws <- as.integer(ndraws)
  data <- as.data.frame(data)
  n <- nrow(data)
  is_hier <- inherits(object, "countreg_hier")
  is_rp <- identical(object$modelType, "countreg.rp") ||
    identical(object$model$modelType, "countreg.rp")
  is_total <- identical(outcome_type, "total")
  dots <- list(...)
  mean_of <- function(obj, z = 0) {
    if (is_hier) return(.fcr_hier_interval_mean(obj, data, z, dots))
    pr <- if (is_total) stats::predict(obj, newdata = data, type = "total")
      else do.call(stats::predict, c(list(obj, newdata = data), dots))
    matrix(as.vector(unname(pr)), nrow = n)
  }
  hd <- NULL
  if (is_rp) {
    K <- .fcr_rp_parts(object, data[0L, , drop = FALSE])$N_rand
    R <- 25L
    # One joint point set for the parameter and coefficient draws
    hd <- .fcr_halton(ndraws * R, K, TRUE, 2L)
  }
  pmf_of <- function(obj, counts, j) {
    P <- if (is_rp) {
      .fcr_rp_pmf(obj, data, counts,
                  hdraws = hd[(j - 1L) * R + seq_len(R), , drop = FALSE])
    } else if (is_hier) {
      do.call(.fcr_hier_pmf, c(list(obj, data, counts, cond_var = TRUE),
                               dots))
    } else if (inherits(obj, "countreg_alloc")) {
      .fcr_alloc_pmf(obj, data, counts,
                     what = if (is_total) "total" else "counts")
    } else {
      .fcr_count_pmf(obj, data, counts)
    }
    if (is.list(P)) P else list(P)
  }
  fit <- mean_of(object)
  J <- ncol(fit)
  g <- .fcr_general_draws(object, ndraws, seed, scrambled)
  objs <- function(j) .margEff_set_theta(object, g$draws[j, ])
  a <- (1 - level) / 2
  if (interval == "confidence") {
    vals <- vapply(seq_len(ndraws), function(j) as.vector(
      suppressWarnings(mean_of(objs(j), g$z[j]))), numeric(n * J))
    vals <- matrix(vals, nrow = n * J)
    qf <- function(pr) matrix(apply(vals, 1, stats::quantile, probs = pr,
                                    na.rm = TRUE, names = FALSE), nrow = n)
    lwr <- qf(a)
    upr <- qf(1 - a)
  } else {
    # Predictive distribution: the PMF averaged over the parameter draws,
    # extended until the upper quantile of every row is reached
    avg <- function(counts) {
      S <- NULL
      acc <- 0
      for (j in seq_len(ndraws)) {
        P <- suppressWarnings(pmf_of(objs(j), counts, j))
        use <- lapply(P, function(m) matrix(is.finite(m), nrow = n))
        P <- Map(function(m, u) {
          m <- matrix(m, nrow = n)
          m[!u] <- 0
          m
        }, P, use)
        S <- if (is.null(S)) P else Map(`+`, S, P)
        acc <- if (is.list(acc)) Map(`+`, acc, use) else use
      }
      Map(`/`, S, acc)
    }
    hi <- max(20, ceiling(3 * max(fit[is.finite(fit)], 1)) + 10)
    hi <- min(hi, max_count)
    S <- avg(0:hi)
    repeat {
      top <- vapply(S, function(m) min(rowSums(m)[!is.na(m[, 1L])], 1),
                    numeric(1))
      if (all(top >= 1 - a) || hi >= max_count) break
      more <- (hi + 1):min(2 * hi, max_count)
      S <- Map(cbind, S, avg(more))
      hi <- max(more)
    }
    q <- function(pr) vapply(S, function(m) {
      apply(m, 1, function(x) {
        k <- which(cumsum(x) >= pr - 1e-12)
        if (length(k)) k[1L] - 1 else NA_real_
      })
    }, numeric(n))
    lwr <- matrix(q(a), nrow = n)
    upr <- matrix(q(1 - a), nrow = n)
    if (any(is.na(upr) & !is.na(lwr))) {
      warning("Some upper prediction limits lie above max_count = ",
              max_count, "; they are NA. Increase `max_count`.",
              call. = FALSE)
    }
  }
  bad <- is.na(fit)
  lwr[bad] <- NA
  upr[bad] <- NA
  info <- list(interval = interval, level = level, uncertainty = "normal",
               draws = ndraws)
  out <- lapply(seq_len(J), function(j) {
    d <- data.frame(fit = fit[, j], lwr = lwr[, j], upr = upr[, j])
    attr(d, "interval") <- info
    d
  })
  if (J == 1L && (is_total || !inherits(object, c("countreg_mv",
                                                  "countreg_alloc")))) {
    return(out[[1L]])
  }
  stats::setNames(out, object$outcomes %||% object$categories)
}

# Expected count of a hierarchical model with the estimated effects shifted
# by z posterior standard deviations (z = 0: predict()).
.fcr_hier_interval_mean <- function(object, data, z, dots) {
  d <- utils::modifyList(list(re.form = NULL, allow.new.levels = TRUE,
                              marginal = FALSE), dots)
  link <- do.call(stats::predict, c(list(object, newdata = data,
                                         type = "link"), dots))
  parts <- .fcr_hier_parts(object, data)
  m <- object$model
  if (z != 0) {
    incl <- .fcr_hier_reform(d$re.form, m$re_terms)
    v <- numeric(nrow(data))
    for (i in which(incl)) {
      v <- v + .fcr_hier_cond_var(m, i, parts$re[[i]])
    }
    link <- link + z * sqrt(v)
  }
  matrix(unname(.fcr_family_mean(m$family, exp(link), parts$alpha,
                                 parts$sigma)), ncol = 1L)
}
