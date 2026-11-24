# Likelihood engines of countreg_hier().
#
# With u ~ N(0, I) the marginal log-likelihood is
#   log L = log integral exp(h(u)) du - (Q/2) log(2 pi),
#   h(u)  = sum_i w_i log f(y_i | eta_i) - u'u / 2,
#   eta   = X beta + offset + A u,   A = Z Lambda.
# The mode u_hat maximizes h by Newton's method; the negative Hessian at u is
#   M = A' D A + I,   D_i = -w_i d^2 log f / d eta_i^2.
#
# Laplace:  log L = h(u_hat) - log|M| / 2.
# Adaptive Gauss-Hermite quadrature (one grouping term): the groups are
# independent, and each group's integral is evaluated at nodes centered at
# its mode and scaled by its curvature,
#   log L_j = -log|R_j| + log sum_k w_k exp(h_j(u_j + R_j^{-1} z_k) + z_k'z_k/2),
# with M_j = R_j' R_j and (z_k, w_k) the Gauss-Hermite rule for N(0, I).
# With one node this is the Laplace approximation.

# Log-density of each row and its first two derivatives with respect to the
# mean's linear predictor eta (all multiplied by the weights w). Analytic for
# POISSON, NB1 and NB2; other families use five-point differences of their
# log-density in eta.
.fcr_hier_logf <- function(family, y, w, dens_args = list()) {
  fam <- .fcr_get_family(family)
  lf <- function(eta, alpha, sigma) {
    out <- tryCatch(
      do.call(fam$dens, c(list(y = y, predicted = exp(eta), alpha = alpha,
                               sigma = sigma), dens_args, list(log = TRUE))),
      error = function(e) rep(NA_real_, length(y)))
    out[is.na(out)] <- -Inf
    out
  }
  lgy <- lgamma(y + 1)
  switch(family,
    POISSON = function(eta, alpha, sigma, deriv = TRUE) {
      mu <- exp(eta)
      list(l = w * (y * eta - mu - lgy), d1 = w * (y - mu), d2 = -w * mu)
    },
    NB2 = function(eta, alpha, sigma, deriv = TRUE) {
      mu <- exp(eta)
      r <- 1 / alpha
      rm <- r + mu
      # dnbinom() is accurate near the Poisson limit (large r), where
      # lgamma(y + r) - lgamma(r) loses its precision
      l <- stats::dnbinom(y, size = r, mu = mu, log = TRUE)
      list(l = w * l, d1 = w * r * (y - mu) / rm,
           d2 = -w * (y + r) * r * mu / rm^2)
    },
    NB1 = function(eta, alpha, sigma, deriv = TRUE) {
      s <- exp(eta) / alpha
      la <- log1p(alpha)
      l <- stats::dnbinom(y, size = s, mu = exp(eta), log = TRUE)
      d1 <- s * (digamma(y + s) - digamma(s) - la)
      d2 <- d1 + s^2 * (trigamma(y + s) - trigamma(s))
      list(l = w * l, d1 = w * d1, d2 = w * d2)
    },
    function(eta, alpha, sigma, deriv = TRUE) {
      h <- 1e-3
      f0 <- lf(eta, alpha, sigma)
      if (!deriv) return(list(l = w * f0))
      fp1 <- lf(eta + h, alpha, sigma)
      fm1 <- lf(eta - h, alpha, sigma)
      fp2 <- lf(eta + 2 * h, alpha, sigma)
      fm2 <- lf(eta - 2 * h, alpha, sigma)
      d1 <- (-fp2 + 8 * fp1 - 8 * fm1 + fm2) / (12 * h)
      d2 <- (-fp2 + 16 * fp1 - 30 * f0 + 16 * fm1 - fm2) / (12 * h^2)
      bad <- !is.finite(d1) | !is.finite(d2)
      d1[bad] <- NA_real_
      d2[bad] <- NA_real_
      list(l = w * f0, d1 = w * d1, d2 = w * d2)
    })
}

# Gauss-Hermite rule for the standard normal (Golub-Welsch): nodes and
# weights summing to 1.
.fcr_gh_rule <- function(n) {
  if (n == 1L) return(list(z = 0, w = 1))
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L))
  J[cbind(seq_len(n - 1L), 2:n)] <- off
  J[cbind(2:n, seq_len(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(z = e$values[o], w = (e$vectors[1, o])^2)
}

# Product rule in q dimensions: a matrix of nodes (one per row) and weights.
.fcr_gh_grid <- function(n, q) {
  r <- .fcr_gh_rule(n)
  Z <- as.matrix(expand.grid(rep(list(r$z), q)))
  W <- as.matrix(expand.grid(rep(list(r$w), q)))
  list(Z = unname(Z), w = apply(W, 1, prod))
}

# Mode of h(u) given eta0 = X beta + offset and A = Z Lambda, by Newton's
# method with step halving from `u`. Returns the mode, the derivatives
# there and the factor of M (with the observed curvature).
.fcr_hier_mode <- function(u, eta0, A, logf, alpha, sigma, tol = 1e-10,
                           max_iter = 100L) {
  Q <- ncol(A)
  I <- Matrix::Diagonal(Q)
  hval <- function(d, u) sum(d$l) - sum(u^2) / 2
  eta <- eta0 + as.vector(A %*% u)
  d <- logf(eta, alpha, sigma)
  if (anyNA(d$d1) || anyNA(d$d2)) return(NULL)
  h0 <- hval(d, u)
  if (!is.finite(h0)) return(NULL)
  for (it in seq_len(max_iter)) {
    g <- as.vector(Matrix::crossprod(A, d$d1)) - u
    D <- pmax(-d$d2, 0)
    if (any(!is.finite(D)) || any(!is.finite(g))) return(NULL)
    M <- Matrix::crossprod(A, A * D) + I
    ch <- tryCatch(suppressWarnings(
      Matrix::Cholesky(Matrix::forceSymmetric(M), perm = TRUE, LDL = FALSE)),
      error = function(e) NULL)
    if (is.null(ch)) return(NULL)
    step <- as.vector(Matrix::solve(ch, g))
    dec <- sum(g * step)
    if (!is.finite(dec)) return(NULL)
    if (dec < tol) {
      # A last full Newton step: log|M| is not stationary at the mode, so
      # its error is first order in the error of u, and the outer
      # derivatives need the mode to near machine precision
      u1 <- u + step
      eta1 <- eta0 + as.vector(A %*% u1)
      d1 <- logf(eta1, alpha, sigma)
      if (!anyNA(d1$d1) && !anyNA(d1$d2)) {
        h1 <- hval(d1, u1)
        if (is.finite(h1)) {
          u <- u1
          eta <- eta1
          d <- d1
          h0 <- h1
        }
      }
      break
    }
    s <- 1
    repeat {
      u1 <- u + s * step
      eta1 <- eta0 + as.vector(A %*% u1)
      d1 <- logf(eta1, alpha, sigma)
      h1 <- if (anyNA(d1$d1) || anyNA(d1$d2)) -Inf else hval(d1, u1)
      if (is.finite(h1) && h1 >= h0 - 1e-12 * abs(h0)) break
      s <- s / 2
      if (s < 1e-10) break
    }
    if (s < 1e-10) break
    u <- u1
    eta <- eta1
    d <- d1
    h0 <- h1
  }
  # Factor of M with the observed curvature
  D <- -d$d2
  M <- Matrix::crossprod(A, A * D) + I
  ch <- if (all(is.finite(D))) tryCatch(suppressWarnings(
    Matrix::Cholesky(Matrix::forceSymmetric(M), perm = TRUE, LDL = FALSE)),
    error = function(e) NULL)
  list(u = u, eta = eta, d = d, h = h0, M = M, ch = ch, D = D)
}

# Laplace approximation of the marginal log-likelihood at the mode.
.fcr_hier_laplace <- function(md) {
  if (is.null(md) || is.null(md$ch)) return(-Inf)
  ld <- as.numeric(Matrix::determinant(md$ch, logarithm = TRUE,
                                       sqrt = FALSE)$modulus)
  md$h - ld / 2
}

# Adaptive Gauss-Hermite quadrature for a single grouping term with q
# effects per group: u is group-major (J x q). `rule` from .fcr_gh_grid().
.fcr_hier_aghq <- function(md, term, eta0, A, logf, alpha, sigma, rule) {
  if (is.null(md) || is.null(md$ch)) return(-Inf)
  q <- term$q
  J <- term$J
  g <- as.integer(term$g)
  # Blocks of M (one per group) and their Cholesky factors
  idx <- function(a, b) cbind((seq_len(J) - 1L) * q + a,
                              (seq_len(J) - 1L) * q + b)
  Mb <- array(0, c(q, q, J))
  for (a in seq_len(q)) for (b in seq_len(q)) {
    Mb[a, b, ] <- as.vector(md$M[idx(a, b)])
  }
  Rinv <- array(0, c(q, q, J))
  logdetR <- numeric(J)
  for (j in seq_len(J)) {
    R <- tryCatch(chol(Mb[, , j]), error = function(e) NULL)
    if (is.null(R)) return(-Inf)
    Rinv[, , j] <- backsolve(R, diag(q))
    logdetR[j] <- sum(log(diag(R)))
  }
  U <- matrix(md$u, J, q, byrow = TRUE)
  K <- nrow(rule$Z)
  lk <- matrix(0, J, K)
  for (k in seq_len(K)) {
    z <- rule$Z[k, ]
    # u_jk = u_j + R_j^{-1} z
    Uk <- U
    for (a in seq_len(q)) {
      Uk[, a] <- U[, a] + colSums(matrix(Rinv[a, , ], q, J) * z)
    }
    uk <- as.vector(t(Uk))
    eta <- eta0 + as.vector(A %*% uk)
    d <- logf(eta, alpha, sigma)
    rl <- rowsum(d$l, g, reorder = TRUE)[, 1]
    lk[, k] <- rl - rowSums(Uk^2) / 2 + sum(z^2) / 2 + log(rule$w[k])
  }
  m <- apply(lk, 1, max)
  lj <- m + log(rowSums(exp(lk - m))) - logdetR
  if (any(!is.finite(lj))) return(-Inf)
  sum(lj)
}

# ---- Simulation engine (one grouping term, any effect distribution) --------
#
# The effects of group j are b_j = L v_j with v_j a vector of independent
# standardized draws (mean 0, variance 1): normal, uniform or symmetric
# triangular. With "gamma" (one effect) exp(b_j) is gamma with mean 1 and
# standard deviation s = |L_11|, the multiplicative heterogeneity of the
# negative binomial random-effects model. The likelihood of each group is
# the average over R scrambled Halton draws of the product of its rows'
# densities.

.fcr_re_dists <- c("normal", "uniform", "triangular", "gamma")

# Standardized draws from uniforms U.
.fcr_re_std <- function(U, dist) {
  switch(dist,
    normal = stats::qnorm(U),
    uniform = sqrt(3) * (2 * U - 1),
    triangular = sqrt(6) * ifelse(U < 0.5, -1 + sqrt(2 * U),
                                  1 - sqrt(2 * (1 - U))),
    stop("Unknown random-effect distribution \"", dist, "\".",
         call. = FALSE))
}

# log of a gamma(1/s^2, 1/s^2) draw (mean 1, sd s) from uniforms U.
.fcr_re_loggamma <- function(U, s) {
  s <- abs(s)
  if (s < 1e-4) return(s * stats::qnorm(U) - s^2 / 2)
  k <- 1 / s^2
  log(stats::qgamma(U, shape = k, rate = k))
}

# Effects of every group at draw r: a J x q matrix.
.fcr_re_draw_effects <- function(S, U1, dist, L, r) {
  q <- nrow(L)
  B <- matrix(0, nrow(S[[1]]), q)
  if (dist[1] == "gamma") {
    B[, 1] <- .fcr_re_loggamma(U1[, r], L[1, 1])
    return(B)
  }
  for (k in seq_len(q)) {
    for (m in seq_len(k)) B[, k] <- B[, k] + L[k, m] * S[[m]][, r]
  }
  B
}

# Simulated log-likelihood; with `posterior = TRUE` also the posterior
# means and covariances of the effects of each group.
#   U         list of q (J x R) matrices of uniforms
#   X         the term's n x q effect design, g its group index (1..J)
#   logf_for  function(c) returning the row log-density for c stacked
#             copies of the data (draws are evaluated c at a time)
.fcr_hier_msl <- function(eta0, X, g, U, dist, L, logf_for, alpha, sigma,
                          posterior = FALSE) {
  J <- nrow(U[[1]])
  R <- ncol(U[[1]])
  q <- ncol(X)
  n <- length(eta0)
  # Effects of every group and draw: q matrices J x R
  if (dist[1] == "gamma") {
    Bk <- list(.fcr_re_loggamma(U[[1]], L[1, 1]))
  } else {
    S <- lapply(seq_len(q), function(k) .fcr_re_std(U[[k]], dist[k]))
    Bk <- lapply(seq_len(q), function(k) {
      out <- matrix(0, J, R)
      for (m in seq_len(k)) out <- out + L[k, m] * S[[m]]
      out
    })
  }
  csize <- max(1L, min(R, floor(2e6 / n)))
  lr <- matrix(0, J, R)
  for (st in seq(1L, R, by = csize)) {
    cols <- st:min(R, st + csize - 1L)
    nc <- length(cols)
    eta <- matrix(eta0, n, nc)
    for (k in seq_len(q)) {
      eta <- eta + X[, k] * Bk[[k]][g, cols, drop = FALSE]
    }
    l <- logf_for(nc)(as.vector(eta), if (!is.null(alpha)) rep(alpha, nc),
                      if (!is.null(sigma)) rep(sigma, nc), deriv = FALSE)$l
    l[is.na(l)] <- -Inf
    lr[, cols] <- rowsum(matrix(l, n, nc), g, reorder = TRUE)
  }
  m <- apply(lr, 1, max)
  if (any(!is.finite(m))) return(list(ll = -Inf))
  w <- exp(lr - m)
  ll <- sum(m + log(rowMeans(w)))
  if (!posterior) return(list(ll = ll))
  w <- w / rowSums(w)
  mean_b <- vapply(Bk, function(B) rowSums(w * B), numeric(J))
  mean_b <- matrix(mean_b, J, q)
  cv <- array(0, c(q, q, J))
  for (a in seq_len(q)) for (b in seq_len(q)) {
    cv[a, b, ] <- rowSums(w * Bk[[a]] * Bk[[b]]) -
      mean_b[, a] * mean_b[, b]
  }
  list(ll = ll, b = mean_b, condvar = cv)
}

# Row log-densities for c stacked copies of the data (cached by c).
.fcr_hier_logf_stack <- function(family, y, w, dens_args) {
  cache <- new.env(parent = emptyenv())
  function(c) {
    key <- as.character(c)
    if (is.null(cache[[key]])) {
      cache[[key]] <- .fcr_hier_logf(family, rep(y, c), rep(w, c),
                                     dens_args)
    }
    cache[[key]]
  }
}

# log E[exp(z'b)] of each row for a term's effects (z the rows of X),
# used for marginal means.
.fcr_re_log_mgf <- function(X, L, dist) {
  if (all(dist == "normal")) {
    S <- L %*% t(L)
    return(rowSums((X %*% S) * X) / 2)
  }
  if (dist[1] == "gamma") {
    s <- abs(L[1, 1])
    if (s < 1e-8) return(rep(0, nrow(X)))
    k <- 1 / s^2
    x <- X[, 1]
    out <- lgamma(k + x) - lgamma(k) - x * log(k)
    out[k + x <= 0] <- Inf
    return(out)
  }
  A <- X %*% L   # (L'z)_m for each row
  out <- 0
  for (m in seq_len(ncol(A))) {
    a <- abs(A[, m])
    out <- out + switch(dist[m],
      normal = a^2 / 2,
      uniform = ifelse(a < 1e-8, 0, .fcr_log_sinhc(sqrt(3) * a)),
      triangular = {
        c6 <- sqrt(6) * a
        ifelse(a < 1e-6, a^2 / 2, log(2 * (cosh(c6) - 1) / c6^2))
      })
  }
  out
}

# Draws of the simulation engine: for each group, R consecutive points of
# the Halton sequence (one dimension per effect), shifted by a uniform
# random amount modulo 1 (Cranley-Patterson rotation), which keeps their
# even spread and makes the average unbiased. The scrambled ("mixed")
# Halton points of randtoolbox are not used: averaged over consecutive
# blocks they left a bias in the log-likelihood that more draws did not
# remove. `seed` fixes the shifts without changing the caller's random
# numbers. Returns a list of q (J x R) matrices.
.fcr_hier_draws <- function(J, R, q, seed = 123) {
  H <- unclass(.fcr_unit_halton(seq_len(J), as.integer(R), q))
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  lapply(H, function(h) (h + stats::runif(J)) %% 1)
}
