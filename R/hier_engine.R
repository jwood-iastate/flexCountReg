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
    POISSON = function(eta, alpha, sigma) {
      mu <- exp(eta)
      list(l = w * (y * eta - mu - lgy), d1 = w * (y - mu), d2 = -w * mu)
    },
    NB2 = function(eta, alpha, sigma) {
      mu <- exp(eta)
      r <- 1 / alpha
      rm <- r + mu
      # dnbinom() is accurate near the Poisson limit (large r), where
      # lgamma(y + r) - lgamma(r) loses its precision
      l <- stats::dnbinom(y, size = r, mu = mu, log = TRUE)
      list(l = w * l, d1 = w * r * (y - mu) / rm,
           d2 = -w * (y + r) * r * mu / rm^2)
    },
    NB1 = function(eta, alpha, sigma) {
      s <- exp(eta) / alpha
      la <- log1p(alpha)
      l <- stats::dnbinom(y, size = s, mu = exp(eta), log = TRUE)
      d1 <- s * (digamma(y + s) - digamma(s) - la)
      d2 <- d1 + s^2 * (trigamma(y + s) - trigamma(s))
      list(l = w * l, d1 = w * d1, d2 = w * d2)
    },
    function(eta, alpha, sigma) {
      h <- 1e-3
      f0 <- lf(eta, alpha, sigma)
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
