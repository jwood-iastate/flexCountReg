# Copulas for countreg_mv(): the joint distribution of the outcomes is
#   P(Y_1 <= y_1, ..., Y_J <= y_J) = C(F_1(y_1), ..., F_J(y_J); theta),
# with each outcome's own count distribution F_j. For counts, the
# probability of the observed combination is the C-volume of the rectangle
# [F_j(y_j - 1), F_j(y_j)] (Genest and Neslehova 2007; Nikoloulopoulos and
# Karlis 2009):
#   P(Y = y) = sum over s in {0, 1}^J of (-1)^|s| C(u_1^s, ..., u_J^s),
# with u_j^s = F_j(y_j) for s_j = 0 and F_j(y_j - 1) for s_j = 1. An outcome
# left out (marginalized) has u_j = 1, and every copula used here then
# reduces to the copula of the other outcomes.
#
# Families and their parameters (theta on the copula's own scale, eta the
# unconstrained parameter that is estimated):
#   gaussian  rho = tanh(eta)                          two outcomes
#   frank     theta = eta (two outcomes, any sign);
#             theta = exp(eta) (more outcomes, positive dependence)
#   clayton   theta = exp(eta)                         lower-tail dependence
#   gumbel    theta = 1 + exp(eta)                     upper-tail dependence

.fcr_copula_families <- c("gaussian", "frank", "clayton", "gumbel")

# Name of the copula parameter on the estimated scale.
.fcr_cop_par_name <- function(family, J) {
  switch(family,
    gaussian = "copula:atanh(rho)",
    frank = if (J == 2L) "copula:theta" else "copula:ln(theta)",
    clayton = "copula:ln(theta)",
    gumbel = "copula:ln(theta - 1)")
}

# The copula parameter from its unconstrained value.
.fcr_cop_theta <- function(eta, family, J) {
  switch(family,
    gaussian = max(min(tanh(eta), 1 - 1e-12), -1 + 1e-12),
    frank = if (J == 2L) eta else exp(eta),
    clayton = exp(eta),
    gumbel = 1 + exp(eta))
}

# Kendall's tau of the copula (of continuous margins: for counts it is the
# copula's dependence, not the rank correlation of the counts).
.fcr_cop_tau <- function(theta, family) {
  switch(family,
    gaussian = 2 / pi * asin(theta),
    clayton = theta / (theta + 2),
    gumbel = 1 - 1 / theta,
    frank = {
      if (abs(theta) < 1e-8) return(0)
      # 1 - 4 / theta + 4 / theta^2 * int_0^theta t / (e^t - 1) dt
      f <- function(t) ifelse(abs(t) < 1e-10, 1 - t / 2, t / expm1(t))
      D <- stats::integrate(f, 0, theta, rel.tol = 1e-10)$value
      1 - 4 / theta + 4 / theta^2 * D
    })
}

# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
.fcr_gl_rule <- function(n) {
  k <- seq_len(n - 1L)
  b <- k / sqrt(4 * k^2 - 1)
  Jm <- matrix(0, n, n)
  Jm[cbind(k, k + 1L)] <- b
  Jm[cbind(k + 1L, k)] <- b
  e <- eigen(Jm, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
}

.fcr_gl_cache <- new.env(parent = emptyenv())
.fcr_gl <- function(n) {
  key <- as.character(n)
  if (is.null(.fcr_gl_cache[[key]])) .fcr_gl_cache[[key]] <- .fcr_gl_rule(n)
  .fcr_gl_cache[[key]]
}

# P(X < a, Y < b) for standard bivariate normal (X, Y) with correlation r,
# vectorized over a and b: Genz's (2004) algorithm (Drezner and Wesolowsky
# 1990) with Gauss-Legendre rules of 6, 12 or 20 points, accurate to about
# 1e-15. Infinite limits are handled exactly.
.fcr_pbvnorm <- function(a, b, r) {
  n <- max(length(a), length(b))
  a <- rep_len(a, n)
  b <- rep_len(b, n)
  out <- numeric(n)
  lo <- a == -Inf | b == -Inf
  ia <- !lo & a == Inf
  ib <- !lo & b == Inf
  out[ia] <- stats::pnorm(b[ia])
  out[ib & !ia] <- stats::pnorm(a[ib & !ia])
  fin <- !lo & !ia & !ib
  if (!any(fin)) return(out)
  # Genz's routine computes the upper probability P(X > h, Y > k)
  h <- -a[fin]
  k <- -b[fin]
  gl <- .fcr_gl(if (abs(r) < 0.3) 6L else if (abs(r) < 0.75) 12L else 20L)
  hk <- h * k
  if (abs(r) < 0.925) {
    hs <- (h * h + k * k) / 2
    asr <- asin(r)
    bvn <- 0
    for (i in seq_along(gl$x)) {
      sn <- sin(asr * (gl$x[i] + 1) / 2)
      bvn <- bvn + gl$w[i] * exp((sn * hk - hs) / (1 - sn * sn))
    }
    bvn <- bvn * asr / (4 * pi) + stats::pnorm(-h) * stats::pnorm(-k)
  } else {
    if (r < 0) {
      k <- -k
      hk <- -hk
    }
    bvn <- 0
    if (abs(r) < 1) {
      as <- (1 - r) * (1 + r)
      A <- sqrt(as)
      bs <- (h - k)^2
      cc <- (4 - hk) / 8
      dd <- (12 - hk) / 16
      asr <- -(bs / as + hk) / 2
      bvn <- ifelse(asr > -100, A * exp(asr) *
                      (1 - cc * (bs - as) * (1 - dd * bs / 5) / 3 +
                         cc * dd * as * as / 5), 0)
      B <- sqrt(bs)
      corr <- exp(-hk / 2) * sqrt(2 * pi) * stats::pnorm(-B / A) * B *
        (1 - cc * bs * (1 - dd * bs / 5) / 3)
      bvn <- bvn - ifelse(-hk < 100, corr, 0)
      A <- A / 2
      for (i in seq_along(gl$x)) {
        xs <- (A * (gl$x[i] + 1))^2
        rs <- sqrt(1 - xs)
        asr <- -(bs / xs + hk) / 2
        term <- A * gl$w[i] * exp(asr) *
          (exp(-hk * (1 - rs) / (2 * (1 + rs))) / rs - (1 + cc * xs * (1 + dd * xs)))
        bvn <- bvn + ifelse(asr > -100, term, 0)
      }
      bvn <- -bvn / (2 * pi)
    }
    if (r > 0) {
      bvn <- bvn + stats::pnorm(-pmax(h, k))
    } else {
      bvn <- -bvn + ifelse(k > h, stats::pnorm(k) - stats::pnorm(h), 0)
    }
  }
  out[fin] <- pmin(pmax(bvn, 0), 1)
  out
}

# The copula C(U) for the rows of U (n x J; columns equal to 1 drop out).
.fcr_cop_C <- function(U, family, theta) {
  U <- pmin(pmax(U, 0), 1)
  switch(family,
    gaussian = .fcr_pbvnorm(stats::qnorm(U[, 1L]), stats::qnorm(U[, 2L]),
                            theta),
    clayton = {
      # u^(-theta) - 1 without cancellation for small theta
      s <- rowSums(expm1(-theta * log(U)))
      out <- exp(-log1p(s) / theta)
      out[rowSums(U == 0) > 0] <- 0
      out
    },
    gumbel = {
      s <- rowSums((-log(U))^theta)
      exp(-s^(1 / theta))
    },
    frank = {
      if (abs(theta) < 1e-10) return(apply(U, 1L, prod))
      J <- ncol(U)
      num <- apply(expm1(-theta * U), 1L, prod) / expm1(-theta)^(J - 1L)
      -log1p(num) / theta
    })
}

# Probability of the rectangles [Flo, Fhi] (n x J): the C-volume, summed
# over the 2^J corners. Columns of Fhi equal to 1 with Flo NA are
# marginalized (not differenced).
.fcr_cop_rect <- function(Fhi, Flo, family, theta) {
  J <- ncol(Fhi)
  act <- which(colSums(!is.na(Flo)) > 0)
  out <- 0
  for (s in 0:(2^length(act) - 1L)) {
    bits <- as.integer(intToBits(s))[seq_along(act)]
    U <- Fhi
    U[is.na(U)] <- 1
    if (any(bits == 1L)) {
      lo <- act[bits == 1L]
      U[, lo] <- Flo[, lo]
    }
    # A corner with some u = 0 contributes 0
    z <- rowSums(U == 0) > 0
    cv <- numeric(nrow(U))
    if (any(!z)) cv[!z] <- .fcr_cop_C(U[!z, , drop = FALSE], family, theta)
    out <- out + (-1)^sum(bits) * cv
  }
  out
}

# The equation's distribution function at q (on the probability scale).
.fcr_mv_cdf <- function(eq, v, q) {
  exp(.fcr_fm_logcdf(eq$dens_args, list(mu = v$mu, alpha = v$alpha,
                                        sigma = v$sigma), eq$family, q))
}

# Rectangle bounds of the counts Y (n x J; NA = marginalized) at the
# parameters p.
.fcr_mv_cop_bounds <- function(p, E, Y) {
  n <- nrow(Y)
  Fhi <- matrix(NA_real_, n, E$J)
  Flo <- Fhi
  for (j in seq_len(E$J)) {
    if (all(is.na(Y[, j]))) next
    eq <- E$eqs[[j]]
    v <- .fcr_mv_eq_values(p, eq)
    if (is.null(v)) return(NULL)
    Fhi[, j] <- .fcr_mv_cdf(eq, v, Y[, j])
    Flo[, j] <- .fcr_mv_cdf(eq, v, Y[, j] - 1)
  }
  list(hi = Fhi, lo = Flo)
}

# Mean and distribution parameters of an equation (NULL where invalid).
.fcr_mv_eq_values <- function(p, eq) {
  S <- eq$S
  n <- length(eq$y)
  mu <- exp(as.vector(eq$X %*% p[S$idx_beta]) + eq$off)
  dp <- .fcr_alloc_dispar(p, S, eq$W1, eq$W2, n)
  if (!all(.fcr_family_valid(eq$family, mu, dp$alpha, dp$sigma))) {
    return(NULL)
  }
  list(mu = mu, alpha = dp$alpha, sigma = dp$sigma)
}

# Log-likelihood of the copula model: one contribution per row.
.fcr_mv_cop_ll <- function(p, E) {
  Y <- vapply(E$eqs, `[[`, numeric(E$n), "y")
  Y <- matrix(Y, nrow = E$n)
  bd <- .fcr_mv_cop_bounds(p, E, Y)
  if (is.null(bd)) return(rep(NA_real_, E$n))
  theta <- .fcr_cop_theta(p[E$cop$idx], E$cop$family, E$J)
  P <- .fcr_cop_rect(bd$hi, bd$lo, E$cop$family, theta)
  ll <- E$w * log(pmax(P, 1e-300))
  ll[E$w == 0] <- 0
  ll
}

# Draws of the copula (n x J matrix of uniforms): conditional inversion for
# two outcomes (Gaussian, Frank of any sign), Marshall-Olkin frailties for
# Archimedean copulas (Clayton: gamma; Gumbel: positive stable, Kanter
# 1975; Frank: logarithmic series, Kemp 1981).
.fcr_rcopula <- function(n, J, family, theta) {
  if (family == "gaussian") {
    z1 <- stats::rnorm(n)
    z2 <- theta * z1 + sqrt(1 - theta^2) * stats::rnorm(n)
    return(cbind(stats::pnorm(z1), stats::pnorm(z2)))
  }
  if (family == "frank" && (J == 2L || theta <= 0)) {
    u1 <- stats::runif(n)
    v <- stats::runif(n)
    if (abs(theta) < 1e-10) return(cbind(u1, v))
    u2 <- -log1p(v * expm1(-theta) / (v + (1 - v) * exp(-theta * u1))) / theta
    return(cbind(u1, u2))
  }
  E <- matrix(stats::rexp(n * J), n, J)
  switch(family,
    clayton = {
      V <- stats::rgamma(n, shape = 1 / theta)
      (1 + E / V)^(-1 / theta)
    },
    gumbel = {
      a <- 1 / theta
      th <- stats::runif(n, 0, pi)
      W <- stats::rexp(n)
      V <- (sin(a * th) / sin(th)^(1 / a)) * (sin((1 - a) * th) / W)^((1 - a) / a)
      exp(-(E / V)^a)
    },
    frank = {
      p <- -expm1(-theta)
      V <- .fcr_rlogseries(n, p)
      -log1p(-p * exp(-E / V)) / theta
    })
}

# Logarithmic series draws, P(V = k) = p^k / (-log(1 - p) k) (Kemp 1981, LK).
.fcr_rlogseries <- function(n, p) {
  h <- log1p(-p)
  u2 <- stats::runif(n)
  out <- rep(1, n)
  todo <- u2 <= p
  if (any(todo)) {
    q <- -expm1(stats::runif(sum(todo)) * h)
    v <- u2[todo]
    k <- ifelse(v < q^2, floor(1 + log(v) / log(q)), ifelse(v > q, 1, 2))
    out[todo] <- k
  }
  out
}

# Smallest count k with F(k) >= u, for each row (the quantile of a
# discrete margin).
.fcr_mv_qcount <- function(eq, v, u, max_count = 1e6) {
  n <- length(u)
  out <- rep(NA_real_, n)
  k <- 0
  act <- seq_len(n)
  while (length(act) && k <= max_count) {
    vk <- list(mu = v$mu[act], alpha = if (!is.null(v$alpha)) v$alpha[act],
               sigma = if (!is.null(v$sigma)) v$sigma[act])
    Fk <- .fcr_mv_cdf(list(dens_args = eq$dens_args, family = eq$family), vk,
                      rep(k, length(act)))
    done <- Fk >= u[act] - 1e-12
    out[act[done]] <- k
    act <- act[!done]
    k <- k + 1
  }
  out
}
