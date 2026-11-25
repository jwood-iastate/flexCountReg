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

# Probability of the rectangles [Flo, Fhi] (n x J): the C-volume (see
# .fcr_cop_volume()). Columns with Flo NA are marginalized (not
# differenced). Ghi, Glo: the upper tails 1 - Fhi and 1 - Flo, accurate
# where F is close to 1 (default: computed from Fhi and Flo).
.fcr_cop_rect <- function(Fhi, Flo, family, theta, Ghi = NULL, Glo = NULL) {
  .fcr_cop_volume(Fhi, Flo, family, theta, Ghi, Glo)$P
}

# The C-volume P of the rectangles [Flo, Fhi] (n x J) and, with grad, its
# derivatives with respect to each Fhi (Ahi) and Flo (Alo), n x J, and to
# theta (Dth).
#
# The volume is the sum over the 2^J corners, accurate to about
# 1e-16 / P relative to P when the corners are of order one. Where an
# outcome is in its upper tail, F(y - 1) and F(y) are both close to 1 and
# the corner sum loses the probability (P(Y = y) falls below the rounding
# of 1). Rows where the corner sum is below 1e-6 of the largest corner
# (a relative error above about 2e-10), or (Gumbel) with an outcome whose
# P(Y > y - 1) is below 1e-3, are computed again: each outcome
# with F(y - 1) > 0.5 is described by its upper tails G = 1 - F
# (Ghi = P(Y > y), Glo = P(Y > y - 1)), and the volume is integrated over
# them with Gauss-Legendre rules on geometric panels, while the other
# outcomes are differenced over their corners (in rows without such an
# outcome, where the corners cancel through strong dependence, every
# outcome with F(y - 1) > 0 is integrated):
#   Archimedean copulas, C(u) = psi(sum_j phi(u_j)): with t_j = phi(u_j),
#   the volume over the outcomes N in their upper tails is the integral of
#   (-1)^k psi^(k)(T + sum_{j in N} t_j) over the box of their t_j (k =
#   |N|, T the sum of phi over the corner of the other outcomes), a
#   function of the sum only: one integral for k = 1 and, with the
#   trapezoid density of t_1 + t_2, for k = 2 (nested for k > 2). The
#   integrand is positive and is evaluated on the log scale.
#   Gaussian copula (two outcomes): the conditional probability of one
#   outcome integrated over the other's interval (.fcr_cop_row_gauss()).
.fcr_cop_volume <- function(Fhi, Flo, family, theta, Ghi = NULL, Glo = NULL,
                            grad = FALSE, Gfun = NULL) {
  n <- nrow(Fhi)
  J <- ncol(Fhi)
  fixed <- is.na(Flo)
  Fhi[fixed] <- 1
  if (is.null(Ghi)) Ghi <- 1 - Fhi
  if (is.null(Glo)) Glo <- 1 - Flo
  Ghi <- matrix(Ghi, n, J)
  Glo <- matrix(Glo, n, J)
  Ghi[fixed] <- 0
  nar <- !fixed & Flo > 0.5
  nar[is.na(nar)] <- FALSE
  slow <- rep(FALSE, n)
  P <- numeric(n)
  Cmax <- numeric(n)
  Ahi <- Alo <- if (grad) matrix(0, n, J)
  Dth <- if (grad) numeric(n)
  fast <- seq_len(n)
  if (length(fast)) {
    act <- which(colSums(!fixed[fast, , drop = FALSE]) > 0)
    Uhi <- Fhi[fast, , drop = FALSE]
    for (s in 0:(2^length(act) - 1L)) {
      bits <- as.integer(intToBits(s))[seq_along(act)]
      U <- Uhi
      lo <- act[bits == 1L]
      if (length(lo)) U[, lo] <- Flo[fast, lo, drop = FALSE]
      # A corner with some u = 0 (or the lower end of an outcome that is
      # left out) contributes 0
      ok <- rowSums(U <= 0 | is.na(U)) == 0
      if (!any(ok)) next
      sg <- (-1)^length(lo)
      rows <- fast[ok]
      if (!grad) {
        cv <- .fcr_cop_C(U[ok, , drop = FALSE], family, theta)
        P[rows] <- P[rows] + sg * cv
        Cmax[rows] <- pmax(Cmax[rows], abs(cv))
        next
      }
      d <- .fcr_cop_derivs(U[ok, , drop = FALSE], family, theta)
      P[rows] <- P[rows] + sg * d$C
      Cmax[rows] <- pmax(Cmax[rows], abs(d$C))
      for (j in act) {
        if (bits[match(j, act)] == 0L) {
          Ahi[rows, j] <- Ahi[rows, j] + sg * d$dU[, j]
        } else {
          Alo[rows, j] <- Alo[rows, j] + sg * d$dU[, j]
        }
      }
      Dth[rows] <- Dth[rows] + sg * d$dth
    }
    # The corner sum loses digits where the probability is far below the
    # corners: such rows are integrated instead. The bivariate normal
    # distribution function is accurate to about 1e-16 in absolute terms
    # only, so improbable rows (such as joint lower tails with a negative
    # correlation) are too.
    # Gumbel corners within 1e-3 of 1 also lose digits of the copula's
    # derivatives (through -log(u)), so rows with an outcome that far in
    # its upper tail are integrated too.
    bad <- P < 1e-6 * Cmax | (family == "gaussian" & P < 1e-9) |
      (family == "gumbel" & rowSums(nar & Glo < 1e-3) > 0)
    bad[is.na(bad)] <- FALSE
    slow <- bad & rowSums(!fixed) > 0
  }
  # Rows without an outcome in its upper tail: every outcome with
  # F(y - 1) > 0 is integrated
  full <- slow & rowSums(nar) == 0
  if (any(slow)) {
    if (!is.null(Gfun)) {
      # Accurate upper tails of the rows integrated (from the caller)
      ws <- which(slow)
      gs <- Gfun(ws)
      Ghi[ws, ] <- gs$hi
      Glo[ws, ] <- gs$lo
      Ghi[fixed] <- 0
    }
    th <- theta
    # Frank at theta = 0 is the independence copula (the limit)
    if (family == "frank" && abs(th) < 1e-9) th <- if (th < 0) -1e-9 else 1e-9
    ar <- if (family != "gaussian") .fcr_arch_pieces(family, th)
    for (i in which(slow)) {
      type <- ifelse(fixed[i, ], "f", ifelse(
        nar[i, ] | (full[i] & Flo[i, ] > 0), "n", "w"))
      r <- if (sum(type != "f") == 1L) {
        # One outcome: its own probability
        j <- which(type != "f")
        a_hi <- a_lo <- numeric(J)
        a_hi[j] <- 1
        a_lo[j] <- -1
        list(P = if (nar[i, j]) Glo[i, j] - Ghi[i, j] else
          Fhi[i, j] - Flo[i, j], Ahi = a_hi, Alo = a_lo, Dth = 0)
      } else if (family == "gaussian") {
        .fcr_cop_row_gauss(th, Fhi[i, ], Ghi[i, ], Flo[i, ], Glo[i, ], type,
                           grad)
      } else {
        .fcr_cop_row_arch(ar, Fhi[i, ], Ghi[i, ], Flo[i, ], Glo[i, ], type,
                          grad)
      }
      # Where the integration fails (generator values that overflow at
      # extreme parameters) the corner sum stays
      if (!is.finite(r$P) || (grad && !all(is.finite(c(r$Ahi, r$Alo,
                                                         r$Dth))))) next
      P[i] <- r$P
      if (grad) {
        Ahi[i, ] <- r$Ahi
        Alo[i, ] <- r$Alo
        Dth[i] <- r$Dth
      }
    }
  }
  list(P = P, Ahi = Ahi, Alo = Alo, Dth = Dth)
}

# Nodes r and weights w of a Gauss-Legendre rule for integrals over
# r in [0, D], on panels geometric in delta + r (delta >= 0 the distance of
# r = 0 from the integrand's singular point) with ratio at most
# exp(lratio), split at `breaks` and no wider than maxw.
.fcr_cop_nodes <- function(D, delta, breaks = numeric(0), maxw = Inf,
                           lratio = log(2), ng = 10L, maxp = 300L) {
  cuts <- sort(unique(c(0, breaks[breaks > 0 & breaks < D], D)))
  st <- numeric(0)
  en <- numeric(0)
  for (i in seq_len(length(cuts) - 1L)) {
    lo <- cuts[i]
    hi <- cuts[i + 1L]
    a <- delta + lo
    b <- delta + hi
    # (the part next to a singular end left out is below 1e-130 of the
    # interval)
    if (!(a > 0)) a <- b * exp(-min(lratio * maxp, 300))
    np <- min(maxp, max(1, ceiling(log(b / a) / lratio)))
    e <- a * exp(log(b / a) * (0:np) / np) - delta
    e[1L] <- lo
    e[np + 1L] <- hi
    e <- cummax(e)
    if (is.finite(maxw)) {
      wd <- diff(e)
      if (any(wd > maxw)) {
        m <- ceiling(wd / maxw)
        p <- rep(seq_along(wd), m)
        j <- sequence(m) - 1
        e <- c(e[p] + wd[p] * j / m[p], hi)
      }
    }
    st <- c(st, e[-length(e)])
    en <- c(en, e[-1L])
  }
  gl <- .fcr_gl(ng)
  half <- (en - st) / 2
  list(r = as.vector(outer(gl$x, half) + rep((st + en) / 2, each = ng)),
       w = as.vector(outer(gl$w, half)))
}

# log(u), from s = 1 - u where u is close to 1.
.fcr_log_u <- function(u, s) ifelse(s < 0.5, log1p(-s), log(u))

# Coefficients of the derivatives of the Gumbel generator inverse
# psi(t) = exp(-t^a): psi^(k)(t) = psi(t) sum_m c_m t^(m a - k), m = 0..k,
# and their derivatives with respect to a.
.fcr_gumbel_coef <- function(k, a) {
  cc <- 1
  dc <- 0
  for (kk in seq_len(k) - 1L) {
    m <- seq_along(cc) - 1
    nc <- numeric(length(cc) + 1L)
    nd <- nc
    top <- seq_along(cc) + 1L
    nc[top] <- nc[top] - a * cc
    nd[top] <- nd[top] - cc - a * dc
    bot <- seq_along(cc)
    nc[bot] <- nc[bot] + (m * a - kk) * cc
    nd[bot] <- nd[bot] + m * cc + (m * a - kk) * dc
    cc <- nc
    dc <- nd
  }
  list(c = cc, dc = dc, m = seq_along(cc) - 1)
}

# Coefficients of the Eulerian polynomial A_n(x) (A_0 = A_1 = 1,
# A_2 = 1 + x, A_3 = 1 + 4 x + x^2, ...), lowest power first.
.fcr_eulerian <- function(n) {
  a <- 1
  if (n <= 1L) return(a)
  for (nn in 2:n) {
    m <- seq_len(nn) - 1
    prev <- c(a, 0)
    prev_m1 <- c(0, a)
    a <- (m + 1) * prev + (nn - m) * prev_m1
  }
  a
}

# Pieces of an Archimedean copula C(u) = psi(sum phi(u_j)) for the
# tail-accurate volumes: the generator phi, its derivative and their
# derivatives with respect to theta (functions of u and s = 1 - u, accurate
# for u close to 1), and lpsi(t, k): the log of (-1)^k psi^(k)(t) (k >= 1,
# positive) and its derivative with respect to theta. `sing` is the
# singular point of psi^(k) on the real line (or a point to its right) and
# maxw the widest panel in t.
.fcr_arch_pieces <- function(family, theta) {
  if (family == "clayton") {
    # phi(u) = u^-theta - 1, psi(t) = (1 + t)^(-1 / theta)
    return(list(
      phi = function(u, s) expm1(-theta * .fcr_log_u(u, s)),
      dphi = function(u, s) -theta * exp((-theta - 1) * .fcr_log_u(u, s)),
      phi_th = function(u, s) {
        lu <- .fcr_log_u(u, s)
        -lu * exp(-theta * lu)
      },
      dphi_th = function(u, s) {
        lu <- .fcr_log_u(u, s)
        exp((-theta - 1) * lu) * (theta * lu - 1)
      },
      lpsi = function(t, k) {
        a <- 1 / theta + seq_len(k) - 1
        l1 <- log1p(t)
        list(lg = sum(log(a)) - (1 / theta + k) * l1,
             dlg = -sum(1 / a) / theta^2 + l1 / theta^2)
      },
      sing = -1, maxw = Inf))
  }
  if (family == "gumbel") {
    # phi(u) = (-log u)^theta, psi(t) = exp(-t^(1 / theta))
    al <- 1 / theta
    return(list(
      phi = function(u, s) (-.fcr_log_u(u, s))^theta,
      dphi = function(u, s) {
        l <- -.fcr_log_u(u, s)
        -theta * l^(theta - 1) / u
      },
      phi_th = function(u, s) {
        l <- -.fcr_log_u(u, s)
        ifelse(l > 0, l^theta * log(l), 0)
      },
      dphi_th = function(u, s) {
        l <- -.fcr_log_u(u, s)
        ifelse(l > 0, -theta * l^(theta - 1) / u * (1 / theta + log(l)), 0)
      },
      lpsi = function(t, k) {
        cf <- .fcr_gumbel_coef(k, al)
        keep <- cf$c != 0
        cm <- cf$c[keep]
        dm <- cf$dc[keep]
        m <- cf$m[keep]
        lt <- log(t)
        ta <- exp(al * lt)
        L <- outer(lt, m * al - k) + rep(log(abs(cm)), each = length(t))
        mx <- apply(L, 1L, max)
        E <- exp(L - mx)
        se <- rowSums(E)
        dl_da <- -ta * lt + rowSums(E * (rep(dm / cm, each = length(t)) +
                                           outer(lt, m))) / se
        list(lg = -ta + mx + log(se), dlg = -dl_da / theta^2)
      },
      sing = 0, maxw = Inf))
  }
  # Frank: phi(u) = -log(expm1(-theta u) / expm1(-theta)),
  # psi(t) = -log(1 - c e^-t) / theta with c = 1 - e^-theta; then
  # (-1)^k psi^(k)(t) = Li_(1-k)(x) / theta, x = c e^-t, and
  # Li_(-n)(x) = x A_n(x) / (1 - x)^(n + 1)
  h <- function(a) ifelse(abs(a) < 1e-8, 1 - a / 2, a / expm1(a))
  dh <- function(a) ifelse(abs(a) < 0.1,
                           -1 / 2 + a / 6 - a^3 / 180 + a^5 / 5040 -
                             a^7 / 151200,
                           (expm1(a) - a * exp(a)) / expm1(a)^2)
  B <- function(a) ifelse(abs(a) < 1e-8, 1 + a / 2, a / -expm1(-a))
  lc <- if (theta > 0) log(-expm1(-theta)) else log(expm1(-theta))
  sc <- sign(theta)
  ex <- 1 / expm1(theta)
  list(
    phi = function(u, s) {
      out <- numeric(length(u))
      up <- s < 0.5
      out[up] <- -log1p(-expm1(theta * s[up]) / expm1(theta))
      out[!up] <- log(expm1(-theta) / expm1(-theta * u[!up]))
      out
    },
    dphi = function(u, s) -h(theta * u) / u,
    phi_th = function(u, s) {
      out <- numeric(length(u))
      up <- s < 0.5
      su <- s[up]
      q <- expm1(theta * su) / expm1(theta)
      omq <- exp(theta * su) * expm1(theta * (1 - su)) / expm1(theta)
      out[up] <- q * (B(theta * su) - B(theta)) / (theta * omq)
      out[!up] <- (h(theta) - h(theta * u[!up])) / theta
      out
    },
    dphi_th = function(u, s) -dh(theta * u),
    lpsi = function(t, k) {
      A <- .fcr_eulerian(k - 1L)
      lx <- lc - t
      x <- sc * exp(lx)
      l1mx <- if (theta > 0) log(-expm1(-t) + exp(-theta - t)) else
        log1p(exp(lx))
      pw <- seq_along(A) - 1
      pa <- as.vector(outer(x, pw, `^`) %*% A)
      dpa <- if (length(A) > 1L) {
        as.vector(outer(x, pw[-1L] - 1, `^`) %*% (A[-1L] * pw[-1L]))
      } else 0
      dx <- x * ex
      list(lg = lx - log(abs(theta)) + log(pa) - k * l1mx,
           dlg = ex - 1 / theta + dpa / pa * dx + k * dx / exp(l1mx))
    },
    sing = 0, maxw = 2)
}

# Integral of (-1)^m psi^(m)(T + sum_i t_i) over the box of t_i in
# [lo_i, lo_i + d_i] (m >= the number of intervals), and of its derivative
# with respect to theta (th).
.fcr_cop_box <- function(ar, m, T, lo, d) {
  k <- length(lo)
  # (generator values that overflow, at extreme parameters)
  if (!all(is.finite(c(T, lo, d)))) return(list(v = NaN, th = NaN))
  if (k == 0L) {
    g <- ar$lpsi(T, m)
    v <- exp(g$lg)
    return(list(v = v, th = if (v > 0) v * g$dlg else 0))
  }
  if (any(d <= 0)) return(list(v = 0, th = 0))
  if (k > 2L) {
    nd <- .fcr_cop_nodes(d[k], T + sum(lo) - ar$sing, maxw = ar$maxw)
    v <- 0
    th <- 0
    for (i in seq_along(nd$r)) {
      b <- .fcr_cop_box(ar, m, T + lo[k] + nd$r[i], lo[-k], d[-k])
      v <- v + nd$w[i] * b$v
      th <- th + nd$w[i] * b$th
    }
    return(list(v = v, th = th))
  }
  if (k == 1L) {
    nd <- .fcr_cop_nodes(d, T + lo - ar$sing, maxw = ar$maxw)
    K <- 1
  } else {
    # t_1 + t_2 - lo_1 - lo_2 = r has the trapezoid density
    # min(r, d_min, d_1 + d_2 - r)
    D <- d[1L] + d[2L]
    nd <- .fcr_cop_nodes(D, T + sum(lo) - ar$sing, breaks = range(d),
                         maxw = ar$maxw)
    K <- pmax(pmin(nd$r, min(d), D - nd$r), 0)
  }
  g <- ar$lpsi(T + sum(lo) + nd$r, m)
  e <- exp(g$lg + log(nd$w * K))
  e[!is.finite(e)] <- 0
  et <- e * g$dlg
  et[!is.finite(et)] <- 0
  list(v = sum(e), th = sum(et))
}

# Tail-accurate volume of one row for an Archimedean copula: u and s = 1 - u
# at F(y) (ub, sb) and F(y - 1) (ua, sa); type "n" (upper tail,
# integrated), "w" (differenced over its corners) or "f" (left out).
# Returns P and its derivatives with respect to F(y) (Ahi), F(y - 1) (Alo)
# and theta (Dth).
.fcr_cop_row_arch <- function(ar, ub, sb, ua, sa, type, grad) {
  act <- type != "f"
  if (anyNA(c(ub[act], sb[act], ua[act], sa[act]))) return(list(P = NaN))
  J <- length(type)
  N <- which(type == "n")
  W <- which(type == "w")
  k <- length(N)
  lo <- ar$phi(ub[N], sb[N])
  hi <- ar$phi(ua[N], sa[N])
  d <- hi - lo
  P <- 0
  Dth <- 0
  Ahi <- Alo <- numeric(J)
  dlo <- dhi <- numeric(k)
  for (cc in 0:(2^length(W) - 1L)) {
    at_a <- as.logical(intToBits(cc))[seq_along(W)]
    uc <- ifelse(at_a, ua[W], ub[W])
    scc <- ifelse(at_a, sa[W], sb[W])
    if (any(uc <= 0)) next
    sg <- (-1)^sum(at_a)
    T <- sum(ar$phi(uc, scc))
    b0 <- .fcr_cop_box(ar, k, T, lo, d)
    P <- P + sg * b0$v
    if (!grad) next
    Dth <- Dth + sg * b0$th
    if (length(W)) {
      # d/dT of the box of order k is minus the box of order k + 1
      b1 <- .fcr_cop_box(ar, k + 1L, T, lo, d)
      Dth <- Dth - sg * b1$v * sum(ar$phi_th(uc, scc))
      du <- -sg * b1$v * ar$dphi(uc, scc)
      Alo[W[at_a]] <- Alo[W[at_a]] + du[at_a]
      Ahi[W[!at_a]] <- Ahi[W[!at_a]] + du[!at_a]
    }
    for (i in seq_len(k)) {
      bl <- .fcr_cop_box(ar, k, T + lo[i], lo[-i], d[-i])
      bh <- .fcr_cop_box(ar, k, T + hi[i], lo[-i], d[-i])
      dlo[i] <- dlo[i] - sg * bl$v
      dhi[i] <- dhi[i] + sg * bh$v
    }
  }
  if (grad) {
    # t = phi(u): the lower end of an interval in t is phi(F(y))
    Ahi[N] <- dlo * ar$dphi(ub[N], sb[N])
    Alo[N] <- dhi * ar$dphi(ua[N], sa[N])
    Dth <- Dth + sum(dlo * ar$phi_th(ub[N], sb[N]) +
                       dhi * ar$phi_th(ua[N], sa[N]))
  }
  list(P = P, Ahi = Ahi, Alo = Alo, Dth = Dth)
}

# Phi(b) - Phi(a) for a <= b, without cancellation in the upper tail.
.fcr_pnorm_diff <- function(a, b) {
  up <- a > 0
  out <- stats::pnorm(b) - stats::pnorm(a)
  out[up] <- stats::pnorm(a[up], lower.tail = FALSE) -
    stats::pnorm(b[up], lower.tail = FALSE)
  out
}

# Tail-accurate volume of one row for the Gaussian copula (two outcomes;
# arguments and value as .fcr_cop_row_arch()): the integral over one
# outcome's interval of the conditional probability of the other's,
#   P = int_a^b [Phi(w(z_b)) - Phi(w(z_a))] du, w(z) = (z - rho x) / sigma,
# x = qnorm(u). The interval is split at u = 0.5 and each part integrated
# over its tail probability (u below 0.5, 1 - u above), so that both
# tails keep their relative accuracy.
.fcr_cop_row_gauss <- function(rho, ub, sb, ua, sa, type, grad) {
  if (anyNA(c(ub, sb, ua, sa))) return(list(P = NaN))
  sig2 <- 1 - rho^2
  sig <- sqrt(sig2)
  zq <- function(u, s) ifelse(s < 0.5, stats::qnorm(s, lower.tail = FALSE),
                              stats::qnorm(u))
  dens <- function(z1, z2) {
    out <- exp(-(rho^2 * (z1^2 + z2^2) - 2 * rho * z1 * z2) / (2 * sig2)) /
      sig
    out[!is.finite(z1) | !is.finite(z2)] <- 0
    out
  }
  # The outcome integrated over: the one with the narrower interval (the
  # other's conditional probability then changes least across it)
  side <- ifelse(ua >= 0.5, 1L, ifelse(ub <= 0.5, -1L, 0L))
  width <- ifelse(side == 1L, sa - sb, ub - ua)
  i <- if (width[2L] < width[1L]) 2L else 1L
  o <- 3L - i
  za <- zq(ua[o], sa[o])
  zb <- zq(ub[o], sb[o])
  wz <- function(z, x) if (is.finite(z)) (z - rho * x) / sig else
    rep(z, length(x))
  D <- function(x) {
    out <- numeric(length(x))
    f <- is.finite(x)
    out[f] <- .fcr_pnorm_diff(wz(za, x[f]), wz(zb, x[f]))
    out
  }
  # Parts of the interval: probabilities s from s0 to s1 in the lower
  # (x = qnorm(s)) or upper (x = qnorm(s, lower = FALSE)) tail
  parts <- switch(as.character(side[i]),
    "1" = list(c(1, sb[i], sa[i])),
    "-1" = list(c(-1, ua[i], ub[i])),
    list(c(-1, ua[i], 0.5), c(1, sb[i], 0.5)))
  x <- numeric(0)
  w <- numeric(0)
  for (pt in parts) {
    if (!(pt[3L] > pt[2L])) next
    x1 <- stats::qnorm(pt[3L], lower.tail = pt[1L] < 0)
    # Panels geometric in s, fine enough that the conditional probability
    # changes little within one: w moves by about rho log(ratio) /
    # (x sigma) across a panel, and the tail probability of w by a factor
    # of about exp(-w dw)
    wmax <- max(1, abs(c(za, zb)[is.finite(c(za, zb))] - rho * x1) / sig)
    lr <- min(log(2), max(1, abs(x1)) * sig / (max(abs(rho), 1e-3) * wmax))
    nd <- .fcr_cop_nodes(pt[3L] - pt[2L], pt[2L], lratio = lr, maxp = 2000L)
    s <- pt[2L] + nd$r
    x <- c(x, stats::qnorm(s, lower.tail = pt[1L] < 0))
    w <- c(w, nd$w)
  }
  Dx <- D(x)
  P <- sum(w * Dx)
  Ahi <- Alo <- numeric(2L)
  Dth <- 0
  if (grad) {
    dw <- function(zo, z) {
      out <- stats::dnorm((zo - rho * z) / sig) * (rho * zo - z) / sig^3
      out[!is.finite(zo) | !is.finite(z)] <- 0
      out
    }
    Dth <- sum(w * (dw(zb, x) - dw(za, x)))
    Ahi[o] <- sum(w * dens(x, zb))
    Alo[o] <- -sum(w * dens(x, za))
    Ahi[i] <- D(zq(ub[i], sb[i]))
    Alo[i] <- -D(zq(ua[i], sa[i]))
  }
  list(P = P, Ahi = Ahi, Alo = Alo, Dth = Dth)
}

# The equation's distribution function at q (on the probability scale).
.fcr_mv_cdf <- function(eq, v, q) {
  exp(.fcr_fm_logcdf(eq$dens_args, list(mu = v$mu, alpha = v$alpha,
                                        sigma = v$sigma), eq$family, q))
}

# Rectangle bounds of the counts Y (n x J; NA = marginalized) at the
# parameters p, with random terms `adds` in the equations' log means:
# F(y) (hi), F(y - 1) (lo) and the upper tails 1 - F (Ghi, Glo).
.fcr_mv_cop_bounds <- function(p, E, Y, adds = NULL) {
  n <- nrow(Y)
  Fhi <- matrix(NA_real_, n, E$J)
  Flo <- Fhi
  Ghi <- Fhi
  Glo <- Fhi
  for (j in seq_len(E$J)) {
    if (all(is.na(Y[, j]))) next
    eq <- E$eqs[[j]]
    v <- .fcr_mv_eq_values(p, eq, adds[[j]])
    if (is.null(v)) return(NULL)
    mg <- .fcr_mv_cop_margin(eq, v, Y[, j])
    Fhi[, j] <- mg$hi$F
    Flo[, j] <- mg$lo$F
    Ghi[, j] <- mg$hi$G
    Glo[, j] <- mg$lo$G
  }
  list(hi = Fhi, lo = Flo, Ghi = Ghi, Glo = Glo)
}

# F(y) (hi) and F(y - 1) (lo) of an equation with their upper tails
# G = 1 - F, accurate where F > 0.5 (see .fcr_mv_upper()), and with grad
# the derivatives of F as .fcr_mv_cdf_grad() (from the upper tail where it
# is used).
.fcr_mv_cop_margin <- function(eq, v, y, grad = FALSE, upper = TRUE) {
  if (grad) {
    cg <- .fcr_mv_cdf_grad(eq, v, y)
  } else {
    cg <- list(hi = list(F = .fcr_mv_cdf(eq, v, y)),
               lo = list(F = .fcr_mv_cdf(eq, v, y - 1)))
  }
  cg$hi$G <- 1 - cg$hi$F
  cg$lo$G <- 1 - cg$lo$F
  if (upper) cg <- .fcr_mv_cop_upper(cg, eq, v, y, seq_along(y), grad)
  cg
}

# The margins cg of .fcr_mv_cop_margin() with accurate upper tails (and,
# with grad, their derivatives) for the rows `rows`.
.fcr_mv_cop_upper <- function(cg, eq, v, y, rows, grad = FALSE) {
  n <- length(y)
  sub <- function(x, i) if (!is.null(x)) rep_len(x, n)[i]
  for (side in c("hi", "lo")) {
    q <- if (side == "hi") y else y - 1
    Fq <- cg[[side]]$F
    idx <- rows[!is.na(Fq[rows]) & !is.na(q[rows]) & Fq[rows] > 0.5]
    if (!length(idx)) next
    up <- .fcr_mv_upper(eq, list(mu = v$mu[idx], alpha = sub(v$alpha, idx),
                                 sigma = sub(v$sigma, idx)),
                        q[idx], Fq[idx], grad)
    cg[[side]]$G[idx] <- up$G
    if (grad && any(up$tail)) {
      t <- idx[up$tail]
      for (nm in c("eta", "a", "s")) {
        cg[[side]][[nm]][t] <- -up[[nm]][up$tail]
      }
    }
  }
  cg
}

# Upper tails G(q) = P(Y > q) of an equation (values v for these rows) for
# rows with F(q) > 0.5. Where 1 - F(q) < 1e-3 it has lost digits: G is
# then the family's own upper tail function when it has one, and otherwise
# (and, with grad, for the derivatives of G with respect to the log mean and
# log distribution parameters) the sum of the probabilities upward from
# q + 1. `tail` marks the rows whose derivatives come from that sum.
.fcr_mv_upper <- function(eq, v, q, Fq, grad = FALSE, max_terms = 20000L) {
  n <- length(q)
  fam <- .fcr_get_family(eq$family)
  G <- 1 - Fq
  small <- G < 1e-3
  out <- list(G = G, eta = numeric(n), a = numeric(n), s = numeric(n),
              tail = rep(FALSE, n))
  if (!any(small)) return(out)
  closed <- !is.null(fam$cdf)
  if (closed) {
    i <- which(small)
    gc <- suppressWarnings(exp(fam$cdf(q[i], v$mu[i], v$alpha[i], v$sigma[i],
                                       FALSE)))
    good <- is.finite(gc)
    out$G[i[good]] <- gc[good]
  }
  if (closed && !grad) return(out)
  idx <- which(small)
  m <- length(idx)
  S <- numeric(m)
  de <- da <- ds <- numeric(m)
  prev <- rep(Inf, m)
  conv <- rep(FALSE, m)
  act <- seq_len(m)
  k <- q[idx] + 1
  sub <- function(x, i) if (!is.null(x)) x[i]
  for (it in seq_len(max_terms)) {
    r <- idx[act]
    kk <- k[act]
    lf <- .fcr_quiet_dens(as.vector(do.call(fam$dens, c(
      list(y = kk, predicted = v$mu[r], alpha = sub(v$alpha, r),
           sigma = sub(v$sigma, r)), eq$dens_args, list(log = TRUE)))),
      log = TRUE)
    f <- exp(lf)
    f[!is.finite(f)] <- 0
    S[act] <- S[act] + f
    if (grad) {
      sc <- .fcr_family_score(eq$family, kk, v$mu[r], sub(v$alpha, r),
                              sub(v$sigma, r), dens_args = eq$dens_args)
      add <- function(x) {
        x <- rep_len(x, length(act)) * f
        x[!is.finite(x)] <- 0
        x
      }
      de[act] <- de[act] + add(sc$eta)
      if (!is.null(sc$alpha)) da[act] <- da[act] + add(sc$alpha)
      if (!is.null(sc$sigma)) ds[act] <- ds[act] + add(sc$sigma)
    }
    done <- (f <= 1e-17 * S[act] & f <= prev[act]) | (S[act] == 0 & it > 1)
    prev[act] <- f
    k[act] <- kk + 1
    conv[act[done]] <- TRUE
    act <- act[!done]
    if (!length(act)) break
  }
  if (!closed) out$G[idx[conv]] <- S[conv]
  if (grad) {
    out$tail[idx[conv]] <- TRUE
    out$eta[idx] <- de
    out$a[idx] <- da
    out$s[idx] <- ds
  }
  out
}

# Mean and distribution parameters of an equation (NULL where invalid).
.fcr_mv_eq_values <- function(p, eq, add_eta = NULL) {
  S <- eq$S
  n <- length(eq$y)
  lp <- as.vector(eq$X %*% p[S$idx_beta]) + eq$off
  if (!is.null(add_eta)) lp <- lp + add_eta
  mu <- exp(lp)
  dp <- .fcr_alloc_dispar(p, S, eq$W1, eq$W2, n)
  if (!all(.fcr_family_valid(eq$family, mu, dp$alpha, dp$sigma))) {
    return(NULL)
  }
  list(mu = mu, alpha = dp$alpha, sigma = dp$sigma)
}

# F(y) and F(y - 1) of an equation and their derivatives with respect to
# the log mean and the log distribution parameters:
#   d F(q) / d theta = sum_{k <= q} f(k) d log f(k) / d theta,
# summed over k = 0, ..., max(y) for the rows that need each k.
.fcr_mv_cdf_grad <- function(eq, v, y) {
  n <- length(y)
  fam <- .fcr_get_family(eq$family)
  z <- function() numeric(n)
  hi <- list(F = z(), eta = z(), a = z(), s = z())
  lo <- list(F = z(), eta = z(), a = z(), s = z())
  for (k in 0:max(y)) {
    act <- which(y >= k)
    sub <- function(x) if (!is.null(x)) x[act]
    yk <- rep(k, length(act))
    lf <- .fcr_quiet_dens(as.vector(do.call(fam$dens, c(
      list(y = yk, predicted = v$mu[act], alpha = sub(v$alpha),
           sigma = sub(v$sigma)), eq$dens_args, list(log = TRUE)))),
      log = TRUE)
    f <- exp(lf)
    sc <- .fcr_family_score(eq$family, yk, v$mu[act], sub(v$alpha),
                            sub(v$sigma), dens_args = eq$dens_args)
    pieces <- list(F = f, eta = f * sc$eta,
                   a = if (!is.null(sc$alpha)) f * sc$alpha else 0,
                   s = if (!is.null(sc$sigma)) f * sc$sigma else 0)
    below <- y[act] > k              # k <= y - 1
    for (nm in names(pieces)) {
      val <- rep_len(pieces[[nm]], length(act))
      val[!is.finite(val)] <- 0
      hi[[nm]][act] <- hi[[nm]][act] + val
      lo[[nm]][act[below]] <- lo[[nm]][act[below]] + val[below]
    }
  }
  hi$F <- pmin(hi$F, 1)
  lo$F <- pmin(lo$F, 1)
  list(hi = hi, lo = lo)
}

# The copula C(U) and its derivatives with respect to each u_j (n x J) and
# to theta, for the rows of U (no u equal to 0).
.fcr_cop_derivs <- function(U, family, theta) {
  U <- pmin(pmax(U, 1e-300), 1)
  J <- ncol(U)
  if (family == "gaussian") {
    a <- stats::qnorm(pmin(U[, 1L], 1 - 1e-16))
    b <- stats::qnorm(pmin(U[, 2L], 1 - 1e-16))
    a[U[, 1L] >= 1] <- Inf
    b[U[, 2L] >= 1] <- Inf
    r <- theta
    sr <- sqrt(1 - r^2)
    C <- .fcr_pbvnorm(a, b, r)
    dU <- cbind(ifelse(is.finite(a), stats::pnorm((b - r * a) / sr), 0),
                ifelse(is.finite(b), stats::pnorm((a - r * b) / sr), 0))
    dU[!is.finite(a) & !is.finite(b), ] <- 0
    dU[is.finite(a) & !is.finite(b), 1L] <- 1
    dU[!is.finite(a) & is.finite(b), 2L] <- 1
    q <- ifelse(is.finite(a) & is.finite(b),
                (a^2 - 2 * r * a * b + b^2) / (1 - r^2), Inf)
    dth <- exp(-q / 2) / (2 * pi * sr)
    return(list(C = C, dU = dU, dth = dth))
  }
  if (family == "clayton") {
    t <- expm1(-theta * log(U))          # u^(-theta) - 1
    s1 <- 1 + rowSums(t)
    C <- exp(-log(s1) / theta)
    dU <- exp((-1 / theta - 1) * log(s1)) * U^(-theta - 1)
    lu <- log(U) * U^(-theta)
    dth <- C * (log(s1) / theta^2 + rowSums(lu) / (theta * s1))
    return(list(C = C, dU = dU, dth = dth))
  }
  if (family == "gumbel") {
    l <- -log(U)
    lt <- l^theta
    s <- rowSums(lt)
    C <- exp(-s^(1 / theta))
    dU <- C * s^(1 / theta - 1) * l^(theta - 1) / U
    dU[!is.finite(dU)] <- 0
    llog <- ifelse(l > 0, lt * log(l), 0)
    dth <- -C * s^(1 / theta) * (-log(s) / theta^2 +
                                   rowSums(llog) / (theta * s))
    dth[!is.finite(dth)] <- 0
    return(list(C = C, dU = dU, dth = dth))
  }
  # frank
  if (abs(theta) < 1e-6) {
    P <- apply(U, 1L, prod)
    dU <- P / U
    dth <- P * ((J - 1 - rowSums(U)) / 2 + P / 2)
    return(list(C = P, dU = dU, dth = dth))
  }
  e <- expm1(-theta * U)
  e0 <- expm1(-theta)
  N <- apply(e, 1L, prod) / e0^(J - 1L)
  C <- -log1p(N) / theta
  dU <- (N / (1 + N)) * exp(-theta * U) / e
  dlogN <- rowSums(-U * exp(-theta * U) / e) + (J - 1) * exp(-theta) / e0
  dth <- log1p(N) / theta^2 - N * dlogN / (theta * (1 + N))
  list(C = C, dU = dU, dth = dth)
}

# d theta / d eta of the copula parameter.
.fcr_cop_dtheta <- function(eta, family, J) {
  switch(family,
    gaussian = 1 - tanh(eta)^2,
    frank = if (J == 2L) 1 else exp(eta),
    clayton = exp(eta),
    gumbel = exp(eta))
}

# Row log-probabilities (unweighted) of the observed counts and, with
# grad, their derivatives with respect to each equation's log mean and log
# distribution parameters (one list per equation, as .fcr_mv_eq_parts())
# and to the copula parameter. `adds` are random terms in the log means.
.fcr_mv_rowparts <- function(p, E, adds = NULL, grad = TRUE) {
  if (is.null(E$cop)) {
    ll <- numeric(E$n)
    parts <- vector("list", E$J)
    for (j in seq_len(E$J)) {
      pt <- .fcr_mv_eq_parts(p, E$eqs[[j]], adds[[j]], grad = grad)
      ll <- ll + pt$ll
      parts[[j]] <- pt
    }
    return(list(ll = ll, parts = parts))
  }
  J <- E$J
  n <- E$n
  fam <- E$cop$family
  eta <- p[E$cop$idx]
  theta <- .fcr_cop_theta(eta, fam, J)
  vals <- lapply(seq_len(J), function(j)
    .fcr_mv_eq_values(p, E$eqs[[j]], adds[[j]]))
  if (any(vapply(vals, is.null, logical(1)))) {
    return(list(ll = rep(NA_real_, n), parts = NULL, s_cop = NULL))
  }
  cg <- lapply(seq_len(J), function(j)
    .fcr_mv_cop_margin(E$eqs[[j]], vals[[j]], E$eqs[[j]]$y, grad = grad,
                       upper = FALSE))
  col <- function(side, nm) matrix(vapply(cg, function(x) x[[side]][[nm]],
                                          numeric(n)), n)
  # Accurate upper tails only for the rows whose volume is integrated
  Gfun <- function(rows) {
    for (j in seq_len(J)) {
      cg[[j]] <<- .fcr_mv_cop_upper(cg[[j]], E$eqs[[j]], vals[[j]],
                                    E$eqs[[j]]$y, rows, grad)
    }
    list(hi = col("hi", "G")[rows, , drop = FALSE],
         lo = col("lo", "G")[rows, , drop = FALSE])
  }
  # P, and with grad the derivatives of P with respect to each F_j(y_j)
  # (Ahi) and F_j(y_j - 1) (Alo) and to theta
  vol <- .fcr_cop_volume(col("hi", "F"), col("lo", "F"), fam, theta,
                         col("hi", "G"), col("lo", "G"), grad = grad,
                         Gfun = Gfun)
  P <- vol$P
  if (!grad) return(list(ll = log(pmax(P, 1e-300))))
  Ahi <- vol$Ahi
  Alo <- vol$Alo
  Dth <- vol$Dth
  pos <- P > 1e-300
  inv <- ifelse(pos, 1 / P, 0)
  parts <- lapply(seq_len(J), function(j) {
    h <- cg[[j]]$hi
    l <- cg[[j]]$lo
    list(s_eta = inv * (Ahi[, j] * h$eta + Alo[, j] * l$eta),
         s_a = inv * (Ahi[, j] * h$a + Alo[, j] * l$a),
         s_s = inv * (Ahi[, j] * h$s + Alo[, j] * l$s))
  })
  list(ll = log(pmax(P, 1e-300)), parts = parts,
       s_cop = inv * Dth * .fcr_cop_dtheta(eta, fam, J))
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
