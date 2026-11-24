# Design of hierarchical (mixed) count models: random-effect terms written
# as in lme4, (1 + x | county) + (1 | county:segment) + (1 | year).
#
# reformulas::findbars() returns the random-effect terms (with nesting,
# a/b, expanded to (1 | a) + (1 | a:b), and double bars, ||, split into
# uncorrelated terms); reformulas::nobars() returns the fixed formula,
# which goes through prepare_model_data() as for countreg().
#
# Each term l has a grouping factor with J_l levels and q_l effects (the
# columns of its own model matrix X_l). Its effects are b_lj = L_l u_lj,
# with u_lj ~ N(0, I) and L_l lower triangular, so Sigma_l = L_l L_l'.
# The diagonal of L_l is unconstrained: a variance can reach zero (a
# singular fit) without a boundary in the parameters, and the sign of a
# column of L_l does not change Sigma_l.
#
# The random effects of all terms are stacked as u (length Q = sum J_l q_l,
# group-major within a term), so the linear predictor is
#   eta = X beta + offset + Z Lambda u,
# with Z the sparse n x Q design and Lambda = blockdiag(I_J (x) L_l).

# Split a formula into its fixed formula and random-effect terms.
.fcr_hier_split <- function(formula) {
  rlang::check_installed("reformulas",
                         reason = "for random-effect terms such as (1 | site).")
  f <- stats::as.formula(formula)
  bars <- reformulas::findbars(f)
  fixed <- reformulas::nobars(f)
  if (!length(bars) && !.fcr_has_smooth(fixed)) {
    stop("`formula` has no random-effect terms, such as (1 | site). Use ",
         "countreg() for a model without them.", call. = FALSE)
  }
  environment(fixed) <- environment(f)
  list(fixed = fixed, bars = bars, env = environment(f))
}

# Label of a grouping expression, e.g. "county:segment".
.fcr_hier_label <- function(gexpr) {
  paste(deparse(gexpr, width.cutoff = 500L), collapse = "")
}

# Grouping factor of a term on data `d`: its variables are made factors,
# so `a:b` is their interaction.
.fcr_hier_group <- function(gexpr, d, env) {
  vars <- all.vars(gexpr)
  miss <- setdiff(vars, names(d))
  if (length(miss)) {
    stop("Grouping variable(s) not found: ", toString(miss), call. = FALSE)
  }
  fl <- lapply(d[vars], function(v) factor(v))
  g <- eval(gexpr, fl, env)
  factor(g)
}

# Model matrix of a term's effects (its left-hand side) on data `d`.
.fcr_hier_term_matrix <- function(t, d) {
  mf <- .fcr_model_frame_xlev(t$terms, d, t$xlevels)
  X <- stats::model.matrix(t$terms, mf, contrasts.arg = t$contrasts)
  if (!identical(colnames(X), t$cnms)) {
    stop("`newdata` gives different random-effect columns for (",
         t$label, ") than the training data.", call. = FALSE)
  }
  unname(X)
}

# Set up the random-effect terms on the estimation rows `d`.
.fcr_hier_terms <- function(bars, d, env) {
  terms <- lapply(bars, function(b) {
    ff <- stats::as.formula(call("~", b[[2]]), env)
    mf <- stats::model.frame(ff, d, na.action = stats::na.pass)
    tt <- stats::terms(mf)
    X <- stats::model.matrix(tt, mf)
    if (!ncol(X)) {
      stop("The random-effect term (", .fcr_hier_label(b), ") has no ",
           "effects.", call. = FALSE)
    }
    if (any(!is.finite(X))) {
      stop("The random-effect term (", .fcr_hier_label(b), ") has ",
           "non-finite values.", call. = FALSE)
    }
    g <- droplevels(.fcr_hier_group(b[[3]], d, env))
    q <- ncol(X)
    list(label = .fcr_hier_label(b),
         group = .fcr_hier_label(b[[3]]), gexpr = b[[3]],
         terms = tt, xlevels = stats::.getXlevels(tt, mf),
         contrasts = attr(X, "contrasts"), cnms = colnames(X),
         X = unname(X), g = g, levels = levels(g), q = q, J = nlevels(g),
         ntheta = as.integer(q * (q + 1L) / 2L))
  })
  .fcr_hier_offsets(terms)
}

# Column offsets of each term in u, and of each term's theta in the vector
# of all covariance parameters.
.fcr_hier_offsets <- function(terms) {
  off <- 0L
  toff <- 0L
  for (i in seq_along(terms)) {
    terms[[i]]$uoff <- off
    terms[[i]]$toff <- toff
    off <- off + as.integer(terms[[i]]$J * terms[[i]]$q)
    toff <- toff + as.integer(terms[[i]]$ntheta)
  }
  terms
}

# Random-effect terms of the penalized smooths of the mean (the mixed-model
# form of mgcv::smooth2random()): each penalty's range space is a block of
# independent N(0, sigma^2) coefficients, sigma^2 = 1 / lambda, and the
# null space becomes fixed columns. Returns the fixed columns, the terms
# and, per smooth, what is needed to return to its original coefficients
# (beta = U (D * c(b_1, ..., b_m, beta_fixed))).
.fcr_hier_smooth_terms <- function(ents, X) {
  fixed <- list()
  terms <- list()
  info <- list()
  for (si in seq_along(ents)) {
    e <- ents[[si]]
    sm <- e$sm
    if (length(e$S) > 1L && !inherits(sm, "t2.smooth")) {
      stop("The smooth ", e$label, " has several penalties; in ",
           "countreg_hier() use s() or t2() (te() and ti() are not ",
           "available).", call. = FALSE)
    }
    sm$X <- X[, e$cols, drop = FALSE]
    sm$S <- e$S
    r <- mgcv::smooth2random(sm, "", type = 2)
    Xf <- r$Xf
    if (is.null(Xf)) Xf <- matrix(0, nrow(X), 0)
    Xf <- as.matrix(Xf)
    if (ncol(Xf)) colnames(Xf) <- paste0(e$label, ".Fx", seq_len(ncol(Xf)))
    fixed[[si]] <- Xf
    ti <- integer(0)
    for (b in seq_along(r$rand)) {
      Zr <- as.matrix(r$rand[[b]])
      terms[[length(terms) + 1L]] <- list(
        type = "smooth", label = e$label, group = e$label, Zs = unname(Zr),
        J = ncol(Zr), q = 1L, ntheta = 1L, cnms = "(smooth)", smooth = si)
      ti <- c(ti, length(terms))
    }
    k <- length(e$cols)
    info[[si]] <- list(label = e$label, cols = e$cols, n_fixed = ncol(Xf),
                       n_rand = vapply(r$rand, ncol, 0L), terms = ti,
                       U = r$trans.U %||% diag(k),
                       D = r$trans.D %||% rep(1, k))
  }
  list(X = do.call(cbind, fixed), terms = terms, info = info)
}

# Sparse design Z (n x Q) of the stacked random effects.
.fcr_hier_Z <- function(terms, n, X_list = NULL, g_list = NULL) {
  blocks <- lapply(seq_along(terms), function(i) {
    t <- terms[[i]]
    if (identical(t$type, "smooth")) {
      return(methods::as(Matrix::Matrix(t$Zs, sparse = TRUE),
                         "CsparseMatrix"))
    }
    X <- if (is.null(X_list)) t$X else X_list[[i]]
    gi <- if (is.null(g_list)) as.integer(t$g) else g_list[[i]]
    q <- t$q
    ok <- !is.na(gi)
    rows <- rep(seq_len(n), q)[rep(ok, q)]
    cols <- (rep((gi - 1L) * q, q) + rep(seq_len(q), each = n))[rep(ok, q)]
    Matrix::sparseMatrix(i = rows, j = cols, x = as.vector(X)[rep(ok, q)],
                         dims = c(n, t$J * q))
  })
  Z <- if (length(blocks) == 1L) blocks[[1]] else do.call(cbind, blocks)
  methods::as(Z, "CsparseMatrix")
}

# Lower-triangular factor L of one term from its theta (column-major lower
# triangle, diagonal included).
.fcr_hier_L <- function(th, q) {
  L <- matrix(0, q, q)
  L[lower.tri(L, diag = TRUE)] <- th
  L
}

# Sparsity template of Lambda: for each term, the positions of the lower
# triangle of L in every group's block.
.fcr_hier_lambda_template <- function(terms) {
  i <- j <- k <- list()
  for (t in terms) {
    q <- t$q
    lt <- which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE)
    base <- t$uoff + (seq_len(t$J) - 1L) * q
    i[[length(i) + 1L]] <- rep(base, each = nrow(lt)) + lt[, 1]
    j[[length(j) + 1L]] <- rep(base, each = nrow(lt)) + lt[, 2]
    k[[length(k) + 1L]] <- rep(t$toff + seq_len(nrow(lt)), t$J)
  }
  Q <- sum(vapply(terms, function(t) t$J * t$q, 0))
  list(i = unlist(i), j = unlist(j), k = unlist(k), Q = Q)
}

.fcr_hier_lambda <- function(tmpl, theta) {
  Matrix::sparseMatrix(i = tmpl$i, j = tmpl$j, x = theta[tmpl$k],
                       dims = c(tmpl$Q, tmpl$Q))
}

# Starting theta: independent effects with standard deviation `sd`.
.fcr_hier_theta_start <- function(terms, sd = 0.5) {
  unlist(lapply(terms, function(t) {
    L <- diag(sd, t$q)
    L[lower.tri(L, diag = TRUE)]
  }))
}

# Covariance matrices of the terms from theta.
.fcr_hier_sigmas <- function(terms, theta) {
  lapply(terms, function(t) {
    L <- .fcr_hier_L(theta[t$toff + seq_len(t$ntheta)], t$q)
    S <- L %*% t(L)
    dimnames(S) <- list(t$cnms, t$cnms)
    S
  })
}

# Standard deviations and correlations of every term, as one named vector
# (lme4's naming: "sd_(Intercept)|site", "cor_lnaadt.(Intercept)|county").
.fcr_hier_vc_vector <- function(terms, theta) {
  S <- .fcr_hier_sigmas(terms, theta)
  out <- list()
  for (i in seq_along(terms)) {
    t <- terms[[i]]
    s <- sqrt(pmax(diag(S[[i]]), 0))
    v <- stats::setNames(s, paste0("sd_", t$cnms, "|", t$group))
    if (t$q > 1L) {
      R <- S[[i]] / outer(s, s)
      lt <- which(lower.tri(R), arr.ind = TRUE)
      v <- c(v, stats::setNames(
        R[lt], paste0("cor_", t$cnms[lt[, 1]], ".", t$cnms[lt[, 2]], "|",
                      t$group)))
    }
    out[[i]] <- v
  }
  unlist(out)
}

# Theta from standard deviations and correlations (the inverse of
# .fcr_hier_vc_vector(), with a non-negative diagonal).
.fcr_hier_theta_from_vc <- function(terms, vc) {
  pos <- 0L
  unlist(lapply(terms, function(t) {
    q <- t$q
    s <- vc[pos + seq_len(q)]
    pos <<- pos + q
    R <- diag(q)
    if (q > 1L) {
      nr <- q * (q - 1L) / 2L
      R[lower.tri(R)] <- vc[pos + seq_len(nr)]
      R[upper.tri(R)] <- t(R)[upper.tri(R)]
      pos <<- pos + nr
    }
    S <- R * outer(s, s)
    L <- tryCatch(t(chol(S)), error = function(e) NULL)
    if (is.null(L)) {
      # Singular: a pivoted factor that reproduces S
      e <- eigen(S, symmetric = TRUE)
      L <- t(qr.R(qr(t(e$vectors %*% diag(sqrt(pmax(e$values, 0)), q)))))
      L <- L %*% diag(sign(diag(L)) + (diag(L) == 0), q)
    }
    L[lower.tri(L, diag = TRUE)]
  }))
}
