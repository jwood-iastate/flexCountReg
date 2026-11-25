# Finite mixtures and latent classes: control settings, family checks and the
# free-parameter map.
#
# Naming: internal helpers of the finite-mixture engine use the prefix
# `.fcr_fm_`. (The prefix `.fcr_mix_` belongs to the numerical integration of
# the mixed-Poisson families, such as PLN, in mixture_quad.R: those families
# mix over a continuous random effect *within* one component.)

#' Control settings for finite-mixture and latent-class count models
#'
#' Settings of the estimation algorithm used by [countreg_mix()] and
#' [countreg_lc()]. `lc_control()` is the same function under the name used
#' with the latent-class interface; both return the same object.
#'
#' @param algorithm `"em"` (default): each start is improved by the
#'   expectation-maximization (EM) algorithm, and the best solutions are then
#'   polished by maximizing the observed-data log-likelihood directly (with
#'   the optimizer chosen by `method`). `"direct"`: every start is maximized
#'   directly, without EM iterations.
#' @param em_maxit maximum number of EM iterations per start.
#' @param em_tol relative convergence tolerance of the EM iterations. EM stops
#'   when the relative change in the observed log-likelihood is below
#'   `em_tol` *and* the largest change in a parameter is below
#'   `sqrt(em_tol)`.
#' @param mstep_maxit maximum number of quasi-Newton iterations in each
#'   numerical M-step (a generalized EM step: it need not reach the
#'   conditional maximum, but it must increase the likelihood).
#' @param n_polish number of distinct EM solutions (the best ones) that are
#'   polished by direct maximization.
#' @param boot_starts number of starts used in each refit of a bootstrap
#'   (`stderr = "boot"`), in addition to the original estimates.
#' @param small_class posterior class share below which a component is
#'   flagged as small in the output. This is a reporting threshold only: it
#'   never changes the fit.
#' @param start_spread spread of the response-based starting partitions (the
#'   standard deviation, in units of the standardized score, of the noise
#'   added to the score that orders the units into starting classes).
#' @param ... arguments passed from `lc_control()` to `mixture_control()`.
#'
#' @returns An object of class `fcr_mixture_control`.
#' @seealso [countreg_mix()], [countreg_lc()]
#' @examples
#' mixture_control(em_maxit = 200)
#' lc_control(algorithm = "direct")
#' @export
mixture_control <- function(algorithm = c("em", "direct"), em_maxit = 200L,
                            em_tol = 1e-6, mstep_maxit = 25L, n_polish = 3L,
                            boot_starts = 5L, small_class = 0.05,
                            start_spread = 0.5) {
  algorithm <- match.arg(algorithm)
  whole <- function(v, nm, lo) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < lo ||
        v != round(v)) {
      stop("`", nm, "` must be a single whole number of at least ", lo, ".",
           call. = FALSE)
    }
    as.integer(v)
  }
  pos <- function(v, nm, hi = Inf) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0 ||
        v >= hi) {
      stop("`", nm, "` must be a single positive number",
           if (is.finite(hi)) paste0(" below ", hi), ".", call. = FALSE)
    }
    as.numeric(v)
  }
  structure(list(
    algorithm = algorithm,
    em_maxit = whole(em_maxit, "em_maxit", 0),
    em_tol = pos(em_tol, "em_tol", 1),
    mstep_maxit = whole(mstep_maxit, "mstep_maxit", 1),
    n_polish = whole(n_polish, "n_polish", 1),
    boot_starts = whole(boot_starts, "boot_starts", 1),
    small_class = pos(small_class, "small_class", 1),
    start_spread = if (identical(start_spread, 0)) 0 else
      pos(start_spread, "start_spread")
  ), class = "fcr_mixture_control")
}

#' @rdname mixture_control
#' @export
lc_control <- function(...) mixture_control(...)

#' @export
print.fcr_mixture_control <- function(x, ...) {
  cat("Finite-mixture control settings\n")
  for (nm in names(x)) cat(" ", nm, ":", format(x[[nm]]), "\n")
  invisible(x)
}

# ---- Families ----------------------------------------------------------------

# Status of each registry family as a mixture component. "validated":
# analytic scores, independent likelihood checks and simulation recovery in
# the test suite. "experimental": the component density, mean, variance,
# distribution function and random generator come from the same registry
# entries as countreg() (one-component agreement is tested), but recovery of
# mixtures of this family has not been studied.
.fcr_fm_validated <- c("POISSON", "NB1", "NB2")

.fcr_fm_family_status <- function(families) {
  ifelse(families %in% .fcr_fm_validated, "validated", "experimental")
}

# Normalize the family argument and the number of components. Returns the
# family vector (length K) and the display labels.
.fcr_fm_families <- function(family, components, components_missing,
                             label_prefix = "Component") {
  # One family, a vector of names, or a list of names and user families
  if (inherits(family, "fcr_count_family")) family <- list(family)
  if (!(is.character(family) || is.list(family)) || !length(family) ||
      (is.character(family) && anyNA(family))) {
    stop("`family` must be family codes or families made with ",
         "new_count_family().", call. = FALSE)
  }
  labels <- names(family)
  fam <- vapply(family, .fcr_resolve_family, character(1), USE.NAMES = FALSE)
  if (!is.numeric(components) || length(components) != 1L ||
      !is.finite(components) || components < 1 ||
      components != round(components)) {
    stop("The number of components (classes) must be a single whole ",
         "number of at least 1.", call. = FALSE)
  }
  K <- as.integer(components)
  if (length(fam) > 1L) {
    if (components_missing) {
      K <- length(fam)
    } else if (length(fam) != K) {
      stop("`family` has ", length(fam), " entries but ", K,
           " components (classes) were requested; give one family, or one ",
           "per component.", call. = FALSE)
    }
  } else {
    fam <- rep(fam, K)
    labels <- NULL
  }
  if (is.null(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    if (!is.null(labels)) {
      stop("Names of `family` (component labels) must be non-empty and ",
           "unique.", call. = FALSE)
    }
    labels <- paste0(label_prefix, seq_len(K))
  }
  specs <- if (is.list(family)) {
    Filter(function(f) inherits(f, "fcr_count_family"), family)
  }
  list(families = unname(fam), labels = unname(labels), K = K,
       specs = unname(specs[!duplicated(vapply(specs, `[[`, "", "key"))]))
}

# ---- Free-parameter map --------------------------------------------------------
#
# All parameters are stored in one named free vector `theta`:
#   Component<k>:Mean:<column>           component-specific mean coefficient
#   Shared:Mean:<column>                 mean coefficient shared by all
#   Component<k>:<label>[:<column>]      distribution parameter (label as in
#                                        countreg(), e.g. ln(alpha))
#   Shared:<label>[:<column>]            shared distribution parameter block
#   Mixing:Component<k>:<column>         multinomial-logit coefficient of
#                                        component k (k < K; component K is
#                                        the reference with coefficients 0)
# Coefficient names always use the immutable component number, never the
# display label. The map gives, for each component, the positions of its
# mean, dis1 and dis2 coefficients in theta (shared positions repeat across
# components) and the positions of the mixing coefficients.

# Logical vector over the columns of X: which belong to shared terms.
.fcr_fm_shared_cols <- function(shared_terms, X, tt) {
  p <- ncol(X)
  out <- rep(FALSE, p)
  if (!length(shared_terms)) return(out)
  if (!is.character(shared_terms) || anyNA(shared_terms)) {
    stop("`shared$mean` must be a character vector of formula terms.",
         call. = FALSE)
  }
  labels <- attr(tt, "term.labels")
  asg <- attr(X, "assign")
  for (s in unique(shared_terms)) {
    if (s == "(Intercept)") {
      if (!"(Intercept)" %in% colnames(X)) {
        stop("`shared$mean` names \"(Intercept)\" but the model has no ",
             "intercept.", call. = FALSE)
      }
      out[colnames(X) == "(Intercept)"] <- TRUE
    } else {
      j <- match(s, labels)
      if (is.na(j)) {
        stop("`shared$mean` term \"", s, "\" is not a term of the formula. ",
             "Terms: ", toString(c("(Intercept)", labels)), ". Name whole ",
             "terms (for a factor, the factor, not one of its dummy ",
             "columns).", call. = FALSE)
      }
      out[asg == j] <- TRUE
    }
  }
  out
}

.fcr_fm_check_shared <- function(shared) {
  if (is.null(shared)) shared <- list()
  if (!is.list(shared)) {
    stop("`shared` must be a list with elements `mean`, `dis1`, `dis2` ",
         "and `rpar`.", call. = FALSE)
  }
  bad <- setdiff(names(shared), c("mean", "dis1", "dis2", "rpar"))
  if (length(bad) || (length(shared) && is.null(names(shared)))) {
    stop("Unknown element(s) of `shared`: ", toString(bad), ". Use `mean` ",
         "(term names), `dis1`, `dis2` and `rpar` (TRUE or FALSE).",
         call. = FALSE)
  }
  flag <- function(v, nm) {
    if (is.null(v)) return(FALSE)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("`shared$", nm, "` must be TRUE or FALSE.", call. = FALSE)
    }
    v
  }
  list(mean = as.character(shared$mean %||% character(0)),
       dis1 = flag(shared$dis1, "dis1"), dis2 = flag(shared$dis2, "dis2"),
       rpar = flag(shared$rpar, "rpar"))
}

# Build the map. D1, D2: design matrices of the distribution parameters (a
# column of ones when the parameter has no formula); d1_formula, d2_formula:
# whether a formula was given (for the coefficient names). Z: mixing design
# at the unit level (NULL for fixed mixing or K = 1). rp: the random
# parameters (NULL without them): their names and whether they are
# correlated; each component gets their scales (`SD:` or, correlated, the
# packed Cholesky factor `Chol:`), unless shared$rpar.
.fcr_fm_map <- function(families, X, tt, D1, D2, d1_formula, d2_formula, Z,
                        shared, rp = NULL) {
  K <- length(families)
  xn <- colnames(X)
  p <- length(xn)
  params <- lapply(families, get_params)
  has1 <- vapply(params, function(pr) !is.null(pr[[1]]), logical(1))
  has2 <- vapply(params, function(pr) !is.null(pr[[2]]), logical(1))
  sh_cols <- if (K > 1L) .fcr_fm_shared_cols(shared$mean, X, tt) else
    rep(FALSE, p)
  for (s in c("dis1", "dis2")) {
    if (!isTRUE(shared[[s]]) || K == 1L) next
    has <- if (s == "dis1") has1 else has2
    if (!all(has)) {
      stop("`shared$", s, "` = TRUE, but not every component has that ",
           "distribution parameter.", call. = FALSE)
    }
    if (length(unique(families)) > 1L) {
      stop("`shared$", s, "` = TRUE needs the same family in every ",
           "component: distribution parameters of different families do ",
           "not have a common meaning, even when their labels agree.",
           call. = FALSE)
    }
  }
  if (K > 1L && length(unique(families)) == 1L && all(sh_cols) &&
      (!has1[1] || isTRUE(shared$dis1)) && (!has2[1] || isTRUE(shared$dis2)) &&
      (is.null(rp) || isTRUE(shared$rpar))) {
    stop("Every parameter is shared, so all components are identical and ",
         "the mixing probabilities are not identified. Leave some ",
         "coefficients component-specific.", call. = FALSE)
  }

  nm <- character(0)
  info <- list()
  add <- function(name, block, comp, col, label, fam) {
    nm <<- c(nm, name)
    info[[length(info) + 1L]] <<- data.frame(
      name = name, block = block, component = comp, column = col,
      label = label, family = fam, stringsAsFactors = FALSE)
    length(nm)
  }
  dis_name <- function(prefix, label, D, j, has_formula) {
    if (has_formula) paste0(prefix, label, ":", colnames(D)[j]) else
      paste0(prefix, label)
  }
  mean_idx <- lapply(seq_len(K), function(k) integer(p))
  d1_idx <- lapply(seq_len(K), function(k) integer(0))
  d2_idx <- d1_idx
  rv_idx <- d1_idx
  rv_names <- if (!is.null(rp)) {
    if (rp$correlated) sub("^Chol:", "Chol:", .fcr_chol_names(rp$names)) else
      paste0("SD:", rp$names)
  }
  rp_shared <- !is.null(rp) && (isTRUE(shared$rpar) || K == 1L)
  comp_specific <- d1_idx          # component-specific positions, in order
  for (k in seq_len(K)) {
    pre <- paste0("Component", k, ":")
    for (j in which(!sh_cols)) {
      mean_idx[[k]][j] <- add(paste0(pre, "Mean:", xn[j]), "mean", k, xn[j],
                              "Mean", families[k])
    }
    if (has1[k] && !(isTRUE(shared$dis1) && K > 1L)) {
      d1_idx[[k]] <- vapply(seq_len(ncol(D1)), function(j)
        add(dis_name(pre, params[[k]][[1]], D1, j, d1_formula), "dis1", k,
            colnames(D1)[j], params[[k]][[1]], families[k]), integer(1))
    }
    if (has2[k] && !(isTRUE(shared$dis2) && K > 1L)) {
      d2_idx[[k]] <- vapply(seq_len(ncol(D2)), function(j)
        add(dis_name(pre, params[[k]][[2]], D2, j, d2_formula), "dis2", k,
            colnames(D2)[j], params[[k]][[2]], families[k]), integer(1))
    }
    if (!is.null(rp) && !rp_shared) {
      rv_idx[[k]] <- vapply(rv_names, function(v)
        add(paste0(pre, v), "rpar", k, sub("^(SD|Chol):", "", v),
            sub(":.*", "", v), families[k]), integer(1), USE.NAMES = FALSE)
    }
    comp_specific[[k]] <- c(mean_idx[[k]][!sh_cols], d1_idx[[k]],
                            d2_idx[[k]], rv_idx[[k]])
  }
  if (rp_shared) {
    i <- vapply(rv_names, function(v)
      add(if (K > 1L) paste0("Shared:", v) else paste0("Component1:", v),
          "rpar", if (K > 1L) NA_integer_ else 1L, sub("^(SD|Chol):", "", v),
          sub(":.*", "", v), families[1]), integer(1), USE.NAMES = FALSE)
    for (k in seq_len(K)) rv_idx[[k]] <- i
  }
  for (j in which(sh_cols)) {
    i <- add(paste0("Shared:Mean:", xn[j]), "mean", NA_integer_, xn[j],
             "Mean", families[1])
    for (k in seq_len(K)) mean_idx[[k]][j] <- i
  }
  if (K > 1L && isTRUE(shared$dis1)) {
    i <- vapply(seq_len(ncol(D1)), function(j)
      add(dis_name("Shared:", params[[1]][[1]], D1, j, d1_formula), "dis1",
          NA_integer_, colnames(D1)[j], params[[1]][[1]], families[1]),
      integer(1))
    for (k in seq_len(K)) d1_idx[[k]] <- i
  }
  if (K > 1L && isTRUE(shared$dis2)) {
    i <- vapply(seq_len(ncol(D2)), function(j)
      add(dis_name("Shared:", params[[1]][[2]], D2, j, d2_formula), "dis2",
          NA_integer_, colnames(D2)[j], params[[1]][[2]], families[1]),
      integer(1))
    for (k in seq_len(K)) d2_idx[[k]] <- i
  }
  mix_idx <- NULL
  if (!is.null(Z) && K > 1L) {
    mix_idx <- matrix(0L, K - 1L, ncol(Z))
    for (k in seq_len(K - 1L)) {
      for (j in seq_len(ncol(Z))) {
        mix_idx[k, j] <- add(paste0("Mixing:Component", k, ":",
                                    colnames(Z)[j]), "mixing", k,
                             colnames(Z)[j], "Mixing", NA_character_)
      }
    }
  }
  info <- do.call(rbind, info)
  rownames(info) <- NULL
  comp_par <- which(info$block != "mixing")
  list(K = K, families = families, names = nm, n_par = length(nm),
       info = info, mean_idx = mean_idx, d1_idx = d1_idx, d2_idx = d2_idx,
       rv_idx = rv_idx, mix_idx = mix_idx, comp_specific = comp_specific,
       comp_par = comp_par, shared_mean = sh_cols,
       shared = shared, has1 = has1, has2 = has2,
       any_shared = any(sh_cols) || (K > 1L &&
                                       (isTRUE(shared$dis1) ||
                                          isTRUE(shared$dis2) ||
                                          (!is.null(rp) && rp_shared))))
}

# Split theta into component and mixing parameters.
.fcr_fm_unpack <- function(theta, map) {
  K <- map$K
  gamma <- NULL
  if (!is.null(map$mix_idx)) {
    gamma <- cbind(t(matrix(theta[map$mix_idx], nrow = K - 1L)), 0)
  }
  list(beta = lapply(map$mean_idx, function(i) theta[i]),
       d1 = lapply(map$d1_idx, function(i) if (length(i)) theta[i]),
       d2 = lapply(map$d2_idx, function(i) if (length(i)) theta[i]),
       rv = lapply(map$rv_idx %||% vector("list", K), function(i)
         if (length(i)) theta[i]),
       gamma = gamma)
}

# Reorder components: new component k is old component perm[k]. Only valid
# between components with the same family (exchangeable blocks of the same
# structure). Mixing coefficients are re-expressed against the new
# reference component (the last one); this is a linear map, so the same
# function transforms estimates, and the covariance follows from its
# matrix (.fcr_fm_perm_matrix()).
.fcr_fm_permute <- function(theta, map, perm) {
  K <- map$K
  if (identical(as.integer(perm), seq_len(K))) return(theta)
  if (!identical(map$families[perm], map$families)) {
    stop("Components can only be relabeled among components of the same ",
         "family.", call. = FALSE)
  }
  out <- theta
  for (k in seq_len(K)) {
    out[map$comp_specific[[k]]] <- theta[map$comp_specific[[perm[k]]]]
  }
  if (!is.null(map$mix_idx)) {
    G <- .fcr_fm_unpack(theta, map)$gamma[, perm, drop = FALSE]
    G <- G - G[, K]
    out[map$mix_idx] <- t(G[, seq_len(K - 1L), drop = FALSE])
  }
  out
}

# Matrix A with permute(theta) = A %*% theta (for transforming covariances)
.fcr_fm_perm_matrix <- function(map, perm) {
  n <- map$n_par
  A <- vapply(seq_len(n), function(j) {
    e <- numeric(n)
    e[j] <- 1
    .fcr_fm_permute(e, map, perm)
  }, numeric(n))
  dimnames(A) <- list(map$names, map$names)
  A
}

# ---- Relabeling ----------------------------------------------------------------
#
# A bootstrap refit is aligned with the original estimates theta by the
# relabeling p that minimizes the scaled squared distance
#   d(p) = sum_i ((permute(est, p)_i - theta_i) / scale_i)^2
# over the permutations that keep each component's family (and fixed
# mixing probability). Enumerating all K! permutations is only feasible for
# a few components, but d(p) separates once the reference component is
# fixed: shared coefficients add a constant, the component-specific
# coefficients of new component k depend only on p[k], and so do its mixing
# coefficients gamma[p[k]] - gamma[p[K]] once p[K] = r is fixed. The
# optimum is therefore the best of at most K linear assignment problems of
# size K - 1 (a single problem of size K for fixed or one-component mixing),
# each solved exactly by the Hungarian algorithm in O(K^3).

# Hungarian algorithm (shortest augmenting paths with potentials) for a
# square matrix of finite costs W: the assignment perm (row i to column
# perm[i]) minimizing sum_i W[i, perm[i]]. Vectors are indexed from 1 with
# position 1 standing for the algorithm's dummy column 0.
.fcr_fm_hungarian <- function(W) {
  n <- nrow(W)
  u <- numeric(n + 1L)
  v <- numeric(n + 1L)
  p <- integer(n + 1L)
  way <- integer(n + 1L)
  for (i in seq_len(n)) {
    p[1L] <- i
    j0 <- 0L
    minv <- rep(Inf, n + 1L)
    used <- rep(FALSE, n + 1L)
    repeat {
      used[j0 + 1L] <- TRUE
      i0 <- p[j0 + 1L]
      free <- which(!used[-1L])
      cur <- W[i0, free] - u[i0 + 1L] - v[free + 1L]
      better <- cur < minv[free + 1L]
      minv[free[better] + 1L] <- cur[better]
      way[free[better] + 1L] <- j0
      m <- which.min(minv[free + 1L])
      delta <- minv[free[m] + 1L]
      j1 <- free[m]
      uj <- which(used)
      u[p[uj] + 1L] <- u[p[uj] + 1L] + delta
      v[uj] <- v[uj] - delta
      minv[!used] <- minv[!used] - delta
      j0 <- j1
      if (p[j0 + 1L] == 0L) break
    }
    repeat {
      j1 <- way[j0 + 1L]
      p[j0 + 1L] <- p[j1 + 1L]
      j0 <- j1
      if (j0 == 0L) break
    }
  }
  perm <- integer(n)
  perm[p[-1L]] <- seq_len(n)
  perm
}

# Optimal linear assignment for a square cost matrix C whose non-finite
# entries are forbidden pairs. Returns the assignment, its cost (Inf when
# no assignment avoids the forbidden pairs) and, with second = TRUE, the
# cost of the second-best assignment: every other assignment leaves out at
# least one pair of the best, so it is the best of the n problems that each
# forbid one of those pairs.
.fcr_fm_lap <- function(C, second = FALSE) {
  n <- nrow(C)
  if (n == 0L) return(list(perm = integer(0), cost = 0, second = Inf))
  fin <- is.finite(C)
  if (!any(fin)) return(list(perm = seq_len(n), cost = Inf, second = Inf))
  # Forbidden pairs get a cost larger than any assignment of allowed pairs,
  # so they are used only when nothing else is possible
  W <- C
  W[!fin] <- 1 + 2 * n * max(abs(C[fin]))
  perm <- .fcr_fm_hungarian(W)
  cost <- sum(C[cbind(seq_len(n), perm)])
  if (is.na(cost)) cost <- Inf
  out <- list(perm = perm, cost = cost, second = Inf)
  if (second && n > 1L && is.finite(cost)) {
    out$second <- min(vapply(seq_len(n), function(i) {
      C2 <- C
      C2[i, perm[i]] <- Inf
      .fcr_fm_lap(C2)$cost
    }, numeric(1)))
  }
  out
}

# Relabeling of the refit `est` closest to `theta` (see above): the
# permutation (for .fcr_fm_permute()), its distance d and the distance of
# the second-best allowed relabeling (Inf when there is none).
.fcr_fm_align <- function(est, theta, map, scale, fixed = NULL) {
  K <- map$K
  ident <- list(perm = seq_len(K), dist = NA_real_, second = Inf)
  if (K == 1L || !all(is.finite(est))) return(ident)
  zsq <- function(v, i) sum(((v - theta[i]) / scale[i])^2)
  # allowed[k, j]: old component j may become new component k
  allowed <- outer(map$families, map$families, "==")
  if (!is.null(fixed)) {
    # Fixed probabilities are positive; equal up to rounding
    allowed <- allowed & abs(outer(fixed, fixed, "-")) <=
      1.5e-8 * outer(fixed, fixed, pmax)
  }
  cs <- map$comp_specific
  C <- matrix(Inf, K, K)
  for (k in seq_len(K)) for (j in which(allowed[k, ])) {
    C[k, j] <- zsq(est[cs[[j]]], cs[[k]])
  }
  shared <- setdiff(map$comp_par, unlist(cs))
  c0 <- zsq(est[shared], shared)
  if (is.null(map$mix_idx)) {
    a <- .fcr_fm_lap(C, second = TRUE)
    return(list(perm = a$perm, dist = c0 + a$cost, second = c0 + a$second))
  }
  G <- .fcr_fm_unpack(est, map)$gamma
  G0 <- .fcr_fm_unpack(theta, map)$gamma
  S <- matrix(scale[map$mix_idx], nrow = K - 1L)
  best <- NULL
  totals <- numeric(0)
  for (r in which(allowed[K, ])) {
    cols <- setdiff(seq_len(K), r)
    # Mixing cost of new component k (< K) taking old component j, with r
    # the new reference
    M <- vapply(cols, function(j) {
      vapply(seq_len(K - 1L), function(k)
        sum(((G[, j] - G[, r] - G0[, k]) / S[k, ])^2), numeric(1))
    }, numeric(K - 1L))
    M <- matrix(M, nrow = K - 1L)
    a <- .fcr_fm_lap(C[-K, cols, drop = FALSE] + M, second = TRUE)
    base <- c0 + C[K, r]
    totals <- c(totals, base + a$cost, base + a$second)
    if (is.null(best) || base + a$cost < best$dist) {
      best <- list(perm = c(cols[a$perm], r), dist = base + a$cost)
    }
  }
  # The second-best relabeling is either the best one with another
  # reference or the second best with the same reference
  best$second <- if (length(totals) > 1L) sort(totals)[2L] else Inf
  best
}
