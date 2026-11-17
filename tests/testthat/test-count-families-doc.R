# ?count_families must match the family registry: the family codes, the
# parameters modeled by dis_param_formula_1/2, the availability of
# underreporting (binomial thinning) and the likelihood-ratio test against
# Poisson.

test_that("reported-count variances follow the binomial-thinning law", {
  # ?count_families: Var(Y*) = p^2 Var(Y) + p (1 - p) E[Y], with the reported
  # counts following the family at mean parameter p * mu and the thinned
  # distribution parameters
  pars <- list(POISSON = c(NA, NA), NB1 = c(0.5, NA), NB2 = c(0.5, NA),
               NBP = c(0.5, 1.5), PLN = c(0.5, NA), PGE = c(2, NA),
               PIG1 = c(0.5, NA), PIG2 = c(0.5, NA), PIG = c(0.3, NA),
               PL = c(1.5, NA), PLG = c(1.5, 0.4), PLL = c(1.5, 0.5),
               PW = c(2, NA), SI = c(exp(0.5), 1), DELAPORTE = c(1.5, 2),
               NBGE = c(1.5, 0.4), HERMITE = c(0.4, NA), PA = c(0.6, NA))
  expect_setequal(names(pars), flexCountReg:::.fcr_thinning_families)
  mu <- 3
  for (fam in names(pars)) {
    a <- if (is.na(pars[[fam]][1])) NULL else pars[[fam]][1]
    s <- if (is.na(pars[[fam]][2])) NULL else pars[[fam]][2]
    v_true <- flexCountReg:::.fcr_family_variance(fam, mu, a, s)
    m_true <- flexCountReg:::.fcr_family_mean(fam, mu, a, s)
    for (p in c(0.2, 0.7)) {
      th <- flexCountReg:::.fcr_thin_params(fam, a, s, p)
      v_rep <- flexCountReg:::.fcr_family_variance(fam, p * mu, th$alpha,
                                                   th$sigma)
      expect_equal(v_rep, p^2 * v_true + p * (1 - p) * m_true,
                   tolerance = 1e-10, label = paste(fam, p))
      # dispersion index of the reported counts: 1 + p (D - 1)
      m_rep <- flexCountReg:::.fcr_family_mean(fam, p * mu, th$alpha,
                                               th$sigma)
      expect_equal(v_rep / m_rep, 1 + p * (v_true / m_true - 1),
                   tolerance = 1e-10, label = paste(fam, p, "index"))
    }
  }
})

count_families_rows <- function() {
  db <- tryCatch(tools::Rd_db("flexCountReg"), error = function(e) NULL)
  rd <- db[["count_families.Rd"]]
  if (is.null(rd)) return(NULL)
  old <- options(useFancyQuotes = FALSE)
  on.exit(options(old))
  txt <- utils::capture.output(tools::Rd2txt(
    rd, options = list(width = 400, underline_titles = FALSE)))
  # The first two tables (one row per family) come before the
  # underreporting section, whose table groups families
  end <- grep("Underreporting and the variance", txt)[1]
  if (!is.na(end)) txt <- txt[seq_len(end - 1)]
  # Table rows start with a quoted family code followed by 2+ spaces
  rows <- grep("^[[:space:]]+'[A-Z0-9]+'[[:space:]]{2,}", txt, value = TRUE)
  # Rd renders "--" as "-" in released R and as an en dash in R-devel
  rows <- gsub("\u2013|\u2014", "-", rows)
  cells <- strsplit(trimws(rows), "[[:space:]]{2,}")
  fam <- gsub("'", "", vapply(cells, `[`, "", 1))
  list(fam = fam, cells = cells)
}

test_that("count_families documents every registered family once per table", {
  r <- count_families_rows()
  skip_if(is.null(r), "help database of the installed package not available")
  fams <- flexCountReg:::.fcr_families
  tab <- table(r$fam)
  expect_setequal(names(tab), fams)
  expect_true(all(tab == 2))   # one row in each of the two tables
})

test_that("count_families matches the registry", {
  r <- count_families_rows()
  skip_if(is.null(r), "help database of the installed package not available")
  first <- !duplicated(r$fam)
  unq <- function(x) gsub("'", "", sub(" \\(.*$", "", x))
  # Table 1: family, distribution, formula 1 parameter, formula 2 parameter
  for (i in which(first)) {
    fam <- r$fam[i]
    cl <- r$cells[[i]]
    pars <- flexCountReg:::get_params(fam)
    lab <- function(p) if (is.null(p)) "-" else p
    expect_identical(unq(cl[3]), lab(pars[[1]]), label = paste(fam, "par 1"))
    expect_identical(unq(cl[4]), lab(pars[[2]]), label = paste(fam, "par 2"))
  }
  # Table 2: family, distribution functions, underreporting, LR test
  nest <- flexCountReg:::.fcr_family_nesting
  for (i in which(!first)) {
    fam <- r$fam[i]
    cl <- r$cells[[i]]
    n <- length(cl)
    thin <- fam %in% flexCountReg:::.fcr_thinning_families
    expect_identical(startsWith(cl[n - 1], "yes"), thin,
                     label = paste(fam, "underreporting"))
    k <- which(nest$restricted == "POISSON" & nest$full == fam)
    expect_identical(cl[n], if (length(k)) nest$type[k] else "-",
                     label = paste(fam, "LR test"))
  }
})
