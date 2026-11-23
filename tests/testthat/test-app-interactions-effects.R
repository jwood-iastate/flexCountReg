# App: interaction terms in the model specification, and marginal effects
# of a fit computed as background jobs.

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("callr")
skip_if_not_installed("jsonlite")
skip_if_not_installed("zip")

wait_job <- function(job, timeout = 180) {
  t0 <- Sys.time()
  while (identical((job <- poll_fit_job(job))$status, "running")) {
    if (difftime(Sys.time(), t0, units = "secs") > timeout) {
      stop("job timed out")
    }
    Sys.sleep(0.2)
  }
  job
}

new_state <- function() {
  st <- new_app_state()
  dir.create(shiny::isolate(st$session_dir))
  st
}

# ---- Interactions ---------------------------------------------------------------

test_that("interactions follow the main effects in the formula", {
  s <- app_spec(recipe = list(response = "y",
                              predictors = c("x 1", "z", "b`q"),
                              interactions = list(c("x 1", "z"),
                                                  c("b`q", "x 1"))))
  expect_identical(compile_formula(s),
                   quote(y ~ `x 1` + z + `b\`q` + `x 1`:z + `b\`q`:`x 1`))
  # The names stay single symbols, whatever they contain
  tl <- attr(stats::terms(eval(compile_formula(s))), "term.labels")
  expect_identical(tl, c("`x 1`", "z", "`b\\`q`", "`x 1`:z", "`x 1`:`b\\`q`"))
  # countreg_rp: an interaction with a random parameter stays fixed
  s_rp <- app_spec(recipe = list(response = "y", predictors = c("x", "z"),
                                 interactions = list(c("x", "z"))),
                   model = list(estimator = "countreg_rp", rpar = "x"))
  expect_identical(compile_formula(s_rp), quote(y ~ z + x:z))
  # No interactions: the formula is unchanged
  s0 <- app_spec(recipe = list(response = "y", predictors = c("x", "z")))
  expect_identical(compile_formula(s0), quote(y ~ x + z))
})

test_that("interactions are validated with readable messages", {
  d <- data.frame(y = c(0, 1, 2), a = 1:3, b = 3:1, c = 1, check.names = FALSE)
  v <- function(inter, preds = c("a", "b")) {
    validate_spec(app_spec(recipe = list(response = "y", predictors = preds,
                                         interactions = inter),
                           model = list(family = "POISSON")), d)
  }
  expect_length(v(list(c("a", "b"))), 0)
  expect_match(v(list(c("a", "a"))), "two different predictors", all = FALSE)
  expect_match(v(list(c("a", "c"))), "\"c\", which is not among the predictors",
               all = FALSE)
  expect_match(v(list("a")), "must name two predictors", all = FALSE)
  expect_match(v(list(c("a", "b"), c("b", "a"))), "more than once",
               all = FALSE)
  # Random parameters are predictors too
  s <- app_spec(recipe = list(response = "y", predictors = c("a", "b"),
                              interactions = list(c("a", "b"))),
                model = list(estimator = "countreg_rp", rpar = "a"))
  expect_length(validate_spec(s, d), 0)
})

test_that("the Build page encodes interaction pairs safely", {
  preds <- c("a b", "c:d", "e\"f")
  ch <- interaction_choices(preds)
  expect_length(ch, 3)
  expect_identical(names(ch)[1], "a b × c:d")
  expect_identical(.pair_decode(ch), list(c("a b", "c:d"), c("a b", "e\"f"),
                                          c("c:d", "e\"f")))
  expect_length(interaction_choices("a"), 0)
  # A pair stays chosen while both predictors do, whatever their order
  kept <- keep_interactions(ch[c(1, 3)], c("c:d", "a b"))
  expect_identical(.pair_decode(kept), list(c("c:d", "a b")))
  expect_length(keep_interactions(ch[3], c("a b", "c:d")), 0)
  # The draft specification receives the pairs
  state <- new_app_state()
  shiny::isolate(app_load_data(state, "bundled", "abalone_rings",
                               response = "rings"))
  shiny::testServer(mod_build_server, args = list(state = state), {
    session$setInputs(estimator = "countreg", family = "POISSON",
                      family_each = NULL, components = 2, nm_family = "NM2",
                      response = "rings",
                      predictors = c("shell_weight", "sex", "length"),
                      factors = "sex", exposure = "", log_offset = "",
                      weights = "", group = "", mixing = NULL,
                      shared_mean = NULL, shared_dis1 = FALSE, rpar = NULL,
                      rpar_intercept = FALSE, rpar_dist = "n",
                      correlated = FALSE, ndraws = 500, n_starts = 3,
                      seed = 1, method = "BFGS", max_iters = 500,
                      stderr = "normal", label = "",
                      interactions = .pair_value(c("shell_weight", "sex")))
    expect_identical(state$spec$recipe$interactions,
                     list(c("shell_weight", "sex")))
    expect_length(validate_spec(state$spec, state$data), 0)
    # A pair whose predictor is no longer chosen leaves the draft
    session$setInputs(predictors = c("shell_weight", "length"))
    expect_length(state$spec$recipe$interactions, 0)
  })
})

test_that("an interaction model fits, predicts and replays like countreg()", {
  skip_on_cran()
  state <- new_state()
  d <- load_bundled_data("washington_roads")
  shiny::isolate({
    app_load_data(state, "bundled", "washington_roads")
    app_add_variable(state, list(type = "log", x = "AADT",
                                 name = "log aadt"))
    s <- app_spec(data = list(source = "bundled", name = "washington_roads",
                              object = "data_1"),
                  recipe = list(response = "Total_crashes",
                                predictors = c("log aadt", "speed50"),
                                factors = "speed50", exposure = "Length",
                                interactions = list(c("log aadt",
                                                      "speed50"))),
                  model = list(family = "NB2"))
    expect_length(validate_spec(s, state$data), 0)
    job <- app_submit_fit(state, s)
    record_finished_job(state, wait_job(job))
    fit <- state$env$fit_1
    sc <- typical_row(s, state$data)
    sc <- rbind(sc, sc)
    sc$speed50[2] <- "1"
    p <- app_run_output(state, list(type = "predict", scenarios = sc,
                                    interval = "none"), job$id)
    pj <- wait_job(app_start_output_job(state, list(
      type = "predict", scenarios = sc, interval = "prediction",
      level = 0.9), job$id))
    record_finished_job(state, pj)
    dir <- tempfile("bundle-")
    dir.create(dir)
    write_bundle(state, dir)
  })
  d$`log aadt` <- log(d$AADT)
  d$speed50 <- factor(d$speed50)
  ref <- countreg(Total_crashes ~ `log aadt` + speed50 + `log aadt`:speed50 +
                    offset(log(Length)), data = d, family = "NB2",
                  method = "BFGS", max.iters = 1000)
  expect_identical(names(coef(fit)), names(coef(ref)))
  expect_equal(coef(fit), coef(ref), tolerance = 1e-8)
  # One interaction coefficient for the non-reference level of speed50
  expect_true("`log aadt`:speed501" %in% names(coef(fit)) ||
                "log aadt:speed501" %in% names(coef(fit)))
  expect_equal(p$value$expected_count,
               unname(predict(ref, newdata = sc)), tolerance = 1e-8)
  expect_identical(shiny::isolate(state$jobs[[pj$id]]$status), "completed")
  # The exported script and settings rebuild the same model
  script <- readLines(file.path(dir, "analysis.R"))
  expect_match(script, "`log aadt`:speed50", fixed = TRUE, all = FALSE)
  out <- callr::r(function(dir) {
    setwd(dir)
    env <- new.env()
    sys.source("analysis.R", envir = env)
    list(coef = stats::coef(env$fit_1), pred = env$pred_1$expected_count)
  }, args = list(dir), libpath = .libPaths())
  expect_equal(out$coef, coef(fit), tolerance = 1e-10)
  expect_equal(out$pred, p$value$expected_count, tolerance = 1e-10)
  sj <- jsonlite::fromJSON(file.path(dir, "settings.json"),
                           simplifyVector = TRUE)
  expect_identical(compile_fit(.spec_from_list(sj$fit_1), "fit_1"),
                   compile_fit(s, "fit_1"))
})

test_that("settings with interactions round-trip; older settings load", {
  s1 <- app_spec(recipe = list(response = "y", predictors = c("a b", "c"),
                               interactions = list(c("a b", "c"))))
  s2 <- app_spec(recipe = list(response = "y",
                               predictors = c("a", "b", "c"),
                               interactions = list(c("a", "b"),
                                                   c("b", "c"))))
  s0 <- app_spec(recipe = list(response = "y", predictors = "a"))
  for (s in list(s0, s1, s2)) {
    r <- spec_from_json(spec_to_json(s))
    expect_identical(r, s)
    expect_identical(compile_formula(r), compile_formula(s))
  }
  # Schema 2 settings (before interactions) are migrated
  old <- jsonlite::fromJSON(spec_to_json(s0), simplifyVector = TRUE)
  old$schema <- 2L
  old$recipe$interactions <- NULL
  r_old <- .spec_from_list(old)
  expect_identical(r_old$schema, .spec_schema)
  expect_identical(r_old$recipe$interactions, list())
  expect_identical(compile_fit(r_old, "fit_1"), compile_fit(s0, "fit_1"))
  # Current settings without the field take the default
  cur <- jsonlite::fromJSON(spec_to_json(s0), simplifyVector = TRUE)
  cur$recipe$interactions <- NULL
  expect_identical(.spec_from_list(cur)$recipe$interactions, list())
  # Other schemas are rejected
  cur$schema <- 1L
  expect_error(.spec_from_list(cur), "schema")
})

# ---- Marginal effects ------------------------------------------------------------

test_that("marginal-effect requests compile to margEffTable() calls", {
  grid <- expand.grid(measure = c("auto", "ame", "elasticity"),
                      indicator = c("pseudo_elasticity", "discrete_change"),
                      se = c("delta", "bootstrap"),
                      stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    req <- margeff_request(g$measure, g$indicator, g$se, bootstraps = 50)
    if (g$se == "bootstrap") req$seed <- 42L
    code <- vapply(compile_output(req, "fit_1", "margeff_1"), .deparse, "")
    call_txt <- paste0(
      "margeff_1 <- flexCountReg::margEffTable(fit_1, measure = \"",
      g$measure, "\", indicator = \"", g$indicator, "\", se = \"", g$se,
      "\"", if (g$se == "bootstrap") ", bootstraps = 50L", ", digits = Inf)")
    expect_identical(gsub("\\s+", " ", code[length(code)]), call_txt,
                     label = paste(g, collapse = "/"))
    if (g$se == "bootstrap") {
      expect_identical(code[1], "set.seed(42L)")
    } else {
      expect_length(code, 1)
    }
  }
  # Restricting the variables
  code <- .deparse(compile_output(margeff_request("ame", "discrete_change",
                                                  "delta",
                                                  vars = c("a b", "c")),
                                  "fit_2", "margeff_3")[[1]])
  expect_match(code, "vars = c(\"a b\", \"c\")", fixed = TRUE)
  expect_error(margeff_request("auto", "pseudo_elasticity", "bootstrap",
                               bootstraps = NA), "at least 1")
  # A bootstrap request gets a frozen seed; the delta method does not
  state <- new_app_state()
  shiny::isolate({
    o1 <- app_prepare_output(state, margeff_request(
      "auto", "pseudo_elasticity", "bootstrap", 20))
    o2 <- app_prepare_output(state, margeff_request(
      "auto", "pseudo_elasticity", "delta"))
  })
  expect_true(is.integer(o1$seed))
  expect_identical(o1$exprs[[1]], call("set.seed", o1$seed))
  expect_null(o2$seed)
})

test_that("marginal effects run as jobs, are shown and are exported", {
  skip_on_cran()
  state <- new_state()
  set.seed(4)
  n <- 250
  d <- data.frame(x = rnorm(n), g = rbinom(n, 1, 0.4))
  d$y <- rpois(n, exp(0.2 + 0.5 * d$x + 0.4 * d$g))
  path <- file.path(shiny::isolate(state$session_dir), "upload-1.csv")
  utils::write.csv(d, path, row.names = FALSE)
  shiny::isolate({
    app_load_data(state, "upload", "me.csv", upload_path = path)
    base <- list(source = "upload", name = "me.csv", object = "data_1")
    s <- app_spec(data = base,
                  recipe = list(response = "y", predictors = c("x", "g"),
                                interactions = list(c("x", "g"))),
                  model = list(family = "POISSON"))
    fj <- app_submit_fit(state, s)
    s_mix <- app_spec(data = base,
                      recipe = list(response = "y", predictors = "x"),
                      model = list(estimator = "countreg_mix",
                                   family = "POISSON", components = 2L,
                                   n_starts = 2L, seed = 3L))
    mj <- app_submit_fit(state, s_mix)
    record_finished_job(state, wait_job(fj))
    record_finished_job(state, wait_job(mj))
    fit <- state$env$fit_1
    state$active <- fj$id
  })
  expect_setequal(margeff_var_choices(s), c("x", "g"))
  shiny::testServer(mod_fit_server, args = list(state = state), {
    session$setInputs(active = fj$id, me_measure = "auto",
                      me_indicator = "pseudo_elasticity", me_se = "delta",
                      me_bootstraps = 200, me_vars = NULL)
    session$setInputs(me_go = 1)
    id <- names(state$jobs)[length(state$jobs)]
    expect_match(output$me_status$html, "Running", fixed = TRUE)
    record_finished_job(state, wait_job(state$jobs[[id]]))
    session$flushReact()
    res <- state$outputs[[id]]
    # The same table as calling margEffTable() on the fit. The call is made
    # in a fresh R process, like the job's worker, so both use the same
    # installed flexCountReg (during development it can differ from the
    # source tree loaded here).
    direct <- callr::r(function(fit) {
      flexCountReg::margEffTable(fit, measure = "auto",
                                 indicator = "pseudo_elasticity",
                                 se = "delta", digits = Inf)
    }, args = list(fit), libpath = .libPaths())
    expect_equal(as.data.frame(res$value), as.data.frame(direct),
                 tolerance = 1e-8)
    expect_identical(res$value$term, c("x", "g"))
    expect_match(output$me_title$html, "margeff_1", fixed = TRUE)
    expect_match(output$me_notes$html, "delta method", fixed = TRUE)
    expect_match(output$me_code, "margeff_1 <- flexCountReg::margEffTable",
                 fixed = TRUE)
    expect_null(output$me_status)
    # Unrounded in the download
    csv <- utils::read.csv(output$me_download)
    expect_equal(csv$estimate, as.data.frame(direct)$estimate,
                 tolerance = 1e-12)

    # An unsupported combination: the error is shown in the panel
    session$setInputs(active = mj$id)
    state$active <- mj$id
    session$setInputs(me_se = "bootstrap", me_bootstraps = 20, me_go = 2)
    id2 <- names(state$jobs)[length(state$jobs)]
    j2 <- wait_job(state$jobs[[id2]])
    expect_identical(j2$status, "failed")
    record_finished_job(state, j2)
    session$flushReact()
    expect_match(output$me_status$html, "alert-danger", fixed = TRUE)
    expect_match(output$me_status$html, "mixture models", fixed = TRUE)
    expect_null(state$outputs[[id2]])
    # A bad setting is reported in the panel before anything runs
    session$setInputs(me_bootstraps = NA, me_go = 3)
    expect_match(output$me_status$html, "at least 1", fixed = TRUE)
  })
  # In the session script, the manifest and a replay
  dir <- tempfile("bundle-")
  dir.create(dir)
  shiny::isolate(write_bundle(state, dir, include_data = TRUE))
  script <- paste(readLines(file.path(dir, "analysis.R")), collapse = "\n")
  expect_match(script, "margeff_1 <- flexCountReg::margEffTable(fit_1",
               fixed = TRUE)
  manifest <- readLines(file.path(dir, "manifest.txt"))
  expect_match(manifest, "^  margeff_1: Marginal effects 1", all = FALSE)
  out <- callr::r(function(dir) {
    setwd(dir)
    env <- new.env()
    sys.source("analysis.R", envir = env)
    as.data.frame(env$margeff_1)
  }, args = list(dir), libpath = .libPaths())
  expect_equal(out$estimate, shiny::isolate(
    as.data.frame(state$env$margeff_1)$estimate), tolerance = 1e-10)
})

test_that("bootstrap marginal effects replay from their frozen seed", {
  skip_on_cran()
  state <- new_state()
  shiny::isolate({
    app_load_data(state, "bundled", "washington_roads")
    s <- app_spec(data = list(source = "bundled", name = "washington_roads",
                              object = "data_1"),
                  recipe = list(response = "Total_crashes",
                                predictors = c("lnaadt", "speed50")),
                  model = list(family = "POISSON"))
    fj <- app_submit_fit(state, s)
    record_finished_job(state, wait_job(fj))
    req <- margeff_request("ame", "discrete_change", "bootstrap", 20)
    j <- wait_job(app_start_output_job(state, req, fj$id))
    expect_identical(j$status, "completed")
    record_finished_job(state, j)
    res <- state$outputs[[j$id]]
    fit <- state$env$fit_1
  })
  expect_true(is.integer(res$seed))
  expect_true(paste0("Seed: ", res$seed) %in% res$details)
  expect_true("Standard errors: bootstrap (20 resamples)" %in% res$details)
  out <- callr::r(function(fit, exprs) {
    set.seed(999)
    env <- new.env()
    env$fit_1 <- fit
    for (e in exprs) eval(e, env)
    as.data.frame(env$margeff_1)
  }, args = list(fit, res$exprs), libpath = .libPaths())
  expect_equal(out, as.data.frame(res$value), tolerance = 1e-12)
})
