## ----include = FALSE----------------------------------------------------------
pkgs <- c("parsnip", "recipes", "workflows", "workflowsets", "rsample",
          "tune", "yardstick", "dials", "hardhat")
have <- all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = have,
                      message = FALSE, warning = FALSE)
options(digits = 4)

## ----setup--------------------------------------------------------------------
library(parsnip)
library(recipes)
library(workflows)
library(workflowsets)
library(rsample)
library(tune)
library(yardstick)
library(flexCountReg)

data("washington_roads")
roads <- tibble::as_tibble(washington_roads)
roads$ID <- as.character(roads$ID)

## ----first--------------------------------------------------------------------
nb_spec <- count_reg(family = "NB2") |>
  set_engine("flexCountReg", method = "BFGS")

nb_fit <- fit(nb_spec, Total_crashes ~ lnaadt + lnlength + speed50,
              data = roads)
tidy(nb_fit)
glance(nb_fit)

## ----first-predict------------------------------------------------------------
new_roads <- roads[1:4, ]
predict(nb_fit, new_roads)
predict(nb_fit, new_roads, type = "conf_int")
predict(nb_fit, new_roads, type = "pred_int", level = 0.9)
augment(nb_fit, new_roads)[, c(".pred", ".resid", "Total_crashes")]

## ----first-raw----------------------------------------------------------------
predict(nb_fit, new_roads, type = "raw",
        opts = list(type = "prob", counts = 0:3, overflow = TRUE))

## ----engines------------------------------------------------------------------
show_engines("count_reg")

## ----workflow-----------------------------------------------------------------
roads$year <- factor(roads$Year)
set.seed(2026)
split <- group_initial_split(roads, group = ID, prop = 0.75)
train <- training(split)
test <- testing(split)

rec <- recipe(Total_crashes ~ lnaadt + lnlength + speed50 + year,
              data = train) |>
  step_normalize(lnaadt, lnlength) |>
  step_dummy(year)

nb_wf <- workflow() |>
  add_recipe(rec) |>
  add_model(nb_spec)
nb_wf_fit <- fit(nb_wf, train)
tidy(nb_wf_fit)

augment(nb_wf_fit, test) |>
  rmse(truth = Total_crashes, estimate = .pred)

## ----offset-formula-----------------------------------------------------------
wf_offset <- workflow() |>
  add_variables(outcomes = Total_crashes,
                predictors = c(lnaadt, lnlength, speed50)) |>
  add_model(count_reg() |>
              set_engine("flexCountReg", method = "BFGS",
                         dis_param_formula_1 = ~ lnlength),
            formula = Total_crashes ~ lnaadt + speed50 + offset(lnlength))
tidy(fit(wf_offset, train))[, 1:4]

## ----offset-recipe, eval = FALSE----------------------------------------------
# workflow() |>
#   add_recipe(rec) |>
#   add_model(nb_spec, formula = Total_crashes ~ s(lnaadt) + speed50)

## ----offset-arg---------------------------------------------------------------
wf_offset2 <- workflow() |>
  add_formula(Total_crashes ~ lnaadt + lnlength) |>
  add_model(count_reg() |> set_engine("flexCountReg", method = "BFGS",
                                      offset = "lnlength"))
tidy(fit(wf_offset2, train))[, 1:4]

## ----workflowset--------------------------------------------------------------
set.seed(2026)
folds <- group_vfold_cv(train, group = ID, v = 3)

specs <- list(
  poisson = count_reg(family = "POISSON") |>
    set_engine("flexCountReg", method = "BFGS"),
  nb2 = nb_spec,
  two_class = count_reg(family = "POISSON") |>
    set_engine("countreg_lc", classes = 2, n_starts = 3, seed = 1,
               group_var = "ID"),
  re_poisson = count_reg() |>
    set_engine("repois", group_var = "ID", method = "BFGS"))

base <- workflow_set(
  preproc = list(base = Total_crashes ~ lnaadt + lnlength + speed50 + ID),
  models = specs[c("two_class", "re_poisson")])
plain <- workflow_set(
  preproc = list(base = Total_crashes ~ lnaadt + lnlength + speed50),
  models = specs[c("poisson", "nb2")])

res <- dplyr::bind_rows(base, plain) |>
  workflow_map("fit_resamples", resamples = folds,
               metrics = metric_set(rmse, mae))
rank_results(res, rank_metric = "rmse")[, c("wflow_id", ".metric", "mean",
                                            "std_err")]

## ----logscore-----------------------------------------------------------------
f <- Total_crashes ~ lnaadt + lnlength + speed50
scores <- lapply(specs, fit_resamples_counts, resamples = folds,
                 formula = f, log_score = TRUE)
do.call(rbind, Map(function(s, nm) {
  cbind(model = nm, s[s$.metric == "log_score", c("mean", "std_err")])
}, scores, names(scores)))

## ----tune-family--------------------------------------------------------------
tune_wf <- workflow() |>
  add_formula(Total_crashes ~ lnaadt + lnlength + speed50) |>
  add_model(count_reg(family = tune()) |>
              set_engine("flexCountReg", method = "BFGS"))

tuned <- tune_grid(tune_wf, resamples = folds,
                   grid = data.frame(family = c("POISSON", "NB1", "NB2",
                                                "PL")),
                   metrics = metric_set(rmse, mae))
show_best(tuned, metric = "rmse")

final_wf <- finalize_workflow(tune_wf, select_best(tuned, metric = "rmse"))
last_fit(final_wf, split, metrics = metric_set(rmse, mae)) |>
  collect_metrics()

## ----tune-classes-------------------------------------------------------------
lc_wf <- workflow() |>
  add_formula(Total_crashes ~ lnaadt + lnlength) |>
  add_model(count_reg(family = "POISSON") |>
              set_engine("countreg_mix", components = tune(), n_starts = 2,
                         seed = 1))
hardhat::extract_parameter_set_dials(lc_wf)
tune_grid(lc_wf, folds, grid = data.frame(components = 2:3),
          metrics = metric_set(rmse)) |>
  collect_metrics()

## ----weights------------------------------------------------------------------
train_w <- train
train_w$w <- hardhat::frequency_weights(rep(1:2, length.out = nrow(train)))
wf_w <- workflow() |>
  add_formula(Total_crashes ~ lnaadt + lnlength) |>
  add_case_weights(w) |>
  add_model(nb_spec)
tidy(fit(wf_w, train_w))[, 1:4]

## ----hier---------------------------------------------------------------------
hier_wf <- workflow() |>
  add_variables(outcomes = Total_crashes,
                predictors = c(lnaadt, lnlength, ID)) |>
  add_model(count_reg(family = "NB2") |>
              set_engine("countreg_hier",
                         predict_args = list(marginal = TRUE)),
            formula = Total_crashes ~ lnaadt + lnlength + (1 | ID))
hier_fit <- fit(hier_wf, train)
predict(hier_fit, test[1:3, ])

## ----rp-----------------------------------------------------------------------
rp_spec <- count_reg(family = "NB2") |>
  set_engine("countreg_rp", rpar_formula = ~ -1 + lnlength,
             ndraws = 100, predict_args = list(method = "Exact"))
rp_wf <- workflow() |>
  add_formula(Total_crashes ~ lnaadt + speed50 + lnlength) |>
  add_model(rp_spec)
rp_fit <- fit(rp_wf, train)
tidy(rp_fit)[, 1:4]
predict(rp_fit, test[1:3, ])

## ----mv-----------------------------------------------------------------------
set.seed(1)
sev <- train
sev$injury <- rbinom(nrow(sev), sev$Total_crashes, 0.35)
sev$pdo <- sev$Total_crashes - sev$injury

mv_wf <- workflow() |>
  add_variables(outcomes = c(injury, pdo),
                predictors = c(lnaadt, lnlength)) |>
  add_model(count_reg(family = "NB2") |>
              set_engine("countreg_mv", copula = "frank"),
            formula = cbind(injury, pdo) ~ lnaadt + lnlength)
mv_fit <- fit(mv_wf, sev)
predict(mv_fit, sev[1:3, ])

## ----mv-int-------------------------------------------------------------------
predict(mv_fit, sev[1:3, ], type = "pred_int", level = 0.9)

## ----mv-validate--------------------------------------------------------------
mv_folds <- group_vfold_cv(sev, group = ID, v = 3)
fit_resamples_counts(mv_wf, mv_folds, log_score = TRUE)

## ----engine, fig.width = 6, fig.height = 4------------------------------------
m <- extract_fit_engine(fit(wf_offset, train))
margEffTable(m)
cureplot(m, indvar = "lnaadt")
sqrt(diag(vcov(m, type = "robust")))[1:3]

