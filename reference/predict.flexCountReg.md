# Predictions for flexCountReg models

Generates predictions for the expected count (lambda) for observations.

For **countreg.rp** (Random Parameters) models, three methods are
available:

- **Simulated**: Uses Halton draws to simulate the random parameters and
  averages the outcomes. This is a simulation-based approximation.

- **Individual**: Estimates observation-specific coefficients
  (conditional on observed outcomes) using Empirical Bayes. Requires the
  outcome variable to be present in `data`.

- **Exact**: Uses the analytical Moment Generating Functions (MGFs) of
  the random parameter distributions to calculate the exact expected
  value. This method is faster and removes simulation error.

For **countreg**, **poisLindRE**, and **RENB** models, the function
calculates the expected value \\\mu = \exp(X\beta)\\ (with appropriate
adjustments for specific families like PLN or underreporting).

## Usage

``` r
# S3 method for class 'flexCountReg'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  a model object estimated using this R package.

- newdata:

  optional dataframe for which to generate predictions.

- ...:

  optional arguments passed to the function. This includes \`method\`.

## Note

optional parameter \`newdata\`: a dataframe that has all of the
variables in the `formula` and `rpar_formula`.

optional parameter \`method\`: Only valid for random parameters models
(\`countreg.rp\`). Options include `Simulated` (default), `Individual`,
or `Exact`.

## References

Wood, J.S., Gayah, V. (2025). Out-of-sample prediction and
interpretation for random parameter generalized linear models. *Accident
Analysis and Prevention*, 220, 108147.

## Examples

``` r
# \donttest{
# Load data and create a dummy variable
data("washington_roads")
washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)

# =========================================================================
# 1. Fixed Parameter Model (countreg)
# =========================================================================
nb2_fixed <- countreg(Total_crashes ~ lnaadt + lnlength + speed50,
                      data = washington_roads, 
                      family = "NB2")
pred_fixed <- predict(nb2_fixed, data = washington_roads)

# =========================================================================
# 2. Random Parameters Model (countreg.rp)
# =========================================================================
rp_nb2 <- countreg.rp(Total_crashes ~ lnaadt + lnlength,
                      rpar_formula = ~ -1 + speed50,
                      data = washington_roads,
                      family = "NB2",
                      rpardists = c(speed50 = "n"),
                      ndraws = 100)

# Method A: Simulated (Default)
pred_sim <- predict(rp_nb2, data = washington_roads, method = "Simulated")

# Method B: Exact (Analytical MGF)
pred_exact <- predict(rp_nb2, data = washington_roads, method = "Exact")

# =========================================================================
# 3. Random Effects Models (poisLindRE / RENB)
# =========================================================================
pl_re <- poisLind.re(Total_crashes ~ lnaadt + lnlength,
                     data = washington_roads,
                     group_var = "ID")
#> Warning: NaNs produced
pred_pl_re <- predict(pl_re, data = washington_roads)
# }
```
