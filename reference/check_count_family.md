# Check a user-defined count family

Runs numerical checks of a family made with
[`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md)
over a grid of means and parameter values: that the probabilities sum to
one, that the mean is \\\mu\\ (when `mean_is_mu = TRUE`) or matches the
`mean` function, that the `variance`, `cdf` and `rand` functions (when
given) agree with the PMF, that inputs are recycled, and that missing
and invalid values give `-Inf` rather than errors.

## Usage

``` r
check_count_family(
  family,
  mu = c(0.05, 1, 20),
  par = list(),
  tol = 1e-06,
  nsim = 2000
)
```

## Arguments

- family:

  a family from
  [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md).

- mu:

  means to check.

- par:

  a named list with values of each parameter (on the parameter's own
  scale); every combination with `mu` is checked.

- tol:

  tolerance for the comparisons.

- nsim:

  draws per grid point to check `rand` (0 to skip).

## Value

A data frame with one row per check (`check`, `ok`, `detail`), printed,
and returned invisibly. A warning lists the checks that failed.

## Examples

``` r
nb_theta <- new_count_family(
  name = "NB-theta",
  params = list(theta = count_param("log")),
  logpmf = function(y, mu, par)
    dnbinom(y, size = par$theta, mu = mu, log = TRUE),
  variance = function(mu, par) mu + mu^2 / par$theta)
check_count_family(nb_theta, mu = c(0.1, 2, 30),
                   par = list(theta = c(0.5, 5)))
#>                        check   ok                                    detail
#>                PMF sums to 1 TRUE               largest |sum - 1| = 1.1e-16
#>               mean equals mu TRUE          largest relative error = 5.2e-11
#>   variance() matches the PMF TRUE          largest relative error = 1.3e-09
#>            vectorized over y TRUE                                        ok
#>              NA mean handled TRUE                                  no error
#>  invalid mean gives no error TRUE no error (the row is given probability 0)
```
