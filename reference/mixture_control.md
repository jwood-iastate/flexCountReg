# Control settings for finite-mixture and latent-class count models

Settings of the estimation algorithm used by
[`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
and
[`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md).
`lc_control()` is the same function under the name used with the
latent-class interface; both return the same object.

## Usage

``` r
mixture_control(
  algorithm = c("em", "direct"),
  em_maxit = 200L,
  em_tol = 1e-06,
  mstep_maxit = 25L,
  n_polish = 3L,
  boot_starts = 5L,
  small_class = 0.05,
  start_spread = 0.5
)

lc_control(...)
```

## Arguments

- algorithm:

  `"em"` (default): each start is improved by the
  expectation-maximization (EM) algorithm, and the best solutions are
  then polished by maximizing the observed-data log-likelihood directly
  (with the optimizer chosen by `method`). `"direct"`: every start is
  maximized directly, without EM iterations.

- em_maxit:

  maximum number of EM iterations per start.

- em_tol:

  relative convergence tolerance of the EM iterations. EM stops when the
  relative change in the observed log-likelihood is below `em_tol` *and*
  the largest change in a parameter is below `sqrt(em_tol)`.

- mstep_maxit:

  maximum number of quasi-Newton iterations in each numerical M-step (a
  generalized EM step: it need not reach the conditional maximum, but it
  must increase the likelihood).

- n_polish:

  number of distinct EM solutions (the best ones) that are polished by
  direct maximization.

- boot_starts:

  number of starts used in each refit of a bootstrap
  (`stderr = "boot"`), in addition to the original estimates.

- small_class:

  posterior class share below which a component is flagged as small in
  the output. This is a reporting threshold only: it never changes the
  fit.

- start_spread:

  spread of the response-based starting partitions (the standard
  deviation, in units of the standardized score, of the noise added to
  the score that orders the units into starting classes).

- ...:

  arguments passed from `lc_control()` to `mixture_control()`.

## Value

An object of class `fcr_mixture_control`.

## See also

[`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md),
[`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)

## Examples

``` r
mixture_control(em_maxit = 200)
#> Finite-mixture control settings
#>   algorithm : em 
#>   em_maxit : 200 
#>   em_tol : 1e-06 
#>   mstep_maxit : 25 
#>   n_polish : 3 
#>   boot_starts : 5 
#>   small_class : 0.05 
#>   start_spread : 0.5 
lc_control(algorithm = "direct")
#> Finite-mixture control settings
#>   algorithm : direct 
#>   em_maxit : 200 
#>   em_tol : 1e-06 
#>   mstep_maxit : 25 
#>   n_polish : 3 
#>   boot_starts : 5 
#>   small_class : 0.05 
#>   start_spread : 0.5 
```
