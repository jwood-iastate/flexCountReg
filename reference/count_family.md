# Count distribution tuning parameter

A dials parameter for the `family` argument of
[`count_reg`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md),
for tuning over count distributions.

## Usage

``` r
count_family(
  values = c("POISSON", "NB1", "NB2", "NBP", "GP1", "GP2", "PIG1", "PIG2", "PL", "PLN",
    "COM")
)
```

## Arguments

- values:

  the candidate families.

## Value

A dials qualitative parameter.

## Examples

``` r
if (requireNamespace("dials", quietly = TRUE)) count_family()
#> Count distribution (qualitative)
#> 11 possible values include:
#> 'POISSON', 'NB1', 'NB2', 'NBP', 'GP1', 'GP2', 'PIG1', 'PIG2', 'PL', 'PLN', and
#> 'COM'
```
