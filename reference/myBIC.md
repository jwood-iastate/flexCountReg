# Calculate Bayesian Information Criterion (BIC)

This function calculates the Bayesian Information Criterion (BIC) for a
given model.

## Usage

``` r
myBIC(LL, nparam, n)
```

## Arguments

- LL:

  Numeric value representing the log-likelihood of the model.

- nparam:

  Numeric value representing the number of parameters in the model.

- n:

  Numeric value representing the number of observations.

## Value

Numeric value representing the BIC.

## Details

The BIC is calculated using the formula: \$\$BIC = -2 \cdot LL + nparam
\cdot \log(n)\$\$ Where \\LL\\ is the log-likelihood of the model,
\\nparam\\ is the number of parameters, and \\n\\ is the number of
observations.

## Examples

``` r
LL <- -120.5
nparam <- 5
n <- 100
myBIC(LL, nparam, n)
#> [1] 264.0259
```
