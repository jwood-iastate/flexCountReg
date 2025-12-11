# Calculate Akaike Information Criterion (AIC)

This function calculates the Akaike Information Criterion (AIC) for a
given model.

## Usage

``` r
myAIC(LL, nparam)
```

## Arguments

- LL:

  Numeric value representing the log-likelihood of the model.

- nparam:

  Numeric value representing the number of parameters in the model.

## Value

Numeric value representing the AIC.

## Details

The AIC is calculated using the formula: \$\$AIC = -2 \cdot LL + 2 \cdot
nparam\$\$ Where \\LL\\ is the log-likelihood of the model and
\\nparam\\ is the number of parameters.

## Examples

``` r
LL <- -120.5
nparam <- 5
myAIC(LL, nparam)
#> [1] 251
```
