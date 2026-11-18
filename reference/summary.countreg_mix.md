# Summaries of finite-mixture and latent-class models

Summaries of finite-mixture and latent-class models

## Usage

``` r
# S3 method for class 'countreg_mix'
summary(object, level = 0.95, ...)

# S3 method for class 'summary.countreg_mix'
print(x, digits = x$digits, ...)
```

## Arguments

- object:

  a model from
  [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  or
  [`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md).

- level:

  confidence level of the coefficient intervals.

- ...:

  further arguments for
  [`summary.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/summary.flexCountReg.md).

- x:

  a `summary.countreg_mix` object.

- digits:

  number of digits to print.

## Value

An object of class `c("summary.countreg_mix", "summary.flexCountReg")`.

## Details

Besides the coefficient table (grouped by component, shared parameters
and mixing model), the summary reports the prior mixing proportions, the
posterior (effective) class sizes, the entropy of the posterior
classification, the starts that reached the selected solution, and
numerical checks of the observed information. No likelihood-ratio test
of the number of components is reported: its usual chi-square reference
distribution does not apply.

Entropy is \\H = -\sum_u w_u \sum_k \tau\_{uk}\log\tau\_{uk}\\ over the
independent units with posterior probabilities \\\tau\\; the relative
entropy \\1 - H / (N \log K)\\ is 1 for a perfect classification and 0
when every unit has probabilities \\1/K\\. ICL is \\BIC + 2H\\ (smaller
is better), a classification-oriented criterion.
