# Likelihood-ratio test of a random-effect term

Refits a
[`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
model without one random-effect term and compares the log-likelihoods.
The null value of a variance is on the boundary of its parameter space,
so for a term with one effect the p-value is that of an equal mixture of
\\\chi^2_0\\ and \\\chi^2_1\\ (Self and Liang 1987): half the naive
p-value. For a term with \\q \> 1\\ effects the mixture of
\\\chi^2\_{k-q}\\ to \\\chi^2_k\\ (\\k = q(q + 1)/2\\) is approximated
by an equal mixture of \\\chi^2\_{k-1}\\ and \\\chi^2_k\\.

## Usage

``` r
re_lrtest(object, term = NULL)
```

## Arguments

- object:

  a
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
  fit.

- term:

  the term to test, as its label (e.g. `"1 | site"`); by default the
  last one.

## Value

A data frame with the log-likelihoods, the statistic, its degrees of
freedom and the p-value.

## References

Self, S.G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82(398),
605-610.
