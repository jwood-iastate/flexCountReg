# Sichel Distribution

Density, distribution function, quantile function, and random generation
for the Sichel distribution.

## Usage

``` r
dsichel(x, mu = 1, sigma = 1, gamma = 1, log = FALSE)

psichel(q, mu = 1, sigma = 1, gamma = 1, lower.tail = TRUE, log.p = FALSE)

qsichel(p, mu = 1, sigma = 1, gamma = 1, lower.tail = TRUE, log.p = FALSE)

rsichel(n, mu = 1, sigma = 1, gamma = 1)
```

## Arguments

- x:

  numeric value or vector of non-negative integer values.

- mu:

  numeric; mean of the distribution (mu \> 0).

- sigma:

  numeric; scale parameter (sigma \> 0).

- gamma:

  numeric; shape parameter (can be any real number).

- log, log.p:

  logical; if TRUE, probabilities are given as log(p).

- q:

  quantile or vector of quantiles.

- lower.tail:

  logical; if TRUE, probabilities are P\[X \<= x\].

- p:

  probability or vector of probabilities.

- n:

  number of random values to generate.

## Details

The Sichel distribution is a three-parameter discrete distribution that
generalizes the Poisson-inverse Gaussian distribution. It is useful for
modeling overdispersed count data.

The PMF is: \$\$f(y\|\mu, \sigma, \gamma) = \frac{(\mu/c)^y
K\_{y+\gamma}(\alpha)}{K\_\gamma(1/\sigma) y!
(\alpha\sigma)^{y+\gamma}}\$\$

## References

Rigby, R. A., Stasinopoulos, D. M., & Akantziliotou, C. (2008). A
framework for modelling overdispersed count data, including the
Poisson-shifted generalized inverse Gaussian distribution. Computational
Statistics & Data Analysis, 53(2), 381-393.

## Examples

``` r
# Basic usage
dsichel(0:10, mu = 5, sigma = 1, gamma = -0.5)
#>  [1] 0.09860584 0.14865390 0.14583707 0.12259787 0.09727854 0.07583610
#>  [7] 0.05907602 0.04630084 0.03659705 0.02918655 0.02347741

# Log-probabilities for numerical stability
dsichel(0:10, mu = 5, sigma = 1, gamma = -0.5, log = TRUE)
#>  [1] -2.316625 -1.906135 -1.925265 -2.098846 -2.330177 -2.579181 -2.828930
#>  [8] -3.072595 -3.307788 -3.534047 -3.751716

# CDF
psichel(5, mu = 5, sigma = 1, gamma = -0.5)
#> [1] 0.6888093
```
