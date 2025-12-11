# Poisson-Lindley Distribution

These functions provide density, distribution function, quantile
function, and random number generation for the Poisson-Lindley (PL)
Distribution

## Usage

``` r
msg1

dplind(x, mean = 1, theta = 1, lambda = NULL, log = FALSE)

pplind(q, mean = 1, theta = 1, lambda = NULL, lower.tail = TRUE, log.p = FALSE)

qplind(p, mean = 1, theta = 1, lambda = NULL)

rplind(n, mean = 1, theta = 1, lambda = NULL)
```

## Format

An object of class `character` of length 1.

## Arguments

- x:

  numeric value or a vector of values.

- mean:

  numeric value or vector of mean values for the distribution (the
  values have to be greater than 0).

- theta:

  single value or vector of values for the theta parameter of the
  distribution (the values have to be greater than 0).

- lambda:

  alternative parameterization (use instead of the mean); numeric value
  or vector of values for lambda parameter of the distribution (the
  values have to be greater than 0).

- log:

  logical; if TRUE, probabilities p are given as log(p).

- q:

  quantile or a vector of quantiles.

- lower.tail:

  logical; if TRUE, probabilities p are \\P\[X\leq x\]\\ otherwise,
  \\P\[X\>x\]\\.

- log.p:

  logical; if TRUE, probabilities p are given as log(p).

- p:

  probability or a vector of probabilities.

- n:

  the number of random numbers to generate.

## Details

The Poisson-Lindley is a 2-parameter count distribution that captures
high densities for small integer values. This makes it ideal for data
that are zero-inflated.

`dplind` computes the density (PDF) of the Poisson-Lindley Distribution.

`pplind` computes the CDF of the Poisson-Lindley Distribution.

`qplind` computes the quantile function of the Poisson-Lindley
Distribution.

`rplind` generates random numbers from the Poisson-Lindley Distribution.

The compound Probability Mass Function (PMF) for the Poisson-Lindley
(PL) distribution is: \$\$f(y\| \theta, \lambda) = \frac{\theta^2
\lambda^y (\theta + \lambda + y + 1)} {(\theta + 1)(\theta +
\lambda)^{y + 2}}\$\$

Where \\\theta\\ and \\\lambda\\ are distribution parameters with the
restrictions that \\\theta \> 0\\ and \\\lambda \> 0\\, and \\y\\ is a
non-negative integer.

The expected value of the distribution is:
\$\$\mu=\frac{\lambda(\theta + 2)}{\theta(\theta + 1)}\$\$

The default is to use the input mean value for the distribution.
However, the lambda parameter can be used as an alternative to the mean
value.

## Examples

``` r
dplind(0, mean = 0.75, theta = 7)
#> [1] 0.57
pplind(c(0, 1, 2, 3, 5, 7, 9, 10), mean = 0.75, theta = 7)
#> [1] 0.5700000 0.8160000 0.9216000 0.9667200 0.9940608 0.9989514 0.9998165
#> [8] 0.9999235
qplind(c(0.1, 0.3, 0.5, 0.9, 0.95), lambda = 4.67, theta = 7)
#> [1] 0 0 0 2 3
rplind(30, mean = 0.75, theta = 7)
#>  [1] 0 0 1 0 0 0 0 0 0 0 1 1 2 0 1 3 1 4 0 1 1 0 1 0 0 0 0 0 1 1
```
