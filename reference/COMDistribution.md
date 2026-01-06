# Conway-Maxwell-Poisson (COM) Distribution

These functions provide the density function, distribution function,
quantile function, and random number generation for the
Conway-Maxwell-Poisson (COM) Distribution

## Usage

``` r
dcom(x, mu = NULL, lambda = 1, nu = 1, log = FALSE)

pcom(q, mu = NULL, lambda = 1, nu = 1, lower.tail = TRUE, log.p = FALSE)

qcom(p, mu = NULL, lambda = 1, nu = 1)

rcom(n, mu = NULL, lambda = 1, nu = 1)
```

## Arguments

- x:

  numeric value or a vector of values.

- mu:

  optional. Numeric value or vector of mean values for the distribution
  (the values have to be greater than 0).

- lambda:

  optional. Numeric value or vector of values for the rate parameter of
  the distribution (the values have to be greater than 0). If \`mu\` is
  provided, \`lambda\` is ignored.

- nu:

  optional. Numeric value or vector of values for the decay parameter of
  the distribution ((the values have to be greater than 0).

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

`dcom` computes the density (PDF) of the COM Distribution.

`pcom` computes the CDF of the COM Distribution.

`qcom` computes the quantile function of the COM Distribution.

`rcom` generates random numbers from the COM Distribution.

The Probability Mass Function (PMF) for the Conway-Maxwell-Poisson
distribution is: \$\$f(x\|\lambda, \nu) = \frac{\lambda^x}{(x!)^\nu
Z(\lambda,\nu)}\$\$

Where \\\lambda\\ and \\\nu\\ are distribution parameters with
\\\lambda\>0\\ and \\\nu\>0\\, and \\Z(\lambda,\nu)\\ is the normalizing
constant.

The normalizing constant is given by:
\$\$Z(\lambda,\nu)=\sum\_{n=0}^{\infty}\frac{\lambda^n}{(n!)^\nu}\$\$

The mean and variance of the distribution are given by:
\$\$E\[x\]=\mu=\lambda \frac{\delta}{\delta \lambda}
\log(Z(\lambda,\nu))\$\$ \$\$Var(x)=\lambda \frac{\delta}{\delta
\lambda} \mu\$\$

When the mean value is given, the rate parameter (\\\lambda\\) is
computed using the mean and the decay parameter (\\\nu\\). This is
useful to allow the calculation of the rate parameter when the mean is
known (e.g., in regression))

dcom gives the density, pcom gives the distribution function, qcom gives
the quantile function, and rcom generates random deviates.

The length of the result is determined by n for rcom, and is the maximum
of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the
result. Only the first elements of the logical arguments are used.

## Examples

``` r
dcom(1, mu=0.75, nu=3)
#> [1] 0.5347848
pcom(c(0,1,2,3,5,7,9,10), lambda=0.75, nu=0.75)
#> [1] 0.4480730 0.7841277 0.9339922 0.9833004 0.9993079 0.9999816 0.9999997
#> [8] 1.0000000
qcom(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, nu=0.75)
#> [1] 0 0 1 2 2
rcom(30, mu=0.75, nu=0.5)
#>  [1] 1 0 0 0 0 0 1 1 2 0 0 0 0 0 0 0 0 3 0 1 1 0 3 1 0 1 1 1 0 0
```
