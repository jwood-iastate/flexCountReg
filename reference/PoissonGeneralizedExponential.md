# Poisson-Generalized-Exponential Distribution

These functions provide density, distribution function, quantile
function, and random number generation for the
Poisson-Generalized-Exponential (PGE) Distribution

## Usage

``` r
dpge(
  x,
  mean = 1,
  shape = 1,
  scale = 1,
  ndraws = 1500,
  log = FALSE,
  haltons = NULL
)

ppge(
  q,
  mean = 1,
  shape = 1,
  scale = 1,
  ndraws = 1500,
  lower.tail = TRUE,
  log.p = FALSE,
  haltons = NULL
)

qpge(p, mean = 1, shape = 1, scale = 1, ndraws = 1500)

rpge(n, mean = 1, shape = 1, scale = 1, ndraws = 1500)
```

## Arguments

- x:

  numeric value or a vector of values.

- mean:

  numeric value or vector of mean values for the distribution (the
  values have to be greater than 0). This is NOT the value of
  \\\lambda\\.

- shape:

  numeric value or vector of shape values for the shape parameter of the
  generalized exponential distribution (the values have to be greater
  than 0).

- scale:

  single value or vector of values for the scale parameter of the
  generalized exponential distribution (the values have to be greater
  than 0).

- ndraws:

  the number of Halton draws to use for the integration.

- log:

  logical; if TRUE, probabilities p are given as log(p).

- haltons:

  an optional vector of Halton draws to use instead of ndraws.

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

## Value

dpge gives the density, ppge gives the distribution function, qpge gives
the quantile function, and rpge generates random deviates.

The length of the result is determined by n for rpge, and is the maximum
of the lengths of the numerical arguments for the other functions.

## Details

`dpge` computes the density (PDF) of the PGE Distribution.

`ppge` computes the CDF of the PGE Distribution.

`qpge` computes the quantile function of the PGE Distribution.

`rpge` generates random numbers from the PGE Distribution.

The Generalized Exponential distribution can be written as a function
with a shape parameter \\\alpha\>0\\ and scale parameter \\\gamma\>0\\.
The distribution has strictly positive continuous values. The PDF of the
distribution is: \$\$f(x\|\alpha,\gamma) = \frac{\alpha}{\gamma}
\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}
e^{-\frac{x}{\gamma}}\$\$

Thus, the compound Probability Mass Function(PMF) for the PGE
distribution is: \$\$f(y\|\lambda,\alpha,\beta) = \int_0^\infty
\frac{\lambda^y x^y e^{-\lambda x}}{y!} \frac{\alpha}{\gamma}
\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}}
dx\$\$

The expected value of the distribution is: \$\$E\[y\]=\mu=\lambda
\left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right)\$\$

Where \\\psi(\cdot)\\ is the digamma function.

The variance is: \$\$\sigma^2 = \lambda
\left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right) +
\left(\frac{-\psi'(\alpha+1)+\psi'(1)}{\gamma^2}\right)\lambda^2\$\$

Where \\\psi'(\cdot)\\ is the trigamma function.

To ensure that \\\mu=e^{X\beta}\\, \\\lambda\\ is replaced with:
\$\$\lambda=\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\$\$

This results in: \$\$ f(y\|\mu,\alpha,\beta) = \int_0^\infty \frac{
\left( \frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)} \right)^y x^y
e^{ -\left( \frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)} \right) x }
}{ y! } \frac{\alpha}{\gamma} \left( 1-e^{-\frac{x}{\gamma}}
\right)^{\alpha-1} e^{-\frac{x}{\gamma}} dx \$\$

Halton draws are used to perform simulation over the lognormal
distribution to solve the integral.

## References

Gupta, R. D., & Kundu, D. (2007). Generalized exponential distribution:
Existing results and some recent developments. Journal of Statistical
planning and inference, 137(11), 3537-3547.

## Examples

``` r
dpge(0, mean=0.75, shape=2, scale=1, ndraws=2000)
#> [1] 0.5338754
ppge(c(0,1,2,3,4,5,6), mean=0.75, shape=2, scale=1, ndraws=500)
#> [1] 0.5351826 0.8199568 0.9357904 0.9782831 0.9929808 0.9978338 0.9993641
qpge(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, shape=2, scale=1, ndraws=500)
#> [1] 0 0 0 2 3
rpge(30, mean=0.75,  shape=2, scale=1, ndraws=500)
#>  [1] 0 1 0 1 0 0 3 0 2 0 1 0 0 0 0 0 2 1 1 0 3 1 0 0 0 0 0 0 0 0
```
