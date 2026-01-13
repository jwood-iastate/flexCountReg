# Poisson-Inverse-Gamma Distribution

These functions provide the density function, distribution function,
quantile function, and random number generation for the
Poisson-Inverse-Gamma (PInvGamma) Distribution

## Usage

``` r
dpinvgamma(x, mu = 1, eta = 1, log = FALSE)

ppinvgamma(q, mu = 1, eta = 1, lower.tail = TRUE, log.p = FALSE)

qpinvgamma(p, mu = 1, eta = 1)

rpinvgamma(n, mu = 1, eta = 1)
```

## Arguments

- x:

  numeric value or a vector of values.

- mu:

  numeric value or vector of mean values for the distribution (the
  values have to be greater than 0).

- eta:

  single value or vector of values for the scale parameter of the
  distribution (the values have to be greater than 0).

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

## Value

dpinvgamma gives the density, ppinvgamma gives the distribution
function, qpinvgamma gives the quantile function, and rcom generates
random deviates.

The length of the result is determined by n for rpinvgamma, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

`dpinvgamma` computes the density (PDF) of the Poisson-Inverse-Gamma
Distribution.

`ppinvgamma` computes the CDF of the Poisson-Inverse-Gama Distribution.

`qpinvgamma` computes the quantile function of the Poisson-Inverse-Gamma
Distribution.

`rpinvgamma` generates random numbers from the Poisson-Inverse-Gamma
Distribution.

The compound Probability Mass Function (PMF) for the
Poisson-Inverse-Gamma distribution is: \$\$
f(x\|\eta,\mu)=\frac{2\left(\mu\left(\frac{1}{\eta}+1\right)\right)^{
\frac{x+\frac{1}{eta}+2}{2}}}{x!\Gamma\left(\frac{1}{\eta}+2\right)}
K\_{x-\frac{1}{\eta}-2}\left(2\sqrt{\mu\left(\frac{1}{\eta}+1\right)}\right)
\$\$

Where \\\eta\\ is a shape parameter with the restriction that
\\\eta\>0\\, \\\mu\>0\\ is the mean value, \\y\\ is a non-negative
integer, and \\K_i(z)\\ is the modified Bessel function of the second
kind. This formulation uses the mean directly.

The variance of the distribution is: \$\$\sigma^2=\mu+\eta\mu^2\$\$

## Examples

``` r
dpinvgamma(1, mu=0.75, eta=1)
#> [1] 0.2935178
ppinvgamma(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3)
#> [1] 0.5623779 0.8373348 0.9352191 0.9701692 0.9905595 0.9957622 0.9976802
#> [8] 0.9981978
qpinvgamma(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5)
#> [1] 0 0 0 2 3
rpinvgamma(30, mu=0.75, eta=1.5)
#>  [1] 0 0 0 0 0 1 0 0 1 1 1 1 0 0 0 0 3 1 1 2 4 1 1 1 0 0 0 3 0 0
```
