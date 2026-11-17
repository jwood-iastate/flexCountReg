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
  haltons = NULL,
  engine = c("quadrature", "halton")
)

ppge(
  q,
  mean = 1,
  shape = 1,
  scale = 1,
  ndraws = 1500,
  lower.tail = TRUE,
  log.p = FALSE,
  haltons = NULL,
  engine = c("quadrature", "halton")
)

qpge(
  p,
  mean = 1,
  shape = 1,
  scale = 1,
  ndraws = 1500,
  engine = c("quadrature", "halton")
)

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

  single value or vector of values for the rate \\b\\ of the generalized
  exponential mixing distribution (greater than 0; despite its name, the
  mixing variable is generated as \\X/b\\ with \\X \sim GE(a, 1)\\).
  Because the distribution is parameterized by its mean, \\b\\ cancels
  and does not change any result; it is only checked for validity.

- ndraws:

  the number of Halton draws (only with `engine = "halton"`).

- log:

  logical; if TRUE, probabilities p are given as log(p).

- haltons:

  an optional vector of Halton draws to use instead of ndraws (only with
  `engine = "halton"`).

- engine:

  the numerical method: `"quadrature"` (default) computes the
  log-probability by adaptive quadrature over the log of the mixing
  variable, accurate for every count including far tails; `"halton"`
  averages the Poisson probabilities over Halton draws (the method of
  flexCountReg 0.1.2 and earlier, kept for reproducibility; badly biased
  far in the tail; the default when `haltons` are given).

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

The mixing variable \\X\\ has a generalized exponential distribution
(Gupta and Kundu, 2007) with shape \\a \> 0\\ and rate \\b \> 0\\:
\$\$f(x \mid a, b) = a b \left(1 - e^{-bx}\right)^{a-1} e^{-bx}, \quad x
\> 0,\$\$ with \$\$E\[X\] = \frac{A}{b}, \quad Var\[X\] = \frac{B}{b^2},
\quad A = \psi(a+1) - \psi(1), \quad B = \psi'(1) - \psi'(a+1),\$\$
where \\\psi\\ and \\\psi'\\ are the digamma and trigamma functions.
Given \\X\\, \\Y \sim Poisson(\lambda X)\\, so \$\$P(Y = y) =
\int_0^\infty \frac{(\lambda x)^y e^{-\lambda x}}{y!} f(x \mid a, b)\\
dx.\$\$ To make \\E\[Y\] = \mu\\ (`mean`, \\e^{X\beta}\\ in a
regression), \\\lambda = \mu b / A\\. Then \$\$E\[Y\] = \mu, \qquad
Var\[Y\] = \mu + \frac{B}{A^2}\mu^2.\$\$ Since \\\lambda X = \mu (bX) /
A\\ and \\bX \sim GE(a, 1)\\, the distribution of \\Y\\ depends on
\\\mu\\ and \\a\\ only: the rate (or scale) cancels. This is why PGE
regression models estimate only the shape (`ln(shape)`); a scale
parameter would not be identified. With \\a = 1\\ the mixing
distribution is exponential and \\Y\\ is geometric with mean \\\mu\\.

The integral has no closed form. By default (`engine = "quadrature"`) it
is computed on the log scale by adaptive quadrature over \\\log X\\,
centered at the mode of the integrand, which is accurate for every
count. `engine = "halton"` averages the Poisson probabilities over
Halton draws of \\X\\, as flexCountReg 0.1.2 and earlier did.

## References

Gupta, R. D., & Kundu, D. (2007). Generalized exponential distribution:
Existing results and some recent developments. Journal of Statistical
planning and inference, 137(11), 3537-3547.

## Examples

``` r
dpge(0, mean=0.75, shape=2, scale=1, ndraws=2000)
#> [1] 0.5333333
ppge(c(0,1,2,3,4,5,6), mean=0.75, shape=2, scale=1, ndraws=500)
#> [1] 0.5333333 0.8177778 0.9339259 0.9769086 0.9920895 0.9973205 0.9990983
qpge(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, shape=2, scale=1, ndraws=500)
#> [1] 0 0 0 2 3
rpge(30, mean=0.75,  shape=2, scale=1, ndraws=500)
#>  [1] 1 1 0 1 1 0 0 0 2 1 0 0 0 0 2 0 2 0 0 0 0 0 0 0 1 1 0 0 1 1
```
