# Poisson-Inverse-Gaussian Distribution

These functions provide the density function, distribution function,
quantile function, and random number generation for the
Poisson-Inverse-Gaussian (PInvGaus) Distribution.

These functions provide the density function, distribution function,
quantile function, and random number generation for the
Poisson-Inverse-Gaussian (PInvGaus) Distribution

## Usage

``` r
dpinvgaus(x, mu = 1, eta = 1, form = "Type 1", log = FALSE)

ppinvgaus(
  q,
  mu = 1,
  eta = 1,
  form = "Type 1",
  lower.tail = TRUE,
  log.p = FALSE
)

qpinvgaus(p, mu = 1, eta = 1, form = "Type 1")

rpinvgaus(n, mu = 1, eta = 1, form = "Type 1")

dpinvgaus(x, mu = 1, eta = 1, form = "Type 1", log = FALSE)

ppinvgaus(
  q,
  mu = 1,
  eta = 1,
  form = "Type 1",
  lower.tail = TRUE,
  log.p = FALSE
)

qpinvgaus(p, mu = 1, eta = 1, form = "Type 1")

rpinvgaus(n, mu = 1, eta = 1, form = "Type 1")
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

- form:

  optional parameter indicating which formulation to use. Options
  include "Type 1" which is the standard form or "Type 2" which follows
  the formulation by Dean et. al. (1987).

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

description dpinvgaus gives the density, ppinvgaus gives the
distribution function, qpinvgaus gives the quantile function, and
rpinvgaus generates random deviates.

The length of the result is determined by n for rpinvgaus, and is the
maximum of the lengths of the numerical arguments for the other
functions.

dpinvgaus gives the density, ppinvgaus gives the distribution function,
qpinvgaus gives the quantile function, and rpinvgaus generates random
deviates.

The length of the result is determined by n for rpinvgaus, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

The Poisson-Inverse-Gaussian distribution is a special case of the
Sichel distribution (Cameron & Trivedi, 2013). It is also known as a
univariate Sichel distribution (Hilbe, 2011).

`dpinvgaus` computes the PDF of the Poisson-Inverse-Gaussian dist.

`ppinvgaus` computes the CDF of the Poisson-Inverse-Gaussian dist.

`qpinvgaus` computes quantiles of the Poisson-Inverse-Gaussian dist.

`rpinvgaus` generates random numbers from the distribution.

The PMF (Type 1) is: \$\$ f(y\|\eta,\mu)=\begin{cases}
f(0)=\exp\left(\frac{\mu}{\eta}(1-\sqrt{1+2\eta})\right)\\
f(y\>0)=f(0)\frac{\mu^y}{y!}(1+2\eta)^{-y/2}
\sum\_{j=0}^{y-1}\frac{\Gamma(y+j)} {\Gamma(y-j)\Gamma(j+1)}
\left(\frac{\eta}{2\mu}\right)^j(1+2\eta)^{-j/2} \end{cases}\$\$

The variance is: \$\$\sigma^2=\mu+\eta\mu\$\$

Type 2 modifies \\\eta\\ → \\\eta\mu\\: \$\$
f(0)=\exp\left(\frac{1}{\eta}(1-\sqrt{1+2\eta\mu})\right) \$\$ \$\$
f(y\>0)=f(0)\frac{\mu^y}{y!}(1+2\eta\mu)^{-y/2}
\sum\_{j=0}^{y-1}\frac{\Gamma(y+j)} {\Gamma(y-j)\Gamma(j+1)}
\left(\frac{\eta}{2}\right)^j(1+2\eta\mu)^{-j/2} \$\$

Resulting variance: \$\$\sigma^2=\mu+\eta\mu^2\$\$

The Poisson-Inverse-Gaussian distribution is a special case of the
Sichel distribution, as noted by Cameron & Trivedi (2013). It is also
known as a univariate Sichel distribution (Hilbe, 2011).

`dpinvgaus` computes the density (PDF) of the Poisson-Inverse-Gaussian
Distribution.

`ppinvgaus` computes the CDF of the Poisson-Inverse-Gaussian
Distribution.

`qpinvgaus` computes the quantile function of the
Poisson-Inverse-Gaussian Distribution.

`rpinvgaus` generates random numbers from the Poisson-Inverse-Gamma
Distribution.

The compound Probability Mass Function (PMF) for the
Poisson-Inverse-Gaussian distribution (Type 1) is (Cameron & Trivedi,
2013): \$\$f(y\|\eta,\mu) = \begin{cases} f(y = 0) = \exp \left(
\frac{\mu}{\eta} \left(1-\sqrt{1 + 2\eta}\right)\right) \\ f(y\|y\>0) =
f(y = 0)\frac{\mu^y}{y!}(1 + 2\eta)^{-y / 2} \cdot \sum\_{j=0}^{y-1}
\frac{\Gamma(y+j)}{\Gamma(y-j)\Gamma(j+1)} \left(
\frac{\eta}{2\mu}\right)^2(1 + 2\eta)^{-j / 2} \end{cases}\$\$

Where \\\eta\\ is a scale parameter with the restriction that
\\\eta\>0\\, \\\mu\\ is the mean value, and \\y\\ is a non-negative
integer.

The variance of the distribution is: \$\$\sigma^2=\mu+\eta\mu\$\$

The alternative parametrization by Dean et. al. (1987) replaces \\\eta\\
with \\\eta\mu\\. This version (Type 2) has the PMF:
\$\$f(y\|\eta,\mu)=\begin{cases} f(y=0) = \exp \left(\frac{1}{\eta}
\left(1-\sqrt{1+2\eta\mu}\right) \right) \\ f(y\|y \> 0) = f(y=0)
\frac{\mu^y}{y!} (1+2\eta\mu)^{-y/2} \cdot \sum\_{j=0}^{y-1}
\frac{\Gamma(y+j)}{\Gamma(y-j)\Gamma(j+1)} \left(\frac{\eta}{2}
\right)^2(1+2\eta\mu)^{-j/2} \end{cases}\$\$

This results in the variance of: \$\$\sigma^2=\mu+\eta\mu^2\$\$

## References

Cameron & Trivedi (2013). Regression Analysis of Count Data. Dean,
Lawless & Willmot (1989). Mixed Poisson–Inverse Gaussian Models. Hilbe
(2011). Negative Binomial Regression.

Cameron, A. C., & Trivedi, P. K. (2013). Regression analysis of count
data, 2nd Edition. Cambridge university press.

Dean, C., Lawless, J. F., & Willmot, G. E. (1989). A mixed
Poisson–Inverse‐Gaussian regression model. Canadian Journal of
Statistics, 17(2), 171-181.

Hilbe, J. M. (2011). Negative binomial regression. Cambridge University
Press.

## Examples

``` r
dpinvgaus(1, mu=0.75, eta=1)
#> [1] 0.250067
ppinvgaus(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3, form="Type 2")
#> [1] 0.6386475 0.8428877 0.9173222 0.9512540 0.9797144 0.9904002 0.9951184
#> [8] 0.9964543
qpinvgaus(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5, form="Type 2")
#> [1] 0 0 0 2 3
rpinvgaus(30, mu=0.75, eta=1.5)
#>  [1] 3 1 0 0 0 1 6 0 0 0 0 0 1 0 4 0 2 2 0 0 0 0 1 0 0 0 2 3 0 0

dpinvgaus(1, mu=0.75, eta=1)
#> [1] 0.250067
ppinvgaus(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3, form="Type 2")
#> [1] 0.6386475 0.8428877 0.9173222 0.9512540 0.9797144 0.9904002 0.9951184
#> [8] 0.9964543
qpinvgaus(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5, form="Type 2")
#> [1] 0 0 0 2 3
rpinvgaus(30, mu=0.75, eta=1.5)
#>  [1]  0  0 11  0  2  0  0  0  0  0  0  0  0  0  0  0  2  0  1  3  0  1  0  0  0
#> [26]  0  1  0  1  0
```
