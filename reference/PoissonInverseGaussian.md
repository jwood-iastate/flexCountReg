# Poisson-Inverse-Gaussian Distribution

These functions provide the density function, distribution function,
quantile function, and random number generation for the
Poisson-Inverse-Gaussian (PInvGaus) Distribution.

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
```

## Arguments

- x:

  numeric value or a vector of values.

- mu:

  numeric value or vector of mean values for the distribution (values
  must be \> 0).

- eta:

  single value or vector of scale parameter values (must be \> 0).

- form:

  optional parameter specifying the formulation. Options: "Type 1"
  (standard) or "Type 2" (Dean et al., 1987).

- log:

  logical; if TRUE, probabilities p are returned as log(p).

- q:

  quantile or a vector of quantiles.

- lower.tail:

  logical; if TRUE, returns \\P\[X\leq x\]\\, otherwise \\P\[X \> x\]\\.

- log.p:

  logical; if TRUE, probabilities p are returned as log(p).

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

## References

Cameron & Trivedi (2013). Regression Analysis of Count Data. Dean,
Lawless & Willmot (1989). Mixed Poisson–Inverse Gaussian Models. Hilbe
(2011). Negative Binomial Regression.

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
#>  [1] 0 9 1 0 0 1 0 0 1 1 1 3 1 1 3 2 0 0 1 0 0 7 0 0 1 1 0 0 0 0
```
