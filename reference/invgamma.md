# Inverse Gamma Distribution

These functions provide the density function, distribution function,
quantile function, and random number generation for the Inverse-Gamma
(IG) Distribution

## Usage

``` r
dinvgamma(x, shape = 2.5, scale = 1, log = FALSE)

pinvgamma(q, shape = 2.5, scale = 1, lower.tail = TRUE, log.p = FALSE)

qinvgamma(p, shape = 2.5, scale = 1, lower.tail = TRUE, log.p = FALSE)

rinvgamma(n, shape = 2.5, scale = 1)
```

## Arguments

- x:

  numeric value or a vector of values.

- shape:

  numeric value or vector of shape values for the distribution (the
  values have to be greater than 0).

- scale:

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

dinvgamma gives the density, pinvgamma gives the distribution function,
qinvgamma gives the quantile function, and rinvgamma generates random
deviates.

The length of the result is determined by n for rinvgamma, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Details

`dinvgamma` computes the density (PDF) of the Inverse-Gamma
Distribution.

`pinvgamma` computes the CDF of the Inverse-Gamma Distribution.

`qinvgamma` computes the quantile function of the Inverse-Gamma
Distribution.

`rinvgamma` generates random numbers from the Inverse-Gamma
Distribution.

The compound Probability Mass Function (PMF) for the Inverse-Gamma
distribution: \$\$f(x \| \alpha, \beta) =
\frac{\beta^\alpha}{\Gamma(\alpha)} \left(\frac{1}{x}\right)^{\alpha+1}
e^{-\frac{\beta}{x}}\$\$

Where \\\alpha\\ is the shape parameter and \\\beta\\ is a scale
parameter with the restrictions that \\\alpha \> 0\\ and \\\eta \> 0\\,
and \\x \> 0\\.

The CDF of the Inverse-Gamma distribution is: \$\$F(x \| \alpha, \beta)
= \frac{\alpha. \Gamma \left(\frac{\beta}{x}\right)}{\Gamma(\alpha)} =
Q\left(\alpha, \frac{\beta}{x} \right)\$\$

Where the numerator is the incomplete gamma function and \\Q(\cdot)\\ is
the regularized gamma function.

The mean of the distribution is (provided \\\alpha\>1\\):
\$\$\mu=\frac{\beta}{\alpha-1}\$\$

The variance of the distribution is (for \\\alpha\>2\\):
\$\$\sigma^2=\frac{\beta^2}{(\alpha-1)^2(\alpha-2)}\$\$

## Examples

``` r
dinvgamma(1, shape = 3, scale = 2)
#> [1] 0.5413411
pinvgamma(c(0.1, 0.5, 1, 3, 5, 10, 30), shape = 3, scale = 2)
#> [1] 4.555150e-07 2.381033e-01 6.766764e-01 9.697879e-01 9.920737e-01
#> [6] 9.988515e-01 9.999530e-01
qinvgamma(c(0.1, 0.3, 0.5, 0.9, 0.95), shape = 3, scale = 2)
#> [1] 0.3757760 0.5531635 0.7479263 1.8147745 2.4459104
rinvgamma(30, shape = 3, scale = 2)
#>  [1] 0.3560428 0.5663151 2.0227915 0.8950786 0.9087604 0.6819983 0.5824978
#>  [8] 0.4527947 0.7087080 3.4054619 0.5423602 0.5749351 1.1928321 3.8493214
#> [15] 0.6759657 1.4620309 0.4670982 0.6110604 0.4027782 1.9430930 0.8898576
#> [22] 0.8698565 0.4380061 1.5957115 0.5638894 0.5808885 0.4808884 0.2358688
#> [29] 0.8996399 1.4472548
```
