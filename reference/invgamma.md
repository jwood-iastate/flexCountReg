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

dinvgamma gives the density, pinvgamma gives the distribution function,
qinvgamma gives the quantile function, and rinvgamma generates random
deviates.

The length of the result is determined by n for rinvgamma, and is the
maximum of the lengths of the numerical arguments for the other
functions.

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
#>  [1] 0.2937822 0.4909691 0.5750880 0.9867140 0.6353773 1.9819641 0.8870580
#>  [8] 1.0057749 0.6611418 0.4277593 2.5997910 0.6948481 0.3511642 0.4924415
#> [15] 0.9255774 6.1444962 1.3036911 1.4501172 0.6066221 0.9371221 0.4320380
#> [22] 1.6987338 2.1125857 0.3442973 0.3036261 1.3134867 0.4868805 1.2093836
#> [29] 1.2170003 0.5214717
```
