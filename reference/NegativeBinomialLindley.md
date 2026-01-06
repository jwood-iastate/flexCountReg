# Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution

These functions provide density, distribution function, quantile
function, and random number generation for the Poisson-Lindley-Gamma
(PLG) Distribution

## Usage

``` r
dplindGamma(x, mean = 1, theta = 1, alpha = 1, log = FALSE)

pplindGamma(
  q,
  mean = 1,
  theta = 1,
  alpha = 1,
  lower.tail = TRUE,
  log.p = FALSE
)

qplindGamma(p, mean = 1, theta = 1, alpha = 1)

rplindGamma(n, mean = 1, theta = 1, alpha = 1)
```

## Arguments

- x:

  numeric value or a vector of values.

- mean:

  numeric value or vector of mean values for the distribution (the
  values have to be greater than 0).

- theta:

  single value or vector of values for the theta parameter of the
  distribution (the values have to be greater than 0).

- alpha:

  single value or vector of values for the \`alpha\` parameter of the
  gamma distribution in the special case that the mean = 1 and the
  variance = \`alpha\` (the values for \`alpha\` have to be greater than
  0).

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

The Poisson-Lindley-Gamma is a count distribution that captures high
densities for small integer values and provides flexibility for heavier
tails.

`dplindGamma` computes the density (PDF) of the Poisson-Lindley-Gamma
Distribution.

`pplindGamma` computes the CDF of the Poisson-Lindley-Gamma
Distribution.

`qplindGamma` computes the quantile function of the
Poisson-Lindley-Gamma Distribution.

`rplindGamma` generates random numbers from the Poisson-Lindley-Gamma
Distribution.

The compound Probability Mass Function (PMF) for the
Poisson-Lindley-Gamma (PLG) distribution is: \$\$
f(x\|\mu,\theta,\alpha)= \frac{ \alpha(\theta+2)^2\Gamma(x+\alpha) }{
\mu^2(\theta+1)^3\Gamma(\alpha) } \left(
\frac{\mu\theta(\theta+1)}{\theta+2} U\left(
x+1,2-\alpha,\frac{\alpha(\theta+2)}{\mu(\theta+1)} \right) +
\alpha(x+1) U\left( x+2,3-\alpha,\frac{\alpha(\theta+2)}{\mu(\theta+1)}
\right) \right) \$\$

Where \\\theta\\ is a distribution parameter from the Poisson-Lindley
distribution with the restrictions that \\\theta\>0\\, \\\alpha\\ is a
parameter for the gamma distribution with the restriction \\\alpha\>0\\,
\\mu\\ is the mean value, and \\x\\ is a non-negative integer, and
\$\$U(a,b,z)\$\$ is the Tricomi's solution to the confluent
hypergeometric function - also known as the confluent hypergeometric
function of the second kind

The expected value of the distribution is: \$\$E\[x\]=\mu\$\$

The variance is: \$\$\sigma^2=\mu+\left(2\alpha+1-\frac{2(1+\alpha)}
{(\theta+2)^2}\right)\mu^2\$\$

While the distribution can be computed using the confluent
hypergeometric function, that function has limitations in value it can
be computed at (along with accuracy, in come cases). For this reason,
the function uses Halton draws to perform simulation over the gamma
distribution to solve the integral. This is sometimes more
computationally efficient as well.

dplindGamma gives the density, pplindGamma gives the distribution
function, qplindGamma gives the quantile function, and rplindGamma
generates random deviates.

The length of the result is determined by n for rplindGamma, and is the
maximum of the lengths of the numerical arguments for the other
functions.

## Examples

``` r
dplindGamma(0, mean=0.75, theta=7, alpha=2)
#> [1] 0.6154674
pplindGamma(c(0,1,2,3,5,7,9,10), mean=0.75, theta=3, alpha=0.5)
#> [1] 0.6965604 0.8477366 0.9106485 0.9430531 0.9733033 0.9859319 0.9920119
#> [8] 0.9938539
qplindGamma(c(0.1,0.3,0.5,0.9,0.95), mean=1.67, theta=0.5, alpha=0.5)
#> [1] 0 0 0 5 8
rplindGamma(30, mean=0.5, theta=0.5, alpha=2)
#>  [1] 2 0 0 0 0 0 0 0 0 3 0 0 1 0 3 0 1 0 0 0 0 1 0 0 0 1 0 0 0 1
```
