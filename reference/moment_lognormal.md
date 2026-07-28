# Raw Moment of a Lognormal Distribution

Computes a raw moment of a lognormal distribution using its closed-form
expression. The lognormal distribution is specified by its log-mean
(\\\mu\\) and log-standard deviation (\\\sigma\\).

## Usage

``` r
moment_lognormal(mu, sigma, n)
```

## Arguments

- mu:

  Numeric vector of means of the log-transformed variable, corresponding
  to \\\mu\\.

- sigma:

  Non-negative numeric vector of standard deviations of the
  log-transformed variable, corresponding to \\\sigma\\.

- n:

  Numeric vector specifying the order of the raw moment.

## Value

A numeric vector containing the raw moments of the specified lognormal
distribution.

## Details

Let \\X \sim N(\mu, \sigma^2)\\ and \\Y = \exp(X)\\. Then \\Y\\ follows
a lognormal distribution, and its raw moment of order \\n\\ is: \$\$
E\[Y^n\] = \exp\left(n\mu + \frac{n^2\sigma^2}{2}\right). \$\$

In particular, the mean is: \$\$ E\[Y\] = \exp\left(\mu +
\frac{\sigma^2}{2}\right), \$\$ and the variance can be calculated as:
\$\$ \mathrm{Var}(Y) = E\[Y^2\] - E\[Y\]^2. \$\$

This function can be used to adjust predictions from generalized linear
mixed models with normally distributed random parameters and a log link.

## Examples

``` r
mu <- 0
sigma <- 1

# Mean of the lognormal distribution
moment_lognormal(mu, sigma, n = 1)
#> [1] 1.648721

# Variance of the lognormal distribution
moment_1 <- moment_lognormal(mu, sigma, n = 1)
moment_2 <- moment_lognormal(mu, sigma, n = 2)
moment_2 - moment_1^2
#> [1] 4.670774
```
