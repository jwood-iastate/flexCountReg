# Moment Generating Function for a Lognormal Distribution

Computes the value of the moment generating function (MGF) for a
lognormal distribution at a given point through numerical integration.
This function is particularly useful for distributions where the MGF
does not have a closed-form solution. The lognormal distribution is
specified by its log-mean (\\\mu\\) and log-standard deviation
(\\\sigma\\).

## Usage

``` r
mgf_lognormal(mu, sigma, n)
```

## Arguments

- mu:

  The mean of the log-transformed variable, corresponding to \\\mu\\ in
  the lognormal distribution's parameters.

- sigma:

  The standard deviation of the log-transformed variable, corresponding
  to \\\sigma\\ in the lognormal distribution's parameters.

- n:

  The point at which to evaluate the MGF, often denoted as \\t\\ in the
  definition of the MGF. This parameter essentially specifies the order
  of the moment generating function.

## Value

The estimated value of the moment generating function (MGF) for the
specified lognormal distribution at the given point.

## Details

The moment generating function (MGF) for the lognormal distribution does
not have a closed form solution. The MGF is defined as: \$\$ M_x(n) =
\int_0^\infty e^{nx}\frac{1}{x\sigma\sqrt{2\pi}}
e^{-\frac{(\ln(x)-\mu)^2}{2\sigma^2}}\\dx \$\$

The MGF for the lognormal distribution is useful for adjusting the
predictions of generalized linear mixed models (GLMMs) that have
parameters that follow a lognormal distribution and use a log link
function. The adjustment for the mean value is the MGF with \\n=1\\ or
\\E\[e^x\]=M_x(n=1)\\. The variance for the lognormal random parameter
is: \$\$Var(e^x)=E\[e^{2x}\]-E\[e^x\]^2=M_x(n=2)-M_x(n=1)^2\$\$

## Examples

``` r
mu <- 0
sigma <- 1
n <- 1
mgf_value <- mgf_lognormal(mu, sigma, n)
print(mgf_value)
#> [1] 2.821428e+36
```
