# Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution

Probability mass, cumulative distribution, quantiles, and random
generation for the mean-parameterized Poisson-Lindley-Gamma (PLG), also
called the negative binomial-Lindley (NB-L), distribution. The PMF uses
a guarded GSL Tricomi-U evaluation with adaptive quadrature as a
fallback.

## Usage

``` r
dplindGamma(
  x,
  mean = 1,
  theta = 1,
  alpha = 1,
  log = FALSE,
  rel.tol = 1e-08,
  subdivisions = 200L,
  method = c("auto", "quadrature")
)

pplindGamma(
  q,
  mean = 1,
  theta = 1,
  alpha = 1,
  lower.tail = TRUE,
  log.p = FALSE,
  rel.tol = 1e-08,
  subdivisions = 200L
)

qplindGamma(p, mean = 1, theta = 1, alpha = 1)

rplindGamma(n, mean = 1, theta = 1, alpha = 1)
```

## Arguments

- x:

  Numeric vector of counts. Negative, infinite, or noninteger counts
  have probability zero; noninteger counts produce a warning.

- mean:

  Marginal mean \\\mu\\, finite and nonnegative. May be a vector for the
  density, CDF, and quantile functions; must be scalar for random
  generation. A zero mean gives a point mass at zero.

- theta:

  Lindley parameter, finite and strictly positive. May be a vector
  except for random generation.

- alpha:

  Gamma dispersion, finite and nonnegative. For positive `alpha`, the
  mean-one gamma factor has shape and rate \\1/\alpha\\, and variance
  \\\alpha\\. Zero selects the exact Poisson-Lindley limit. May be a
  vector except for random generation.

- log:

  Logical; return log probabilities from `dplindGamma`.

- rel.tol:

  Relative tolerance, a scalar between `1e-12` and `1e-3`; default
  `1e-8`. Controls adaptive quadrature and, for the PMF, screening of
  GSL error estimates and log-formula roundoff. It is a numerical
  target, not a certified bound on the final error.

- subdivisions:

  Maximum number of subintervals for each side of each component
  integral; an integer of at least 50, default 200.

- method:

  PMF evaluation method. `"auto"` (default) attempts GSL outside the
  guarded regions and uses quadrature if any check fails. `"quadrature"`
  always uses adaptive quadrature and is useful for validation. There is
  deliberately no option to force an unchecked GSL value. This argument
  applies only to `dplindGamma`.

- q:

  Numeric vector of quantiles. Finite values are rounded down.

- lower.tail:

  Logical; if `TRUE`, return \\P(Y\le q)\\; otherwise return
  \\P(Y\>q)\\.

- log.p:

  Logical; return log probabilities from `pplindGamma`.

- p:

  Numeric vector of probabilities in \\\[0,1\]\\.

- n:

  Number of observations to generate: a nonnegative integer scalar.

## Value

`dplindGamma` returns probability masses or log masses; `pplindGamma`
returns cumulative or survival probabilities, possibly on the log scale;
`qplindGamma` returns quantiles; and `rplindGamma` returns `n` random
counts.

## Details

**Parameterization and probability mass.** Let \\L\\ have Lindley
density \$\$f_L(l)=\frac{\theta^2}{1+\theta}(1+l)e^{-\theta l},\quad
l\>0,\$\$ and let \\W=L/E(L)\\, where
\\E(L)=(\theta+2)/\[\theta(\theta+1)\]\\. Conditional on \\W\\, \\Y\\ is
negative binomial with size \\1/\alpha\\ and mean \\\mu W\\.
Equivalently, \\Y\mid W,G\\ is Poisson with mean \\\mu WG\\, with
independent \\G\sim\mathrm{Gamma}(1/\alpha,1/\alpha)\\, using the
shape-rate parameterization. Thus `mean` is the marginal mean, not the
unnormalized Lindley scale.

The compound probability mass function (PMF) for the
Poisson-Lindley-Gamma distribution is \$\$ f(x\mid\mu,\theta,\alpha)=
\frac{ (\theta+2)^2\Gamma(x+1/\alpha) }{
\alpha\mu^2(\theta+1)^3\Gamma(1/\alpha) } \left(
\frac{\mu\theta(\theta+1)}{\theta+2} U\left(
x+1,2-1/\alpha,\frac{\theta+2}{\alpha\mu(\theta+1)} \right)
+\frac{x+1}{\alpha} U\left(
x+2,3-1/\alpha,\frac{\theta+2}{\alpha\mu(\theta+1)} \right) \right).
\$\$ This expression applies for \\\mu\>0\\, \\\theta\>0\\,
\\\alpha\>0\\, and nonnegative integer \\x\\. Here \\U\\ is Tricomi's
confluent hypergeometric function of the second kind, with
positive-integrand representation
\$\$U(a,b,z)=\frac{1}{\Gamma(a)}\int_0^\infty
e^{-zt}t^{a-1}(1+t)^{b-a-1}\\dt,\quad a\>0,\\ z\>0.\$\$ This integral
exists for every real \\b\\; integer or nonpositive \\b\\ is not, by
itself, a mathematical singularity of \\U\\. The expected value of the
distribution is \$\$E\[x\]=\mu.\$\$ The variance is
\$\$\sigma^2=\mu+\left((1+\alpha)
\left(2-\frac{2}{(\theta+2)^2}\right)-1\right)\mu^2.\$\$

**GSL instability and automatic fallback.** In the PMF above, both
Tricomi terms have third argument
\$\$z=\frac{\theta+2}{\alpha\mu(\theta+1)}.\$\$ The ROOT discussion
cited below reports GSL series-convergence failures for \\U(a,b,z)\\
when \\a\>8\\ and \\z\\ is close to zero; an example is
\\U(9.50606,0.5,0.000160903)\\. The final correction in that thread says
\\a\>8\\, superseding the earlier statement \\a\<10\\. These are
historical implementation failures, not singularities of the function.
The discussion also concerns Kummer's first-kind function at negative
arguments; those examples are not the basis of this PMF's guard, since
its Tricomi argument \\z\\ is positive.

Additional checks against positive-integrand quadrature with GSL 2.7.1
found underestimated error for large negative second arguments, e.g.
`x = 0, mean = 1000, theta = 1, alpha = 0.001` (\\1/\alpha=1000\\,
\\z=1.5\\), where the relative PMF discrepancy was about
\\4\times10^{-8}\\ despite successful GSL status. The guard
\\0\<\alpha\<0.01\\ conservatively excludes this large-negative-b
regime; it too is an implementation policy, not a mathematical boundary.

In `method = "auto"`, both U calls are bypassed whenever either first
argument exceeds 8 and \\z\<1\\. Here the first arguments are \\x+1\\
and \\x+2\\, so this guard applies to \\x\ge7\\. The cutoff `1` is a
conservative implementation choice around the reported
small-positive-argument failures, not a published universal boundary.
Tests with GSL 2.7.1 also found inaccurate, successful-status results
outside \\z\<0.01\\ (including \\z\\ near 0.12), motivating the wider
\\z\<1\\ guard. Behavior depends on the linked GSL version and on all
three arguments. For numerical and computational safeguards, quadrature
is also used directly when \\0\<\alpha\<0.01\\, either first argument
exceeds \\10^4\\, or \\z\\ cannot be represented as a positive normal
double.

Otherwise, the two U terms are evaluated together with
`gsl::hyperg_U(..., give = TRUE, strict = TRUE)`. The complete PMF is
recomputed by quadrature if GSL throws an error or warning; returns
`NULL`, missing/malformed fields, `NA`, `NaN`, or infinite values;
reports nonzero status; returns a nonpositive or subnormal U value; or
supplies an invalid error estimate. Each estimated relative U error must
be at most `rel.tol / 4`. A conservative estimate of roundoff in the log
formula must also be at most `rel.tol / 4`. Invalid component
probabilities or a final probability outside \\(0,1\]\\ trigger the same
fallback. GSL errors are not exposed as warnings when quadrature
subsequently succeeds. An unavailable GSL namespace also causes fallback
when this file is sourced independently.

Gamma ratios, powers, and the sum of the two positive terms are
evaluated on the log scale. A zero U result is treated as possible
underflow, never as proof of a zero probability. Neither success status
nor an error estimate proves accuracy for every argument. Use
`method = "quadrature"` for independent checks in the parameter range of
an application, especially after changing the linked GSL version.

**Adaptive quadrature and cumulative probabilities.** Write \\X=\theta
L\\. With probability \\\theta/(1+\theta)\\,
\\X\sim\mathrm{Gamma}(1,1)\\; otherwise, \\X\sim\mathrm{Gamma}(2,1)\\.
With \\c=(\theta+2)/(\theta+1)\\, \\W=X/c\\. Quadrature integrates the
conditional NB PMF separately against these two positive gamma
densities. Each integral uses \\t=\log X\\, is scaled at its mode, and
is integrated on both sides of the mode with separate width scales.
Component probabilities are combined using log-sum-exp. The fallback
does not use a difference of hypergeometric functions, truncate the
count support, or substitute an arbitrary positive probability floor.

`pplindGamma` always uses quadrature of the conditional NB CDF or
survival function. Upper tails are integrated directly; probabilities
near one are obtained from the opposite small tail. At `alpha = 0`, the
conditional routines are Poisson and GSL is bypassed. Small positive
`alpha` is never silently replaced by zero. If quadrature itself fails,
the affected result is `NaN` with a warning; extreme inputs need not be
resolvable in double precision. On the ordinary probability scale,
sufficiently small probabilities can still underflow to zero; use
`log = TRUE` or `log.p = TRUE` for such tails.

`qplindGamma` uses sequential inversion of `pplindGamma` with its
default integration controls, returning the smallest nonnegative integer
with CDF at least `p`. The endpoints are zero for `p = 0` and infinity
for `p = 1` when `mean > 0`; for `mean = 0`, every valid quantile is
zero. `rplindGamma` retains inverse-CDF generation using uniform draws
and this quantile function. Quantiles and simulation therefore do not
benefit from the GSL PMF path and can be slow for large means or
probabilities near one.

Numeric arguments of `dplindGamma` and `pplindGamma` are recycled to
their maximum length; any zero-length numeric argument produces
`numeric(0)`. Their logical flags and numerical controls must be scalar.
Missing values propagate. Invalid distribution parameters produce `NaN`
with a warning. The quantile function uses `Vectorize` recycling and
simplification; random generation requires scalar parameters. Names on
`x` or `q` are retained when their length equals the output length.

## References

ROOT Forum, *Problem in Hypergeometric functions with GSL*, including
the April 19, 2010 correction and October 3, 2014 follow-up:
<https://root-forum.cern.ch/t/problem-in-hypergeometric-functions-with-gsl/9507>.

GNU Scientific Library Reference Manual, Hypergeometric Functions:
<https://www.gnu.org/software/gsl/doc/html/specfunc.html>.

## See also

[`hyperg_U`](https://rdrr.io/pkg/gsl/man/Hyperg.html),
[`integrate`](https://rdrr.io/r/stats/integrate.html),
[`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html)

## Examples

``` r
dplindGamma(0:5, mean = 0.75, theta = 7, alpha = 2)
#> [1] 0.69879293 0.15001417 0.06221417 0.03203763 0.01849751 0.01148442

# Compare the hybrid evaluator with the quadrature reference.
x <- 0:10
fast <- dplindGamma(x, mean = 2, theta = 1, alpha = 0.5)
reference <- dplindGamma(x, mean = 2, theta = 1, alpha = 0.5,
                        method = "quadrature")
max(abs(fast - reference))
#> [1] 6.27276e-15

# x >= 7 and z = (theta+2)/(alpha*mean*(theta+1)) < 1:
# the small-z guard sends this PMF directly to quadrature.
dplindGamma(8, mean = 1e4, theta = 1, alpha = 2/3, log = TRUE)
#> [1] -8.444742

# Log probabilities remain useful when ordinary probabilities underflow.
dplindGamma(10000, mean = 0.1, theta = 1, alpha = 0.5, log = TRUE)
#> [1] -1067.923
pplindGamma(10, mean = 0.75, theta = 3, alpha = 0.5,
            lower.tail = FALSE, log.p = TRUE)
#> [1] -6.722337
dplindGamma(0:3, mean = 0.75, theta = 2, alpha = 0)
#> [1] 0.56320000 0.25190400 0.10838016 0.04538696
qplindGamma(c(0, 0.5, 0.9, 1), mean = 1.67, theta = 0.5, alpha = 0.5)
#> [1]   0   1   4 Inf
rplindGamma(5, mean = 0.5, theta = 0.5, alpha = 2)
#> [1] 0 2 2 0 0
```
