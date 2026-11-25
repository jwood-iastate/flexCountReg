# Count Families for countreg() and countreg_rp()

Overview of the count distributions that can be used as the `family` of
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
and
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md):
their distribution parameters (the parameters modeled by
`dis_param_formula_1` and `dis_param_formula_2`), variances,
distribution functions, and the model features available for each.

## Details

In every family the mean model is \\\mu = \exp(X\beta)\\ (times the
exposure given by an offset, and the reporting probability in an
underreporting model). The distribution parameters are
\\\exp(W\delta)\\, where \\W\\ is the design matrix of
`dis_param_formula_1` or `dis_param_formula_2` (an intercept when the
formula is not given); the coefficient labels in the fitted model show
the scale of \\W\delta\\. For example `ln(alpha)` means \\\alpha =
\exp(W\delta)\\, `logit(rho)` means \\\rho = 1 / (1 + \exp(-W\delta))\\,
and the Sichel `gamma` is \\W\delta\\ itself.

In the variance column, \\\mu\\ is the mean parameter \\\exp(X\beta)\\;
it is the mean of the distribution except where noted.
`A = digamma(k + 1) - digamma(1)` and
`B = trigamma(1) - trigamma(k + 1)` for a generalized exponential shape
`k`, and \\M_2 = 1 + B / A^2\\. Family names are not case-sensitive and
non-alphanumeric characters are ignored (`"nb-2"` is `"NB2"`).

### Parameters and moments

|  |  |  |  |  |
|----|----|----|----|----|
| `family` | Distribution | `dis_param_formula_1` | `dis_param_formula_2` | Variance |
| `POISSON` | Poisson | – | – | `mu` |
| `NB1` | Negative binomial 1 | `ln(alpha)` | – | `mu * (1 + alpha)` |
| `NB2` | Negative binomial 2 | `ln(alpha)` | – | `mu + alpha * mu^2` |
| `NBP` | Negative binomial P | `ln(alpha)` | `ln(p)` | `mu + alpha * mu^p` |
| `GP1` | Generalized Poisson 1 | `ln(phi+1)` | – | `(1 + phi)^2 * mu` for `phi >= 0` (1) |
| `GP2` | Generalized Poisson 2 | `ln(alpha+1)` | – | `(1 + alpha * mu)^2 * mu` for `alpha >= 0` (1) |
| `COM` | Conway-Maxwell-Poisson | `ln(nu)` | – | computed numerically |
| `PLN` | Poisson-lognormal | `ln(sigma)` | – | `m + (exp(sigma^2) - 1) * m^2`, mean `m = mu * exp(sigma^2 / 2)` |
| `PGE` | Poisson-generalized exponential | `ln(shape)` | – | `mu + mu^2 * B / A^2` |
| `PIG1` | Poisson-inverse Gaussian 1 | `ln(eta)` | – | `mu * (1 + eta)` |
| `PIG2` | Poisson-inverse Gaussian 2 | `ln(eta)` | – | `mu + eta * mu^2` |
| `PIG` | Poisson-inverse gamma | `ln(eta)` | – | `mu + eta * mu^2` |
| `PL` | Poisson-Lindley | `ln(theta)` | – | `mu + mu^2 * (1 - 2 / (theta + 2)^2)` |
| `PLG` | Poisson-Lindley-gamma (negative binomial-Lindley) | `ln(theta)` | `ln(alpha)` | `mu + ((1 + alpha) * (2 - 2 / (theta + 2)^2) - 1) * mu^2` |
| `PLL` | Poisson-Lindley-lognormal | `ln(theta)` | `ln(sigma)` | `m + ((2 - 2 / (theta + 2)^2) * exp(sigma^2) - 1) * m^2`, mean `m = mu * exp(sigma^2 / 2)` |
| `PW` | Poisson-Weibull | `ln(alpha)` | – | `mu + (gamma(1 + 2 / alpha) / gamma(1 + 1 / alpha)^2 - 1) * mu^2` |
| `SI` | Sichel | `gamma` (identity link) | `ln(sigma)` | see [`dsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md) |
| `GW` | Generalized Waring | `ln(k)` | `ln(rho-1)` | finite only for `rho > 2`; see [`dgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md) |
| `DW` | Discrete Weibull | `ln(shape)` | – | computed numerically (2) |
| `DPO` | Double Poisson | `ln(phi)` | – | computed numerically, about `mu / phi` (2) |
| `DELAPORTE` | Delaporte | `logit(rho)` | `ln(k)` | `mu + (rho * mu)^2 / k` |
| `HERMITE` | Hermite (order 2) | `logit(rho)` | – | `mu * (1 + rho)` |
| `HERMITE3` | Generalized Hermite, order 3 | `logit(rho)` | – | `mu * (1 + 2 * rho)` |
| `HERMITE4` | Generalized Hermite, order 4 | `logit(rho)` | – | `mu * (1 + 3 * rho)` |
| `HYPERPOIS` | Hyper-Poisson | `ln(gamma)` | – | computed numerically (2) |
| `GAMMACOUNT` | Gamma renewal count | `ln(shape)` | – | computed numerically (2) |
| `PA` | Polya-Aeppli (Poisson-Tweedie, power 1.5) | `ln(phi)` | – | `mu + phi * mu^1.5` |
| `NBGE` | Negative binomial-generalized exponential | `ln(shape)` | `ln(alpha)` | `mu + ((1 + alpha) * M2 - 1) * mu^2` |

\(1\) With a negative dispersion parameter the generalized Poisson
distributions are truncated and renormalized, \\\mu\\ is a location
parameter, and the moments are computed from the renormalized
probabilities (see
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)).

\(2\) The natural parameter of each observation is found by solving the
mean equation numerically, so \\\mu\\ is exactly the mean. These
families are slower to fit, especially in
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md);
their shape parameter must lie in \\\[10^{-6}, 10^6\]\\ and the mean be
at most \\10^9\\ (outside this numerical domain the likelihood is
treated as 0). For `DPO` the probabilities are normalized by summing at
most two million terms, so in addition \\\mu / \phi\\ (about the
variance) must be below about \\10^{10}\\; very small \\\phi\\ (below
about \\10^{-4}\\) also makes each evaluation slow.

### Distribution functions and model features

The distribution functions use each distribution's own parameterization,
which can differ from the regression parameterization (see each help
page). "Underreporting" is the binomial-thinning model of
`underreport_formula`; it requires the reported counts to follow the
same family. "LR vs Poisson" is the likelihood-ratio test reported by
[`summary.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/summary.flexCountReg.md):
"boundary" when Poisson is on the edge of the parameter space (50:50
chi-square mixture), "interior" when it is an interior special case, and
"–" when no regular test exists.

|  |  |  |  |
|----|----|----|----|
| `family` | Distribution functions | Underreporting | LR vs Poisson |
| `POISSON` | [`stats::dpois()`](https://rdrr.io/r/stats/Poisson.html) | yes | – |
| `NB1` | [`stats::dnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) with `size = mu / alpha` | yes | boundary |
| `NB2` | [`stats::dnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) with `size = 1 / alpha` | yes | boundary |
| `NBP` | [`stats::dnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) with `size = mu^(2 - p) / alpha` | yes | – |
| `GP1` | [`dgp1()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson.md) | – | interior |
| `GP2` | [`dgp2()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson2.md) | – | interior |
| `COM` | [`dcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md) | – | interior |
| `PLN` | [`dpLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLognormal.md) | yes | boundary |
| `PGE` | [`dpge()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md) | yes | – |
| `PIG1` | [`dpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md) | yes | boundary |
| `PIG2` | [`dpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md) with `form = "Type 2"` | yes | boundary |
| `PIG` | [`dpinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md) | yes | boundary |
| `PL` | [`dplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md) | yes | – |
| `PLG` | [`dplindGamma()`](https://jwood-iastate.github.io/flexCountReg/reference/NegativeBinomialLindley.md) | yes | – |
| `PLL` | [`dplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md) | yes | – |
| `PW` | [`dpoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md) | yes | – |
| `SI` | [`dsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md) | yes | – |
| `GW` | [`dgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md) | – | – |
| `DW` | [`ddweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md) | – | – |
| `DPO` | [`ddpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md) | – | interior |
| `DELAPORTE` | [`ddelaporte()`](https://jwood-iastate.github.io/flexCountReg/reference/Delaporte.md) | yes | – |
| `HERMITE` | [`dghermite()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md) with `m = 2` | yes (`rho` becomes `p * rho`) | – |
| `HERMITE3` | [`dghermite()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md) with `m = 3` | – | – |
| `HERMITE4` | [`dghermite()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md) with `m = 4` | – | – |
| `HYPERPOIS` | [`dhpois()`](https://jwood-iastate.github.io/flexCountReg/reference/HyperPoisson.md) | – | interior |
| `GAMMACOUNT` | [`dgcount()`](https://jwood-iastate.github.io/flexCountReg/reference/GammaCount.md) | – | interior |
| `PA` | [`dpaeppli()`](https://jwood-iastate.github.io/flexCountReg/reference/PolyaAeppli.md) | yes (`phi` becomes `sqrt(p) * phi`) | – |
| `NBGE` | [`dnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md) | yes | – |

Zero inflation and hurdles (`zero_formula`) are available for every
family in
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
and every family can be used with random parameters in
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
(GP2 only with a nonnegative dispersion).

### Finite mixtures and latent classes

Every family can be a component of a finite mixture or latent-class
model
([`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md),
[`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)),
with the same component density, mean, variance, distribution function
and random generator as in
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
(one-component mixtures reproduce
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)).
Mixtures of `POISSON`, `NB1` and `NB2` components are validated
(analytic scores, independent likelihood checks and simulation recovery
in the tests). Mixtures of the other families are experimental: they are
computed correctly, but how well their components and mixing proportions
can be recovered has not been studied, and adding outer components to a
family that already mixes over a continuous random effect (PLN, PGE,
PLL, NBGE, and others) can make the model weakly identified. A finite
mixture selects one component and draws one count from it; it is not a
convolution (a sum of counts) and not the within-family mixing of these
families. A Generalized Waring component with \\\rho \le 2\\ has an
infinite variance, so the mixture variance is infinite and Pearson
residuals are not available.

### Underreporting and the variance of the reported counts

In an underreporting model the true count \\Y\\ has mean \\\mu =
\exp(X\beta)\\ and variance \\V\\, as in the first table, and each event
is reported independently with probability \\p = F(Z\gamma)\\ (\\F\\ the
logistic or normal distribution function of `underreport_family`, \\Z\\
the design matrix of `underreport_formula`). The reported count is a
binomial thinning of the true count, \\Y^\* \mid Y \sim
\mathrm{Binomial}(Y, p)\\, so by the laws of total expectation and
variance \$\$E\[Y^\*\] = p\mu, \qquad \mathrm{Var}(Y^\*) = p^2 V + p(1 -
p)\mu.\$\$ The second term is the binomial reporting noise. In terms of
the dispersion index (variance-to-mean ratio) \\D = V / \mu\\ of the
true counts, the reported counts have \$\$D^\* =
\frac{\mathrm{Var}(Y^\*)}{E\[Y^\*\]} = 1 + p(D - 1).\$\$ Underreporting
therefore shrinks any departure from the Poisson variance by the factor
\\p\\: the reported counts are less overdispersed (or less
underdispersed) than the true counts, and a Poisson true count stays
Poisson. The same \\p\\ lowers the mean and shrinks the excess variance;
a heavily underreported outcome can look close to Poisson even when the
true counts are strongly overdispersed.

For the families that allow underreporting, the reported counts follow
the same family, with mean parameter \\p\mu\\ (written `m` below; it is
the reported mean except for PLN and PLL, whose mean is
\\p\mu\exp(\sigma^2 / 2)\\) and the distribution parameters shown below
(the estimated distribution parameters are those of the true counts):

|  |  |  |
|----|----|----|
| `family` | Variance of the reported counts (`m = p * mu`) | Reported-count parameters |
| `POISSON` | `m` | – |
| `NB1` | `m * (1 + p * alpha)` | `p * alpha` |
| `NB2` | `m + alpha * m^2` | unchanged |
| `NBP` | `m + alpha * p^(2 - P) * m^P` | `alpha * p^(2 - P)` |
| `PIG1` | `m * (1 + p * eta)` | `p * eta` |
| `PIG2`, `PIG` | `m + eta * m^2` | unchanged |
| `PLN`, `PGE`, `PL`, `PLG`, `PLL`, `PW`, `SI`, `NBGE` | the variance formula of the first table, evaluated at the mean parameter `p * mu` | unchanged |
| `DELAPORTE` | `m + (rho * m)^2 / k` | unchanged |
| `HERMITE` | `m * (1 + p * rho)` | `p * rho` |
| `PA` | `m + sqrt(p) * phi * m^1.5` | `sqrt(p) * phi` |

Two consequences for interpretation:

- For NB2, PIG2, PIG, DELAPORTE and the mixed Poisson families whose
  mixing distribution does not depend on the mean (PLN, PGE, PL, PLG,
  PLL, PW, SI and NBGE), the excess variance is quadratic in the mean
  and the distribution parameters are the same for true and reported
  counts: the variance of the reported counts is the family variance at
  the reported mean parameter. For NB1, PIG1, NBP, HERMITE and PA, the
  dispersion of the reported counts depends on \\p\\, so observations
  with different reporting probabilities have different reported-count
  dispersion even with a constant dispersion parameter.

- The reporting probability is identified mainly through variables that
  affect reporting but not the true counts (see `underreport_formula` in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)).
  For NB1, PIG1, NBP, HERMITE and PA the dependence of the reported
  dispersion on \\p\\ adds some information, but it is confounded with
  the dispersion parameter itself and should not be relied on to
  separate \\\mu\\ from \\p\\.

Fitted values and diagnostics use the reported counts:
`predict(type = "response")` is \\p\mu\\, `predict(type = "true")` is
\\\mu\\ and `predict(type = "prob_reported")` is \\p\\; Pearson
residuals divide by the reported-count variance \\\mathrm{Var}(Y^\*)\\;
and quantile residuals, rootograms and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) use the
reported-count distribution. The reporting probability is a fixed
function of the covariates: extra variation in reporting between
observations with the same covariates is not part of the model.

## See also

[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md),
[`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md),
[`summary.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/summary.flexCountReg.md)

## Examples

``` r
# The same mean model with two dispersion-flexible families
data("washington_roads")
f <- Total_crashes ~ lnaadt + lnlength + speed50
fit_nb2 <- countreg(f, data = washington_roads, family = "NB2",
                    method = "BFGS")
fit_pa <- countreg(f, data = washington_roads, family = "PA",
                   method = "BFGS")
c(NB2 = AIC(fit_nb2), PA = AIC(fit_pa))
#>      NB2       PA 
#> 2179.884 2179.167 
```
