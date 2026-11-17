# Abalone Shell-Ring Counts

Counts of shell growth rings and physical measurements of 4,177 blacklip
abalone (*Haliotis rubra*) from Tasmania, Australia, from the UCI
Machine Learning Repository Abalone dataset. The number of rings is an
integer count that is roughly equidispersed overall but conditionally
underdispersed given the size measurements, which makes it a useful
example for the dispersion-flexible count distributions in this package.

## Usage

``` r
abalone_rings
```

## Format

A data frame with 4,177 rows (one per abalone) and 9 columns:

- sex:

  Factor with levels `F` (female), `I` (infant) and `M` (male). Source
  name `Sex`.

- length:

  Longest shell measurement, in millimeters divided by 200. Source name
  `Length`.

- diameter:

  Shell diameter perpendicular to the length, in millimeters divided
  by 200. Source name `Diameter`.

- height:

  Height with meat in the shell, in millimeters divided by 200. Source
  name `Height`.

- whole_weight:

  Weight of the whole abalone, in grams divided by 200. Source name
  `Whole weight`.

- shucked_weight:

  Weight of the meat, in grams divided by 200. Source name
  `Shucked weight`.

- viscera_weight:

  Gut weight after bleeding, in grams divided by 200. Source name
  `Viscera weight`.

- shell_weight:

  Shell weight after drying, in grams divided by 200. Source name
  `Shell weight`.

- rings:

  Number of shell rings (integer, 1 to 29). The count response. Source
  name `Rings`.

## Source

Nash, W., Sellers, T., Talbot, S., Cawthorn, A. and Ford, W. (1994).
Abalone (dataset). UCI Machine Learning Repository.
[doi:10.24432/C55C7W](https://doi.org/10.24432/C55C7W) . Archive
<https://archive.ics.uci.edu/static/public/1/abalone.zip> (retrieved
2026-11-17).

## Details

The observation unit is one abalone. In the source, examples with
missing values were removed (most of them missing the ring count) and
the continuous measurements were divided by 200; the stored values keep
that scaling. The data were collected for a 1994 population-biology
study and should not be treated as representative of current
populations.

Preparation for this package (script `data-raw/abalone_rings.R` in the
source repository, which checks the SHA-256 hashes of the publisher's
archive and data file): columns renamed to lower-case snake_case, `sex`
stored as a factor and `rings` as an integer. No rows were removed and
no values were changed.

Data-quality notes: two abalone have `height` equal to 0 and one has
`height` 1.13 (226 mm), far above all others; check these before using
logarithms of the measurements. The age in years is often approximated
as `rings + 1.5`; that derived age is not a count and should not be used
as a predictor of `rings`. The absence of zero counts reflects the
biology (every abalone has at least one ring), not a zero-truncated
sampling design.

The data are licensed under the Creative Commons Attribution 4.0
International license (CC BY 4.0,
<https://creativecommons.org/licenses/by/4.0/>); the MIT license of the
package code does not apply to them. Credit the creators when using or
redistributing them. See also the file `DATA-LICENSES.md` installed with
the package.

## References

Nash, W. J., Sellers, T. L., Talbot, S. R., Cawthorn, A. J. and Ford, W.
B. (1994). The population biology of abalone (*Haliotis* species) in
Tasmania. I. Blacklip abalone (*H. rubra*) from the north coast and
islands of Bass Strait. Sea Fisheries Division, Technical Report No. 48.

## Examples

``` r
data("abalone_rings", package = "flexCountReg")
str(abalone_rings)
#> 'data.frame':    4177 obs. of  9 variables:
#>  $ sex           : Factor w/ 3 levels "F","I","M": 3 3 1 3 2 2 1 1 3 1 ...
#>  $ length        : num  0.455 0.35 0.53 0.44 0.33 0.425 0.53 0.545 0.475 0.55 ...
#>  $ diameter      : num  0.365 0.265 0.42 0.365 0.255 0.3 0.415 0.425 0.37 0.44 ...
#>  $ height        : num  0.095 0.09 0.135 0.125 0.08 0.095 0.15 0.125 0.125 0.15 ...
#>  $ whole_weight  : num  0.514 0.226 0.677 0.516 0.205 ...
#>  $ shucked_weight: num  0.2245 0.0995 0.2565 0.2155 0.0895 ...
#>  $ viscera_weight: num  0.101 0.0485 0.1415 0.114 0.0395 ...
#>  $ shell_weight  : num  0.15 0.07 0.21 0.155 0.055 0.12 0.33 0.26 0.165 0.32 ...
#>  $ rings         : int  15 7 9 10 7 8 20 16 9 19 ...
c(mean = mean(abalone_rings$rings), var = var(abalone_rings$rings))
#>      mean       var 
#>  9.933684 10.395266 

# \donttest{
fit <- countreg(rings ~ sex + length + diameter + shell_weight,
                data = abalone_rings, family = "POISSON", method = "BFGS")
summary(fit)
#> Call:
#> countreg(formula = rings ~ sex + length + diameter + shell_weight, 
#>     data = abalone_rings, family = "POISSON", method = "BFGS")
#> 
#> Model: countreg | family: POISSON | optimizer: BFGS
#> Iterations: 3 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    1.9801     0.0390   50.81  < 2e-16 ***
#> sexI          -0.1197     0.0151   -7.93  2.2e-15 ***
#> sexM          -0.0140     0.0114   -1.22     0.22    
#> length        -1.0483     0.2461   -4.26  2.1e-05 ***
#> diameter       1.5894     0.3100    5.13  3.0e-07 ***
#> shell_weight   1.0014     0.0755   13.26  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -9751.611 | AIC: 19515.22 | BIC: 19553.25 
#> N: 4177
# Conditional dispersion: the Pearson statistic per degree of freedom is
# well below 1, i.e. the counts are underdispersed given the predictors
mu <- predict(fit, type = "response")
sum((abalone_rings$rings - mu)^2 / mu) / (nrow(abalone_rings) - 6)
#> [1] 0.6119697
# }
```
