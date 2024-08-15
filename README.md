


<!-- README.md is generated from README.Rmd. Please edit that file -->

# flexCountReg

<!-- badges: start -->

[![codecov](https://codecov.io/gh/jwood-iastate/flexCountReg/branch/main/graph/badge.svg)](https://codecov.io/gh/jwood-iastate/flexCountReg)
<!-- badges: end -->

The goal of flexCountReg is to provide functions that allow the analyst
to estimate count regression models that can handle multiple analysis
issues including excess zeros, overdispersion as a function of variables
(i.e., generalized count models), random parameters, etc.

## Installation

You can install the development version of flexCountReg like using:

``` r
# install.packages("devtools")
devtools::install_github("jwood-iastate/flexCountReg")
```

## Functions and Data

The following functions are included in the `flexCountReg` package,
grouped by continuous and count distributions.

**Distribution Functions**

*Continuous Distributions*

- Inverse Gamma Distribution

  - `dinvgamma` for the density function
  - `pinvgamma` for the cumulative density function
  - `qinvgamma` for the quantile function
  - `rinvgamma` for random number generation

- Triangle Distribution

  - `dtri` for the density function
  - `ptri` for the cumulative density function
  - `qtri` for the quantile function
  - `rtri` for random number generation

- Lognormal Distribution

  - `mgf_lognormal` for estimating the moment generating function

*Count Distributions*

- Generalized Waring Distribution

  - `dgwar` for the density function
  - `pgwar` for the cumulative density function
  - `qgwar` for the quantile function
  - `rgwar` for random number generation

- Poisson-Generalized-Exponential Distribution

  - `dpge` for the density function
  - `ppge` for the cumulative density function
  - `qpge` for the quantile function
  - `rpge` for random number generation

- Poisson-Inverse-Gaussian Distribution

  - `dpinvgaus` for the density function
  - `ppinvgaus` for the cumulative density function
  - `qpinvgaus` for the quantile function
  - `rpinvgaus` for random number generation

- Poisson-Lindley Distribution

  - `dplind` for the density function
  - `pplind` for the cumulative density function
  - `qplind` for the quantile function
  - `rplind` for random number generation

- Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution

  - `dplindGamma` for the density function
  - `pplindGamma` for the cumulative density function
  - `qplindGamma` for the quantile function
  - `rplindGamma` for random number generation

- Poisson-Lindley-Lognormal Distribution

  - `dplindLnorm` for the density function
  - `pplindLnorm` for the cumulative density function
  - `qplindLnorm` for the quantile function
  - `rplindLnorm` for random number generation

- Poisson-Lognormal Distribution

  - `dpLnorm` for the density function
  - `ppLnorm` for the cumulative density function
  - `qtpLnorm` for the quantile function
  - `rpLnorm` for random number generation

- Poisson-Weibull Distribution

  - `dpoisweibull` for the density function
  - `ppoisweibull` for the cumulative density function
  - `qpoisweibull` for the quantile function
  - `rpoisweibull` for random number generation

- Sichel Distribution

  - `dsichel` for the density function
  - `psichel` for the cumulative density function
  - `qsichel` for the quantile function
  - `rsichel` for random number generation

**Model Estimation Functions**

- `flexCountReg` is a general function for estimating any of the
  regression models
- `nbg` estimates negative binomial regression (NB-1, NB-2, or NB-P) and
  allows the overdispersion parameter to be specified as a function of
  predictors.
- `poisGE` estimates the Poisson-Generalized-Exponential regression
  model. It allows the scale parameter to be specified as a function of
  predictors.
- `poisInvGaus` estimates the Poisson-Inverse-Gaussian regression model.
- `poisLind` estimates the Poisson-Lindley regression model.
- `poisLindGamma` estimates the Poisson-Lindley-Gamma (i.e., Negative
  Binomial-Lindley) regression model.
- `poisLindLnorm` estimates the Poisson-Lindley-Lognormal regression
  model.
- `poisLogn` estimates the Poisson-Lognormal regression model. It allows
  the standard deviation parameter ($\sigma$) to be specified as a
  function of predictors.
- `pwiebreg` estimates the Poisson-Weibull regression model. It allows
  the shape and scale parameters to be specified as functions of
  predictors.
- `rpnb` estimates the random parameters negative binomial regression
  (NB-1, NB-2, or NB-P).
- `sichel` estimates the Sichel regression model. It allows the scale
  parameter to be specified as a function of predictors.

**Model Evaluation, Comparison, and Convenience Functions**

- `cureplot` generates a CURE plot for the specified model, based on the
  [cureplots package](https://gbasulto.github.io/cureplots/).
- `mae` computes the Mean Absolute Error (MAE).
- `myAIC` computes the Akaike Information Criterion (AIC) value.
- `myBIC` computes the Bayesian Information Criterion (BIC) value.
- `regCompTable` creates a publication-ready table comparing multiple
  models. This can include the regression estimate results, AIC, BIC,
  and Pseudo R-Square values.
- `regCompTest` compares any given model with a base model. This can be
  used to perform a likelihood ratio test between models.
- `rmse` computes the Root Mean Squared Error (RMSE).
- `predict` allows the predict function to be used for out-of-sample
  predictions for any of the flexCountReg models.
- `summary` allows the use of the summary function to get a model
  summary from a flexCountReg regression object.

**Data** A dataset, `washington_roads`, is included. It is based on a
sample of Washington primary 2-lane roads from the years 2016-2018. Data
for the roads, traffic volumes (AADT) and associated crashes were
obtained from the [Highway Safety Information System
(HSIS)](https://highways.dot.gov/research/safety/hsis).

### Probability Distributions

As noted in the list of functions, the probability distributions below
are included in the flexCountReg package. Details of the distributions
are provided in the documentation (help files).

**Continuous Distributions**

- Inverse Gamma Distribution
- Triangle Distribution

**Count Distributions**

- Generalized Waring Distribution
- Negative Binomial in various forms (NB-1, NB-2, and NB-P)
- Poisson-Generalized-Exponential Distribution
- Poisson-Inverse-Gaussian Distribution
- Poisson-Lindley Distribution
- Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution
- Poisson-Lindley-Lognormal Distribution
- Poisson-Lognormal Distribution
- Poisson-Weibull Distribution
- Sichel Distribution

## Example

The following is an example of using flexCountReg to estimate a negative
binomial (NB-2) regression model with the overdispersion parameter as a
function of predictor variables:

``` r
library(gt) # used to format summary tables here
library(flexCountReg)
#> Registered S3 method overwritten by 'flexCountReg':
#>   method         from  
#>   summary.maxLik maxLik
```

``` r
library(knitr)

data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
gen.nb2 <- flexCountReg(Total_crashes ~ lnaadt + lnlength + speed50 +
                                ShouldWidth04 + AADTover10k,
                                ln.alpha.formula = ~ 1+lnlength,
                                data=washington_roads,
                                dist="NB2",
                                method = 'NM')
```

``` r
kable(summary(gen.nb2), caption = "NB-2 Model Summary")
#> Call:
#>  Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 +      AADTover10k 
#> 
#>  Method:  Nelder-Mead maximization 
#> Iterations:  201 
#> Convergence:  iteration limit exceeded  
#> Log-likelihood:  -1073.026 
#> 
#> Parameter Estimates:
#> # A tibble: 8 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -8.72        0.049 -178.        0         -8.81      -8.62 
#> 2 lnaadt              1.06        0.006  185.        0          1.05       1.07 
#> 3 lnlength            0.934       0.042   22.2       0          0.852      1.02 
#> 4 speed50            -0.427       0.101   -4.21      0         -0.625     -0.228
#> 5 ShouldWidth04       0.285       0.07     4.07      0          0.147      0.422
#> 6 AADTover10k         0.747       0.116    6.42      0          0.519      0.975
#> 7 ln(alpha):  (Inte…  0.581       0.172    3.37      0.001      0.243      0.919
#> 8 ln(alpha):  lnlen…  0.058       0.151    0.383     0.701     -0.239      0.355
```

| parameter              |  coeff | Std. Err. |   t-stat | p-value | lower CI | upper CI |
|:-----------------------|-------:|----------:|---------:|--------:|---------:|---------:|
| (Intercept)            | -8.716 |     0.049 | -177.505 |   0.000 |   -8.812 |   -8.620 |
| lnaadt                 |  1.063 |     0.006 |  185.240 |   0.000 |    1.052 |    1.074 |
| lnlength               |  0.934 |     0.042 |   22.248 |   0.000 |    0.852 |    1.017 |
| speed50                | -0.427 |     0.101 |   -4.209 |   0.000 |   -0.625 |   -0.228 |
| ShouldWidth04          |  0.285 |     0.070 |    4.066 |   0.000 |    0.147 |    0.422 |
| AADTover10k            |  0.747 |     0.116 |    6.417 |   0.000 |    0.519 |    0.975 |
| ln(alpha): (Intercept) |  0.581 |     0.172 |    3.371 |   0.001 |    0.243 |    0.919 |
| ln(alpha): lnlength    |  0.058 |     0.151 |    0.383 |   0.701 |   -0.239 |    0.355 |

NB-2 Model Summary

``` r
teststats <- regCompTest(gen.nb2)
kable(teststats$statistics)
```

| Statistic             |     Model | BaseModel |
|:----------------------|----------:|----------:|
| AIC                   | 2162.0527 |  3049.659 |
| BIC                   | 2204.5638 |  3054.973 |
| LR Test Statistic     |  901.6065 |        NA |
| LR degrees of freedom |    7.0000 |        NA |
| LR p-value            |    0.0000 |        NA |
| McFadden’s Pseudo R^2 |    0.2958 |        NA |

Checking the CURE plot:

``` r
cureplot(gen.nb2, indvar  ="lnaadt")
#> Covariate: indvar_values
#> CURE data frame was provided. Its first column, lnaadt, will be used.
```

<img src="man/figures/README-cureplot-initial-1.png" width="100%" />

Modifying the model to fit better:

``` r
gen.nb2 <- flexCountReg(Total_crashes ~ lnaadt  + lnlength + speed50 +
                                ShouldWidth04 + AADTover10k + I(AADTover10k/lnaadt),
                                ln.alpha.formula = ~ 1+lnlength,
                                data=washington_roads,
                                dist="NB2",
                                method = 'NM')
kable(summary(gen.nb2), caption = "Modified NB-2 Model Summary")
#> Call:
#>  Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 +      AADTover10k + I(AADTover10k/lnaadt) 
#> 
#>  Method:  Nelder-Mead maximization 
#> Iterations:  202 
#> Convergence:  iteration limit exceeded  
#> Log-likelihood:  -1062.504 
#> 
#> Parameter Estimates:
#> # A tibble: 9 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -7.36        0.043 -173.        0         -7.45      -7.28 
#> 2 lnaadt              0.894       0.005  181.        0          0.884      0.904
#> 3 lnlength            0.862       0.037   23.4       0          0.79       0.935
#> 4 speed50            -0.409       0.092   -4.46      0         -0.589     -0.229
#> 5 ShouldWidth04       0.251       0.059    4.22      0          0.134      0.367
#> 6 AADTover10k         5.09        0.09    56.8       0          4.91       5.26 
#> 7 I(AADTover10k/ln… -41.2         0.856  -48.1       0        -42.9      -39.5  
#> 8 ln(alpha):  (Int…   1.56        0.35     4.46      0          0.873      2.24 
#> 9 ln(alpha):  lnle…   0.022       0.345    0.065     0.948     -0.653      0.698
```

| parameter              |   coeff | Std. Err. |   t-stat | p-value | lower CI | upper CI |
|:-----------------------|--------:|----------:|---------:|--------:|---------:|---------:|
| (Intercept)            |  -7.364 |     0.043 | -173.202 |   0.000 |   -7.447 |   -7.281 |
| lnaadt                 |   0.894 |     0.005 |  181.133 |   0.000 |    0.884 |    0.904 |
| lnlength               |   0.862 |     0.037 |   23.373 |   0.000 |    0.790 |    0.935 |
| speed50                |  -0.409 |     0.092 |   -4.458 |   0.000 |   -0.589 |   -0.229 |
| ShouldWidth04          |   0.251 |     0.059 |    4.224 |   0.000 |    0.134 |    0.367 |
| AADTover10k            |   5.086 |     0.090 |   56.822 |   0.000 |    4.910 |    5.261 |
| I(AADTover10k/lnaadt)  | -41.179 |     0.856 |  -48.098 |   0.000 |  -42.857 |  -39.501 |
| ln(alpha): (Intercept) |   1.558 |     0.350 |    4.459 |   0.000 |    0.873 |    2.244 |
| ln(alpha): lnlength    |   0.022 |     0.345 |    0.065 |   0.948 |   -0.653 |    0.698 |

Modified NB-2 Model Summary

``` r
teststats <- regCompTest(gen.nb2)
kable(teststats$statistics)
```

| Statistic             |     Model | BaseModel |
|:----------------------|----------:|----------:|
| AIC                   | 2143.0082 |  3049.659 |
| BIC                   | 2190.8331 |  3054.973 |
| LR Test Statistic     |  922.6510 |        NA |
| LR degrees of freedom |    8.0000 |        NA |
| LR p-value            |    0.0000 |        NA |
| McFadden’s Pseudo R^2 |    0.3027 |        NA |

``` r
cureplot(gen.nb2, indvar  ="lnaadt")
#> Covariate: indvar_values
#> CURE data frame was provided. Its first column, lnaadt, will be used.
```

<img src="man/figures/README-cureplot-updt-1.png" width="100%" />

Estimating another model (NB-P) - without the interaction:

``` r
gen.nbp <- flexCountReg(Total_crashes ~ lnaadt + lnlength + speed50 +
                                ShouldWidth04 + AADTover10k,
                                ln.alpha.formula = ~ 1+lnlength,
                                data=washington_roads,
                                dist="NBP",
                                method = 'NM')
kable(summary(gen.nbp), caption = "NB-P Model Summary")
#> Call:
#>  Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 +      AADTover10k 
#> 
#>  Method:  Nelder-Mead maximization 
#> Iterations:  202 
#> Convergence:  iteration limit exceeded  
#> Log-likelihood:  -1062.206 
#> 
#> Parameter Estimates:
#> # A tibble: 9 × 7
#>   parameter           coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>               <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)        -7.73        0.043 -180.        0         -7.81      -7.64 
#> 2 lnaadt              0.935       0.005  189.        0          0.925      0.944
#> 3 lnlength            0.839       0.037   22.5       0          0.766      0.913
#> 4 speed50            -0.39        0.093   -4.18      0         -0.573     -0.207
#> 5 ShouldWidth04       0.255       0.059    4.29      0          0.139      0.372
#> 6 AADTover10k         0.693       0.088    7.92      0          0.522      0.865
#> 7 ln(alpha):  (Inte… -1.45        0.297   -4.88      0         -2.03      -0.865
#> 8 ln(alpha):  lnlen… -0.116       0.254   -0.458     0.647     -0.615      0.382
#> 9 P                   1.69        0.291    5.80      0          1.12       2.26
```

| parameter              |  coeff | Std. Err. |   t-stat | p-value | lower CI | upper CI |
|:-----------------------|-------:|----------:|---------:|--------:|---------:|---------:|
| (Intercept)            | -7.727 |     0.043 | -180.434 |   0.000 |   -7.811 |   -7.643 |
| lnaadt                 |  0.935 |     0.005 |  188.818 |   0.000 |    0.925 |    0.944 |
| lnlength               |  0.839 |     0.037 |   22.456 |   0.000 |    0.766 |    0.913 |
| speed50                | -0.390 |     0.093 |   -4.183 |   0.000 |   -0.573 |   -0.207 |
| ShouldWidth04          |  0.255 |     0.059 |    4.293 |   0.000 |    0.139 |    0.372 |
| AADTover10k            |  0.693 |     0.088 |    7.918 |   0.000 |    0.522 |    0.865 |
| ln(alpha): (Intercept) | -1.447 |     0.297 |   -4.878 |   0.000 |   -2.028 |   -0.865 |
| ln(alpha): lnlength    | -0.116 |     0.254 |   -0.458 |   0.647 |   -0.615 |    0.382 |
| P                      |  1.686 |     0.291 |    5.796 |   0.000 |    1.116 |    2.256 |

NB-P Model Summary

``` r
teststats <- regCompTest(gen.nbp)
kable(teststats$statistics)
```

| Statistic             |     Model | BaseModel |
|:----------------------|----------:|----------:|
| AIC                   | 2142.4125 |  3049.659 |
| BIC                   | 2190.2375 |  3054.973 |
| LR Test Statistic     |  923.2467 |        NA |
| LR degrees of freedom |    8.0000 |        NA |
| LR p-value            |    0.0000 |        NA |
| McFadden’s Pseudo R^2 |    0.3029 |        NA |

Checking the CURE plot (notice that the CURE plot is MUCH better in this
case than the NB-2 without the interaction and still better than the
modified NB-2):

``` r
cureplot(gen.nbp, indvar  ="lnaadt")
#> Covariate: indvar_values
#> CURE data frame was provided. Its first column, lnaadt, will be used.
```

<img src="man/figures/README-cureplot-lnaadt-1.png" width="100%" />

Creating a table to compare the models:

``` r
regCompTable(list("Generalized NB-2"=gen.nb2, "Generalized NB-P"=gen.nbp), tableType="tibble") |> 
  kable()
```

| Parameter              | Generalized NB-2      | Generalized NB-P     |
|:-----------------------|:----------------------|:---------------------|
| (Intercept)            | -7.364 (0.043)\*\*\*  | -7.727 (0.043)\*\*\* |
| lnaadt                 | 0.894 (0.005)\*\*\*   | 0.935 (0.005)\*\*\*  |
| lnlength               | 0.862 (0.037)\*\*\*   | 0.839 (0.037)\*\*\*  |
| speed50                | -0.409 (0.092)\*\*\*  | -0.39 (0.093)\*\*\*  |
| ShouldWidth04          | 0.251 (0.059)\*\*\*   | 0.255 (0.059)\*\*\*  |
| AADTover10k            | 5.086 (0.09)\*\*\*    | 0.693 (0.088)\*\*\*  |
| I(AADTover10k/lnaadt)  | -41.179 (0.856)\*\*\* | —                    |
| ln(alpha): (Intercept) | 1.558 (0.35)\*\*\*    | -1.447 (0.297)\*\*\* |
| ln(alpha): lnlength    | 0.022 (0.345)         | -0.116 (0.254)       |
| P                      | —                     | 1.686 (0.291)\*\*\*  |
| N Obs.                 | 1501                  | 1501                 |
| LL                     | -1062.504             | -1062.206            |
| AIC                    | 2143.008              | 2142.412             |
| BIC                    | 2190.833              | 2190.237             |
| Pseudo-R-Sq.           | 0.303                 | 0.303                |

Note that the metrics for comparison are similar. While the models both
have the same number of parameters, the NB-P was able to get better
performance without requiring the interaction terms (which leads to
strange relationships between the exposure metric and the outcome).
