


<!-- README.md is generated from README.Rmd. Please edit that file -->

# flexCountReg <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![codecov](https://codecov.io/gh/jwood-iastate/flexCountReg/graph/badge.svg?token=BX2FJQNPK2)](https://codecov.io/gh/jwood-iastate/flexCountReg)
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

- Conway-Maxwell-Poisson Distribution

  - `dcom` for the density function
  - `pcom` for the cumulative density function
  - `qcom` for the quantile function
  - `rcom` for random number generation

**Model Estimation Functions**

- `countreg` is a general function for estimating the non-panel,
  non-random parameters count regression models
- `rpnb` estimates the random parameters negative binomial regression
  (NB-1, NB-2, or NB-P).
- `rppLind` estimates the random parameters Poisson-Lindley regression
  model.
- `poislind.re` estimates the random effects Poisson-Lindley model
- `renb` estimates the random effects negative binomial regression
  model.

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
- Conway-Maxwell-Poisson (COM) Distribution

## Example

The following is an example of using flexCountReg to estimate a negative
binomial (NB-2) regression model with the overdispersion parameter as a
function of predictor variables:

``` r
library(gt) # used to format summary tables here
library(flexCountReg)
library(knitr)

data("washington_roads")
washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
gen.nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
               data = washington_roads, family = "NB2",
               dis_param_formula_1 = ~ speed50,  method='BFGS')
```

``` r
kable(summary(gen.nb2), caption = "NB-2 Model Summary")
#> Call:
#>  Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus 
#> 
#>  Method:  BFGS maximization 
#> Iterations:  53 
#> Convergence:  successful convergence  
#> Log-likelihood:  -1064.876 
#> 
#> Parameter Estimates:
#> # A tibble: 7 × 7
#>   parameter              coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>                  <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)           -7.40        0.043  -171.       0         -7.49      -7.32 
#> 2 lnaadt                 0.912       0.005   182.       0          0.902      0.921
#> 3 lnlength               0.843       0.037    22.9      0          0.771      0.915
#> 4 speed50               -0.47        0.102    -4.62     0         -0.669     -0.27 
#> 5 AADT10kplus            0.77        0.089     8.61     0          0.594      0.945
#> 6 ln(alpha):(Intercept)  1.62        0.294     5.51     0          1.04       2.19 
#> 7 ln(alpha):speed50     -1.31        0.469    -2.78     0.005     -2.22      -0.387
```

| parameter             |  coeff | Std. Err. |   t-stat | p-value | lower CI | upper CI |
|:----------------------|-------:|----------:|---------:|--------:|---------:|---------:|
| (Intercept)           | -7.401 |     0.043 | -171.489 |   0.000 |   -7.486 |   -7.317 |
| lnaadt                |  0.912 |     0.005 |  182.454 |   0.000 |    0.902 |    0.921 |
| lnlength              |  0.843 |     0.037 |   22.878 |   0.000 |    0.771 |    0.915 |
| speed50               | -0.470 |     0.102 |   -4.619 |   0.000 |   -0.669 |   -0.270 |
| AADT10kplus           |  0.770 |     0.089 |    8.607 |   0.000 |    0.594 |    0.945 |
| ln(alpha):(Intercept) |  1.619 |     0.294 |    5.514 |   0.000 |    1.044 |    2.195 |
| ln(alpha):speed50     | -1.306 |     0.469 |   -2.785 |   0.005 |   -2.225 |   -0.387 |

NB-2 Model Summary

``` r
teststats <- regCompTest(gen.nb2)
kable(teststats$statistics)
```

| Statistic             |     Model | BaseModel |
|:----------------------|----------:|----------:|
| AIC                   | 2143.7522 |  3049.659 |
| BIC                   | 2180.9494 |  3054.973 |
| LR Test Statistic     |  917.9070 |        NA |
| LR degrees of freedom |    6.0000 |        NA |
| LR p-value            |    0.0000 |        NA |
| McFadden’s Pseudo R^2 |    0.3012 |        NA |

Checking the CURE plot:

``` r
cureplot(gen.nb2, indvar  ="lnaadt")
#>    [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>   [77] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [153] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [229] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [305] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [381] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [457] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [533] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [609] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [685] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [761] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [837] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [913] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [989] 1 1 1 1 1 1 1 1 1 1 1 1
#>  [ reached 'max' / getOption("max.print") -- omitted 501 entries ]
#> Covariate: indvar_values
#> CURE data frame was provided. Its first column, lnaadt, will be used.
```

<img src="man/figures/README-cureplot-initial-1.png" width="100%" />

Modifying the model to fit better:

``` r


gen.nb2 <- countreg(Total_crashes ~ lnaadt  + lnlength + speed50 +
                                ShouldWidth04 + AADT10kplus + 
                                I(AADT10kplus/lnaadt),
                                data = washington_roads, family = "NB2",
                                dis_param_formula_1 = ~ lnlength,  method='BFGS')

kable(summary(gen.nb2), caption = "Modified NB-2 Model Summary")
#> Call:
#>  Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 +      AADT10kplus + I(AADT10kplus/lnaadt) 
#> 
#>  Method:  BFGS maximization 
#> Iterations:  65 
#> Convergence:  successful convergence  
#> Log-likelihood:  -1061.935 
#> 
#> Parameter Estimates:
#> # A tibble: 9 × 7
#>   parameter               coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>                   <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)            -7.70        0.043  -180.       0         -7.78      -7.61 
#> 2 lnaadt                  0.932       0.005   188.       0          0.922      0.942
#> 3 lnlength                0.852       0.038    22.6      0          0.778      0.926
#> 4 speed50                -0.4         0.091    -4.39     0         -0.579     -0.222
#> 5 ShouldWidth04           0.26        0.06      4.34     0          0.142      0.377
#> 6 AADT10kplus             5.12        0.092    55.6      0          4.94       5.30 
#> 7 I(AADT10kplus/lnaadt) -42.1         0.938   -44.9      0        -44.0      -40.3  
#> 8 ln(alpha):(Intercept)   1.90        0.328     5.79     0          1.25       2.54 
#> 9 ln(alpha):lnlength      0.417       0.242     1.72     0.085     -0.057      0.892
```

| parameter             |   coeff | Std. Err. |   t-stat | p-value | lower CI | upper CI |
|:----------------------|--------:|----------:|---------:|--------:|---------:|---------:|
| (Intercept)           |  -7.695 |     0.043 | -180.089 |   0.000 |   -7.779 |   -7.611 |
| lnaadt                |   0.932 |     0.005 |  188.011 |   0.000 |    0.922 |    0.942 |
| lnlength              |   0.852 |     0.038 |   22.569 |   0.000 |    0.778 |    0.926 |
| speed50               |  -0.400 |     0.091 |   -4.388 |   0.000 |   -0.579 |   -0.222 |
| ShouldWidth04         |   0.260 |     0.060 |    4.336 |   0.000 |    0.142 |    0.377 |
| AADT10kplus           |   5.123 |     0.092 |   55.599 |   0.000 |    4.942 |    5.304 |
| I(AADT10kplus/lnaadt) | -42.146 |     0.938 |  -44.937 |   0.000 |  -43.984 |  -40.307 |
| ln(alpha):(Intercept) |   1.895 |     0.328 |    5.787 |   0.000 |    1.253 |    2.537 |
| ln(alpha):lnlength    |   0.417 |     0.242 |    1.723 |   0.085 |   -0.057 |    0.892 |

Modified NB-2 Model Summary

``` r
teststats <- regCompTest(gen.nb2)
kable(teststats$statistics)
```

| Statistic             |     Model | BaseModel |
|:----------------------|----------:|----------:|
| AIC                   | 2141.8694 |  3049.659 |
| BIC                   | 2189.6944 |  3054.973 |
| LR Test Statistic     |  923.7898 |        NA |
| LR degrees of freedom |    8.0000 |        NA |
| LR p-value            |    0.0000 |        NA |
| McFadden’s Pseudo R^2 |    0.3031 |        NA |

``` r
cureplot(gen.nb2, indvar  ="lnaadt")
#>    [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>   [77] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [153] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [229] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [305] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [381] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [457] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [533] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [609] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [685] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [761] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [837] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [913] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [989] 1 1 1 1 1 1 1 1 1 1 1 1
#>  [ reached 'max' / getOption("max.print") -- omitted 501 entries ]
#> Covariate: indvar_values
#> CURE data frame was provided. Its first column, lnaadt, will be used.
```

<img src="man/figures/README-cureplot-updt-1.png" width="100%" />

Estimating another model (NB-P) - without the interaction:

``` r
gen.nbp <- countreg(Total_crashes ~ lnaadt  + lnlength + speed50 +
                                ShouldWidth04 + AADT10kplus,
                                data = washington_roads, family = "NBp",
                                dis_param_formula_1 = ~ lnlength,  method='BFGS')
kable(summary(gen.nbp), caption = "NB-P Model Summary")
#> Call:
#>  Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 +      AADT10kplus 
#> 
#>  Method:  BFGS maximization 
#> Iterations:  65 
#> Convergence:  successful convergence  
#> Log-likelihood:  -1062.195 
#> 
#> Parameter Estimates:
#> # A tibble: 9 × 7
#>   parameter              coeff `Std. Err.` `t-stat` `p-value` `lower CI` `upper CI`
#>   <chr>                  <dbl>       <dbl>    <dbl>     <dbl>      <dbl>      <dbl>
#> 1 (Intercept)           -7.76        0.043 -181.        0         -7.85      -7.68 
#> 2 lnaadt                 0.938       0.005  189.        0          0.928      0.948
#> 3 lnlength               0.836       0.037   22.3       0          0.763      0.91 
#> 4 speed50               -0.384       0.093   -4.13      0         -0.567     -0.202
#> 5 ShouldWidth04          0.258       0.06     4.33      0          0.141      0.374
#> 6 AADT10kplus            0.689       0.088    7.85      0          0.517      0.861
#> 7 ln(alpha):(Intercept) -1.50        0.297   -5.05      0         -2.08      -0.915
#> 8 ln(alpha):lnlength    -0.168       0.244   -0.688     0.491     -0.646      0.31 
#> 9 ln(p)                  0.525       0.173    3.03      0.002      0.186      0.864
```

| parameter             |  coeff | Std. Err. |   t-stat | p-value | lower CI | upper CI |
|:----------------------|-------:|----------:|---------:|--------:|---------:|---------:|
| (Intercept)           | -7.764 |     0.043 | -181.133 |   0.000 |   -7.848 |   -7.680 |
| lnaadt                |  0.938 |     0.005 |  189.457 |   0.000 |    0.928 |    0.948 |
| lnlength              |  0.836 |     0.037 |   22.307 |   0.000 |    0.763 |    0.910 |
| speed50               | -0.384 |     0.093 |   -4.126 |   0.000 |   -0.567 |   -0.202 |
| ShouldWidth04         |  0.258 |     0.060 |    4.333 |   0.000 |    0.141 |    0.374 |
| AADT10kplus           |  0.689 |     0.088 |    7.854 |   0.000 |    0.517 |    0.861 |
| ln(alpha):(Intercept) | -1.496 |     0.297 |   -5.046 |   0.000 |   -2.078 |   -0.915 |
| ln(alpha):lnlength    | -0.168 |     0.244 |   -0.688 |   0.491 |   -0.646 |    0.310 |
| ln(p)                 |  0.525 |     0.173 |    3.034 |   0.002 |    0.186 |    0.864 |

NB-P Model Summary

``` r
teststats <- regCompTest(gen.nbp)
kable(teststats$statistics)
```

| Statistic             |     Model | BaseModel |
|:----------------------|----------:|----------:|
| AIC                   | 2142.3895 |  3049.659 |
| BIC                   | 2190.2144 |  3054.973 |
| LR Test Statistic     |  923.2697 |        NA |
| LR degrees of freedom |    8.0000 |        NA |
| LR p-value            |    0.0000 |        NA |
| McFadden’s Pseudo R^2 |    0.3029 |        NA |

Checking the CURE plot (notice that the CURE plot is MUCH better in this
case than the NB-2 without the interaction and still better than the
modified NB-2):

``` r
cureplot(gen.nbp, indvar  ="lnaadt")
#>    [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>   [77] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [153] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [229] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [305] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [381] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [457] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [533] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [609] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [685] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [761] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [837] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [913] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [989] 1 1 1 1 1 1 1 1 1 1 1 1
#>  [ reached 'max' / getOption("max.print") -- omitted 501 entries ]
#> Covariate: indvar_values
#> CURE data frame was provided. Its first column, lnaadt, will be used.
```

<img src="man/figures/README-cureplot-lnaadt-1.png" width="100%" />

Creating a table to compare the models:

``` r
regCompTable(list("Generalized NB-2"=gen.nb2, "Generalized NB-P"=gen.nbp), tableType="tibble") |> 
  kable()
```

| Parameter             | Generalized NB-2      | Generalized NB-P     |
|:----------------------|:----------------------|:---------------------|
| (Intercept)           | -7.695 (0.043)\*\*\*  | -7.764 (0.043)\*\*\* |
| lnaadt                | 0.932 (0.005)\*\*\*   | 0.938 (0.005)\*\*\*  |
| lnlength              | 0.852 (0.038)\*\*\*   | 0.836 (0.037)\*\*\*  |
| speed50               | -0.4 (0.091)\*\*\*    | -0.384 (0.093)\*\*\* |
| ShouldWidth04         | 0.26 (0.06)\*\*\*     | 0.258 (0.06)\*\*\*   |
| AADT10kplus           | 5.123 (0.092)\*\*\*   | 0.689 (0.088)\*\*\*  |
| I(AADT10kplus/lnaadt) | -42.146 (0.938)\*\*\* | —                    |
| ln(alpha):(Intercept) | 1.895 (0.328)\*\*\*   | -1.496 (0.297)\*\*\* |
| ln(alpha):lnlength    | 0.417 (0.242)         | -0.168 (0.244)       |
| ln(p)                 | —                     | 0.525 (0.173)\*\*    |
| N Obs.                | 1501                  | 1501                 |
| LL                    | -1061.935             | -1062.195            |
| AIC                   | 2141.869              | 2142.389             |
| BIC                   | 2189.694              | 2190.214             |
| Pseudo-R-Sq.          | 0.303                 | 0.303                |

Note that the metrics for comparison are similar. While the models both
have the same number of parameters, the NB-P was able to get better
performance without requiring the interaction terms (which leads to
strange relationships between the exposure metric and the outcome).
