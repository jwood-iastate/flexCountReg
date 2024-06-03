
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flexCountReg

<!-- badges: start -->
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

## Functions

| Function                        | Category                 | Description                                                                                                                        |
|---------------------------------|--------------------------|------------------------------------------------------------------------------------------------------------------------------------|
| dinvgamma                       | Probability Distribution | Inverse Gamma Distribution                                                                                                         |
| dpge                            | Probability Distribution | Poisson-Generalized-Exponential Distribution                                                                                       |
| dpinvgaus                       | Probability Distribution | Poisson-Inverse-Gaussian Distribution                                                                                              |
| dplind                          | Probability Distribution | Poisson-Lindley Distribution                                                                                                       |
| dplindGamma                     | Probability Distribution | Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution                                                                     |
| dplindLnorm                     | Probability Distribution | Poisson-Lindley-Lognormal Distribution                                                                                             |
| dpLnorm                         | Probability Distribution | Poisson-Lognormal Distribution                                                                                                     |
| dsichel                         | Probability Distribution | Sichel Distribution                                                                                                                |
| dtri                            | Probability Distribution | Triangle Distribution                                                                                                              |
| flexCountReg                    | Regression Fit           | Estimate flexible count regression models with and without random parameters                                                       |
| flexCountReg.predict            |                          | Function for generating predictions based on the random parameters negative binomial with multiple optional methods                |
| genWaring                       |                          | Function for estimating a Generalized Waring regression model                                                                      |
| invgamma                        |                          | Inverse Gamma Distribution                                                                                                         |
| mgf_lognormal                   |                          | Moment Generating Function for a Lognormal Distribution                                                                            |
| mle.predict                     |                          | Predictions for non-random parameters count models                                                                                 |
| nbg                             |                          | Function for estimating a variety of negative binomial models (NB-1, NB-2, NB-P and generalized versions of each)                  |
| Negative-Binomial-Lindley       |                          | Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution                                                                     |
| pinvgamma                       |                          | Inverse Gamma Distribution                                                                                                         |
| poisGE                          |                          | Poisson-Generalized-Exponential Regression                                                                                         |
| poisInvGaus                     |                          | Function for estimating a Poisson-Inverse-Gaussian regression model                                                                |
| poisLind                        |                          | Function for estimating a Poisson-Lindley regression model                                                                         |
| poisLindGamma                   |                          | Function for estimating a Poisson-Lindley-Gamma (i.e., Negative Binomial-Lindley) regression model                                 |
| poisLindLnorm                   |                          | Function for estimating a Poisson-Lindley-Lognormal regression model                                                               |
| poisLogn                        |                          | Poisson-Lognormal Regression                                                                                                       |
| Poisson-Generalized-Exponential |                          | Poisson-Generalized-Exponential Distribution                                                                                       |
| Poisson-Inverse-Gaussian        |                          | Poisson-Inverse-Gaussian Distribution                                                                                              |
| Poisson-Lindley                 |                          | Poisson-Lindley Distribution                                                                                                       |
| Poisson-Lindley-Lognormal       |                          | Poisson-Lindley-Lognormal Distribution                                                                                             |
| Poisson-Lognormal               |                          | Poisson-Lognormal Distribution                                                                                                     |
| ppge                            |                          | Poisson-Generalized-Exponential Distribution                                                                                       |
| ppinvgaus                       |                          | Poisson-Inverse-Gaussian Distribution                                                                                              |
| pplind                          |                          | Poisson-Lindley Distribution                                                                                                       |
| pplindGamma                     |                          | Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution                                                                     |
| pplindLnorm                     |                          | Poisson-Lindley-Lognormal Distribution                                                                                             |
| ppLnorm                         |                          | Poisson-Lognormal Distribution                                                                                                     |
| psichel                         |                          | Sichel Distribution                                                                                                                |
| ptri                            |                          | Triangle Distribution                                                                                                              |
| qinvgamma                       |                          | Inverse Gamma Distribution                                                                                                         |
| qpge                            |                          | Poisson-Generalized-Exponential Distribution                                                                                       |
| qpinvgaus                       |                          | Poisson-Inverse-Gaussian Distribution                                                                                              |
| qplind                          |                          | Poisson-Lindley Distribution                                                                                                       |
| qplindGamma                     |                          | Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution                                                                     |
| qplindLnorm                     |                          | Poisson-Lindley-Lognormal Distribution                                                                                             |
| qpLnorm                         |                          | Poisson-Lognormal Distribution                                                                                                     |
| qsichel                         |                          | Sichel Distribution                                                                                                                |
| qtri                            |                          | Triangle Distribution                                                                                                              |
| rinvgamma                       |                          | Inverse Gamma Distribution                                                                                                         |
| rpge                            |                          | Poisson-Generalized-Exponential Distribution                                                                                       |
| rpinvgaus                       |                          | Poisson-Inverse-Gaussian Distribution                                                                                              |
| rplind                          |                          | Poisson-Lindley Distribution                                                                                                       |
| rplindGamma                     |                          | Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution                                                                     |
| rplindLnorm                     |                          | Poisson-Lindley-Lognormal Distribution                                                                                             |
| rpLnorm                         |                          | Poisson-Lognormal Distribution                                                                                                     |
| rpnb                            |                          | Function for estimating a random parameter negative binomial with the ability to specify if the NB-1, NB-2, or NB-P should be used |
| rpnb.predict                    |                          | Function for generating predictions based on the random parameters negative binomial with multiple optional methods                |
| rsichel                         |                          | Sichel Distribution                                                                                                                |
| rtri                            |                          | Triangle Distribution                                                                                                              |
| sichel                          |                          | Function for estimating a Sichel regression model                                                                                  |
| Sichel-Distribution             |                          | Sichel Distribution                                                                                                                |
| Triangular                      |                          | Triangle Distribution                                                                                                              |
| washington_roads                |                          | Washington Road Crashes                                                                                                            |

## Example

The following is an example of using flexCountReg to estimate a negative
binomial regression model with the overdispersion parameter as a
function of predictor variables:

``` r
library(flexCountReg)

data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
gen.nb2 <- flexCountReg(Total_crashes ~ lnaadt + lnlength + speed50 +
                                ShouldWidth04 + AADTover10k,
                                ln.alpha.formula = ~ 1+lnlength,
                                data=washington_roads,
                                dist="NB2",
                                method = 'NM')
#> [1] "The Likelihood Ratio (LR) Test for H0: NB is No Better than Poisson"
#> [1] "LR =  -3.5182"
#> [1] "LR degrees of freedom =  2"
#> [1] "LR p-value =  1"
#> [1] "Macfadden's Pseudo R^2 =  0.2958"
```

``` r
summary(gen.nb2)
#> --------------------------------------------
#> Maximum Likelihood estimation
#> Nelder-Mead maximization, 201 iterations
#> Return code 1: iteration limit exceeded 
#> Log-Likelihood: -1073.026 
#> 8  free parameters
#> Estimates:
#>                         Estimate Std. error t value  Pr(> t)    
#> (Intercept)             -8.71622    0.58075 -15.009  < 2e-16 ***
#> lnaadt                   1.06305    0.06864  15.488  < 2e-16 ***
#> lnlength                 0.93446    0.07981  11.709  < 2e-16 ***
#> speed50                 -0.42663    0.12115  -3.521 0.000429 ***
#> ShouldWidth04            0.28456    0.10472   2.717 0.006581 ** 
#> AADTover10k              0.74651    0.16500   4.524 6.06e-06 ***
#> ln(alpha):  (Intercept)  0.58113    0.30326   1.916 0.055335 .  
#> ln(alpha):  lnlength     0.05803    0.26243   0.221 0.824979    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> --------------------------------------------
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
