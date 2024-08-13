


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
gt(summary(gen.nb2), caption = "NB-2 Model Summary")
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

<div id="ndrnpakokm" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ndrnpakokm table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#ndrnpakokm thead, #ndrnpakokm tbody, #ndrnpakokm tfoot, #ndrnpakokm tr, #ndrnpakokm td, #ndrnpakokm th {
  border-style: none;
}
&#10;#ndrnpakokm p {
  margin: 0;
  padding: 0;
}
&#10;#ndrnpakokm .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#ndrnpakokm .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#ndrnpakokm .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#ndrnpakokm .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#ndrnpakokm .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#ndrnpakokm .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#ndrnpakokm .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#ndrnpakokm .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#ndrnpakokm .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#ndrnpakokm .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#ndrnpakokm .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#ndrnpakokm .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#ndrnpakokm .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#ndrnpakokm .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#ndrnpakokm .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ndrnpakokm .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#ndrnpakokm .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#ndrnpakokm .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#ndrnpakokm .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ndrnpakokm .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#ndrnpakokm .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ndrnpakokm .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#ndrnpakokm .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ndrnpakokm .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#ndrnpakokm .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ndrnpakokm .gt_left {
  text-align: left;
}
&#10;#ndrnpakokm .gt_center {
  text-align: center;
}
&#10;#ndrnpakokm .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#ndrnpakokm .gt_font_normal {
  font-weight: normal;
}
&#10;#ndrnpakokm .gt_font_bold {
  font-weight: bold;
}
&#10;#ndrnpakokm .gt_font_italic {
  font-style: italic;
}
&#10;#ndrnpakokm .gt_super {
  font-size: 65%;
}
&#10;#ndrnpakokm .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#ndrnpakokm .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#ndrnpakokm .gt_indent_1 {
  text-indent: 5px;
}
&#10;#ndrnpakokm .gt_indent_2 {
  text-indent: 10px;
}
&#10;#ndrnpakokm .gt_indent_3 {
  text-indent: 15px;
}
&#10;#ndrnpakokm .gt_indent_4 {
  text-indent: 20px;
}
&#10;#ndrnpakokm .gt_indent_5 {
  text-indent: 25px;
}
&#10;#ndrnpakokm .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#ndrnpakokm div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <caption>NB-2 Model Summary</caption>
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="parameter">parameter</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="coeff">coeff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Std. Err.">Std. Err.</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="t-stat">t-stat</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="p-value">p-value</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lower CI">lower CI</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="upper CI">upper CI</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="parameter" class="gt_row gt_left">(Intercept)</td>
<td headers="coeff" class="gt_row gt_right">-8.716</td>
<td headers="Std. Err." class="gt_row gt_right">0.049</td>
<td headers="t-stat" class="gt_row gt_right">-177.505</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">-8.812</td>
<td headers="upper CI" class="gt_row gt_right">-8.620</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">lnaadt</td>
<td headers="coeff" class="gt_row gt_right">1.063</td>
<td headers="Std. Err." class="gt_row gt_right">0.006</td>
<td headers="t-stat" class="gt_row gt_right">185.240</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">1.052</td>
<td headers="upper CI" class="gt_row gt_right">1.074</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">lnlength</td>
<td headers="coeff" class="gt_row gt_right">0.934</td>
<td headers="Std. Err." class="gt_row gt_right">0.042</td>
<td headers="t-stat" class="gt_row gt_right">22.248</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.852</td>
<td headers="upper CI" class="gt_row gt_right">1.017</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">speed50</td>
<td headers="coeff" class="gt_row gt_right">-0.427</td>
<td headers="Std. Err." class="gt_row gt_right">0.101</td>
<td headers="t-stat" class="gt_row gt_right">-4.209</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">-0.625</td>
<td headers="upper CI" class="gt_row gt_right">-0.228</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ShouldWidth04</td>
<td headers="coeff" class="gt_row gt_right">0.285</td>
<td headers="Std. Err." class="gt_row gt_right">0.070</td>
<td headers="t-stat" class="gt_row gt_right">4.066</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.147</td>
<td headers="upper CI" class="gt_row gt_right">0.422</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">AADTover10k</td>
<td headers="coeff" class="gt_row gt_right">0.747</td>
<td headers="Std. Err." class="gt_row gt_right">0.116</td>
<td headers="t-stat" class="gt_row gt_right">6.417</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.519</td>
<td headers="upper CI" class="gt_row gt_right">0.975</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ln(alpha):  (Intercept)</td>
<td headers="coeff" class="gt_row gt_right">0.581</td>
<td headers="Std. Err." class="gt_row gt_right">0.172</td>
<td headers="t-stat" class="gt_row gt_right">3.371</td>
<td headers="p-value" class="gt_row gt_right">0.001</td>
<td headers="lower CI" class="gt_row gt_right">0.243</td>
<td headers="upper CI" class="gt_row gt_right">0.919</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ln(alpha):  lnlength</td>
<td headers="coeff" class="gt_row gt_right">0.058</td>
<td headers="Std. Err." class="gt_row gt_right">0.151</td>
<td headers="t-stat" class="gt_row gt_right">0.383</td>
<td headers="p-value" class="gt_row gt_right">0.701</td>
<td headers="lower CI" class="gt_row gt_right">-0.239</td>
<td headers="upper CI" class="gt_row gt_right">0.355</td></tr>
  </tbody>
  &#10;  
</table>
</div>

``` r
tesstats <- regCompTest(gen.nb2)
tesstats$gtTable
```

<div id="mobiljqngs" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#mobiljqngs table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#mobiljqngs thead, #mobiljqngs tbody, #mobiljqngs tfoot, #mobiljqngs tr, #mobiljqngs td, #mobiljqngs th {
  border-style: none;
}
&#10;#mobiljqngs p {
  margin: 0;
  padding: 0;
}
&#10;#mobiljqngs .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#mobiljqngs .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#mobiljqngs .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#mobiljqngs .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#mobiljqngs .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#mobiljqngs .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#mobiljqngs .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#mobiljqngs .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#mobiljqngs .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#mobiljqngs .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#mobiljqngs .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#mobiljqngs .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#mobiljqngs .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#mobiljqngs .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#mobiljqngs .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#mobiljqngs .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#mobiljqngs .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#mobiljqngs .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#mobiljqngs .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#mobiljqngs .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#mobiljqngs .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#mobiljqngs .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#mobiljqngs .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#mobiljqngs .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#mobiljqngs .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#mobiljqngs .gt_left {
  text-align: left;
}
&#10;#mobiljqngs .gt_center {
  text-align: center;
}
&#10;#mobiljqngs .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#mobiljqngs .gt_font_normal {
  font-weight: normal;
}
&#10;#mobiljqngs .gt_font_bold {
  font-weight: bold;
}
&#10;#mobiljqngs .gt_font_italic {
  font-style: italic;
}
&#10;#mobiljqngs .gt_super {
  font-size: 65%;
}
&#10;#mobiljqngs .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#mobiljqngs .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#mobiljqngs .gt_indent_1 {
  text-indent: 5px;
}
&#10;#mobiljqngs .gt_indent_2 {
  text-indent: 10px;
}
&#10;#mobiljqngs .gt_indent_3 {
  text-indent: 15px;
}
&#10;#mobiljqngs .gt_indent_4 {
  text-indent: 20px;
}
&#10;#mobiljqngs .gt_indent_5 {
  text-indent: 25px;
}
&#10;#mobiljqngs .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#mobiljqngs div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Model Comparison Statistics</td>
    </tr>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Statistic">Statistic</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Model">Model</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="BaseModel">BaseModel</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
<td headers="Model" class="gt_row gt_right">2,162.0527</td>
<td headers="BaseModel" class="gt_row gt_right">3,049.6592</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
<td headers="Model" class="gt_row gt_right">2,204.5638</td>
<td headers="BaseModel" class="gt_row gt_right">3,054.9731</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
<td headers="Model" class="gt_row gt_right">901.6065</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
<td headers="Model" class="gt_row gt_right">7.0000</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
<td headers="Model" class="gt_row gt_right">0.0000</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
<td headers="Model" class="gt_row gt_right">0.2958</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
  </tbody>
  &#10;  
</table>
</div>

Checking the CURE plot:

``` r
cureplot(gen.nb2, indvar  ="lnaadt")
#> [1] 1501   14
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
gt(summary(gen.nb2), caption = "Modified NB-2 Model Summary")
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
#> 4 speed50            -0.409       0.092   -4.46      0         -0.588     -0.229
#> 5 ShouldWidth04       0.251       0.059    4.23      0          0.134      0.367
#> 6 AADTover10k         5.09        0.09    56.8       0          4.91       5.26 
#> 7 I(AADTover10k/ln… -41.2         0.938  -43.9       0        -43.0      -39.3  
#> 8 ln(alpha):  (Int…   1.56        0.34     4.58      0          0.892      2.22 
#> 9 ln(alpha):  lnle…   0.022       0.35     0.064     0.949     -0.663      0.708
```

<div id="twqidtaqbo" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#twqidtaqbo table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#twqidtaqbo thead, #twqidtaqbo tbody, #twqidtaqbo tfoot, #twqidtaqbo tr, #twqidtaqbo td, #twqidtaqbo th {
  border-style: none;
}
&#10;#twqidtaqbo p {
  margin: 0;
  padding: 0;
}
&#10;#twqidtaqbo .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#twqidtaqbo .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#twqidtaqbo .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#twqidtaqbo .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#twqidtaqbo .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#twqidtaqbo .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#twqidtaqbo .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#twqidtaqbo .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#twqidtaqbo .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#twqidtaqbo .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#twqidtaqbo .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#twqidtaqbo .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#twqidtaqbo .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#twqidtaqbo .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#twqidtaqbo .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#twqidtaqbo .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#twqidtaqbo .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#twqidtaqbo .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#twqidtaqbo .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#twqidtaqbo .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#twqidtaqbo .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#twqidtaqbo .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#twqidtaqbo .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#twqidtaqbo .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#twqidtaqbo .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#twqidtaqbo .gt_left {
  text-align: left;
}
&#10;#twqidtaqbo .gt_center {
  text-align: center;
}
&#10;#twqidtaqbo .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#twqidtaqbo .gt_font_normal {
  font-weight: normal;
}
&#10;#twqidtaqbo .gt_font_bold {
  font-weight: bold;
}
&#10;#twqidtaqbo .gt_font_italic {
  font-style: italic;
}
&#10;#twqidtaqbo .gt_super {
  font-size: 65%;
}
&#10;#twqidtaqbo .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#twqidtaqbo .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#twqidtaqbo .gt_indent_1 {
  text-indent: 5px;
}
&#10;#twqidtaqbo .gt_indent_2 {
  text-indent: 10px;
}
&#10;#twqidtaqbo .gt_indent_3 {
  text-indent: 15px;
}
&#10;#twqidtaqbo .gt_indent_4 {
  text-indent: 20px;
}
&#10;#twqidtaqbo .gt_indent_5 {
  text-indent: 25px;
}
&#10;#twqidtaqbo .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#twqidtaqbo div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <caption>Modified NB-2 Model Summary</caption>
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="parameter">parameter</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="coeff">coeff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Std. Err.">Std. Err.</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="t-stat">t-stat</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="p-value">p-value</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lower CI">lower CI</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="upper CI">upper CI</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="parameter" class="gt_row gt_left">(Intercept)</td>
<td headers="coeff" class="gt_row gt_right">-7.364</td>
<td headers="Std. Err." class="gt_row gt_right">0.043</td>
<td headers="t-stat" class="gt_row gt_right">-173.167</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">-7.447</td>
<td headers="upper CI" class="gt_row gt_right">-7.281</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">lnaadt</td>
<td headers="coeff" class="gt_row gt_right">0.894</td>
<td headers="Std. Err." class="gt_row gt_right">0.005</td>
<td headers="t-stat" class="gt_row gt_right">181.134</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.884</td>
<td headers="upper CI" class="gt_row gt_right">0.904</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">lnlength</td>
<td headers="coeff" class="gt_row gt_right">0.862</td>
<td headers="Std. Err." class="gt_row gt_right">0.037</td>
<td headers="t-stat" class="gt_row gt_right">23.373</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.790</td>
<td headers="upper CI" class="gt_row gt_right">0.935</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">speed50</td>
<td headers="coeff" class="gt_row gt_right">-0.409</td>
<td headers="Std. Err." class="gt_row gt_right">0.092</td>
<td headers="t-stat" class="gt_row gt_right">-4.462</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">-0.588</td>
<td headers="upper CI" class="gt_row gt_right">-0.229</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ShouldWidth04</td>
<td headers="coeff" class="gt_row gt_right">0.251</td>
<td headers="Std. Err." class="gt_row gt_right">0.059</td>
<td headers="t-stat" class="gt_row gt_right">4.226</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.134</td>
<td headers="upper CI" class="gt_row gt_right">0.367</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">AADTover10k</td>
<td headers="coeff" class="gt_row gt_right">5.086</td>
<td headers="Std. Err." class="gt_row gt_right">0.090</td>
<td headers="t-stat" class="gt_row gt_right">56.822</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">4.910</td>
<td headers="upper CI" class="gt_row gt_right">5.261</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">I(AADTover10k/lnaadt)</td>
<td headers="coeff" class="gt_row gt_right">-41.179</td>
<td headers="Std. Err." class="gt_row gt_right">0.938</td>
<td headers="t-stat" class="gt_row gt_right">-43.907</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">-43.017</td>
<td headers="upper CI" class="gt_row gt_right">-39.341</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ln(alpha):  (Intercept)</td>
<td headers="coeff" class="gt_row gt_right">1.558</td>
<td headers="Std. Err." class="gt_row gt_right">0.340</td>
<td headers="t-stat" class="gt_row gt_right">4.581</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.892</td>
<td headers="upper CI" class="gt_row gt_right">2.225</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ln(alpha):  lnlength</td>
<td headers="coeff" class="gt_row gt_right">0.022</td>
<td headers="Std. Err." class="gt_row gt_right">0.350</td>
<td headers="t-stat" class="gt_row gt_right">0.064</td>
<td headers="p-value" class="gt_row gt_right">0.949</td>
<td headers="lower CI" class="gt_row gt_right">-0.663</td>
<td headers="upper CI" class="gt_row gt_right">0.708</td></tr>
  </tbody>
  &#10;  
</table>
</div>

``` r

tesstats <- regCompTest(gen.nb2)
tesstats$gtTable
```

<div id="pcbkgiknjo" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#pcbkgiknjo table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#pcbkgiknjo thead, #pcbkgiknjo tbody, #pcbkgiknjo tfoot, #pcbkgiknjo tr, #pcbkgiknjo td, #pcbkgiknjo th {
  border-style: none;
}
&#10;#pcbkgiknjo p {
  margin: 0;
  padding: 0;
}
&#10;#pcbkgiknjo .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#pcbkgiknjo .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#pcbkgiknjo .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#pcbkgiknjo .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#pcbkgiknjo .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#pcbkgiknjo .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#pcbkgiknjo .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#pcbkgiknjo .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#pcbkgiknjo .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#pcbkgiknjo .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#pcbkgiknjo .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#pcbkgiknjo .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#pcbkgiknjo .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#pcbkgiknjo .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#pcbkgiknjo .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pcbkgiknjo .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#pcbkgiknjo .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#pcbkgiknjo .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#pcbkgiknjo .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pcbkgiknjo .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#pcbkgiknjo .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pcbkgiknjo .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#pcbkgiknjo .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pcbkgiknjo .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#pcbkgiknjo .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pcbkgiknjo .gt_left {
  text-align: left;
}
&#10;#pcbkgiknjo .gt_center {
  text-align: center;
}
&#10;#pcbkgiknjo .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#pcbkgiknjo .gt_font_normal {
  font-weight: normal;
}
&#10;#pcbkgiknjo .gt_font_bold {
  font-weight: bold;
}
&#10;#pcbkgiknjo .gt_font_italic {
  font-style: italic;
}
&#10;#pcbkgiknjo .gt_super {
  font-size: 65%;
}
&#10;#pcbkgiknjo .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#pcbkgiknjo .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#pcbkgiknjo .gt_indent_1 {
  text-indent: 5px;
}
&#10;#pcbkgiknjo .gt_indent_2 {
  text-indent: 10px;
}
&#10;#pcbkgiknjo .gt_indent_3 {
  text-indent: 15px;
}
&#10;#pcbkgiknjo .gt_indent_4 {
  text-indent: 20px;
}
&#10;#pcbkgiknjo .gt_indent_5 {
  text-indent: 25px;
}
&#10;#pcbkgiknjo .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#pcbkgiknjo div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Model Comparison Statistics</td>
    </tr>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Statistic">Statistic</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Model">Model</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="BaseModel">BaseModel</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
<td headers="Model" class="gt_row gt_right">2,143.0082</td>
<td headers="BaseModel" class="gt_row gt_right">3,049.6592</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
<td headers="Model" class="gt_row gt_right">2,190.8331</td>
<td headers="BaseModel" class="gt_row gt_right">3,054.9731</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
<td headers="Model" class="gt_row gt_right">922.6510</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
<td headers="Model" class="gt_row gt_right">8.0000</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
<td headers="Model" class="gt_row gt_right">0.0000</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
<td headers="Model" class="gt_row gt_right">0.3027</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
  </tbody>
  &#10;  
</table>
</div>

``` r
cureplot(gen.nb2, indvar  ="lnaadt")
#> [1] 1501   14
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
gt(summary(gen.nbp), caption = "NB-P Model Summary")
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
#> 8 ln(alpha):  lnlen… -0.116       0.256   -0.455     0.649     -0.619      0.386
#> 9 P                   1.69        0.291    5.80      0          1.12       2.26
```

<div id="rblteafarb" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#rblteafarb table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#rblteafarb thead, #rblteafarb tbody, #rblteafarb tfoot, #rblteafarb tr, #rblteafarb td, #rblteafarb th {
  border-style: none;
}
&#10;#rblteafarb p {
  margin: 0;
  padding: 0;
}
&#10;#rblteafarb .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#rblteafarb .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#rblteafarb .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#rblteafarb .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#rblteafarb .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#rblteafarb .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rblteafarb .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#rblteafarb .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#rblteafarb .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#rblteafarb .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#rblteafarb .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#rblteafarb .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#rblteafarb .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#rblteafarb .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#rblteafarb .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#rblteafarb .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#rblteafarb .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#rblteafarb .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#rblteafarb .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rblteafarb .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#rblteafarb .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#rblteafarb .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#rblteafarb .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rblteafarb .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#rblteafarb .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#rblteafarb .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rblteafarb .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rblteafarb .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#rblteafarb .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#rblteafarb .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#rblteafarb .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#rblteafarb .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#rblteafarb .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rblteafarb .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#rblteafarb .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#rblteafarb .gt_left {
  text-align: left;
}
&#10;#rblteafarb .gt_center {
  text-align: center;
}
&#10;#rblteafarb .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#rblteafarb .gt_font_normal {
  font-weight: normal;
}
&#10;#rblteafarb .gt_font_bold {
  font-weight: bold;
}
&#10;#rblteafarb .gt_font_italic {
  font-style: italic;
}
&#10;#rblteafarb .gt_super {
  font-size: 65%;
}
&#10;#rblteafarb .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#rblteafarb .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#rblteafarb .gt_indent_1 {
  text-indent: 5px;
}
&#10;#rblteafarb .gt_indent_2 {
  text-indent: 10px;
}
&#10;#rblteafarb .gt_indent_3 {
  text-indent: 15px;
}
&#10;#rblteafarb .gt_indent_4 {
  text-indent: 20px;
}
&#10;#rblteafarb .gt_indent_5 {
  text-indent: 25px;
}
&#10;#rblteafarb .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#rblteafarb div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <caption>NB-P Model Summary</caption>
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="parameter">parameter</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="coeff">coeff</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Std. Err.">Std. Err.</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="t-stat">t-stat</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="p-value">p-value</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lower CI">lower CI</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="upper CI">upper CI</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="parameter" class="gt_row gt_left">(Intercept)</td>
<td headers="coeff" class="gt_row gt_right">-7.727</td>
<td headers="Std. Err." class="gt_row gt_right">0.043</td>
<td headers="t-stat" class="gt_row gt_right">-180.434</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">-7.811</td>
<td headers="upper CI" class="gt_row gt_right">-7.643</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">lnaadt</td>
<td headers="coeff" class="gt_row gt_right">0.935</td>
<td headers="Std. Err." class="gt_row gt_right">0.005</td>
<td headers="t-stat" class="gt_row gt_right">188.818</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.925</td>
<td headers="upper CI" class="gt_row gt_right">0.944</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">lnlength</td>
<td headers="coeff" class="gt_row gt_right">0.839</td>
<td headers="Std. Err." class="gt_row gt_right">0.037</td>
<td headers="t-stat" class="gt_row gt_right">22.456</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.766</td>
<td headers="upper CI" class="gt_row gt_right">0.913</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">speed50</td>
<td headers="coeff" class="gt_row gt_right">-0.390</td>
<td headers="Std. Err." class="gt_row gt_right">0.093</td>
<td headers="t-stat" class="gt_row gt_right">-4.183</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">-0.573</td>
<td headers="upper CI" class="gt_row gt_right">-0.207</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ShouldWidth04</td>
<td headers="coeff" class="gt_row gt_right">0.255</td>
<td headers="Std. Err." class="gt_row gt_right">0.059</td>
<td headers="t-stat" class="gt_row gt_right">4.293</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.139</td>
<td headers="upper CI" class="gt_row gt_right">0.372</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">AADTover10k</td>
<td headers="coeff" class="gt_row gt_right">0.693</td>
<td headers="Std. Err." class="gt_row gt_right">0.088</td>
<td headers="t-stat" class="gt_row gt_right">7.918</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">0.522</td>
<td headers="upper CI" class="gt_row gt_right">0.865</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ln(alpha):  (Intercept)</td>
<td headers="coeff" class="gt_row gt_right">-1.447</td>
<td headers="Std. Err." class="gt_row gt_right">0.297</td>
<td headers="t-stat" class="gt_row gt_right">-4.878</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">-2.028</td>
<td headers="upper CI" class="gt_row gt_right">-0.865</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">ln(alpha):  lnlength</td>
<td headers="coeff" class="gt_row gt_right">-0.116</td>
<td headers="Std. Err." class="gt_row gt_right">0.256</td>
<td headers="t-stat" class="gt_row gt_right">-0.455</td>
<td headers="p-value" class="gt_row gt_right">0.649</td>
<td headers="lower CI" class="gt_row gt_right">-0.619</td>
<td headers="upper CI" class="gt_row gt_right">0.386</td></tr>
    <tr><td headers="parameter" class="gt_row gt_left">P</td>
<td headers="coeff" class="gt_row gt_right">1.686</td>
<td headers="Std. Err." class="gt_row gt_right">0.291</td>
<td headers="t-stat" class="gt_row gt_right">5.796</td>
<td headers="p-value" class="gt_row gt_right">0.000</td>
<td headers="lower CI" class="gt_row gt_right">1.116</td>
<td headers="upper CI" class="gt_row gt_right">2.256</td></tr>
  </tbody>
  &#10;  
</table>
</div>

``` r

tesstats <- regCompTest(gen.nbp)
tesstats$gtTable
```

<div id="swxyxinati" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#swxyxinati table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#swxyxinati thead, #swxyxinati tbody, #swxyxinati tfoot, #swxyxinati tr, #swxyxinati td, #swxyxinati th {
  border-style: none;
}
&#10;#swxyxinati p {
  margin: 0;
  padding: 0;
}
&#10;#swxyxinati .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#swxyxinati .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#swxyxinati .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#swxyxinati .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#swxyxinati .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#swxyxinati .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#swxyxinati .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#swxyxinati .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#swxyxinati .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#swxyxinati .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#swxyxinati .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#swxyxinati .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#swxyxinati .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#swxyxinati .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#swxyxinati .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#swxyxinati .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#swxyxinati .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#swxyxinati .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#swxyxinati .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#swxyxinati .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#swxyxinati .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#swxyxinati .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#swxyxinati .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#swxyxinati .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#swxyxinati .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#swxyxinati .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#swxyxinati .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#swxyxinati .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#swxyxinati .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#swxyxinati .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#swxyxinati .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#swxyxinati .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#swxyxinati .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#swxyxinati .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#swxyxinati .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#swxyxinati .gt_left {
  text-align: left;
}
&#10;#swxyxinati .gt_center {
  text-align: center;
}
&#10;#swxyxinati .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#swxyxinati .gt_font_normal {
  font-weight: normal;
}
&#10;#swxyxinati .gt_font_bold {
  font-weight: bold;
}
&#10;#swxyxinati .gt_font_italic {
  font-style: italic;
}
&#10;#swxyxinati .gt_super {
  font-size: 65%;
}
&#10;#swxyxinati .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#swxyxinati .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#swxyxinati .gt_indent_1 {
  text-indent: 5px;
}
&#10;#swxyxinati .gt_indent_2 {
  text-indent: 10px;
}
&#10;#swxyxinati .gt_indent_3 {
  text-indent: 15px;
}
&#10;#swxyxinati .gt_indent_4 {
  text-indent: 20px;
}
&#10;#swxyxinati .gt_indent_5 {
  text-indent: 25px;
}
&#10;#swxyxinati .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#swxyxinati div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Model Comparison Statistics</td>
    </tr>
    &#10;    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Statistic">Statistic</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Model">Model</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="BaseModel">BaseModel</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
<td headers="Model" class="gt_row gt_right">2,142.4125</td>
<td headers="BaseModel" class="gt_row gt_right">3,049.6592</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
<td headers="Model" class="gt_row gt_right">2,190.2375</td>
<td headers="BaseModel" class="gt_row gt_right">3,054.9731</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
<td headers="Model" class="gt_row gt_right">923.2467</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
<td headers="Model" class="gt_row gt_right">8.0000</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
<td headers="Model" class="gt_row gt_right">0.0000</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
    <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
<td headers="Model" class="gt_row gt_right">0.3029</td>
<td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
  </tbody>
  &#10;  
</table>
</div>

Checking the CURE plot (notice that the CURE plot is MUCH better in this
case than the NB-2 without the interaction and still better than the
modified NB-2):

``` r
cureplot(gen.nbp, indvar  ="lnaadt")
#> [1] 1501   14
#> Covariate: indvar_values
#> CURE data frame was provided. Its first column, lnaadt, will be used.
```

<img src="man/figures/README-cureplot-lnaadt-1.png" width="100%" />

Creating a table to compare the models:

``` r
regCompTable(list("Generalized NB-2"=gen.nb2, "Generalized NB-P"=gen.nbp), tableType="gt")
```

<div id="pxhcpudcts" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#pxhcpudcts table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#pxhcpudcts thead, #pxhcpudcts tbody, #pxhcpudcts tfoot, #pxhcpudcts tr, #pxhcpudcts td, #pxhcpudcts th {
  border-style: none;
}
&#10;#pxhcpudcts p {
  margin: 0;
  padding: 0;
}
&#10;#pxhcpudcts .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#pxhcpudcts .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#pxhcpudcts .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#pxhcpudcts .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#pxhcpudcts .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#pxhcpudcts .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#pxhcpudcts .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#pxhcpudcts .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#pxhcpudcts .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#pxhcpudcts .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#pxhcpudcts .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#pxhcpudcts .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#pxhcpudcts .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#pxhcpudcts .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#pxhcpudcts .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pxhcpudcts .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#pxhcpudcts .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#pxhcpudcts .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#pxhcpudcts .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pxhcpudcts .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#pxhcpudcts .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pxhcpudcts .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#pxhcpudcts .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pxhcpudcts .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#pxhcpudcts .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#pxhcpudcts .gt_left {
  text-align: left;
}
&#10;#pxhcpudcts .gt_center {
  text-align: center;
}
&#10;#pxhcpudcts .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#pxhcpudcts .gt_font_normal {
  font-weight: normal;
}
&#10;#pxhcpudcts .gt_font_bold {
  font-weight: bold;
}
&#10;#pxhcpudcts .gt_font_italic {
  font-style: italic;
}
&#10;#pxhcpudcts .gt_super {
  font-size: 65%;
}
&#10;#pxhcpudcts .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#pxhcpudcts .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#pxhcpudcts .gt_indent_1 {
  text-indent: 5px;
}
&#10;#pxhcpudcts .gt_indent_2 {
  text-indent: 10px;
}
&#10;#pxhcpudcts .gt_indent_3 {
  text-indent: 15px;
}
&#10;#pxhcpudcts .gt_indent_4 {
  text-indent: 20px;
}
&#10;#pxhcpudcts .gt_indent_5 {
  text-indent: 25px;
}
&#10;#pxhcpudcts .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#pxhcpudcts div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_title gt_font_normal" style>Model Comparisons</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="3" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>Mean (Standard Error)</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Parameter">Parameter</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Generalized NB-2">Generalized NB-2</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Generalized NB-P">Generalized NB-P</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Parameter" class="gt_row gt_left">(Intercept)</td>
<td headers="Generalized NB-2" class="gt_row gt_right">-7.364 (0.043)***</td>
<td headers="Generalized NB-P" class="gt_row gt_right">-7.727 (0.043)***</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">lnaadt</td>
<td headers="Generalized NB-2" class="gt_row gt_right">0.894 (0.005)***</td>
<td headers="Generalized NB-P" class="gt_row gt_right">0.935 (0.005)***</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">lnlength</td>
<td headers="Generalized NB-2" class="gt_row gt_right">0.862 (0.037)***</td>
<td headers="Generalized NB-P" class="gt_row gt_right">0.839 (0.037)***</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">speed50</td>
<td headers="Generalized NB-2" class="gt_row gt_right">-0.409 (0.092)***</td>
<td headers="Generalized NB-P" class="gt_row gt_right">-0.39 (0.093)***</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">ShouldWidth04</td>
<td headers="Generalized NB-2" class="gt_row gt_right">0.251 (0.059)***</td>
<td headers="Generalized NB-P" class="gt_row gt_right">0.255 (0.059)***</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">AADTover10k</td>
<td headers="Generalized NB-2" class="gt_row gt_right">5.086 (0.09)***</td>
<td headers="Generalized NB-P" class="gt_row gt_right">0.693 (0.088)***</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">I(AADTover10k/lnaadt)</td>
<td headers="Generalized NB-2" class="gt_row gt_right">-41.179 (0.938)***</td>
<td headers="Generalized NB-P" class="gt_row gt_right">---</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">ln(alpha):  (Intercept)</td>
<td headers="Generalized NB-2" class="gt_row gt_right">1.558 (0.34)***</td>
<td headers="Generalized NB-P" class="gt_row gt_right">-1.447 (0.297)***</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">ln(alpha):  lnlength</td>
<td headers="Generalized NB-2" class="gt_row gt_right">0.022 (0.35)</td>
<td headers="Generalized NB-P" class="gt_row gt_right">-0.116 (0.256)</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">P</td>
<td headers="Generalized NB-2" class="gt_row gt_right">---</td>
<td headers="Generalized NB-P" class="gt_row gt_right">1.686 (0.291)***</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">N Obs.</td>
<td headers="Generalized NB-2" class="gt_row gt_right">1501</td>
<td headers="Generalized NB-P" class="gt_row gt_right">1501</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">LL</td>
<td headers="Generalized NB-2" class="gt_row gt_right">-1062.504</td>
<td headers="Generalized NB-P" class="gt_row gt_right">-1062.206</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">AIC</td>
<td headers="Generalized NB-2" class="gt_row gt_right">2143.008</td>
<td headers="Generalized NB-P" class="gt_row gt_right">2142.412</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">BIC</td>
<td headers="Generalized NB-2" class="gt_row gt_right">2190.833</td>
<td headers="Generalized NB-P" class="gt_row gt_right">2190.237</td></tr>
    <tr><td headers="Parameter" class="gt_row gt_left">Pseudo-R-Sq.</td>
<td headers="Generalized NB-2" class="gt_row gt_right">0.303</td>
<td headers="Generalized NB-P" class="gt_row gt_right">0.303</td></tr>
  </tbody>
  <tfoot class="gt_sourcenotes">
    <tr>
      <td class="gt_sourcenote" colspan="3">p-value codes: *=(p&lt;=0.05), **=(p&lt;=0.01), ***=(p&lt;=0.001)</td>
    </tr>
  </tfoot>
  &#10;</table>
</div>

Note that the metrics for comparison are similar. While the models both
have the same number of parameters, the NB-P was able to get better
performance without requiring the interaction terms (which leads to
strange relationships between the exposure metric and the outcome).
