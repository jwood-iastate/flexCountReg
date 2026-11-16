# Compare Regression Models with Likelihood Ratio Test, AIC, and BIC

This function compares fitted regression model objects using the
Likelihood Ratio (LR) test, Akaike Information Criterion (AIC), and
Bayesian Information Criterion (BIC).

## Usage

``` r
regCompTest(
  model,
  data = NULL,
  basemodel = "Poisson",
  variables = FALSE,
  print = FALSE,
  ...,
  model2 = NULL
)
```

## Arguments

- model:

  A fitted regression model object.

- data:

  An optional data frame containing the variables in the model. If not
  supplied, the original data used to estimate the model will be used
  when available.

- basemodel:

  A character string specifying the family of base model to compare
  against (options include the family from
  [`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  or "Poisson"). Default is "Poisson".

- variables:

  Logical. If `TRUE`, the base model will include the same variables as
  the provided model. If `FALSE`, the base model will be an
  intercept-only model. Default is `FALSE`.

- print:

  Logical. If `TRUE`, a table of the results will be shown. If `FALSE`,
  the table of results will not be printed to the console.

- ...:

  Additional arguments to be passed to the base model fitting function.
  These arguments are used only when a single fitted model is supplied
  and the function needs to estimate the base model.

- model2:

  An optional second fitted regression model object. If supplied, the
  function compares `model` and `model2` directly and does not estimate
  a separate base model.

## Value

A list containing the comparison results.

For the single-model workflow, the returned list includes `LL`,
`LLbase`, `LR`, `LRdof`, `AIC`, `AICbase`, `BIC`, `BICbase`,
`LR_pvalue`, `PseudoR2`, `statistics`, `gtTable`, `latexTable`, and
`htmlTable`.

For the two-model workflow, the returned list includes `model1`,
`model2`, `base_model`, `comparison_model`, `base_index`,
`comparison_index`, `LL_base`, `LL_comparison`, `LR`, `LRdof`,
`LR_pvalue`, `statistics`, `gtTable`, `latexTable`, and `htmlTable`.

## Details

When only one fitted model is supplied, the function estimates a base
model on the same rows, offset and weights and compares the fitted model
against it. When two fitted models are supplied, the function compares
them directly; the model with fewer parameters is the restricted (base)
model of the LR test.

The function performs the following steps for the single-model workflow:

1.  Fits the base model, either a Poisson regression or another
    specified model.

2.  Computes the log-likelihoods of the provided model and the base
    model.

3.  Calculates the AIC and BIC for both models.

4.  Conducts a Likelihood Ratio test to compare the fitted model against
    the estimated base model when the fitted model has more parameters
    than the base model.

5.  Computes McFadden's Pseudo R^2, \\1 - LL/LL_0\\, where \\LL_0\\ is
    the log-likelihood of an intercept-only model of the same family for
    [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
    models (an intercept-only Poisson model for the other estimators).

For the two-model workflow, the model with fewer parameters is the
restricted (base) model and the other is the comparison (full) model.
The models must be fit to the same observations (an error is raised
otherwise) and should be nested; a warning is given when the restricted
model has the higher log-likelihood (a sign that the models are not
nested or that one did not converge), and no LR test is reported when
they have the same number of parameters. The LR test assumes the
restricted model is in the interior of the parameter space; when it is
on the boundary (for example Poisson within NB2), the p-value is
conservative.

The Likelihood-Ratio test is computed as \$\$LR = 2 (LL\_{comparison} -
LL\_{base})\$\$ and is chi-squared with degrees of freedom equal to the
difference in the number of parameters. AIC and BIC are those of
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html); for panel models, BIC uses
the number of panels.

## Examples

``` r

# Example 1: Compare one fitted model against an estimated Poisson base model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)

nbp.base <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
                     ShouldWidth04 + AADTover10k,
                     data = washington_roads, family = "NBP", method = "NM",
                     max.iters = 3000)

regCompTest(nbp.base, data = washington_roads, basemodel = "Poisson",
            variables = TRUE, print = TRUE)
#> <div id="ujyntsapwa" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#ujyntsapwa table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #ujyntsapwa thead, #ujyntsapwa tbody, #ujyntsapwa tfoot, #ujyntsapwa tr, #ujyntsapwa td, #ujyntsapwa th {
#>   border-style: none;
#> }
#> 
#> #ujyntsapwa p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #ujyntsapwa .gt_table {
#>   display: table;
#>   border-collapse: collapse;
#>   line-height: normal;
#>   margin-left: auto;
#>   margin-right: auto;
#>   color: #333333;
#>   font-size: 16px;
#>   font-weight: normal;
#>   font-style: normal;
#>   background-color: #FFFFFF;
#>   width: auto;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #A8A8A8;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #A8A8A8;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #ujyntsapwa .gt_title {
#>   color: #333333;
#>   font-size: 125%;
#>   font-weight: initial;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-color: #FFFFFF;
#>   border-bottom-width: 0;
#> }
#> 
#> #ujyntsapwa .gt_subtitle {
#>   color: #333333;
#>   font-size: 85%;
#>   font-weight: initial;
#>   padding-top: 3px;
#>   padding-bottom: 5px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-color: #FFFFFF;
#>   border-top-width: 0;
#> }
#> 
#> #ujyntsapwa .gt_heading {
#>   background-color: #FFFFFF;
#>   text-align: center;
#>   border-bottom-color: #FFFFFF;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_col_headings {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_col_heading {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: normal;
#>   text-transform: inherit;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: bottom;
#>   padding-top: 5px;
#>   padding-bottom: 6px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   overflow-x: hidden;
#> }
#> 
#> #ujyntsapwa .gt_column_spanner_outer {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: normal;
#>   text-transform: inherit;
#>   padding-top: 0;
#>   padding-bottom: 0;
#>   padding-left: 4px;
#>   padding-right: 4px;
#> }
#> 
#> #ujyntsapwa .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #ujyntsapwa .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #ujyntsapwa .gt_column_spanner {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   vertical-align: bottom;
#>   padding-top: 5px;
#>   padding-bottom: 5px;
#>   overflow-x: hidden;
#>   display: inline-block;
#>   width: 100%;
#> }
#> 
#> #ujyntsapwa .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #ujyntsapwa .gt_group_heading {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: middle;
#>   text-align: left;
#> }
#> 
#> #ujyntsapwa .gt_empty_group_heading {
#>   padding: 0.5px;
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   vertical-align: middle;
#> }
#> 
#> #ujyntsapwa .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #ujyntsapwa .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #ujyntsapwa .gt_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   margin: 10px;
#>   border-top-style: solid;
#>   border-top-width: 1px;
#>   border-top-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: middle;
#>   overflow-x: hidden;
#> }
#> 
#> #ujyntsapwa .gt_stub {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-right-style: solid;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #ujyntsapwa .gt_stub_row_group {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-right-style: solid;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   vertical-align: top;
#> }
#> 
#> #ujyntsapwa .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #ujyntsapwa .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #ujyntsapwa .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #ujyntsapwa .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #ujyntsapwa .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #ujyntsapwa .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #ujyntsapwa .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_footnotes {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   border-bottom-style: none;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #ujyntsapwa .gt_sourcenotes {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   border-bottom-style: none;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #ujyntsapwa .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #ujyntsapwa .gt_left {
#>   text-align: left;
#> }
#> 
#> #ujyntsapwa .gt_center {
#>   text-align: center;
#> }
#> 
#> #ujyntsapwa .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #ujyntsapwa .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #ujyntsapwa .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #ujyntsapwa .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #ujyntsapwa .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #ujyntsapwa .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #ujyntsapwa .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #ujyntsapwa .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #ujyntsapwa .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #ujyntsapwa .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #ujyntsapwa .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #ujyntsapwa .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #ujyntsapwa .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #ujyntsapwa div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
#>   height: 0px !important;
#> }
#> </style>
#>   <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
#>   <thead>
#>     <tr class="gt_heading">
#>       <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Model Comparison Statistics</td>
#>     </tr>
#>     
#>     <tr class="gt_col_headings">
#>       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Statistic">Statistic</th>
#>       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Model">Model</th>
#>       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="BaseModel">BaseModel</th>
#>     </tr>
#>   </thead>
#>   <tbody class="gt_table_body">
#>     <tr><td headers="Statistic" class="gt_row gt_left">Log-likelihood</td>
#> <td headers="Model" class="gt_row gt_right">−1,062.2659</td>
#> <td headers="BaseModel" class="gt_row gt_right">−1,071.2672</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">Number of parameters</td>
#> <td headers="Model" class="gt_row gt_right">8.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">6.0000</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
#> <td headers="Model" class="gt_row gt_right">2,140.5317</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,154.5345</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
#> <td headers="Model" class="gt_row gt_right">2,183.0428</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,186.4178</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
#> <td headers="Model" class="gt_row gt_right">18.0027</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">2.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">0.0001</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.2083</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>   </tbody>
#>   
#> </table>
#> </div>
#> $model
#> Maximum Likelihood estimation
#> Nelder-Mead maximization, 917 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1062.266 (8 free parameter(s))
#> Estimate(s): -7.796102 0.9417162 0.8342337 -0.3870385 0.2592118 0.6860279 -1.326803 0.5023309 
#> 
#> $y
#>    [1]  0  2  2  0  0  1  2  0  1  0  0  0  0  1  0  0  1  0  0  0  0  2  0  0
#>   [25]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#>   [49]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  1  0
#>   [73]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0
#>   [97]  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1
#>  [121]  0  0  1  0  0  0  0  0  2  0  0  0  0  0  1  0  0  2  0  0  0  0  2  0
#>  [145]  0  1  1  0  0  0  0  0  0  0  2  2  1  1  3  0  1  2  1  1  0  0  0  0
#>  [169]  0  1  0  1  3  4  1  3  4  2  1  1  3  0  1  2  1  0  2  0  0  1  1  0
#>  [193]  8  1  2  2  1  3  4  5  0  1  6  4  0  0  1  4  0  0  0  0  0  0  1  1
#>  [217]  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2
#>  [241]  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0
#>  [265]  0  0  0  0  1  1  0  0  0  0  0  1  0  0  0  2  0  0  0  0  1  1  2  0
#>  [289]  1  2  2  2  1  1  2  0  2  0  1  4  0  0  0  2  1  1  4 10  3  0  1  0
#>  [313]  0  0  1  3  0  2  2  1  0  1  2  3  0  0  1  1  2  2  0  2  3  1  0  0
#>  [337]  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1
#>  [361]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#>  [385]  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  2  2  0  0  1  0
#>  [409]  1  1  0  0  0  1  1  1  0  0  1  1  0  0  0  0  0  0  1  2  0  1  0  0
#>  [433]  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  1
#>  [457]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  1
#>  [481]  0  0  1  1  0  0  0  0  3  0  1  0  0  1  0  0  2  0  0  0  7  0  0  0
#>  [505]  2  1  0  1  0  0  0  1  0  0  0  0  1  3  0  0  0  0  0  1  0  0  0  0
#>  [529]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0
#>  [553]  0  1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  1  0
#>  [577]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0
#>  [601]  0  0  0  2  0  0  1  0  0  0  0  0  0  1  1  2  0  0  0  0  0  0  0  1
#>  [625]  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  2  0  0  0  0  0  0  0
#>  [649]  1  1  0  0  1  3  0  3  4  2  3  1  0  1  1  0  0  0  0  1  0  1  0  2
#>  [673]  0  1  2  1  2  4  2  2  3  2  3  2  1  0  0  0  0  0  1  1  1  5  0  2
#>  [697]  5  3  1  3  2  5  2  2  2  1  4  0  0  0  0  0  0  0  0  0  1  0  0  0
#>  [721]  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0
#>  [745]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  1
#>  [769]  0  1  0  0  0  0  0  1  0  1  0  0  0  0  1  0  1  0  0  0  2  2  1  1
#>  [793]  1  2  1  1  0  0  0  0  1  1  1  2  0  2  1  4  0  2  1  3  0  1  2  2
#>  [817]  1  0  4  1  0  0  1  1  1  1  0  0  0  0  0  0  0  1  0  0  1  0  0  0
#>  [841]  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#>  [865]  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
#>  [889]  0  1  0  0  0  1  0  1  0  0  0  0  2  0  1  2  1  0  1  1  1  1  0  0
#>  [913]  0  0  3  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#>  [937]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#>  [961]  1  0  0  0  0  1  0  1  0  0  1  0  0  0  1  0  0  0  0  2  0  0  1  0
#>  [985]  0  0  0  0  0  0  0  0  0  0  0  0  3  3  2  2  8  1  3  0  1  0  1  0
#> [1009]  0  0  1  1  0  0  3  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> [1033]  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0
#> [1057]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0
#> [1081]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0  0  0  0  0  0  0
#> [1105]  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0  0  0  1  0  2  1  4  1  0
#> [1129]  0  0  0  0  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  1  1  0  0  0
#> [1153]  0  1  1  1  7  0  4  3  0  1  1  0  0  0  3  1  1  1  1  0  0  1  3  0
#> [1177]  4  2  1  4  2  2  1  1  0  0  0  0  0  0  4  2  0  4  0  3  7  2  4  2
#> [1201]  4  2  6  0  1  2  0  0  0  0  0  0  0  1  1  0  0  0  1  0  0  0  0  1
#> [1225]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0
#> [1249]  0  0  0  0  0  1  0  0  0  0  0  0  0  1  2  0  0  0  0  2  0  0  0  0
#> [1273]  0  0  0  1  0  0  0  1  0  0  1  2  1  2  1  0  2  0  2  1  1  2  1  1
#> [1297]  2  0  0  0  0  0  1  0  0  2  0  4  1  3  0  0  0  0  0  0  4  1  5  0
#> [1321]  0  1  3  2  0  1  2  0  0  0  0  0  2  1  0  0  0  0  0  1  0  0  0  0
#> [1345]  0  0  0  0  0  0  0  0  1  1  1  0  1  0  0  0  0  0  0  0  0  0  0  0
#> [1369]  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
#> [1393]  0  0  0  0  0  0  0  0  4  0  0  0  0  1  0  1  1  0  2  1  0  0  2  0
#> [1417]  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> [1441]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0  0
#> [1465]  0  1  1  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  2  0  1  0  0
#> [1489]  1  0  0  0  1  0  0  0  0  1  2  1  5
#> 
#> $LL
#> [1] -1062.266
#> 
#> $LLbase
#> [1] -1071.267
#> 
#> $LR
#> [1] 18.00274
#> 
#> $LRdof
#> [1] 2
#> 
#> $AIC
#> [1] 2140.532
#> 
#> $AICbase
#> [1] 2154.534
#> 
#> $BIC
#> [1] 2183.043
#> 
#> $BICbase
#> [1] 2186.418
#> 
#> $LR_pvalue
#> [1] 0.0001232406
#> 
#> $PseudoR2
#> [1] 0.2083299
#> 
#> $statistics
#> # A tibble: 8 × 3
#>   Statistic                  Model BaseModel
#>   <chr>                      <dbl>     <dbl>
#> 1 Log-likelihood        -1062.        -1071.
#> 2 Number of parameters      8             6 
#> 3 AIC                    2141.         2155.
#> 4 BIC                    2183.         2186.
#> 5 LR Test Statistic        18.0          NA 
#> 6 LR degrees of freedom     2            NA 
#> 7 LR p-value                0.0001       NA 
#> 8 McFadden's Pseudo R^2     0.208        NA 
#> 
#> $gtTable
#> <div id="zzmwnfmrog" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#zzmwnfmrog table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #zzmwnfmrog thead, #zzmwnfmrog tbody, #zzmwnfmrog tfoot, #zzmwnfmrog tr, #zzmwnfmrog td, #zzmwnfmrog th {
#>   border-style: none;
#> }
#> 
#> #zzmwnfmrog p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #zzmwnfmrog .gt_table {
#>   display: table;
#>   border-collapse: collapse;
#>   line-height: normal;
#>   margin-left: auto;
#>   margin-right: auto;
#>   color: #333333;
#>   font-size: 16px;
#>   font-weight: normal;
#>   font-style: normal;
#>   background-color: #FFFFFF;
#>   width: auto;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #A8A8A8;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #A8A8A8;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #zzmwnfmrog .gt_title {
#>   color: #333333;
#>   font-size: 125%;
#>   font-weight: initial;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-color: #FFFFFF;
#>   border-bottom-width: 0;
#> }
#> 
#> #zzmwnfmrog .gt_subtitle {
#>   color: #333333;
#>   font-size: 85%;
#>   font-weight: initial;
#>   padding-top: 3px;
#>   padding-bottom: 5px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-color: #FFFFFF;
#>   border-top-width: 0;
#> }
#> 
#> #zzmwnfmrog .gt_heading {
#>   background-color: #FFFFFF;
#>   text-align: center;
#>   border-bottom-color: #FFFFFF;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_col_headings {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_col_heading {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: normal;
#>   text-transform: inherit;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: bottom;
#>   padding-top: 5px;
#>   padding-bottom: 6px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   overflow-x: hidden;
#> }
#> 
#> #zzmwnfmrog .gt_column_spanner_outer {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: normal;
#>   text-transform: inherit;
#>   padding-top: 0;
#>   padding-bottom: 0;
#>   padding-left: 4px;
#>   padding-right: 4px;
#> }
#> 
#> #zzmwnfmrog .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #zzmwnfmrog .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #zzmwnfmrog .gt_column_spanner {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   vertical-align: bottom;
#>   padding-top: 5px;
#>   padding-bottom: 5px;
#>   overflow-x: hidden;
#>   display: inline-block;
#>   width: 100%;
#> }
#> 
#> #zzmwnfmrog .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #zzmwnfmrog .gt_group_heading {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: middle;
#>   text-align: left;
#> }
#> 
#> #zzmwnfmrog .gt_empty_group_heading {
#>   padding: 0.5px;
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   vertical-align: middle;
#> }
#> 
#> #zzmwnfmrog .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #zzmwnfmrog .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #zzmwnfmrog .gt_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   margin: 10px;
#>   border-top-style: solid;
#>   border-top-width: 1px;
#>   border-top-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: middle;
#>   overflow-x: hidden;
#> }
#> 
#> #zzmwnfmrog .gt_stub {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-right-style: solid;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #zzmwnfmrog .gt_stub_row_group {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-right-style: solid;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   vertical-align: top;
#> }
#> 
#> #zzmwnfmrog .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #zzmwnfmrog .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #zzmwnfmrog .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #zzmwnfmrog .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #zzmwnfmrog .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #zzmwnfmrog .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #zzmwnfmrog .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_footnotes {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   border-bottom-style: none;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #zzmwnfmrog .gt_sourcenotes {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   border-bottom-style: none;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #zzmwnfmrog .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #zzmwnfmrog .gt_left {
#>   text-align: left;
#> }
#> 
#> #zzmwnfmrog .gt_center {
#>   text-align: center;
#> }
#> 
#> #zzmwnfmrog .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #zzmwnfmrog .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #zzmwnfmrog .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #zzmwnfmrog .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #zzmwnfmrog .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #zzmwnfmrog .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #zzmwnfmrog .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #zzmwnfmrog .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #zzmwnfmrog .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #zzmwnfmrog .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #zzmwnfmrog .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #zzmwnfmrog .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #zzmwnfmrog .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #zzmwnfmrog div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
#>   height: 0px !important;
#> }
#> </style>
#>   <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
#>   <thead>
#>     <tr class="gt_heading">
#>       <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Model Comparison Statistics</td>
#>     </tr>
#>     
#>     <tr class="gt_col_headings">
#>       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Statistic">Statistic</th>
#>       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Model">Model</th>
#>       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="BaseModel">BaseModel</th>
#>     </tr>
#>   </thead>
#>   <tbody class="gt_table_body">
#>     <tr><td headers="Statistic" class="gt_row gt_left">Log-likelihood</td>
#> <td headers="Model" class="gt_row gt_right">−1,062.2659</td>
#> <td headers="BaseModel" class="gt_row gt_right">−1,071.2672</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">Number of parameters</td>
#> <td headers="Model" class="gt_row gt_right">8.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">6.0000</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
#> <td headers="Model" class="gt_row gt_right">2,140.5317</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,154.5345</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
#> <td headers="Model" class="gt_row gt_right">2,183.0428</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,186.4178</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
#> <td headers="Model" class="gt_row gt_right">18.0027</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">2.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">0.0001</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.2083</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>   </tbody>
#>   
#> </table>
#> </div>
#> 
#> $latexTable
#> \begin{table}
#> 
#> \caption{Model Comparison Statistics}
#> \centering
#> \begin{tabular}[t]{lrr}
#> \toprule
#> Statistic & Model & BaseModel\\
#> \midrule
#> Log-likelihood & -1062.2659 & -1071.267\\
#> Number of parameters & 8.0000 & 6.000\\
#> AIC & 2140.5317 & 2154.535\\
#> BIC & 2183.0428 & 2186.418\\
#> LR Test Statistic & 18.0027 & NA\\
#> \addlinespace
#> LR degrees of freedom & 2.0000 & NA\\
#> LR p-value & 0.0001 & NA\\
#> McFadden's Pseudo R\textasciicircum{}2 & 0.2083 & NA\\
#> \bottomrule
#> \end{tabular}
#> \end{table}
#> 
#> $htmlTable
#> <table class='table table-striped'>
#> <caption>Model Comparison Statistics</caption>
#>  <thead>
#>   <tr>
#>    <th style="text-align:left;"> Statistic </th>
#>    <th style="text-align:right;"> Model </th>
#>    <th style="text-align:right;"> BaseModel </th>
#>   </tr>
#>  </thead>
#> <tbody>
#>   <tr>
#>    <td style="text-align:left;"> Log-likelihood </td>
#>    <td style="text-align:right;"> -1062.2659 </td>
#>    <td style="text-align:right;"> -1071.267 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> Number of parameters </td>
#>    <td style="text-align:right;"> 8.0000 </td>
#>    <td style="text-align:right;"> 6.000 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> AIC </td>
#>    <td style="text-align:right;"> 2140.5317 </td>
#>    <td style="text-align:right;"> 2154.535 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> BIC </td>
#>    <td style="text-align:right;"> 2183.0428 </td>
#>    <td style="text-align:right;"> 2186.418 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR Test Statistic </td>
#>    <td style="text-align:right;"> 18.0027 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR degrees of freedom </td>
#>    <td style="text-align:right;"> 2.0000 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR p-value </td>
#>    <td style="text-align:right;"> 0.0001 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> McFadden's Pseudo R^2 </td>
#>    <td style="text-align:right;"> 0.2083 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#> </tbody>
#> </table>
#> 

# Example 2: Compare two fitted models directly
nb2.base <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
                     ShouldWidth04 + AADTover10k,
                     data = washington_roads, family = "NB2", method = "NM",
                     max.iters = 3000)
#> Warning: The optimizer did not report convergence (iteration limit exceeded). Consider a different `method` or more `max.iters`.

regCompTest(nbp.base, data = washington_roads, model2 = nb2.base,
            print = TRUE)
#> <div id="vxkvczrocj" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#vxkvczrocj table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #vxkvczrocj thead, #vxkvczrocj tbody, #vxkvczrocj tfoot, #vxkvczrocj tr, #vxkvczrocj td, #vxkvczrocj th {
#>   border-style: none;
#> }
#> 
#> #vxkvczrocj p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #vxkvczrocj .gt_table {
#>   display: table;
#>   border-collapse: collapse;
#>   line-height: normal;
#>   margin-left: auto;
#>   margin-right: auto;
#>   color: #333333;
#>   font-size: 16px;
#>   font-weight: normal;
#>   font-style: normal;
#>   background-color: #FFFFFF;
#>   width: auto;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #A8A8A8;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #A8A8A8;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #vxkvczrocj .gt_title {
#>   color: #333333;
#>   font-size: 125%;
#>   font-weight: initial;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-color: #FFFFFF;
#>   border-bottom-width: 0;
#> }
#> 
#> #vxkvczrocj .gt_subtitle {
#>   color: #333333;
#>   font-size: 85%;
#>   font-weight: initial;
#>   padding-top: 3px;
#>   padding-bottom: 5px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-color: #FFFFFF;
#>   border-top-width: 0;
#> }
#> 
#> #vxkvczrocj .gt_heading {
#>   background-color: #FFFFFF;
#>   text-align: center;
#>   border-bottom-color: #FFFFFF;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_col_headings {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_col_heading {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: normal;
#>   text-transform: inherit;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: bottom;
#>   padding-top: 5px;
#>   padding-bottom: 6px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   overflow-x: hidden;
#> }
#> 
#> #vxkvczrocj .gt_column_spanner_outer {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: normal;
#>   text-transform: inherit;
#>   padding-top: 0;
#>   padding-bottom: 0;
#>   padding-left: 4px;
#>   padding-right: 4px;
#> }
#> 
#> #vxkvczrocj .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #vxkvczrocj .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #vxkvczrocj .gt_column_spanner {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   vertical-align: bottom;
#>   padding-top: 5px;
#>   padding-bottom: 5px;
#>   overflow-x: hidden;
#>   display: inline-block;
#>   width: 100%;
#> }
#> 
#> #vxkvczrocj .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #vxkvczrocj .gt_group_heading {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: middle;
#>   text-align: left;
#> }
#> 
#> #vxkvczrocj .gt_empty_group_heading {
#>   padding: 0.5px;
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   vertical-align: middle;
#> }
#> 
#> #vxkvczrocj .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #vxkvczrocj .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #vxkvczrocj .gt_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   margin: 10px;
#>   border-top-style: solid;
#>   border-top-width: 1px;
#>   border-top-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: middle;
#>   overflow-x: hidden;
#> }
#> 
#> #vxkvczrocj .gt_stub {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-right-style: solid;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #vxkvczrocj .gt_stub_row_group {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-right-style: solid;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   vertical-align: top;
#> }
#> 
#> #vxkvczrocj .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #vxkvczrocj .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #vxkvczrocj .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #vxkvczrocj .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #vxkvczrocj .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #vxkvczrocj .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #vxkvczrocj .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_footnotes {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   border-bottom-style: none;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #vxkvczrocj .gt_sourcenotes {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   border-bottom-style: none;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #vxkvczrocj .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #vxkvczrocj .gt_left {
#>   text-align: left;
#> }
#> 
#> #vxkvczrocj .gt_center {
#>   text-align: center;
#> }
#> 
#> #vxkvczrocj .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #vxkvczrocj .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #vxkvczrocj .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #vxkvczrocj .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #vxkvczrocj .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #vxkvczrocj .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #vxkvczrocj .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #vxkvczrocj .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #vxkvczrocj .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #vxkvczrocj .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #vxkvczrocj .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #vxkvczrocj .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #vxkvczrocj .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #vxkvczrocj div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
#>   height: 0px !important;
#> }
#> </style>
#>   <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
#>   <thead>
#>     <tr class="gt_heading">
#>       <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Model Comparison Statistics</td>
#>     </tr>
#>     
#>     <tr class="gt_col_headings">
#>       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Statistic">Statistic</th>
#>       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Model">Model</th>
#>       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="BaseModel">BaseModel</th>
#>     </tr>
#>   </thead>
#>   <tbody class="gt_table_body">
#>     <tr><td headers="Statistic" class="gt_row gt_left">Log-likelihood</td>
#> <td headers="Model" class="gt_row gt_right">−1,062.2659</td>
#> <td headers="BaseModel" class="gt_row gt_right">−1,063.4206</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">Number of parameters</td>
#> <td headers="Model" class="gt_row gt_right">8.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">7.0000</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
#> <td headers="Model" class="gt_row gt_right">2,140.5317</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,140.8411</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
#> <td headers="Model" class="gt_row gt_right">2,183.0428</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,178.0384</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
#> <td headers="Model" class="gt_row gt_right">2.3094</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">1.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">0.1286</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.2083</td>
#> <td headers="BaseModel" class="gt_row gt_right">0.2075</td></tr>
#>   </tbody>
#>   
#> </table>
#> </div>
#> $base_model
#> Maximum Likelihood estimation
#> Nelder-Mead maximization, 3002 iterations
#> Return code 1: iteration limit exceeded 
#> Log-Likelihood: -1063.421 (7 free parameter(s))
#> Estimate(s): -7.350311 0.8911364 0.8436279 -0.370684 0.254564 0.7938497 -1.320848 
#> 
#> $comparison_model
#> Maximum Likelihood estimation
#> Nelder-Mead maximization, 917 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1062.266 (8 free parameter(s))
#> Estimate(s): -7.796102 0.9417162 0.8342337 -0.3870385 0.2592118 0.6860279 -1.326803 0.5023309 
#> 
#> $base_index
#> [1] 2
#> 
#> $comparison_index
#> [1] 1
#> 
#> $LL_base
#> [1] -1063.421
#> 
#> $LL_comparison
#> [1] -1062.266
#> 
#> $LR
#> [1] 2.309402
#> 
#> $LRdof
#> [1] 1
#> 
#> $LR_pvalue
#> [1] 0.1285935
#> 
#> $PseudoR2_base
#> [1] 0.2074693
#> 
#> $PseudoR2_comparison
#> [1] 0.2083299
#> 
#> $statistics
#> # A tibble: 8 × 3
#>   Statistic                 Model BaseModel
#>   <chr>                     <dbl>     <dbl>
#> 1 Log-likelihood        -1062.    -1063.   
#> 2 Number of parameters      8         7    
#> 3 AIC                    2141.     2141.   
#> 4 BIC                    2183.     2178.   
#> 5 LR Test Statistic         2.31     NA    
#> 6 LR degrees of freedom     1        NA    
#> 7 LR p-value                0.129    NA    
#> 8 McFadden's Pseudo R^2     0.208     0.208
#> 
#> $gtTable
#> <div id="pvteseuwbf" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#pvteseuwbf table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #pvteseuwbf thead, #pvteseuwbf tbody, #pvteseuwbf tfoot, #pvteseuwbf tr, #pvteseuwbf td, #pvteseuwbf th {
#>   border-style: none;
#> }
#> 
#> #pvteseuwbf p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #pvteseuwbf .gt_table {
#>   display: table;
#>   border-collapse: collapse;
#>   line-height: normal;
#>   margin-left: auto;
#>   margin-right: auto;
#>   color: #333333;
#>   font-size: 16px;
#>   font-weight: normal;
#>   font-style: normal;
#>   background-color: #FFFFFF;
#>   width: auto;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #A8A8A8;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #A8A8A8;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #pvteseuwbf .gt_title {
#>   color: #333333;
#>   font-size: 125%;
#>   font-weight: initial;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-color: #FFFFFF;
#>   border-bottom-width: 0;
#> }
#> 
#> #pvteseuwbf .gt_subtitle {
#>   color: #333333;
#>   font-size: 85%;
#>   font-weight: initial;
#>   padding-top: 3px;
#>   padding-bottom: 5px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-color: #FFFFFF;
#>   border-top-width: 0;
#> }
#> 
#> #pvteseuwbf .gt_heading {
#>   background-color: #FFFFFF;
#>   text-align: center;
#>   border-bottom-color: #FFFFFF;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_col_headings {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_col_heading {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: normal;
#>   text-transform: inherit;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: bottom;
#>   padding-top: 5px;
#>   padding-bottom: 6px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   overflow-x: hidden;
#> }
#> 
#> #pvteseuwbf .gt_column_spanner_outer {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: normal;
#>   text-transform: inherit;
#>   padding-top: 0;
#>   padding-bottom: 0;
#>   padding-left: 4px;
#>   padding-right: 4px;
#> }
#> 
#> #pvteseuwbf .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #pvteseuwbf .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #pvteseuwbf .gt_column_spanner {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   vertical-align: bottom;
#>   padding-top: 5px;
#>   padding-bottom: 5px;
#>   overflow-x: hidden;
#>   display: inline-block;
#>   width: 100%;
#> }
#> 
#> #pvteseuwbf .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #pvteseuwbf .gt_group_heading {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: middle;
#>   text-align: left;
#> }
#> 
#> #pvteseuwbf .gt_empty_group_heading {
#>   padding: 0.5px;
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   vertical-align: middle;
#> }
#> 
#> #pvteseuwbf .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #pvteseuwbf .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #pvteseuwbf .gt_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   margin: 10px;
#>   border-top-style: solid;
#>   border-top-width: 1px;
#>   border-top-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 1px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 1px;
#>   border-right-color: #D3D3D3;
#>   vertical-align: middle;
#>   overflow-x: hidden;
#> }
#> 
#> #pvteseuwbf .gt_stub {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-right-style: solid;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #pvteseuwbf .gt_stub_row_group {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   font-size: 100%;
#>   font-weight: initial;
#>   text-transform: inherit;
#>   border-right-style: solid;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   vertical-align: top;
#> }
#> 
#> #pvteseuwbf .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #pvteseuwbf .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #pvteseuwbf .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #pvteseuwbf .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #pvteseuwbf .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #pvteseuwbf .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #pvteseuwbf .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_footnotes {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   border-bottom-style: none;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #pvteseuwbf .gt_sourcenotes {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   border-bottom-style: none;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#>   border-left-style: none;
#>   border-left-width: 2px;
#>   border-left-color: #D3D3D3;
#>   border-right-style: none;
#>   border-right-width: 2px;
#>   border-right-color: #D3D3D3;
#> }
#> 
#> #pvteseuwbf .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #pvteseuwbf .gt_left {
#>   text-align: left;
#> }
#> 
#> #pvteseuwbf .gt_center {
#>   text-align: center;
#> }
#> 
#> #pvteseuwbf .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #pvteseuwbf .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #pvteseuwbf .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #pvteseuwbf .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #pvteseuwbf .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #pvteseuwbf .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #pvteseuwbf .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #pvteseuwbf .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #pvteseuwbf .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #pvteseuwbf .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #pvteseuwbf .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #pvteseuwbf .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #pvteseuwbf .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #pvteseuwbf div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
#>   height: 0px !important;
#> }
#> </style>
#>   <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
#>   <thead>
#>     <tr class="gt_heading">
#>       <td colspan="3" class="gt_heading gt_title gt_font_normal gt_bottom_border" style>Model Comparison Statistics</td>
#>     </tr>
#>     
#>     <tr class="gt_col_headings">
#>       <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Statistic">Statistic</th>
#>       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Model">Model</th>
#>       <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="BaseModel">BaseModel</th>
#>     </tr>
#>   </thead>
#>   <tbody class="gt_table_body">
#>     <tr><td headers="Statistic" class="gt_row gt_left">Log-likelihood</td>
#> <td headers="Model" class="gt_row gt_right">−1,062.2659</td>
#> <td headers="BaseModel" class="gt_row gt_right">−1,063.4206</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">Number of parameters</td>
#> <td headers="Model" class="gt_row gt_right">8.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">7.0000</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
#> <td headers="Model" class="gt_row gt_right">2,140.5317</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,140.8411</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
#> <td headers="Model" class="gt_row gt_right">2,183.0428</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,178.0384</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
#> <td headers="Model" class="gt_row gt_right">2.3094</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">1.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">0.1286</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.2083</td>
#> <td headers="BaseModel" class="gt_row gt_right">0.2075</td></tr>
#>   </tbody>
#>   
#> </table>
#> </div>
#> 
#> $latexTable
#> \begin{table}
#> 
#> \caption{Model Comparison Statistics}
#> \centering
#> \begin{tabular}[t]{lrr}
#> \toprule
#> Statistic & Model & BaseModel\\
#> \midrule
#> Log-likelihood & -1062.2659 & -1063.4206\\
#> Number of parameters & 8.0000 & 7.0000\\
#> AIC & 2140.5317 & 2140.8411\\
#> BIC & 2183.0428 & 2178.0384\\
#> LR Test Statistic & 2.3094 & NA\\
#> \addlinespace
#> LR degrees of freedom & 1.0000 & NA\\
#> LR p-value & 0.1286 & NA\\
#> McFadden's Pseudo R\textasciicircum{}2 & 0.2083 & 0.2075\\
#> \bottomrule
#> \end{tabular}
#> \end{table}
#> 
#> $htmlTable
#> <table class='table table-striped'>
#> <caption>Model Comparison Statistics</caption>
#>  <thead>
#>   <tr>
#>    <th style="text-align:left;"> Statistic </th>
#>    <th style="text-align:right;"> Model </th>
#>    <th style="text-align:right;"> BaseModel </th>
#>   </tr>
#>  </thead>
#> <tbody>
#>   <tr>
#>    <td style="text-align:left;"> Log-likelihood </td>
#>    <td style="text-align:right;"> -1062.2659 </td>
#>    <td style="text-align:right;"> -1063.4206 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> Number of parameters </td>
#>    <td style="text-align:right;"> 8.0000 </td>
#>    <td style="text-align:right;"> 7.0000 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> AIC </td>
#>    <td style="text-align:right;"> 2140.5317 </td>
#>    <td style="text-align:right;"> 2140.8411 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> BIC </td>
#>    <td style="text-align:right;"> 2183.0428 </td>
#>    <td style="text-align:right;"> 2178.0384 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR Test Statistic </td>
#>    <td style="text-align:right;"> 2.3094 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR degrees of freedom </td>
#>    <td style="text-align:right;"> 1.0000 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR p-value </td>
#>    <td style="text-align:right;"> 0.1286 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> McFadden's Pseudo R^2 </td>
#>    <td style="text-align:right;"> 0.2083 </td>
#>    <td style="text-align:right;"> 0.2075 </td>
#>   </tr>
#> </tbody>
#> </table>
#> 
```
