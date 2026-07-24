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
model from the supplied data and compares the fitted model against that
base model. When two fitted models are supplied, the function compares
the two fitted models directly and uses the model with the larger
log-likelihood as the base model for the LR test.

The function performs the following steps for the single-model workflow:

1.  Fits the base model, either a Poisson regression or another
    specified model.

2.  Computes the log-likelihoods of the provided model and the base
    model.

3.  Calculates the AIC and BIC for both models.

4.  Conducts a Likelihood Ratio test to compare the fitted model against
    the estimated base model when the fitted model has more parameters
    than the base model.

5.  Computes McFadden's Pseudo R^2.

For the two-model workflow, the function compares the two fitted models
directly. The model with the larger log-likelihood is treated as the
base model for the LR test, and the other model is treated as the
comparison model.

The Likelihood-Ratio test is computed as \$\$LR = -2 (LL\_{base} -
LL\_{comparison})\$\$. The test is chi-squared with degrees of freedom
\$\$dof = \|N\_{base \\ params} - N\_{comparison \\ params}\|\$\$. The
AIC is calculated as \$\$AIC = -2 \cdot LL + 2 \cdot nparam\$\$, and the
BIC is calculated as \$\$BIC = -2 \cdot LL + nparam \cdot \log(n)\$\$.

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
#> <div id="nprbjsesge" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#nprbjsesge table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #nprbjsesge thead, #nprbjsesge tbody, #nprbjsesge tfoot, #nprbjsesge tr, #nprbjsesge td, #nprbjsesge th {
#>   border-style: none;
#> }
#> 
#> #nprbjsesge p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #nprbjsesge .gt_table {
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
#> #nprbjsesge .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #nprbjsesge .gt_title {
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
#> #nprbjsesge .gt_subtitle {
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
#> #nprbjsesge .gt_heading {
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
#> #nprbjsesge .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #nprbjsesge .gt_col_headings {
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
#> #nprbjsesge .gt_col_heading {
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
#> #nprbjsesge .gt_column_spanner_outer {
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
#> #nprbjsesge .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #nprbjsesge .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #nprbjsesge .gt_column_spanner {
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
#> #nprbjsesge .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #nprbjsesge .gt_group_heading {
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
#> #nprbjsesge .gt_empty_group_heading {
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
#> #nprbjsesge .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #nprbjsesge .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #nprbjsesge .gt_row {
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
#> #nprbjsesge .gt_stub {
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
#> #nprbjsesge .gt_stub_row_group {
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
#> #nprbjsesge .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #nprbjsesge .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #nprbjsesge .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #nprbjsesge .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #nprbjsesge .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #nprbjsesge .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #nprbjsesge .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #nprbjsesge .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #nprbjsesge .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #nprbjsesge .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #nprbjsesge .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #nprbjsesge .gt_footnotes {
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
#> #nprbjsesge .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #nprbjsesge .gt_sourcenotes {
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
#> #nprbjsesge .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #nprbjsesge .gt_left {
#>   text-align: left;
#> }
#> 
#> #nprbjsesge .gt_center {
#>   text-align: center;
#> }
#> 
#> #nprbjsesge .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #nprbjsesge .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #nprbjsesge .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #nprbjsesge .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #nprbjsesge .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #nprbjsesge .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #nprbjsesge .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #nprbjsesge .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #nprbjsesge .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #nprbjsesge .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #nprbjsesge .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #nprbjsesge .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #nprbjsesge .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #nprbjsesge div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
#> <td headers="Model" class="gt_row gt_right">−1,062.2664</td>
#> <td headers="BaseModel" class="gt_row gt_right">−1,071.2672</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">Number of parameters</td>
#> <td headers="Model" class="gt_row gt_right">8.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">6.0000</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
#> <td headers="Model" class="gt_row gt_right">2,140.5327</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,154.5345</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
#> <td headers="Model" class="gt_row gt_right">2,183.0438</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,186.4178</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
#> <td headers="Model" class="gt_row gt_right">18.0018</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">2.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">0.0001</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.3029</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>   </tbody>
#>   
#> </table>
#> </div>
#> $model
#> Maximum Likelihood estimation
#> Nelder-Mead maximization, 703 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1062.266 (8 free parameter(s))
#> Estimate(s): -7.79761 0.9418564 0.8335215 -0.3860324 0.2577976 0.6863257 -1.327278 0.4993697 
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
#> [1] 18.00176
#> 
#> $LRdof
#> [1] 2
#> 
#> $AIC
#> [1] 2140.533
#> 
#> $AICbase
#> [1] 2154.534
#> 
#> $BIC
#> [1] 2183.044
#> 
#> $BICbase
#> [1] 2186.418
#> 
#> $LR_pvalue
#> [1] 0.000123301
#> 
#> $PseudoR2
#> [1] 0.3028969
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
#> 8 McFadden's Pseudo R^2     0.303        NA 
#> 
#> $gtTable
#> <div id="sovsoqmgyn" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#sovsoqmgyn table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #sovsoqmgyn thead, #sovsoqmgyn tbody, #sovsoqmgyn tfoot, #sovsoqmgyn tr, #sovsoqmgyn td, #sovsoqmgyn th {
#>   border-style: none;
#> }
#> 
#> #sovsoqmgyn p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #sovsoqmgyn .gt_table {
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
#> #sovsoqmgyn .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #sovsoqmgyn .gt_title {
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
#> #sovsoqmgyn .gt_subtitle {
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
#> #sovsoqmgyn .gt_heading {
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
#> #sovsoqmgyn .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #sovsoqmgyn .gt_col_headings {
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
#> #sovsoqmgyn .gt_col_heading {
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
#> #sovsoqmgyn .gt_column_spanner_outer {
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
#> #sovsoqmgyn .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #sovsoqmgyn .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #sovsoqmgyn .gt_column_spanner {
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
#> #sovsoqmgyn .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #sovsoqmgyn .gt_group_heading {
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
#> #sovsoqmgyn .gt_empty_group_heading {
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
#> #sovsoqmgyn .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #sovsoqmgyn .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #sovsoqmgyn .gt_row {
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
#> #sovsoqmgyn .gt_stub {
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
#> #sovsoqmgyn .gt_stub_row_group {
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
#> #sovsoqmgyn .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #sovsoqmgyn .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #sovsoqmgyn .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #sovsoqmgyn .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #sovsoqmgyn .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #sovsoqmgyn .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #sovsoqmgyn .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #sovsoqmgyn .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #sovsoqmgyn .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #sovsoqmgyn .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #sovsoqmgyn .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #sovsoqmgyn .gt_footnotes {
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
#> #sovsoqmgyn .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #sovsoqmgyn .gt_sourcenotes {
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
#> #sovsoqmgyn .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #sovsoqmgyn .gt_left {
#>   text-align: left;
#> }
#> 
#> #sovsoqmgyn .gt_center {
#>   text-align: center;
#> }
#> 
#> #sovsoqmgyn .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #sovsoqmgyn .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #sovsoqmgyn .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #sovsoqmgyn .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #sovsoqmgyn .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #sovsoqmgyn .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #sovsoqmgyn .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #sovsoqmgyn .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #sovsoqmgyn .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #sovsoqmgyn .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #sovsoqmgyn .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #sovsoqmgyn .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #sovsoqmgyn .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #sovsoqmgyn div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
#> <td headers="Model" class="gt_row gt_right">−1,062.2664</td>
#> <td headers="BaseModel" class="gt_row gt_right">−1,071.2672</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">Number of parameters</td>
#> <td headers="Model" class="gt_row gt_right">8.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">6.0000</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
#> <td headers="Model" class="gt_row gt_right">2,140.5327</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,154.5345</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
#> <td headers="Model" class="gt_row gt_right">2,183.0438</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,186.4178</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
#> <td headers="Model" class="gt_row gt_right">18.0018</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">2.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">0.0001</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.3029</td>
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
#> Log-likelihood & -1062.2664 & -1071.267\\
#> Number of parameters & 8.0000 & 6.000\\
#> AIC & 2140.5327 & 2154.535\\
#> BIC & 2183.0438 & 2186.418\\
#> LR Test Statistic & 18.0018 & NA\\
#> \addlinespace
#> LR degrees of freedom & 2.0000 & NA\\
#> LR p-value & 0.0001 & NA\\
#> McFadden's Pseudo R\textasciicircum{}2 & 0.3029 & NA\\
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
#>    <td style="text-align:right;"> -1062.2664 </td>
#>    <td style="text-align:right;"> -1071.267 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> Number of parameters </td>
#>    <td style="text-align:right;"> 8.0000 </td>
#>    <td style="text-align:right;"> 6.000 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> AIC </td>
#>    <td style="text-align:right;"> 2140.5327 </td>
#>    <td style="text-align:right;"> 2154.535 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> BIC </td>
#>    <td style="text-align:right;"> 2183.0438 </td>
#>    <td style="text-align:right;"> 2186.418 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR Test Statistic </td>
#>    <td style="text-align:right;"> 18.0018 </td>
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
#>    <td style="text-align:right;"> 0.3029 </td>
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

regCompTest(nbp.base, data = washington_roads, model2 = nb2.base,
            print = TRUE)
#> <div id="gtrogrnzut" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#gtrogrnzut table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #gtrogrnzut thead, #gtrogrnzut tbody, #gtrogrnzut tfoot, #gtrogrnzut tr, #gtrogrnzut td, #gtrogrnzut th {
#>   border-style: none;
#> }
#> 
#> #gtrogrnzut p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #gtrogrnzut .gt_table {
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
#> #gtrogrnzut .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #gtrogrnzut .gt_title {
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
#> #gtrogrnzut .gt_subtitle {
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
#> #gtrogrnzut .gt_heading {
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
#> #gtrogrnzut .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #gtrogrnzut .gt_col_headings {
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
#> #gtrogrnzut .gt_col_heading {
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
#> #gtrogrnzut .gt_column_spanner_outer {
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
#> #gtrogrnzut .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #gtrogrnzut .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #gtrogrnzut .gt_column_spanner {
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
#> #gtrogrnzut .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #gtrogrnzut .gt_group_heading {
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
#> #gtrogrnzut .gt_empty_group_heading {
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
#> #gtrogrnzut .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #gtrogrnzut .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #gtrogrnzut .gt_row {
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
#> #gtrogrnzut .gt_stub {
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
#> #gtrogrnzut .gt_stub_row_group {
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
#> #gtrogrnzut .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #gtrogrnzut .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #gtrogrnzut .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #gtrogrnzut .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #gtrogrnzut .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #gtrogrnzut .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #gtrogrnzut .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #gtrogrnzut .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #gtrogrnzut .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #gtrogrnzut .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #gtrogrnzut .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #gtrogrnzut .gt_footnotes {
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
#> #gtrogrnzut .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #gtrogrnzut .gt_sourcenotes {
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
#> #gtrogrnzut .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #gtrogrnzut .gt_left {
#>   text-align: left;
#> }
#> 
#> #gtrogrnzut .gt_center {
#>   text-align: center;
#> }
#> 
#> #gtrogrnzut .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #gtrogrnzut .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #gtrogrnzut .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #gtrogrnzut .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #gtrogrnzut .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #gtrogrnzut .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #gtrogrnzut .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #gtrogrnzut .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #gtrogrnzut .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #gtrogrnzut .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #gtrogrnzut .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #gtrogrnzut .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #gtrogrnzut .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #gtrogrnzut div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
#> <td headers="Model" class="gt_row gt_right">−1,062.2664</td>
#> <td headers="BaseModel" class="gt_row gt_right">−1,063.6726</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">Number of parameters</td>
#> <td headers="Model" class="gt_row gt_right">8.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">7.0000</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
#> <td headers="Model" class="gt_row gt_right">2,140.5327</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,141.3452</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
#> <td headers="Model" class="gt_row gt_right">2,183.0438</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,178.5424</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
#> <td headers="Model" class="gt_row gt_right">2.8125</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">1.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">0.0935</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.3029</td>
#> <td headers="BaseModel" class="gt_row gt_right">0.3020</td></tr>
#>   </tbody>
#>   
#> </table>
#> </div>
#> $base_model
#> Maximum Likelihood estimation
#> Nelder-Mead maximization, 703 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1062.266 (8 free parameter(s))
#> Estimate(s): -7.79761 0.9418564 0.8335215 -0.3860324 0.2577976 0.6863257 -1.327278 0.4993697 
#> 
#> $comparison_model
#> Maximum Likelihood estimation
#> Nelder-Mead maximization, 890 iterations
#> Return code 0: successful convergence 
#> Log-Likelihood: -1063.673 (7 free parameter(s))
#> Estimate(s): -7.150729 0.8663599 0.8358287 -0.4012173 0.230993 0.8275592 -1.32798 
#> 
#> $base_index
#> [1] 1
#> 
#> $comparison_index
#> [1] 2
#> 
#> $LL_base
#> [1] -1062.266
#> 
#> $LL_comparison
#> [1] -1063.673
#> 
#> $LR
#> [1] 2.812464
#> 
#> $LRdof
#> [1] 1
#> 
#> $LR_pvalue
#> [1] 0.0935346
#> 
#> $PseudoR2_base
#> [1] 0.3028969
#> 
#> $PseudoR2_comparison
#> [1] 0.301974
#> 
#> $statistics
#> # A tibble: 8 × 3
#>   Statistic                  Model BaseModel
#>   <chr>                      <dbl>     <dbl>
#> 1 Log-likelihood        -1062.     -1064.   
#> 2 Number of parameters      8          7    
#> 3 AIC                    2141.      2141.   
#> 4 BIC                    2183.      2179.   
#> 5 LR Test Statistic         2.81      NA    
#> 6 LR degrees of freedom     1         NA    
#> 7 LR p-value                0.0935    NA    
#> 8 McFadden's Pseudo R^2     0.303      0.302
#> 
#> $gtTable
#> <div id="utdkelafqq" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#utdkelafqq table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #utdkelafqq thead, #utdkelafqq tbody, #utdkelafqq tfoot, #utdkelafqq tr, #utdkelafqq td, #utdkelafqq th {
#>   border-style: none;
#> }
#> 
#> #utdkelafqq p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #utdkelafqq .gt_table {
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
#> #utdkelafqq .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #utdkelafqq .gt_title {
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
#> #utdkelafqq .gt_subtitle {
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
#> #utdkelafqq .gt_heading {
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
#> #utdkelafqq .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #utdkelafqq .gt_col_headings {
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
#> #utdkelafqq .gt_col_heading {
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
#> #utdkelafqq .gt_column_spanner_outer {
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
#> #utdkelafqq .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #utdkelafqq .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #utdkelafqq .gt_column_spanner {
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
#> #utdkelafqq .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #utdkelafqq .gt_group_heading {
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
#> #utdkelafqq .gt_empty_group_heading {
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
#> #utdkelafqq .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #utdkelafqq .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #utdkelafqq .gt_row {
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
#> #utdkelafqq .gt_stub {
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
#> #utdkelafqq .gt_stub_row_group {
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
#> #utdkelafqq .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #utdkelafqq .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #utdkelafqq .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #utdkelafqq .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #utdkelafqq .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #utdkelafqq .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #utdkelafqq .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #utdkelafqq .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #utdkelafqq .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #utdkelafqq .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #utdkelafqq .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #utdkelafqq .gt_footnotes {
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
#> #utdkelafqq .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #utdkelafqq .gt_sourcenotes {
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
#> #utdkelafqq .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #utdkelafqq .gt_left {
#>   text-align: left;
#> }
#> 
#> #utdkelafqq .gt_center {
#>   text-align: center;
#> }
#> 
#> #utdkelafqq .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #utdkelafqq .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #utdkelafqq .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #utdkelafqq .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #utdkelafqq .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #utdkelafqq .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #utdkelafqq .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #utdkelafqq .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #utdkelafqq .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #utdkelafqq .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #utdkelafqq .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #utdkelafqq .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #utdkelafqq .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #utdkelafqq div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
#> <td headers="Model" class="gt_row gt_right">−1,062.2664</td>
#> <td headers="BaseModel" class="gt_row gt_right">−1,063.6726</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">Number of parameters</td>
#> <td headers="Model" class="gt_row gt_right">8.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">7.0000</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">AIC</td>
#> <td headers="Model" class="gt_row gt_right">2,140.5327</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,141.3452</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">BIC</td>
#> <td headers="Model" class="gt_row gt_right">2,183.0438</td>
#> <td headers="BaseModel" class="gt_row gt_right">2,178.5424</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR Test Statistic</td>
#> <td headers="Model" class="gt_row gt_right">2.8125</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">1.0000</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">0.0935</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.3029</td>
#> <td headers="BaseModel" class="gt_row gt_right">0.3020</td></tr>
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
#> Log-likelihood & -1062.2664 & -1063.673\\
#> Number of parameters & 8.0000 & 7.000\\
#> AIC & 2140.5327 & 2141.345\\
#> BIC & 2183.0438 & 2178.542\\
#> LR Test Statistic & 2.8125 & NA\\
#> \addlinespace
#> LR degrees of freedom & 1.0000 & NA\\
#> LR p-value & 0.0935 & NA\\
#> McFadden's Pseudo R\textasciicircum{}2 & 0.3029 & 0.302\\
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
#>    <td style="text-align:right;"> -1062.2664 </td>
#>    <td style="text-align:right;"> -1063.673 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> Number of parameters </td>
#>    <td style="text-align:right;"> 8.0000 </td>
#>    <td style="text-align:right;"> 7.000 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> AIC </td>
#>    <td style="text-align:right;"> 2140.5327 </td>
#>    <td style="text-align:right;"> 2141.345 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> BIC </td>
#>    <td style="text-align:right;"> 2183.0438 </td>
#>    <td style="text-align:right;"> 2178.542 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR Test Statistic </td>
#>    <td style="text-align:right;"> 2.8125 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR degrees of freedom </td>
#>    <td style="text-align:right;"> 1.0000 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR p-value </td>
#>    <td style="text-align:right;"> 0.0935 </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> McFadden's Pseudo R^2 </td>
#>    <td style="text-align:right;"> 0.3029 </td>
#>    <td style="text-align:right;"> 0.302 </td>
#>   </tr>
#> </tbody>
#> </table>
#> 
```
