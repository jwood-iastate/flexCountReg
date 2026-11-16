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
  or "Poisson"). Default is "Poisson". It is used for the LR test of a
  single model and, in both workflows, as the baseline of McFadden's
  Pseudo R^2.

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

The LR test is reported only when the restricted model is nested in the
larger one with a regular restriction: the same estimator, zero and
underreporting components and random parameters, the restricted model's
parameters a subset of the larger model's, and the same family or a
nested family (Poisson in NB1, NB2, PLN, PIG1, PIG2, PIG, GP1, GP2 and
COM; NB1 and NB2 in NBP; PL in PLG and PLL). When the restriction puts
one parameter on the boundary of its space (e.g. Poisson in NB2), the
p-value uses a 50:50 mixture of chi-square distributions with df - 1 and
df degrees of freedom. Two fitted models must use the same observations,
outcome, offsets and weights. Otherwise the LR statistic and p-value are
`NA` (`LR_type` gives the kind of restriction); AIC and BIC are always
reported.

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
    the log-likelihood of the base model (by default an intercept-only
    Poisson model, the base count model, so that the pseudo R^2 is
    comparable across models and distributions).

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
#> No LR test: the POISSON base model is not nested in the fitted model with a regular restriction (the POISSON family is not nested in the NBP family with a regular restriction). Compare the models with AIC and BIC.
#> <div id="apeqfnqalm" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#apeqfnqalm table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #apeqfnqalm thead, #apeqfnqalm tbody, #apeqfnqalm tfoot, #apeqfnqalm tr, #apeqfnqalm td, #apeqfnqalm th {
#>   border-style: none;
#> }
#> 
#> #apeqfnqalm p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #apeqfnqalm .gt_table {
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
#> #apeqfnqalm .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #apeqfnqalm .gt_title {
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
#> #apeqfnqalm .gt_subtitle {
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
#> #apeqfnqalm .gt_heading {
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
#> #apeqfnqalm .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #apeqfnqalm .gt_col_headings {
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
#> #apeqfnqalm .gt_col_heading {
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
#> #apeqfnqalm .gt_column_spanner_outer {
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
#> #apeqfnqalm .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #apeqfnqalm .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #apeqfnqalm .gt_column_spanner {
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
#> #apeqfnqalm .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #apeqfnqalm .gt_group_heading {
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
#> #apeqfnqalm .gt_empty_group_heading {
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
#> #apeqfnqalm .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #apeqfnqalm .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #apeqfnqalm .gt_row {
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
#> #apeqfnqalm .gt_stub {
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
#> #apeqfnqalm .gt_stub_row_group {
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
#> #apeqfnqalm .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #apeqfnqalm .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #apeqfnqalm .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #apeqfnqalm .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #apeqfnqalm .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #apeqfnqalm .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #apeqfnqalm .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #apeqfnqalm .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #apeqfnqalm .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #apeqfnqalm .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #apeqfnqalm .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #apeqfnqalm .gt_footnotes {
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
#> #apeqfnqalm .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #apeqfnqalm .gt_sourcenotes {
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
#> #apeqfnqalm .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #apeqfnqalm .gt_left {
#>   text-align: left;
#> }
#> 
#> #apeqfnqalm .gt_center {
#>   text-align: center;
#> }
#> 
#> #apeqfnqalm .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #apeqfnqalm .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #apeqfnqalm .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #apeqfnqalm .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #apeqfnqalm .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #apeqfnqalm .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #apeqfnqalm .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #apeqfnqalm .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #apeqfnqalm .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #apeqfnqalm .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #apeqfnqalm .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #apeqfnqalm .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #apeqfnqalm .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #apeqfnqalm div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
#> <td headers="Model" class="gt_row gt_right">NA</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">NA</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">NA</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.0084</td>
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
#> [1] NA
#> 
#> $LRdof
#> [1] NA
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
#> [1] NA
#> 
#> $LR_type
#> [1] NA
#> 
#> $PseudoR2
#> [1] 0.008402545
#> 
#> $statistics
#> # A tibble: 8 × 3
#>   Statistic                  Model BaseModel
#>   <chr>                      <dbl>     <dbl>
#> 1 Log-likelihood        -1062.        -1071.
#> 2 Number of parameters      8             6 
#> 3 AIC                    2141.         2155.
#> 4 BIC                    2183.         2186.
#> 5 LR Test Statistic        NA            NA 
#> 6 LR degrees of freedom    NA            NA 
#> 7 LR p-value               NA            NA 
#> 8 McFadden's Pseudo R^2     0.0084       NA 
#> 
#> $gtTable
#> <div id="gqshkgaaeg" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#gqshkgaaeg table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #gqshkgaaeg thead, #gqshkgaaeg tbody, #gqshkgaaeg tfoot, #gqshkgaaeg tr, #gqshkgaaeg td, #gqshkgaaeg th {
#>   border-style: none;
#> }
#> 
#> #gqshkgaaeg p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #gqshkgaaeg .gt_table {
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
#> #gqshkgaaeg .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #gqshkgaaeg .gt_title {
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
#> #gqshkgaaeg .gt_subtitle {
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
#> #gqshkgaaeg .gt_heading {
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
#> #gqshkgaaeg .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #gqshkgaaeg .gt_col_headings {
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
#> #gqshkgaaeg .gt_col_heading {
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
#> #gqshkgaaeg .gt_column_spanner_outer {
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
#> #gqshkgaaeg .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #gqshkgaaeg .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #gqshkgaaeg .gt_column_spanner {
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
#> #gqshkgaaeg .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #gqshkgaaeg .gt_group_heading {
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
#> #gqshkgaaeg .gt_empty_group_heading {
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
#> #gqshkgaaeg .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #gqshkgaaeg .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #gqshkgaaeg .gt_row {
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
#> #gqshkgaaeg .gt_stub {
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
#> #gqshkgaaeg .gt_stub_row_group {
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
#> #gqshkgaaeg .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #gqshkgaaeg .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #gqshkgaaeg .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #gqshkgaaeg .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #gqshkgaaeg .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #gqshkgaaeg .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #gqshkgaaeg .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #gqshkgaaeg .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #gqshkgaaeg .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #gqshkgaaeg .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #gqshkgaaeg .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #gqshkgaaeg .gt_footnotes {
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
#> #gqshkgaaeg .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #gqshkgaaeg .gt_sourcenotes {
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
#> #gqshkgaaeg .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #gqshkgaaeg .gt_left {
#>   text-align: left;
#> }
#> 
#> #gqshkgaaeg .gt_center {
#>   text-align: center;
#> }
#> 
#> #gqshkgaaeg .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #gqshkgaaeg .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #gqshkgaaeg .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #gqshkgaaeg .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #gqshkgaaeg .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #gqshkgaaeg .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #gqshkgaaeg .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #gqshkgaaeg .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #gqshkgaaeg .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #gqshkgaaeg .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #gqshkgaaeg .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #gqshkgaaeg .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #gqshkgaaeg .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #gqshkgaaeg div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
#> <td headers="Model" class="gt_row gt_right">NA</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR degrees of freedom</td>
#> <td headers="Model" class="gt_row gt_right">NA</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">LR p-value</td>
#> <td headers="Model" class="gt_row gt_right">NA</td>
#> <td headers="BaseModel" class="gt_row gt_right">NA</td></tr>
#>     <tr><td headers="Statistic" class="gt_row gt_left">McFadden's Pseudo R^2</td>
#> <td headers="Model" class="gt_row gt_right">0.0084</td>
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
#> LR Test Statistic & NA & NA\\
#> \addlinespace
#> LR degrees of freedom & NA & NA\\
#> LR p-value & NA & NA\\
#> McFadden's Pseudo R\textasciicircum{}2 & 0.0084 & NA\\
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
#>    <td style="text-align:right;"> NA </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR degrees of freedom </td>
#>    <td style="text-align:right;"> NA </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> LR p-value </td>
#>    <td style="text-align:right;"> NA </td>
#>    <td style="text-align:right;"> NA </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;"> McFadden's Pseudo R^2 </td>
#>    <td style="text-align:right;"> 0.0084 </td>
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
#> <div id="xxnbvxjcbz" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#xxnbvxjcbz table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #xxnbvxjcbz thead, #xxnbvxjcbz tbody, #xxnbvxjcbz tfoot, #xxnbvxjcbz tr, #xxnbvxjcbz td, #xxnbvxjcbz th {
#>   border-style: none;
#> }
#> 
#> #xxnbvxjcbz p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #xxnbvxjcbz .gt_table {
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
#> #xxnbvxjcbz .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #xxnbvxjcbz .gt_title {
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
#> #xxnbvxjcbz .gt_subtitle {
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
#> #xxnbvxjcbz .gt_heading {
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
#> #xxnbvxjcbz .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #xxnbvxjcbz .gt_col_headings {
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
#> #xxnbvxjcbz .gt_col_heading {
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
#> #xxnbvxjcbz .gt_column_spanner_outer {
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
#> #xxnbvxjcbz .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #xxnbvxjcbz .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #xxnbvxjcbz .gt_column_spanner {
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
#> #xxnbvxjcbz .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #xxnbvxjcbz .gt_group_heading {
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
#> #xxnbvxjcbz .gt_empty_group_heading {
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
#> #xxnbvxjcbz .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #xxnbvxjcbz .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #xxnbvxjcbz .gt_row {
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
#> #xxnbvxjcbz .gt_stub {
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
#> #xxnbvxjcbz .gt_stub_row_group {
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
#> #xxnbvxjcbz .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #xxnbvxjcbz .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #xxnbvxjcbz .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #xxnbvxjcbz .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #xxnbvxjcbz .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #xxnbvxjcbz .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #xxnbvxjcbz .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #xxnbvxjcbz .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #xxnbvxjcbz .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #xxnbvxjcbz .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #xxnbvxjcbz .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #xxnbvxjcbz .gt_footnotes {
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
#> #xxnbvxjcbz .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #xxnbvxjcbz .gt_sourcenotes {
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
#> #xxnbvxjcbz .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #xxnbvxjcbz .gt_left {
#>   text-align: left;
#> }
#> 
#> #xxnbvxjcbz .gt_center {
#>   text-align: center;
#> }
#> 
#> #xxnbvxjcbz .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #xxnbvxjcbz .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #xxnbvxjcbz .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #xxnbvxjcbz .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #xxnbvxjcbz .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #xxnbvxjcbz .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #xxnbvxjcbz .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #xxnbvxjcbz .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #xxnbvxjcbz .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #xxnbvxjcbz .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #xxnbvxjcbz .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #xxnbvxjcbz .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #xxnbvxjcbz .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #xxnbvxjcbz div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
#> <td headers="Model" class="gt_row gt_right">0.3029</td>
#> <td headers="BaseModel" class="gt_row gt_right">0.3021</td></tr>
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
#> $LR_type
#> [1] "interior"
#> 
#> $PseudoR2_base
#> [1] 0.3021394
#> 
#> $PseudoR2_comparison
#> [1] 0.3028972
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
#> 8 McFadden's Pseudo R^2     0.303     0.302
#> 
#> $gtTable
#> <div id="utvaysqzap" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
#>   <style>#utvaysqzap table {
#>   font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
#>   -webkit-font-smoothing: antialiased;
#>   -moz-osx-font-smoothing: grayscale;
#> }
#> 
#> #utvaysqzap thead, #utvaysqzap tbody, #utvaysqzap tfoot, #utvaysqzap tr, #utvaysqzap td, #utvaysqzap th {
#>   border-style: none;
#> }
#> 
#> #utvaysqzap p {
#>   margin: 0;
#>   padding: 0;
#> }
#> 
#> #utvaysqzap .gt_table {
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
#> #utvaysqzap .gt_caption {
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#> }
#> 
#> #utvaysqzap .gt_title {
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
#> #utvaysqzap .gt_subtitle {
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
#> #utvaysqzap .gt_heading {
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
#> #utvaysqzap .gt_bottom_border {
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #utvaysqzap .gt_col_headings {
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
#> #utvaysqzap .gt_col_heading {
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
#> #utvaysqzap .gt_column_spanner_outer {
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
#> #utvaysqzap .gt_column_spanner_outer:first-child {
#>   padding-left: 0;
#> }
#> 
#> #utvaysqzap .gt_column_spanner_outer:last-child {
#>   padding-right: 0;
#> }
#> 
#> #utvaysqzap .gt_column_spanner {
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
#> #utvaysqzap .gt_spanner_row {
#>   border-bottom-style: hidden;
#> }
#> 
#> #utvaysqzap .gt_group_heading {
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
#> #utvaysqzap .gt_empty_group_heading {
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
#> #utvaysqzap .gt_from_md > :first-child {
#>   margin-top: 0;
#> }
#> 
#> #utvaysqzap .gt_from_md > :last-child {
#>   margin-bottom: 0;
#> }
#> 
#> #utvaysqzap .gt_row {
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
#> #utvaysqzap .gt_stub {
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
#> #utvaysqzap .gt_stub_row_group {
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
#> #utvaysqzap .gt_row_group_first td {
#>   border-top-width: 2px;
#> }
#> 
#> #utvaysqzap .gt_row_group_first th {
#>   border-top-width: 2px;
#> }
#> 
#> #utvaysqzap .gt_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #utvaysqzap .gt_first_summary_row {
#>   border-top-style: solid;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #utvaysqzap .gt_first_summary_row.thick {
#>   border-top-width: 2px;
#> }
#> 
#> #utvaysqzap .gt_last_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #utvaysqzap .gt_grand_summary_row {
#>   color: #333333;
#>   background-color: #FFFFFF;
#>   text-transform: inherit;
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #utvaysqzap .gt_first_grand_summary_row {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-top-style: double;
#>   border-top-width: 6px;
#>   border-top-color: #D3D3D3;
#> }
#> 
#> #utvaysqzap .gt_last_grand_summary_row_top {
#>   padding-top: 8px;
#>   padding-bottom: 8px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#>   border-bottom-style: double;
#>   border-bottom-width: 6px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #utvaysqzap .gt_striped {
#>   background-color: rgba(128, 128, 128, 0.05);
#> }
#> 
#> #utvaysqzap .gt_table_body {
#>   border-top-style: solid;
#>   border-top-width: 2px;
#>   border-top-color: #D3D3D3;
#>   border-bottom-style: solid;
#>   border-bottom-width: 2px;
#>   border-bottom-color: #D3D3D3;
#> }
#> 
#> #utvaysqzap .gt_footnotes {
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
#> #utvaysqzap .gt_footnote {
#>   margin: 0px;
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #utvaysqzap .gt_sourcenotes {
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
#> #utvaysqzap .gt_sourcenote {
#>   font-size: 90%;
#>   padding-top: 4px;
#>   padding-bottom: 4px;
#>   padding-left: 5px;
#>   padding-right: 5px;
#> }
#> 
#> #utvaysqzap .gt_left {
#>   text-align: left;
#> }
#> 
#> #utvaysqzap .gt_center {
#>   text-align: center;
#> }
#> 
#> #utvaysqzap .gt_right {
#>   text-align: right;
#>   font-variant-numeric: tabular-nums;
#> }
#> 
#> #utvaysqzap .gt_font_normal {
#>   font-weight: normal;
#> }
#> 
#> #utvaysqzap .gt_font_bold {
#>   font-weight: bold;
#> }
#> 
#> #utvaysqzap .gt_font_italic {
#>   font-style: italic;
#> }
#> 
#> #utvaysqzap .gt_super {
#>   font-size: 65%;
#> }
#> 
#> #utvaysqzap .gt_footnote_marks {
#>   font-size: 75%;
#>   vertical-align: 0.4em;
#>   position: initial;
#> }
#> 
#> #utvaysqzap .gt_asterisk {
#>   font-size: 100%;
#>   vertical-align: 0;
#> }
#> 
#> #utvaysqzap .gt_indent_1 {
#>   text-indent: 5px;
#> }
#> 
#> #utvaysqzap .gt_indent_2 {
#>   text-indent: 10px;
#> }
#> 
#> #utvaysqzap .gt_indent_3 {
#>   text-indent: 15px;
#> }
#> 
#> #utvaysqzap .gt_indent_4 {
#>   text-indent: 20px;
#> }
#> 
#> #utvaysqzap .gt_indent_5 {
#>   text-indent: 25px;
#> }
#> 
#> #utvaysqzap .katex-display {
#>   display: inline-flex !important;
#>   margin-bottom: 0.75em !important;
#> }
#> 
#> #utvaysqzap div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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
#> <td headers="Model" class="gt_row gt_right">0.3029</td>
#> <td headers="BaseModel" class="gt_row gt_right">0.3021</td></tr>
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
#> McFadden's Pseudo R\textasciicircum{}2 & 0.3029 & 0.3021\\
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
#>    <td style="text-align:right;"> 0.3029 </td>
#>    <td style="text-align:right;"> 0.3021 </td>
#>   </tr>
#> </tbody>
#> </table>
#> 
```
