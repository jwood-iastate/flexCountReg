# Interactive analysis app

`launch_app()` opens an interactive interface to flexCountReg in a web
browser, with a tab for each step: load data (a bundled example or a CSV
file) and create new variables (logs, powers, standardized values,
indicators, products, ratios, categories, factors); specify and estimate
count regression
([`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)),
random-parameter
([`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md)),
random-effects panel
([`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
[`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md),
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)),
finite-mixture and latent-class
([`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md),
[`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md))
models; check their numerical status; diagnose them (rootograms,
randomized quantile residuals, CURE plots); compare fits; compute
predictions, confidence or prediction intervals and tolerance intervals
for a table of scenarios; and download the analysis.

## Usage

``` r
launch_app(..., max_upload_mb = 50)

app_object(max_upload_mb = 50)
```

## Arguments

- ...:

  arguments passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html), such
  as `port` or `launch.browser`.

- max_upload_mb:

  largest CSV file that can be uploaded, in megabytes.

## Value

`launch_app()` runs the app (it returns when the app is closed).
`app_object()` returns the app as a `shiny.appobj`, for deployment or
testing.

## Details

Everything the app runs is compiled into R code that calls the package's
functions, evaluated from that code, and appended to a session script:
data loading, new variables, fits, diagnostics, predictions and
comparisons, in order. The "R code" panel at the bottom shows the
session script (and the code of the model being built), and the Export
tab downloads it as `analysis.R`, which reproduces the whole session
without the app. Column names and typed values enter the code only as
names, numbers or strings, never as code. Models are estimated in
background R processes, so the interface stays responsive and long fits
can be canceled. Data, fits and jobs belong to one browser session and
are removed when it ends; uploaded data are read as CSV only and stay on
your computer.

The app needs the packages shiny, bslib, DT, callr, jsonlite and zip
(suggested, not required, by flexCountReg); `launch_app()` offers to
install any that are missing.

## Examples

``` r
if (interactive()) {
  launch_app()
}
```
