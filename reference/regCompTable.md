# Create a Table Comparing Regression Models with AIC, BIC, and McFadden's Pseudo-R-Squared

This function creates tables comparing the flexCountReg package models
supplied to the function.

## Usage

``` r
regCompTable(
  models,
  coefs = TRUE,
  AIC = TRUE,
  BIC = TRUE,
  RSquare = TRUE,
  tableType = "tibble",
  digits = 3,
  basemodel = "Poisson",
  variables = FALSE
)
```

## Arguments

- models:

  A named list of fitted flexCountReg model objects. This must include 2
  or more models.

- coefs:

  A logical. The default value `TRUE` indicates that the coefficients
  from the models should be included in the table of comparisons.

- AIC:

  A logical. The default value `TRUE` indicates that AIC values for the
  models should be included.

- BIC:

  A logical. The default value `TRUE` indicates that BIC values for the
  models should be included.

- RSquare:

  A logical. The default value `TRUE` indicates that the McFadden's
  Pseudo-R-Squared statistic should be included. By default it compares
  every model against the same intercept-only Poisson model (the base
  count model), so the values are comparable across models; see
  `basemodel` and `variables`.

- tableType:

  The type of table format to return. Options include "tibble" for
  returning the table as a tibble, "gt" for a
  [gt](https://gt.rstudio.com/reference/gt.html) table object, or
  "latex" for a latex table. The default is "tibble".

- digits:

  An integer value indicating the number of decimals to round the table
  values to.

- basemodel:

  the family of the baseline model for McFadden's Pseudo-R-Squared:
  "Poisson" (default) or a family of
  [`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

- variables:

  logical; if `TRUE`, the baseline model includes the mean-model
  covariates of each model; if `FALSE` (default), it is an
  intercept-only model (with the model's offset).

## Value

A table comparing the models supplied to the function in the format
specified The table includes coefficients, standard errors, statistical
significance, AIC, BIC, and McFadden's Pseudo-R-Squared. Table formats
include: "tibble", "gt", and "latex".

## Examples

``` r

# Comparing the NBP model with the NB2 and NB1 models
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)

nb.1 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
                    ShouldWidth04 + AADTover10k,
                    data=washington_roads, family = 'NB1', method = 'NM')
                    
nb.2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
                    ShouldWidth04 + AADTover10k,
                    data=washington_roads, family = 'NB2', method = 'NM')
                    
                    
nb.p <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
                    ShouldWidth04 + AADTover10k,
                    data=washington_roads, family = 'NBP', method = 'NM')
                    
                    
comptable <- 
 regCompTable(list("NB-1"=nb.1, "NB-2"=nb.2, "NB-P"=nb.p), tableType="latex")
print(comptable)
#> \begin{table}
#> 
#> \caption{Model Comparison Statistics}
#> \centering
#> \begin{tabular}[t]{llll}
#> \toprule
#> Parameter & NB-1 & NB-2 & NB-P\\
#> \midrule
#> (Intercept) & -7.661 (0.493)*** & -7.742 (0.489)*** & -7.796 (0.5)***\\
#> lnaadt & 0.924 (0.057)*** & 0.936 (0.057)*** & 0.942 (0.058)***\\
#> lnlength & 0.815 (0.066)*** & 0.84 (0.069)*** & 0.834 (0.069)***\\
#> speed50 & -0.37 (0.11)*** & -0.385 (0.109)*** & -0.387 (0.112)***\\
#> ShouldWidth04 & 0.256 (0.089)** & 0.254 (0.092)** & 0.259 (0.092)**\\
#> \addlinespace
#> AADTover10k & 0.672 (0.122)*** & 0.713 (0.134)*** & 0.686 (0.132)***\\
#> ln(alpha) & -1.707 (0.347)*** & -1.46 (0.32)*** & -1.327 (0.293)***\\
#> ln(p) & --- & --- & 0.502 (0.176)**\\
#> N Obs. & 1501 & 1501 & 1501\\
#> LL & -1065.03 & -1062.794 & -1062.266\\
#> \addlinespace
#> AIC & 2144.06 & 2139.589 & 2140.532\\
#> BIC & 2181.257 & 2176.786 & 2183.043\\
#> Pseudo-R-Sq. & 0.301 & 0.303 & 0.303\\
#> \bottomrule
#> \end{tabular}
#> \end{table}
```
