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
  digits = 3
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
  Pseudo-R-Squared statistic (comparing against a Poisson regression
  model) should be included.

- tableType:

  The type of table format to return. Options include "tibble" for
  returning the table as a tibble, "gt" for a
  [gt](https://gt.rstudio.com/reference/gt.html) table object, or
  "latex" for a latex table. The default is "tibble".

- digits:

  An integer value indicating the number of decimals to round the table
  values to.

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
#> (Intercept) & -7.811 (0.483)*** & -7.153 (0.5)*** & -7.796 (0.504)***\\
#> lnaadt & 0.942 (0.056)*** & 0.867 (0.059)*** & 0.942 (0.059)***\\
#> lnlength & 0.819 (0.066)*** & 0.836 (0.07)*** & 0.834 (0.07)***\\
#> speed50 & -0.368 (0.11)*** & -0.402 (0.11)*** & -0.387 (0.111)***\\
#> ShouldWidth04 & 0.257 (0.089)** & 0.23 (0.093)* & 0.259 (0.092)**\\
#> \addlinespace
#> AADTover10k & 0.653 (0.121)*** & 0.829 (0.14)*** & 0.686 (0.133)***\\
#> ln(alpha) & -1.707 (0.342)*** & -1.326 (0.289)*** & -1.327 (0.293)***\\
#> ln(p) & --- & --- & 0.502 (0.176)**\\
#> N Obs. & 1501 & 1501 & 1501\\
#> LL & -1065.084 & -1063.672 & -1062.266\\
#> \addlinespace
#> AIC & 2144.167 & 2141.345 & 2140.532\\
#> BIC & 2181.364 & 2178.542 & 2183.043\\
#> Pseudo-R-Sq. & 0.301 & 0.302 & 0.303\\
#> \bottomrule
#> \end{tabular}
#> \end{table}
```
