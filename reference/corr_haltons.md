# Generate Correlated Random Variables Using Halton or Scrambled Halton Draws

This function generates `N` correlated random variables using Halton or
scrambled Halton draws. The function supports normal and truncated
normal distributions.

## Usage

``` r
corr_haltons(
  means,
  cholesky = NULL,
  stdev = NULL,
  correlations = NULL,
  hdraws = NULL,
  ndraws = 500,
  scrambled = FALSE,
  dist = "normal",
  lower = -Inf,
  upper = Inf
)
```

## Arguments

- means:

  A numeric vector of means for each variable.

- cholesky:

  A Cholesky decomposition matrix to introduce correlation.

- stdev:

  A numeric vector of standard deviations for each variable. If
  provided, the function will use these values instead of the Cholesky
  decomposition matrix (must also provide a correlation matrix if
  providing standard deviations). Default is NULL.

- correlations:

  A correlation matrix to introduce correlation. If provided, the
  function will use these values instead of the Cholesky decomposition
  matrix (must also provide standard deviations). Default is NULL.

- hdraws:

  A matrix of Halton or scrambled Halton draws. If provided, the
  function will use these draws instead of generating new ones. Default
  is NULL.

- ndraws:

  An integer specifying the number of values to simulate for each
  variable. Default is 500.

- scrambled:

  A logical value indicating whether to use scrambled Halton draws.
  Default is FALSE.

- dist:

  A character string specifying the distribution type. Options are
  "normal" and "truncated_normal". Default is "normal".

- lower:

  A numeric value specifying the lower bound for truncated normal
  distribution. Default is -Inf.

- upper:

  A numeric value specifying the upper bound for truncated normal
  distribution. Default is Inf.

## Value

A matrix with `N` columns and `ndraws` rows containing the simulated
values for the correlated random variables.

## Examples

``` r
# Define mean, correlation, and standard deviations
means <- c(3, 2, 0.9)
sdevs <- c(0.25,1.5,0.8)
CORR <- matrix(c(1, -0.3, 0.5, -0.3, 1, -0.2, 0.5, -0.2, 1), 3, 3)

# Create the Cholesky decomposition matrix and set values for ndraws, etc.
ndraws <- 5000
scrambled <- TRUE
dist <- "normal"

# simulated the data
simulated_data <- corr_haltons(means, stdev=sdevs, correlations=CORR,
                                ndraws=ndraws, scrambled=scrambled,
                                dist=dist)

# look at the mean, standard deviation, and correlation of the simulated data
apply(simulated_data, 2, mean)
#> [1] 2.9991697 2.0247787 0.9048085
apply(simulated_data, 2, sd)
#> [1] 0.2473552 1.4834057 0.7913246
cor(simulated_data)
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.2942589  0.5091853
#> [2,] -0.2942589  1.0000000 -0.1922064
#> [3,]  0.5091853 -0.1922064  1.0000000

# providing a cholesky decomposition matrix
dist <- "normal"
cholesky <- chol(cor2cov(CORR, sdevs))
simulated_data <- corr_haltons(means, cholesky=cholesky, ndraws=ndraws,
                                scrambled=scrambled, dist=dist)
apply(simulated_data, 2, mean)
#> [1] 2.9991697 2.0247787 0.9048085
apply(simulated_data, 2, sd)
#> [1] 0.2473552 1.4834057 0.7913246
cor(simulated_data)
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.2942589  0.5091853
#> [2,] -0.2942589  1.0000000 -0.1922064
#> [3,]  0.5091853 -0.1922064  1.0000000

# Truncated normal
dist <- "truncated_normal"
lower <- 0
upper <- 30
simulated_data <- corr_haltons(means, cholesky=cholesky, ndraws=ndraws,
                                scrambled=scrambled, dist=dist,
                                lower=lower, upper=upper)
apply(simulated_data, 2, mean)
#> [1] 2.069626 3.317236 2.333870
apply(simulated_data, 2, sd)
#> [1] 0.8399966 5.3108743 1.6579825
cor(simulated_data)
#>           [,1]        [,2]        [,3]
#> [1,] 1.0000000  0.13032493  0.39069072
#> [2,] 0.1303249  1.00000000 -0.03192843
#> [3,] 0.3906907 -0.03192843  1.00000000
```
