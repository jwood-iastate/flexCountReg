# Generate a covariance matrix using a correlation matrix and vector of standard deviations

Generate a covariance matrix using a correlation matrix and vector of
standard deviations

## Usage

``` r
cor2cov(C, S)
```

## Arguments

- C:

  A correlation matrix.

- S:

  A vector of standard deviations.

## Value

A covariance matrix

## Examples

``` r
C <- matrix(c(1,-0.3,0.7,-0.3,1,-0.2,0.7,-0.2,1), 3, 3)
S <- c(0.5, 2, 1.25)
cor2cov(C,S)
#>         [,1] [,2]    [,3]
#> [1,]  0.2500 -0.3  0.4375
#> [2,] -0.3000  4.0 -0.5000
#> [3,]  0.4375 -0.5  1.5625
```
