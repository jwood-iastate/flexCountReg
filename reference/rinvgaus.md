# Generate inverse-Gaussian random numbers

Generates random variables from an inverse-Gaussian distribution using
the Michael, Schucany, and Haas (1976) method.

## Usage

``` r
rinvgaus(n, mean = 1, shape = 1)
```

## Arguments

- n:

  Number of observations.

- mean:

  Mean parameter (\> 0).

- shape:

  Shape parameter (\> 0).

## Value

A numeric vector of length n.

## Details

Parameterization: IG(mean = mu, shape = lambda)
