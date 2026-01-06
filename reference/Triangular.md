# Triangle Distribution

These functions provide density, distribution function, quantile
function, and random number generation for the Triangle Distribution,
specified by its mean, standard deviation, and optional lower and upper
bounds.

## Usage

``` r
dtri(x, mode = 0, sigma = 1, upper = NA, lower = NA, log = FALSE)

ptri(
  q,
  mode = 0,
  sigma = 1,
  upper = NA,
  lower = NA,
  lower.tail = TRUE,
  log.p = FALSE
)

qtri(p, mode = 0, sigma = 1, upper = NA, lower = NA)

rtri(n, mode = 0, sigma = 1, upper = NA, lower = NA)
```

## Arguments

- x:

  numeric value or a vector of values.

- mode:

  numeric value or vector of mode values for the distribution.

- sigma:

  single value or vector indicating both the positive and negative max
  differences from the mean (if the difference is the same).

- upper:

  single value or vector for the upper limit of the distribution (must
  be used with \`lower\`).

- lower:

  single value or vector for the lower limit of the distribution (must
  be used with \`upper\`).

- log:

  logical; if TRUE, probabilities p are given as log(p).

- q:

  quantile or a vector of quantiles.

- lower.tail:

  logical; if TRUE, probabilities p are \\P\[X\leq x\]\\ otherwise,
  \\P\[X\>x\]\\.

- log.p:

  logical; if TRUE, probabilities p are given as log(p).

- p:

  probability or a vector of probabilities.

- n:

  the number of random numbers to generate.

## Details

The Triangle Distribution is defined by three points: a (minimum), b
(maximum), and c (mode), where the density is zero outside the interval
\[a, b\], increases linearly from a to c, and decreases linearly from c
to b.

`dtri` computes the density (PDF) of the Triangle Distribution.

`ptri` computes the CDF of the Triangle Distribution.

`qtri` computes the quantile function of the Triangle Distribution.

`rtri` generates random numbers from the Triangle Distribution.

The mode and standard deviation parameters define the distribution's
location and scale, respectively, while the lower and upper bounds
explicitly set the minimum and maximum values of the distribution.

dtri gives the density, ptri gives the distribution function, qtri gives
the quantile function, and rtri generates random deviates.

The length of the result is determined by n for rtri, and is the maximum
of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the
result. Only the first elements of the logical arguments are used.

## Examples

``` r
dtri(4, mode=8, upper=13, lower=1)
#> [1] 0.07142857
ptri(c(0, 1, 2, 3, 5, 7, 9, 10), mode = 3, upper=9, lower = 1)
#> [1] 0.0000000 0.0000000 0.0625000 0.2500000 0.6666667 0.9166667 1.0000000
#> [8] 1.0000000
qtri(c(0.1, 0.3, 0.5, 0.9, 0.95), mode = 3, upper = 9, lower = 1)
#> [1] 2.264911 3.203449 4.101021 6.809110 7.450807
rtri(30, mode = 5, sigma = 3)
#>  [1] 3.346247 7.517553 6.051029 6.340004 7.138958 6.173909 5.818189 5.710741
#>  [9] 6.565418 6.103555 4.177710 7.140990 4.860410 3.535905 5.287116 5.806178
#> [17] 3.279292 5.607881 4.385498 4.784582 4.841973 5.938685 4.365143 4.641275
#> [25] 3.164350 6.724128 5.803815 5.024051 7.093128 5.451739
```
