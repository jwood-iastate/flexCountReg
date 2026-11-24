# Plot a two-variable smooth as a surface

Contours (over a color image) of a smooth of two variables, such as
`te(lnaadt, lnlength)`, masked where there are no data nearby.

## Usage

``` r
plot_surface(
  object,
  term = NULL,
  parameter = NULL,
  n = 400,
  too_far = 0.1,
  est = NULL,
  ...
)
```

## Arguments

- object:

  a fitted model with a smooth of two variables.

- term:

  the smooth's label (default: the first two-variable smooth).

- parameter:

  optional parameter label, when the same smooth appears in several
  components.

- n:

  grid size (as in
  [`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md)).

- too_far:

  grid points farther than this (on the unit square of the two
  variables) from every observation are not shown
  ([`mgcv::exclude.too.far()`](https://rdrr.io/pkg/mgcv/man/exclude.too.far.html));
  `0` shows the whole grid.

- est:

  optional output of
  [`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md)
  for the term.

- ...:

  further arguments passed to
  [`image()`](https://rdrr.io/r/graphics/image.html).

## Value

The plotted grid, invisibly.
