# Predictions from multivariate count models

Predictions from multivariate count models

## Usage

``` r
# S3 method for class 'countreg_mv'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "variance", "prob", "joint", "conditional"),
  max_count = NULL,
  y = NULL,
  outcome = NULL,
  given = NULL,
  random = c("marginal", "posterior"),
  ...
)
```

## Arguments

- object:

  a
  [`countreg_mv()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mv.md)
  fit.

- newdata:

  optional new data; the estimation data by default.

- type:

  the prediction:

  `"response"` (default)

  :   a matrix of the expected counts, one column per outcome.

  `"link"`

  :   the linear predictors (log of the mean parameter, at the mean of
      the random parameters).

  `"variance"`

  :   the variances (including the variance from the random parameters).

  `"prob"`

  :   a list with one matrix per outcome of the probabilities of the
      counts `0:max_count`.

  `"joint"`

  :   the joint probability of the counts `y` of all outcomes, one value
      per row.

  `"conditional"`

  :   a matrix of the probabilities of the counts `0:max_count` of
      `outcome` given the counts `given` of other outcomes, such as the
      probability of a fatal crash at a site given its injury crashes.

- max_count:

  for `"prob"` and `"conditional"`: the largest count (default: the
  largest observed count of the outcome).

- y:

  for `"joint"`: the counts of every outcome, a vector named by outcome
  (the same for every row) or a matrix or data frame with one column per
  outcome.

- outcome, given:

  for `"conditional"`: the name of the outcome, and the counts of the
  outcomes conditioned on, named as `y` (they can be some of the other
  outcomes).

- random:

  for fits with random parameters: `"marginal"` (default) averages over
  their distribution (the prediction for a new unit); `"posterior"`
  averages over each unit's posterior distribution given its counts in
  the estimation data (Empirical Bayes; units not in the estimation data
  get the marginal prediction). `newdata` then needs the `group_var`
  column(s).

- ...:

  unused.

## Value

A matrix with one column per outcome (`"response"`, `"link"`,
`"variance"`), a list of matrices (`"prob"`), a vector (`"joint"`) or a
matrix with one column per count (`"conditional"`).

## Details

Without random parameters the outcomes are independent given the
covariates, so joint probabilities are products of the outcomes'
probabilities and conditional probabilities equal the marginal ones.
With random parameters shared by the equations they are not: the
probabilities are averages over the random parameters of products of the
outcomes' probabilities.
