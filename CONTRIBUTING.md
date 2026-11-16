# Contributing to flexCountReg

Bug reports and pull requests are welcome at
<https://github.com/jwood-iastate/flexCountReg/issues>. Please include a
minimal reproducible example (`reprex::reprex()`).

## Development workflow

- Code style: the tidyverse style guide, 80-character lines.
- Documentation: roxygen2 (run `devtools::document()`; the version in
  `DESCRIPTION` under `Config/roxygen2/version` is the one used to
  generate the `man/` files).
- Tests: testthat edition 3. Every bug fix gets a regression test that
  fails before the fix. Prefer tests that check results (closed forms,
  reference implementations, simulation recovery, invariances) over
  tests that only check that code runs. Mark slow stochastic tests with
  `skip_on_cran()`.
- Model families are defined once, in `R/families.R`; data preparation
  and prediction blueprints in `R/model_data.R`.

## Release checklist

1.  `devtools::document()` and `devtools::build_readme()`.
2.  `devtools::test()` with no failures and no unexpected warnings.
3.  `devtools::check(remote = TRUE, manual = TRUE)` with no errors,
    warnings or notes; `urlchecker::url_check()`;
    `devtools::spell_check()`.
4.  Check on the minimum R version in `DESCRIPTION`, on R-devel
    (`rhub::rhub_check()`), and with only the hard dependencies
    installed (the `no-suggests` CI job).
5.  Update `NEWS.md` (user-visible changes, including changed results)
    and the version number in `DESCRIPTION`.
6.  Build the pkgdown site and check the vignettes.
7.  `devtools::submit_cran()`; tag the release on GitHub.
