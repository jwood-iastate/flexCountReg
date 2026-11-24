# Package index

## Count regression models

Estimate count models with flexible distributions, underreporting,
random parameters, random effects and fixed effects.

- [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  : Count regression models
- [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  : Random Parameters Count Regression Models
- [`countreg_fe()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_fe.md)
  : Fixed-effects count regression
- [`cre_test()`](https://jwood-iastate.github.io/flexCountReg/reference/cre_test.md)
  : Test of correlated random effects (random against fixed effects)
- [`count_families`](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)
  [`count-families`](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)
  : Count Families for countreg() and countreg_rp()
- [`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md)
  : Estimate a Random Effects Negative Binomial regression model
- [`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md)
  : Function for estimating a Random Effects Poisson-Lindley regression
  model
- [`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md)
  : Estimate a Random Effects Poisson Panel Model
- [`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)
  : Negative Multinomial Family Panel Regression (NM1, NM2, NM-P)

## Finite mixtures and latent classes

Mixtures of count distributions and latent-class count regression.

- [`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  [`countreg_lc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md)
  : Finite-mixture and latent-class count regression
- [`mixture_control()`](https://jwood-iastate.github.io/flexCountReg/reference/mixture_control.md)
  [`lc_control()`](https://jwood-iastate.github.io/flexCountReg/reference/mixture_control.md)
  : Control settings for finite-mixture and latent-class count models
- [`predict(`*`<countreg_mix>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_mix.md)
  : Predictions for finite-mixture and latent-class models
- [`summary(`*`<countreg_mix>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/summary.countreg_mix.md)
  [`print(`*`<summary.countreg_mix>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/summary.countreg_mix.md)
  : Summaries of finite-mixture and latent-class models
- [`component_probs()`](https://jwood-iastate.github.io/flexCountReg/reference/component_probs.md)
  [`class_probs()`](https://jwood-iastate.github.io/flexCountReg/reference/component_probs.md)
  : Component (class) membership probabilities

## Interactive app

- [`launch_app()`](https://jwood-iastate.github.io/flexCountReg/reference/launch_app.md)
  [`app_object()`](https://jwood-iastate.github.io/flexCountReg/reference/launch_app.md)
  : Interactive analysis app

## Working with fitted models

- [`summary(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/summary.flexCountReg.md)
  [`print(`*`<summary.flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/summary.flexCountReg.md)
  [`as_tibble(`*`<summary.flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/summary.flexCountReg.md)
  : Summarize a flexCountReg model
- [`coef(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`nobs(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`logLik(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`formula(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`terms(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`confint(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`fitted(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`residuals(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`simulate(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  [`print(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)
  : Standard methods for flexCountReg models
- [`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md)
  : Estimated smooth functions
- [`plot_smooths()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_smooths.md)
  : Plot estimated smooth functions
- [`vcov(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/vcov.flexCountReg.md)
  : Covariance matrix of the estimates
- [`estfun(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/estfun.flexCountReg.md)
  [`bread(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/estfun.flexCountReg.md)
  : Scores and bread for the sandwich package
- [`tidy(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-tidiers.md)
  [`glance(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-tidiers.md)
  [`augment(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-tidiers.md)
  : Tidy, glance and augment flexCountReg models
- [`predict(`*`<flexCountReg>`*`)`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md)
  : Predictions for flexCountReg models

## Prediction, effects and diagnostics

- [`margEffTable()`](https://jwood-iastate.github.io/flexCountReg/reference/margEffTable.md)
  : Marginal Effects, Elasticities, and Pseudo-Elasticities for
  flexCountReg Models
- [`tolerance_interval()`](https://jwood-iastate.github.io/flexCountReg/reference/tolerance_interval.md)
  : Tolerance Intervals for Count Regression Models
- [`cureplot()`](https://jwood-iastate.github.io/flexCountReg/reference/cureplot.md)
  : Cumulative Residuals (CURE) Plot for Count Models
- [`rootogram()`](https://jwood-iastate.github.io/flexCountReg/reference/rootogram.md)
  : Rootogram for a count regression model

## Model comparison and fit metrics

- [`regCompTest()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md)
  : Compare Regression Models with Likelihood Ratio Test, AIC, and BIC
- [`regCompTable()`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTable.md)
  : Create a Table Comparing Regression Models with AIC, BIC, and
  McFadden's Pseudo-R-Squared
- [`mae()`](https://jwood-iastate.github.io/flexCountReg/reference/mae.md)
  : Calculate Mean Absolute Error (MAE)
- [`rmse()`](https://jwood-iastate.github.io/flexCountReg/reference/rmse.md)
  : Calculate Root Mean Squared Error (RMSE)
- [`myAIC()`](https://jwood-iastate.github.io/flexCountReg/reference/myAIC.md)
  : Calculate Akaike Information Criterion (AIC)
- [`myBIC()`](https://jwood-iastate.github.io/flexCountReg/reference/myBIC.md)
  : Calculate Bayesian Information Criterion (BIC)

## tidymodels

A parsnip model specification with the flexCountReg engine.

- [`count_reg()`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md)
  : Count regression model specification for tidymodels
- [`count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/count_family.md)
  : Count distribution tuning parameter

## User-defined count distributions

Define a count distribution from its log probability function and use it
in any of the models.

- [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md)
  : Define a count distribution for the regression functions
- [`count_param()`](https://jwood-iastate.github.io/flexCountReg/reference/count_param.md)
  : A distribution parameter of a user-defined count family
- [`check_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/check_count_family.md)
  : Check a user-defined count family

## Count distributions

- [`dcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md)
  [`pcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md)
  [`qcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md)
  [`rcom()`](https://jwood-iastate.github.io/flexCountReg/reference/COMDistribution.md)
  : Conway-Maxwell-Poisson (COM) Distribution
- [`ddelaporte()`](https://jwood-iastate.github.io/flexCountReg/reference/Delaporte.md)
  [`pdelaporte()`](https://jwood-iastate.github.io/flexCountReg/reference/Delaporte.md)
  [`qdelaporte()`](https://jwood-iastate.github.io/flexCountReg/reference/Delaporte.md)
  [`rdelaporte()`](https://jwood-iastate.github.io/flexCountReg/reference/Delaporte.md)
  : Delaporte Distribution
- [`ddweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md)
  [`pdweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md)
  [`qdweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md)
  [`rdweib()`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md)
  : Discrete Weibull Distribution
- [`ddpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md)
  [`pdpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md)
  [`qdpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md)
  [`rdpo()`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md)
  : Double Poisson Distribution
- [`dgcount()`](https://jwood-iastate.github.io/flexCountReg/reference/GammaCount.md)
  [`pgcount()`](https://jwood-iastate.github.io/flexCountReg/reference/GammaCount.md)
  [`qgcount()`](https://jwood-iastate.github.io/flexCountReg/reference/GammaCount.md)
  [`rgcount()`](https://jwood-iastate.github.io/flexCountReg/reference/GammaCount.md)
  : Gamma Count Distribution
- [`dghermite()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md)
  [`pghermite()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md)
  [`qghermite()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md)
  [`rghermite()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md)
  : Generalized Hermite Distribution
- [`dhpois()`](https://jwood-iastate.github.io/flexCountReg/reference/HyperPoisson.md)
  [`phpois()`](https://jwood-iastate.github.io/flexCountReg/reference/HyperPoisson.md)
  [`qhpois()`](https://jwood-iastate.github.io/flexCountReg/reference/HyperPoisson.md)
  [`rhpois()`](https://jwood-iastate.github.io/flexCountReg/reference/HyperPoisson.md)
  : Hyper-Poisson Distribution
- [`dnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md)
  [`pnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md)
  [`qnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md)
  [`rnbge()`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md)
  : Negative Binomial-Generalized Exponential Distribution
- [`dpaeppli()`](https://jwood-iastate.github.io/flexCountReg/reference/PolyaAeppli.md)
  [`ppaeppli()`](https://jwood-iastate.github.io/flexCountReg/reference/PolyaAeppli.md)
  [`qpaeppli()`](https://jwood-iastate.github.io/flexCountReg/reference/PolyaAeppli.md)
  [`rpaeppli()`](https://jwood-iastate.github.io/flexCountReg/reference/PolyaAeppli.md)
  : Polya-Aeppli Distribution
- [`dgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md)
  [`pgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md)
  [`qgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md)
  [`rgwar()`](https://jwood-iastate.github.io/flexCountReg/reference/Generalized-Waring.md)
  : Generalized Waring Distribution
- [`dgp1()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson.md)
  [`pgp1()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson.md)
  [`qgp1()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson.md)
  [`rgp1()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson.md)
  : Generalized Poisson Version 1 Distribution
- [`dgp2()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson2.md)
  [`pgp2()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson2.md)
  [`qgp2()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson2.md)
  [`rgp2()`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedPoisson2.md)
  : Generalized Poisson Version 2 Distribution
- [`dplindGamma()`](https://jwood-iastate.github.io/flexCountReg/reference/NegativeBinomialLindley.md)
  [`pplindGamma()`](https://jwood-iastate.github.io/flexCountReg/reference/NegativeBinomialLindley.md)
  [`qplindGamma()`](https://jwood-iastate.github.io/flexCountReg/reference/NegativeBinomialLindley.md)
  [`rplindGamma()`](https://jwood-iastate.github.io/flexCountReg/reference/NegativeBinomialLindley.md)
  : Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution
- [`dpge()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md)
  [`ppge()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md)
  [`qpge()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md)
  [`rpge()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonGeneralizedExponential.md)
  : Poisson-Generalized-Exponential Distribution
- [`dpinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md)
  [`ppinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md)
  [`qpinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md)
  [`rpinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGamma.md)
  : Poisson-Inverse-Gamma Distribution
- [`dpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md)
  [`ppinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md)
  [`qpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md)
  [`rpinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonInverseGaussian.md)
  : Poisson-Inverse-Gaussian Distribution
- [`dplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md)
  [`pplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md)
  [`qplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md)
  [`rplind()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindley.md)
  : Poisson-Lindley Distribution
- [`dplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  [`pplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  [`qplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  [`rplindLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md)
  : Poisson-Lindley-Lognormal Distribution
- [`dpLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLognormal.md)
  [`ppLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLognormal.md)
  [`qpLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLognormal.md)
  [`rpLnorm()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLognormal.md)
  : Poisson-Lognormal Distribution
- [`dpoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md)
  [`ppoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md)
  [`qpoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md)
  [`rpoisweibull()`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonWeibull.md)
  : Poisson-Weibull Distribution Functions
- [`dsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  [`psichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  [`qsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  [`rsichel()`](https://jwood-iastate.github.io/flexCountReg/reference/SichelDistribution.md)
  : Sichel Distribution

## Continuous distributions and simulation tools

- [`dlindley()`](https://jwood-iastate.github.io/flexCountReg/reference/Lindley.md)
  [`plindley()`](https://jwood-iastate.github.io/flexCountReg/reference/Lindley.md)
  [`qlindley()`](https://jwood-iastate.github.io/flexCountReg/reference/Lindley.md)
  [`rlindley()`](https://jwood-iastate.github.io/flexCountReg/reference/Lindley.md)
  : One-Parameter Lindley Distribution
- [`dinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/invgamma.md)
  [`pinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/invgamma.md)
  [`qinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/invgamma.md)
  [`rinvgamma()`](https://jwood-iastate.github.io/flexCountReg/reference/invgamma.md)
  : Inverse Gamma Distribution
- [`rinvgaus()`](https://jwood-iastate.github.io/flexCountReg/reference/rinvgaus.md)
  : Generate inverse-Gaussian random numbers
- [`dtri()`](https://jwood-iastate.github.io/flexCountReg/reference/Triangular.md)
  [`ptri()`](https://jwood-iastate.github.io/flexCountReg/reference/Triangular.md)
  [`qtri()`](https://jwood-iastate.github.io/flexCountReg/reference/Triangular.md)
  [`rtri()`](https://jwood-iastate.github.io/flexCountReg/reference/Triangular.md)
  : Triangle Distribution
- [`moment_lognormal()`](https://jwood-iastate.github.io/flexCountReg/reference/moment_lognormal.md)
  : Raw Moment of a Lognormal Distribution
- [`halton_dists()`](https://jwood-iastate.github.io/flexCountReg/reference/halton_dists.md)
  : Generate pseudo-random draws from specified distributions using
  Halton draws
- [`corr_haltons()`](https://jwood-iastate.github.io/flexCountReg/reference/corr_haltons.md)
  : Generate Correlated Random Variables Using Halton or Scrambled
  Halton Draws
- [`cor2cov()`](https://jwood-iastate.github.io/flexCountReg/reference/cor2cov.md)
  : Generate a covariance matrix using a correlation matrix and vector
  of standard deviations

## Data

- [`washington_roads`](https://jwood-iastate.github.io/flexCountReg/reference/washington_roads.md)
  : Washington Road Crashes
- [`abalone_rings`](https://jwood-iastate.github.io/flexCountReg/reference/abalone_rings.md)
  : Abalone Shell-Ring Counts
- [`stream_salamanders`](https://jwood-iastate.github.io/flexCountReg/reference/stream_salamanders.md)
  : Stream Salamander Counts
