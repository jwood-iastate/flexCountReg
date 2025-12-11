# One-Parameter Lindley Distribution

Distribution function for the one-parameter Lindley distribution with
parameter theta.

## Usage

``` r
dlindley(x, theta = 1, log = FALSE)

plindley(q, theta = 1, lower.tail = TRUE, log.p = FALSE)

qlindley(p, theta = 1, lower.tail = TRUE, log.p = FALSE)

rlindley(n, theta = 1)
```

## Arguments

- x:

  a single value or vector of positive values.

- theta:

  distribution parameter value. Default is 1.

- log, log.p:

  logical; If TRUE, probabilities p are given as log(p). If FALSE,
  probabilities p are given directly. Default is FALSE.

- q:

  a single value or vector of quantiles.

- lower.tail:

  logical; If TRUE, (default), \\P(X \leq x)\\ are returned, otherwise
  \\P(X \> x)\\ is returned. Default is TRUE.

- p:

  a single value or vector of probabilities.

- n:

  number of random values to generate.

## Details

Probability density function (PDF) \$\$f(x\mid \theta )=\frac{\theta
^{2}}{(1+\theta )}(1+x)e^{-\theta x}\$\$

Cumulative distribution function (CDF) \$\$F(x\mid \theta ) = 1 -
\left(1+ \frac{\theta x}{1+\theta }\right)e^{-\theta x}\$\$

Quantile function (Inverse CDF) \$\$ Q(p\mid\theta) = -1 -
\frac{1}{\theta} - \frac{1}{\theta}
W\_{-1}\\\left((1+\theta)(p-1)e^{-(1+\theta)}\right) \$\$

where \\W\_{-1}()\\ is the negative branch of the Lambert W function.

The moment generating function (MGF) is:
\$\$M_X(t)=\frac{\theta^2(\theta-t+1)}{(\theta+1)(\theta-t)^2}\$\$

The distribution mean and variance are:
\$\$\mu=\frac{\theta+2}{\theta(1+\theta)}\$\$
\$\$\sigma^2=\frac{\mu}{\theta+2}\left(\frac{6}{\theta}-4\right)-\mu^2\$\$

## Examples

``` r
x <- seq(0, 5, by = 0.1)
p <- seq(0.1, 0.9, by = 0.1)
q <- c(0.2, 3, 0.2)
dlindley(x, theta = 1.5)
#>  [1] 0.000000000 0.852100897 0.800083678 0.746024937 0.691502661 0.637694846
#>  [7] 0.585460310 0.535404756 0.487934623 0.443300846 0.401634288 0.362974327
#> [13] 0.327291799 0.294507328 0.264505885 0.237148255 0.212280011 0.189738448
#> [19] 0.169357892 0.150973677 0.134425085 0.119557434 0.106223522 0.094284540
#> [25] 0.083610591 0.074080899 0.065583793 0.058016508 0.051284873 0.045302912
#> [31] 0.039992388 0.035282311 0.031108444 0.027412793 0.024143116 0.021252450
#> [37] 0.018698645 0.016443944 0.014454571 0.012700355 0.011154385 0.009792681
#> [43] 0.008593906 0.007539091 0.006611389 0.005795854 0.005079239 0.004449808
#> [49] 0.003897178 0.003412165 0.002986656
dlindley(x, theta=0.5, log=TRUE)
#>  [1]      -Inf -1.746449 -1.709438 -1.679395 -1.655287 -1.636294 -1.621756
#>  [8] -1.611131 -1.603973 -1.599906 -1.598612 -1.599822 -1.603302 -1.608850
#> [15] -1.616291 -1.625469 -1.636248 -1.648508 -1.662140 -1.677049 -1.693147
#> [22] -1.710357 -1.728609 -1.747837 -1.767984 -1.788997 -1.810826 -1.833427
#> [29] -1.856758 -1.880783 -1.905465 -1.930772 -1.956675 -1.983144 -2.010155
#> [36] -2.037682 -2.065703 -2.094197 -2.123144 -2.152524 -2.182322 -2.212519
#> [43] -2.243101 -2.274053 -2.305361 -2.337011 -2.368993 -2.401293 -2.433902
#> [50] -2.466807 -2.500000
plindley(q, theta = 1.5)
#> [1] 0.1702836 0.9688948 0.1702836
plindley(q, theta = 0.5, lower.tail = FALSE)
#> [1] 0.9651599 0.4462603 0.9651599
qlindley(p, theta = 1.5)
#> [1] 0.1145570 0.2376153 0.3721440 0.5222804 0.6942468 0.8982645 1.1532321
#> [8] 1.5010908 2.0740189
qlindley(p, theta = 0.5)
#> [1] 0.5440196 1.0430722 1.5435368 2.0718301 2.6536848 3.3240888 4.1429816
#> [8] 5.2395404 7.0163914

set.seed(123154)
rlindley(5, theta = 1.5)
#> Warning: 'n' must be a non-negative integer
#> [1] 1.40543787 0.97518942 0.42922042 0.06554152 0.24698987
rlindley(5, theta = 0.5)
#> Warning: 'n' must be a non-negative integer
#> [1] 2.5089347 2.5340820 1.4101438 0.1856275 3.6330299
```
