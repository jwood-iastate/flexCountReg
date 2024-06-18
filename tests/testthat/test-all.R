library(flexCountReg)

### Name: flexCountReg.predict
### Title: Function for generating predictions based on the random
###   parameters negative binomial with multiple optional methods
### Aliases: flexCountReg.predict

### ** Examples


## Poisson-Lindley Model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
poislind.mod <- flexCountReg(Animal ~ lnaadt + lnlength + speed50 +
                               ShouldWidth04 + AADTover10k,
                             data=washington_roads, dist="Poisson Lindley",
                             method="BHHH")
summary(poislind.mod)

pred <- flexCountReg.predict(poislind.mod, washington_roads)

hist(pred)




### Name: flexCountReg
### Title: Estimate flexible count regression models with and without
###   random parameters
### Aliases: flexCountReg

### ** Examples

## Poisson-Lindley Model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
poislind.mod <- flexCountReg(Animal ~ lnaadt + lnlength + speed50 +
                               ShouldWidth04 + AADTover10k,
                             data=washington_roads, dist="Poisson Lindley",
                             method="sann")
summary(poislind.mod)




### Name: Generalized-Waring
### Title: Generalized Waring Distribution
### Aliases: Generalized-Waring dgwar pgwar qgwar rgwar

### ** Examples

dgwar(0, mu=1, alpha=1, rho=3)
pgwar(c(0,1,2,3), mu=1, alpha=2, rho=3)
qgwar(c(0.1, 0.5, 0.9), mu=1, alpha=2, rho=3)
rgwar(10, mu=1, alpha=2, rho=3)




### Name: genWaring
### Title: Function for estimating a Generalized Waring regression model
### Aliases: genWaring

### ** Examples


# Generalized Waring Model
data("washington_roads")
genwaring.mod <- genWaring(Total_crashes ~ lnaadt + lnlength,
                           data=washington_roads,
                           method='BHHH')
summary(genwaring.mod)



### Name: invgamma
### Title: Inverse Gamma Distribution
### Aliases: invgamma dinvgamma pinvgamma qinvgamma rinvgamma

### ** Examples

dinvgamma(1, shape=3, scale=2)
pinvgamma(c(0.1, 0.5, 1, 3, 5, 10, 30), shape=3, scale=2)
qinvgamma(c(0.1,0.3,0.5,0.9,0.95), shape=3, scale=2)
rinvgamma(30, shape=3, scale=2)




### Name: mae
### Title: Calculate Mean Absolute Error (MAE)
### Aliases: mae

### ** Examples

y <- c(1, 2, 3)
mu <- c(1.1, 1.9, 3.2)
mae(y, mu)




### Name: mgf_lognormal
### Title: Moment Generating Function for a Lognormal Distribution
### Aliases: mgf_lognormal

### ** Examples

mu <- 0
sigma <- 1
n <- 1
mgf_value <- mgf_lognormal(mu, sigma, n)
print(mgf_value)


### Name: myAIC
### Title: Calculate Akaike Information Criterion (AIC)
### Aliases: myAIC

### ** Examples

LL <- -120.5
nparam <- 5
myAIC(LL, nparam)




### Name: myBIC
### Title: Calculate Bayesian Information Criterion (BIC)
### Aliases: myBIC

### ** Examples

LL <- -120.5
nparam <- 5
n <- 100
myBIC(LL, nparam, n)




### Name: nbg
### Title: Function for estimating a variety of negative binomial models
###   (NB-1, NB-2, NB-P and generalized versions of each)
### Aliases: nbg

### ** Examples


## NB-P model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)

nbp.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
                  ShouldWidth04 + AADTover10k,
                data=washington_roads, form = 'nbp', method = 'BHHH',
                max.iters=3000)
summary(nbp.base)

## Generalized NB-P model

nbp.overdispersion <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
                            ShouldWidth04 + AADTover10k,
                          data=washington_roads,
                          form = 'nbp',
                          method = 'NM',
                          max.iters=3000,
                          ln.alpha.formula = ~ 1+lnlength)
summary(nbp.overdispersion)

## NB-1 Model
nb1.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
                  ShouldWidth04 + AADTover10k,
                data=washington_roads, form = 'nb1',
                method = 'NM',
                max.iters=3000)
summary(nb1.base)

## Generalize NB-1 Model
nb1.overdispersion <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
                            ShouldWidth04 + AADTover10k,
                          data=washington_roads, form = 'nb1',
                          method = 'NM',
                          max.iters=3000, ln.alpha.formula = ~ 1+lnlength)
summary(nb1.overdispersion)

## NB-2 Model
nb2.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
                  ShouldWidth04 + AADTover10k,
                data=washington_roads, form = 'nb2',
                method = 'BFGS',
                max.iters=3000)
summary(nb2.base)

## Generalize NB-2 Model
nb2.overdispersion <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
                            ShouldWidth04 + AADTover10k,
                          data=washington_roads, form = 'nb2',
                          method = 'NM',
                          max.iters=3000, ln.alpha.formula = ~ 1+lnlength)
summary(nb2.overdispersion)




### Name: Negative-Binomial-Lindley
### Title: Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution
### Aliases: Negative-Binomial-Lindley dplindGamma pplindGamma qplindGamma
###   rplindGamma

### ** Examples

dplindGamma(0, mean=0.75, theta=7, alpha=2, ndraws=100)
pplindGamma(c(0,1,2,3,5,7,9,10), mean=0.75, theta=7, alpha=2, ndraws=100)
qplindGamma(c(0.1,0.3,0.5,0.9,0.95), lambda=4.67, theta=7, alpha=2,
            ndraws=100)
rplindGamma(1, lambda=4.67, theta=7, alpha=2, ndraws=100)




### Name: poisGE
### Title: Poisson-Generalized-Exponential Regression
### Aliases: poisGE

### ** Examples

# Generalized Poisson-Generalized-Exponential
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)

poisge.mod <- poisGE(Total_crashes ~ lnaadt + lnlength + speed50 +
                       ShouldWidth04 + AADTover10k,
                     ln.scale.formula = ~ lnaadt,
                     data=washington_roads[1:500,], 
                     ndraws = 10,
                     max.iters = 500)
summary(poisge.mod)




### Name: poisInvGaus
### Title: Function for estimating a Poisson-Inverse-Gaussian regression
###   model
### Aliases: poisInvGaus

### ** Examples

# Poisson-Inverse-Gaussian Model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable

poisinvgaus.mod <- poisInvGaus(Total_crashes ~ lnaadt + lnlength + speed50 + 
                                 ShouldWidth04 + AADTover10k,
                               data=washington_roads,
                               method="nm",
                               max.iters = 1000)
summary(poisinvgaus.mod)



### Name: poisLind
### Title: Function for estimating a Poisson-Lindley regression model
### Aliases: poisLind

### ** Examples

# Poisson-Lindley Model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
poislind.mod <- poisLind(Animal ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
                         data=washington_roads,
                         method="nm",
                         max.iters = 1000)
summary(poislind.mod)



### Name: poisLindGamma
### Title: Function for estimating a Poisson-Lindley-Gamma (i.e., Negative
###   Binomial-Lindley) regression model
### Aliases: poisLindGamma

### ** Examples


## Poisson-Lindley-Gamma Model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
poislindgamma.mod <- poisLindGamma(Animal ~ lnaadt + lnlength + speed50 +
                                     ShouldWidth04 + AADTover10k,
                                   data=washington_roads,
                                   ndraws=10,
                                   max.iters = 1000)
summary(poislindgamma.mod)



### Name: poisLindLnorm
### Title: Function for estimating a Poisson-Lindley-Lognormal regression
###   model
### Aliases: poisLindLnorm

### ** Examples


## Poisson-Lindley-Lognormal Model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
poislindlnorm.mod <- poisLindLnorm(Animal ~ lnaadt + lnlength + speed50 
                                   + ShouldWidth04 + AADTover10k,
                                   data=washington_roads,
                                   method="nm",
                                   ndraws=10,
                                   max.iters = 1000)
summary(poislindlnorm.mod)



### Name: poisLogn
### Title: Poisson-Lognormal Regression
### Aliases: poisLogn

### ** Examples


# Generalized Poisson-Lognormal
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)

poslogn.mod <- poisLogn(Total_crashes ~ lnaadt + lnlength + speed50 +
                          ShouldWidth04 + AADTover10k,
                        ln.sigma.formula = ~ lnaadt + AADTover10k,
                        data=washington_roads, 
                        ndraws = 10, 
                        method = 'BHHH')
summary(poslogn.mod)




### Name: Poisson-Generalized-Exponential
### Title: Poisson-Generalized-Exponential Distribution
### Aliases: Poisson-Generalized-Exponential dpge ppge qpge rpge

### ** Examples

dpge(0, mean=0.75, shape=2, scale=1, ndraws=100)
ppge(c(0,1,2,3,4,5,6), mean=0.75, shape=2, scale=1, ndraws=100)
qpge(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, shape=2, scale=1, ndraws=100)
rpge(30, mean=0.75,  shape=2, scale=1, ndraws=100)




### Name: Poisson-Inverse-Gaussian
### Title: Poisson-Inverse-Gaussian Distribution
### Aliases: Poisson-Inverse-Gaussian dpinvgaus ppinvgaus qpinvgaus
###   rpinvgaus

### ** Examples

dpinvgaus(1, mu=0.75, eta=1)
ppinvgaus(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3, form="Type 2")
qpinvgaus(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5, form="Type 2")
rpinvgaus(30, mu=0.75, eta=1.5)

dpinvgaus(1, mu=0.75, eta=1)
ppinvgaus(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3, form="Type 2")
qpinvgaus(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5, form="Type 2")
rpinvgaus(30, mu=0.75, eta=1.5)




### Name: Poisson-Lindley-Lognormal
### Title: Poisson-Lindley-Lognormal Distribution
### Aliases: Poisson-Lindley-Lognormal dplindLnorm pplindLnorm qplindLnorm
###   rplindLnorm

### ** Examples

dplindLnorm(0, mean=0.75, theta=7, sigma=2, ndraws=100)
pplindLnorm(c(0,1,2,3,5,7,9,10), mean=0.75, theta=7, sigma=2, ndraws=100)
qplindLnorm(c(0.1,0.3,0.5,0.9,0.95), lambda=4.67, theta=7, sigma=2, ndraws=100)
rplindLnorm(30, mean=0.75, theta=7, sigma=2, ndraws=100)




### Name: Poisson-Lindley
### Title: Poisson-Lindley Distribution
### Aliases: Poisson-Lindley dplind pplind qplind rplind

### ** Examples

dplind(0, mean=0.75, theta=7)
pplind(c(0,1,2,3,5,7,9,10), mean=0.75, theta=7)
qplind(c(0.1,0.3,0.5,0.9,0.95), lambda=4.67, theta=7)
rplind(30, mean=0.75, theta=7)




### Name: Poisson-Lognormal
### Title: Poisson-Lognormal Distribution
### Aliases: Poisson-Lognormal dpLnorm ppLnorm qpLnorm rpLnorm

### ** Examples

dpLnorm(0, mean=0.75, sigma=2, ndraws=100)
ppLnorm(c(0,1,2,3,5,7,9,10), mean=0.75, sigma=2, ndraws=100)
qpLnorm(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, sigma=2, ndraws=100)
rpLnorm(30, mean=0.75,  sigma=2, ndraws=100)




### Name: PoissonWeibull
### Title: Poisson-Weibull Distribution Functions
### Aliases: PoissonWeibull dpoisweibull ppoisweibull qpoisweibull
###   rpoisweibull

### ** Examples

dpoisweibull(4, alpha = 1.5, beta = 0.5)
ppoisweibull(4, alpha = 1.5, beta = 0.5)
qpoisweibull(0.95, alpha = 1.5, beta = 0.5)
rpoisweibull(10, alpha = 1.5, beta = 0.5)




### Name: pwiebreg
### Title: Poisson-Weibull Regression with Optional Random Parameters
### Aliases: pwiebreg

### ** Examples

data("washington_roads")
pw_rp <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
                  rpar_formula = ~ speed50,
                  alpha_formula = ~ lnlength,
                  beta_formula = ~ lnaadt,
                  data = washington_roads,
                  ndraws = 10,
                  correlated = FALSE)
print(summary(pw_rp))




### Name: regCompTest
### Title: Compare Regression Models with Likelihood Ratio Test, AIC, and
###   BIC
### Aliases: regCompTest

### ** Examples


# Comparing the NBP model with the NB2 model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)

nbp.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
                  ShouldWidth04 + AADTover10k,
                data=washington_roads, form = 'nbp', method = 'NM',
                max.iters=3000)
comptests <- regCompTest(nbp.base, washington_roads, basemodel="NB2", print=TRUE)




### Name: rmse
### Title: Calculate Root Mean Squared Error (RMSE)
### Aliases: rmse

### ** Examples

y <- c(1, 2, 3)
mu <- c(1.1, 1.9, 3.2)
rmse(y, mu)




### Name: rpnb.predict
### Title: Function for generating predictions based on the random
###   parameters negative binomial with multiple optional methods
### Aliases: rpnb.predict

### ** Examples

## No test: 

## Random Parameters Negative Binomial model (NB-2)

data("washington_roads")

nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
               rpar_formula = ~ speed50,
               data = washington_roads,
               ndraws = 10,
               correlated = TRUE,
               form = 'nb2',
               method = "bfgs",
               print.level = 1)

## Exact Prediction
hist(rpnb.predict(nb2.rp, washington_roads))

## Simulated Prediction
hist(rpnb.predict(nb2.rp, washington_roads, method="Simulated"))

## Individual-Specific Coefficient Based Prediction
hist(rpnb.predict(nb2.rp, washington_roads, method="Individual"))
## End(No test)



### Name: rpnb
### Title: Function for estimating a random parameter negative binomial
###   with the ability to specify if the NB-1, NB-2, or NB-P should be used
### Aliases: rpnb

### ** Examples

## Not run: 
##D ## Random Parameters Negative Binomial model (NB-2)
##D data("washington_roads")
##D nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
##D                rpar_formula = ~ speed50,
##D                data = washington_roads,
##D                ndraws = 10,
##D                correlated = TRUE,
##D                form = 'nb2',
##D                method = "nr",
##D                print.level = 1)
##D 
##D summary(nb2.rp)
## End(Not run)



### Name: rppoisson
### Title: Random Parameters Poisson Model
### Aliases: rppoisson

### ** Examples

## Not run: 
##D data("washington_roads")
##D poisson_rp <- rppoisson(Total_crashes ~ lnlength + lnaadt,
##D                         rpar_formula = ~ speed50,
##D                         data = washington_roads,
##D                         ndraws = 10,
##D                         correlated = FALSE)
##D print(summary(poisson_rp))
## End(Not run)



### Name: Sichel-Distribution
### Title: Sichel Distribution
### Aliases: Sichel-Distribution dsichel psichel qsichel rsichel

### ** Examples

dsichel(1, mu=0.75, sigma=1, gamma=-3)
psichel(c(0,1,2,3,5,7,9,10), mu=0.75, sigma=2, gamma=3)
qsichel(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, sigma=1, gamma=15)
rsichel(30, mu=0.75, sigma=0.5, gamma=1)




### Name: sichel
### Title: Function for estimating a Sichel regression model
### Aliases: sichel

### ** Examples


# Sichel Model
data("washington_roads")

sichel.mod <- sichel(Total_crashes ~ lnaadt + lnlength,
                     data=washington_roads,
                     method="NM",
                     max.iters = 1000)
summary(sichel.mod)



### Name: summary.boot
### Title: Custom summary method for models with bootstrapped standard
###   errors
### Aliases: summary.boot

### ** Examples

# Poisson-Weibull
pw_rp <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
                  data = washington_roads,
                  ndraws = 10,
                  bootstraps = 10)
summary.boot(pw_rp)




### Name: summary.maxLik
### Title: Custom summary method for maxLik objects
### Aliases: summary.maxLik

### ** Examples

# NB2 Model
data("washington_roads")
washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)

nb2.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
                  ShouldWidth04 + AADTover10k,
                data=washington_roads, form = 'nb2',
                method = 'NR',
                max.iters=3000)

summary.maxLik(nb2.base)           




### Name: Triangular
### Title: Triangle Distribution
### Aliases: Triangular dtri ptri qtri rtri

### ** Examples

dtri(4, mode=8, upper=13, lower=1)
ptri(c(0, 1, 2, 3, 5, 7, 9, 10), mode = 3, upper=9, lower = 1)
qtri(c(0.1, 0.3, 0.5, 0.9, 0.95), mode = 3, upper = 9, lower = 1)
rtri(30, mode = 5, sigma = 3)
