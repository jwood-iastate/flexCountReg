# Stream Salamander Counts

Repeated counts of stream salamanders in 23 headwater streams in
south-eastern Kentucky, USA, 11 of them affected by mountaintop removal
mining and valley filling, from Price et al. (2016). Counts are
zero-heavy and strongly overdispersed, and the repeated samples of each
site form groups, which makes the data a useful example for
zero-inflated, hurdle and grouped count models.

## Usage

``` r
stream_salamanders
```

## Format

A data frame with 644 rows (23 sites x 4 samples x 7 species groups) and
9 columns:

- site:

  Factor: stream site, `R-1` to `R-12` (reference streams) and `VF-1` to
  `VF-11` (valley-filled streams). Each site was sampled four times; use
  it as the grouping variable.

- mined:

  Factor with levels `no` and `yes`: whether the stream is affected by
  mountaintop removal mining and valley filling.

- cover:

  Amount of cover objects in the stream reach (site level; centred and
  scaled in the source).

- sample:

  Sample (visit) number, 1 to 4. Integer.

- DOP:

  Days since the last precipitation (centred and scaled).

- Wtemp:

  Water temperature (centred and scaled).

- DOY:

  Day of the year of the sample (centred and scaled).

- spp:

  Factor: species and life stage. `GP`: *Gyrinophilus porphyriticus*;
  `PR`: *Pseudotriton ruber*; `DM`: *Desmognathus monticola*; `EC-A`:
  *Eurycea cirrigera* adults; `EC-L`: *E. cirrigera* larvae; `DES-L`:
  *Desmognathus* spp. larvae; `DF`: *D. fuscus*.

- count:

  Number of salamanders of the group counted in the sample. Integer. The
  count response.

## Source

Price, S. J., Muncy, B. L., Bonner, S. J., Drayer, A. N. and Barton, C.
D. (2016). Data from: Effects of mountaintop removal mining and valley
filling on the occupancy and abundance of stream salamanders. Dryad.
[doi:10.5061/dryad.5m8f6](https://doi.org/10.5061/dryad.5m8f6) (CC0 1.0;
downloaded 2026-11-17).

## Details

Salamanders were counted in a 10-m reach of each stream on four visits.
The counts are observed numbers of animals, not detection-corrected
abundance: the original study modeled occupancy and abundance with an
explicit detection model, which ordinary zero-inflated or hurdle count
models do not reproduce. Results from the samples of one site are not
independent.

Preparation for this package (script `data-raw/stream_salamanders.R` in
the source repository) started from the original Dryad file
(`Price et al. JAPPL-2015-count-data.csv`, one row per site), checked
its published MD5 digest, and reshaped it to one row per site, sample
and group. Site labels were normalized (e.g. `"VF - 4"` to `"VF-4"`),
the 0/1 `Mined` indicator was stored as the factor `mined` and the count
columns were mapped to the group codes above. The covariates were
already centred and scaled in the source and are unchanged. The result
matches the `Salamanders` data of the glmmTMB package (version 1.1.14)
on every site, sample and group, except that here `mined` has `no` as
its first (reference) level and `site` is an unordered factor.

The Dryad deposit is released under the Creative Commons CC0 1.0
public-domain dedication. Please cite the study and the data deposit.

## References

Price, S. J., Muncy, B. L., Bonner, S. J., Drayer, A. N. and Barton, C.
D. (2016). Effects of mountaintop removal mining and valley filling on
the occupancy and abundance of stream salamanders. *Journal of Applied
Ecology*, 53(2), 459-468.
[doi:10.1111/1365-2664.12585](https://doi.org/10.1111/1365-2664.12585)

## Examples

``` r
data("stream_salamanders", package = "flexCountReg")
str(stream_salamanders)
#> 'data.frame':    644 obs. of  9 variables:
#>  $ site  : Factor w/ 23 levels "R-1","R-2","R-3",..: 13 14 15 1 2 3 4 5 6 7 ...
#>  $ mined : Factor w/ 2 levels "no","yes": 2 2 2 1 1 1 1 1 1 1 ...
#>  $ cover : num  -1.442 0.298 0.398 -0.448 0.597 ...
#>  $ sample: int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ DOP   : num  -0.596 -0.596 -1.191 0 0.596 ...
#>  $ Wtemp : num  -1.2294 0.0848 1.0142 -3.0234 -0.1443 ...
#>  $ DOY   : num  -1.497 -1.497 -1.294 -2.712 -0.687 ...
#>  $ spp   : Factor w/ 7 levels "GP","PR","DM",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ count : int  0 0 0 2 2 1 1 2 4 1 ...
mean(stream_salamanders$count == 0)          # share of zeros
#> [1] 0.6009317
with(stream_salamanders, tapply(count, mined, mean))
#>        no       yes 
#> 2.2648810 0.2954545 

# \donttest{
fit <- countreg(count ~ spp + mined, data = stream_salamanders,
                family = "NB2", method = "BFGS")
summary(fit)
#> Call:
#> countreg(formula = count ~ spp + mined, data = stream_salamanders, 
#>     family = "NB2", method = "BFGS")
#> 
#> Model: countreg | family: NB2 | optimizer: BFGS
#> Iterations: 42 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    0.578      0.166    3.49  0.00049 ***
#> sppPR         -1.228      0.291   -4.23  2.4e-05 ***
#> sppDM          0.404      0.230    1.76  0.07813 .  
#> sppEC-A       -0.671      0.257   -2.61  0.00900 ** 
#> sppEC-L        0.639      0.220    2.90  0.00369 ** 
#> sppDES-L       0.822      0.222    3.71  0.00021 ***
#> sppDF          0.360      0.235    1.53  0.12588    
#> minedyes      -2.038      0.148  -13.79  < 2e-16 ***
#> ln(alpha)      0.217      0.128    1.69  0.09079 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -837.653 | AIC: 1693.306 | BIC: 1733.515 
#> N: 644
#> LR test vs Poisson: chi-sq = 358.299 on 1 df, p-value = <2e-16 (boundary-corrected) 
# }
```
