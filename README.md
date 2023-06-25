
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesianPlatformDesignTimeTrend

<!-- badges: start -->
<!-- badges: end -->

The goal of BayesianPlatformDesignTimeTrend is to simulates the
multi-arm multi-stage or platform trial with Bayesian approach using the
‘rstan’ package, which provides the R interface for to the stan. The
package uses Thall’s and Trippa’s randomisation approach for Bayesian
adaptive randomisation. In addition, the time trend problem of platform
trial can be studied in this package. There is a demo for multi-arm
multi-stage trial for two different null scenario in this package.

## Installation

You can install the ‘BayesianPlatformDesignTimeTrend’ package 1.1.1 like
so:

``` r
# install.packages("BayesianPlatformDesignTimeTrend")
# devtools::install_github("ZXW834/PlatFormDesignTime", build_vignettes = TRUE)
```

## Demo

-   `Demo_CutoffScreening()` is a demo function performing cutoff
    screening process
-   `Demo_multiplescrenariotrialsimulation()` is a demo function
    performing MAMS trial simulation

## Tutorials

-   `MAMS-CutoffScreening-tutorial` is a tutorial document of how to do
    cutoff screening under Bayesian MAMS trial
-   `MAMS-trial-simulation-tutorial` is a tutorial document of how to do
    Bayesian MAMS trial simulation with or without time trend

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# library(BayesianPlatformDesignTimeTrend)
## basic example code
```

``` r
output=Trial.simulation(ntrials = 10000,
trial.fun = simulatetrial,
 input.info = list(
   response.probs = c(0.4, 0.4),
   ns = c(60, 120, 180, 240, 300),
   max.ar = 0.75,
   rand.algo = "Urn",
   max.deviation = 3,
   model.inf = list(
     model = "tlr",
     ibb.inf = list(
       pi.star = 0.5,
       pess = 2,
       betabinomialmodel = ibetabinomial.post
     ),
     tlr.inf = list(
       beta0_prior_mu = 0,
       beta1_prior_mu = 0,
       beta0_prior_sigma = 2.5,
       beta1_prior_sigma = 2.5,
       beta0_df = 7,
       beta1_df = 7,
       reg.inf =  "main",
       variable.inf = "Fixeffect"
     )
   ),
   Stopbound.inf = Stopboundinf(
     Stop.type = "Early-Pocock",
       Boundary.type = "Symmetric",
         cutoff = c(0.9925, 0.0075)
         ),
   Random.inf = list(
     Fixratio = FALSE,
     Fixratiocontrol = NA,
     BARmethod = "Thall",
     Thall.tuning.inf = list(tuningparameter = "Fixed",  fixvalue = 1)
   ),
   trend.inf = list(
     trend.type = "step",
     trend.effect = c(0, 0),
     trend_add_or_multip = "mult"
   )
 ),
 cl = 2)
```

Here is the operational characteristics table for previous single null
scenario simulation.

``` r
output$OPC
```

``` r
#>  $OPC
#>                          Type.I.Error.or.Power       Bias        rMSE     N.per.arm.1
#>0404TimeTrend00stage5main                0.0444  0.0007538   0.3390904        146.4978
#>                          N.per.arm.2    Survive.per.arm.1    Survive.per.arm.2          N
#>0404TimeTrend00stage5main    146.7282               58.552              58.6241    293.226
```
