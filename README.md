
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesianPlatformDesignTimeTrend

<!-- badges: start -->
<!-- badges: end -->

The goal of BayesianPlatformDesignTimeTrend is to simulates the multi-arm
multi-stage or platform trial with Bayesian approach using the ‘rstan’
package, which provides the R interface for to the stan. The package
uses Thall’s and Trippa’s randomisation approach for Bayesian adaptive
randomisation. In addition, the time trend problem of platform trial can
be studied in this package. There is a demo for multi-arm multi-stage
trial for two different null scenario in this package.

## Installation

You can install the ‘BayesianPlatformDesignTimeTrend’ package 1.0.0 like so:

``` r
# install.packages("BayesianPlatformDesignTimeTrend")
devtools::install_github("ZXW834/PlatFormDesignTime", build_vignettes = TRUE)
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
library(BayesianPlatformDesignTimeTrend)
#> Loading required package: rstan
#> Loading required package: StanHeaders
#> Loading required package: ggplot2
#> rstan (Version 2.21.8, GitRev: 2e1f913d3ca3)
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores()).
#> To avoid recompilation of unchanged Stan programs, we recommend calling
#> rstan_options(auto_write = TRUE)
## basic example code
Trial.simulation (ntrials = 10,
trial.fun = simulatetrial,
 input.info = list(
   response.probs = c(0.4, 0.4),
   ns = c(30, 60, 90, 120, 150),
   max.ar = 0.75,
   rand.type = "Urn",
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
         cutoff = c(0.99, 0.01)
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
#> [1] "Start trial information initialisation"
#> $result
#> $result$`0404TimeTrend00stage5main`
#> $result$`0404TimeTrend00stage5main`[[1]]
#>     PP1C nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.9376 16  6  14   9        0    -0.481     1.071    0.525
#> 2 0.9792 24  9  36  24        0    -0.483     1.162    0.300
#> 3 0.9904 30 10  60  36        1    -0.699     1.099    0.212
#> 4     NA NA NA  NA  NA       NA        NA        NA       NA
#> 5     NA NA NA  NA  NA       NA        NA        NA       NA
#> 
#> $result$`0404TimeTrend00stage5main`[[2]]
#>    PP1C nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.002 14  8  16   1        1     0.203    -2.733      1.1
#> 2    NA NA NA  NA  NA       NA        NA        NA       NA
#> 3    NA NA NA  NA  NA       NA        NA        NA       NA
#> 4    NA NA NA  NA  NA       NA        NA        NA       NA
#> 5    NA NA NA  NA  NA       NA        NA        NA       NA
#> 
#> $result$`0404TimeTrend00stage5main`[[3]]
#>     PP1C  nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.0152  15  9  15   3        0     0.359    -1.688    0.669
#> 2 0.0868  38 17  22   6        0    -0.234    -0.765    0.316
#> 3 0.1012  60 22  30   7        0    -0.550    -0.652    0.254
#> 4 0.1320  82 30  38  10        0    -0.557    -0.481    0.181
#> 5 0.1060 104 40  46  13        0    -0.471    -0.467    0.145
#> 
#> $result$`0404TimeTrend00stage5main`[[4]]
#>     PP1C  nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.1184  15  8  15   5        0     0.121    -0.813    0.517
#> 2 0.0620  38 16  22   5        0    -0.339    -0.891    0.357
#> 3 0.0468  61 28  29   8        0    -0.171    -0.792    0.247
#> 4 0.0592  83 35  37  10        0    -0.329    -0.662    0.178
#> 5 0.0612 106 45  44  13        0    -0.304    -0.569    0.145
#> 
#> $result$`0404TimeTrend00stage5main`[[5]]
#>     PP1C nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.8824 14  4  16   8        0    -0.930     0.908    0.566
#> 2 0.4800 22  9  38  15        0    -0.398    -0.034    0.285
#> 3 0.7272 37 15  53  25        0    -0.392     0.262    0.192
#> 4 0.3760 45 20  75  31        0    -0.231    -0.126    0.144
#> 5 0.4952 64 27  86  36        0    -0.320    -0.007    0.111
#> 
#> $result$`0404TimeTrend00stage5main`[[6]]
#>     PP1C  nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.0192  15  7  15   2        0    -0.179    -1.665    0.811
#> 2 0.0320  37 17  23   5        0    -0.188    -1.107    0.366
#> 3 0.0292  59 25  31   7        0    -0.322    -0.922    0.248
#> 4 0.0292  81 37  39  11        0    -0.175    -0.758    0.172
#> 5 0.0276 103 47  47  14        0    -0.181    -0.683    0.136
#> 
#> $result$`0404TimeTrend00stage5main`[[7]]
#>     PP1C nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.4004 14  4  16   4        0    -0.954    -0.189    0.640
#> 2 0.3624 31 10  29   8        0    -0.768    -0.217    0.314
#> 3 0.6588 49 16  41  15        0    -0.730     0.171    0.189
#> 4 0.9252 59 17  61  25        0    -0.911     0.538    0.150
#> 5 0.8340 67 21  83  32        0    -0.800     0.335    0.117
#> 
#> $result$`0404TimeTrend00stage5main`[[8]]
#>     PP1C nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.2660 16 10  14   7        0     0.497    -0.471    0.511
#> 2 0.1836 37 17  23   8        0    -0.174    -0.470    0.281
#> 3 0.3804 60 26  30  12        0    -0.275    -0.142    0.207
#> 4 0.7688 78 30  42  19        0    -0.470     0.272    0.141
#> 5 0.4984 85 33  65  25        0    -0.460    -0.005    0.113
#> 
#> $result$`0404TimeTrend00stage5main`[[9]]
#>     PP1C nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.1852 15  4  15   2        0    -1.084    -0.832    0.879
#> 2 0.3760 38 12  22   6        0    -0.790    -0.200    0.356
#> 3 0.4276 57 20  33  11        0    -0.633    -0.075    0.206
#> 4 0.1996 74 28  46  14        0    -0.504    -0.337    0.156
#> 5 0.4492 96 35  54  19        0    -0.559    -0.049    0.125
#> 
#> $result$`0404TimeTrend00stage5main`[[10]]
#>     PP1C nC yC nE1 yE1 H1^1tpIE Intercept Trt1_Mean Trt1_Var
#> 1 0.9856 14  2  16   8        0    -1.781     1.724    0.766
#> 2 0.8308 21  5  39  14        0    -1.183     0.588    0.389
#> 3 0.9816 29  5  61  23        0    -1.575     1.062    0.280
#> 4 0.8676 37 10  83  31        0    -1.011     0.489    0.187
#> 5 0.9308 45 12 105  41        0    -1.027     0.574    0.154
#> 
#> 
#> 
#> $OPC
#>                           Type.I.Error.or.Power    Bias    rMSE N.per.arm.1
#> 0404TimeTrend00stage5main                   0.2 -0.2505 1.05761        71.4
#>                           N.per.arm.2 Survive.per.arm.1 Survive.per.arm.2   N
#> 0404TimeTrend00stage5main        60.6              27.8                23 132
#> 
#> $Nameofsaveddata
#> $Nameofsaveddata$nameTable
#> [1] "TABLENOTRENDEARLY-POCOCKSYMMETRICTHALL.RData"
#> 
#> $Nameofsaveddata$nameData
#> [1] "DATANOTRENDEARLY-POCOCKSYMMETRICTHALL.RData"
```
